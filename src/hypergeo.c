// SPDX-License-Identifier: AGPL-3.0-or-later
/**
 * This file is part of kycg.
 *
 * Copyright (C) 2026-present Wanding Zhou
 *
 * kycg is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * kycg is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with kycg.  If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * Log-space hypergeometric tail probabilities and Benjamini-Hochberg
 * correction.
 *
 * GOAL
 *   Reproduce, in C and to near machine precision, the two numerical
 *   operations knowYourCG performs in R on top of YAME's counts:
 *
 *     log10.p.value <- phyper(...) / log(10)      [testEnrichment.R:252]
 *     FDR           <- p.adjust(p.value, "fdr")   [testEnrichment.R:179]
 *
 * WHY NOT REUSE htslib
 *   YAME vendors htslib, which carries kt_fisher_exact() in htslib/kfunc.c.
 *   It is unusable here on two counts. It takes `int` counts, and our universe
 *   is ~29M rows with products that overflow well before that; and it
 *   accumulates tail mass in *linear* space, looping "until underflow" by
 *   construction. The results users care about are exactly the ones whose
 *   p-values underflow to zero in double precision, so a linear accumulator
 *   returns 0 for every interesting row. Everything here stays logarithmic.
 *
 * WHY THE SADDLE-POINT EXPANSION
 *   The obvious implementation of the point mass is a difference of three
 *   lgamma() calls. For a 29M-row universe those terms are each ~5e8, so
 *   cancellation leaves roughly 1e-7 absolute error in the natural log — which
 *   is far coarser than the ~1e-10 agreement with knowYourCG that is the
 *   acceptance bar for `kycg test`.
 *
 *   So we use Catherine Loader's saddle-point formulation (Loader 2000,
 *   "Fast and Accurate Computation of Binomial Probabilities"), which is what
 *   R's own dhyper()/dbinom() use internally. It reorganizes the point mass
 *   into a Stirling-series remainder (stirlerr) plus a relative-entropy term
 *   (bd0), neither of which suffers the cancellation. Matching R's algorithm
 *   rather than merely R's formula is what makes agreement to ~1e-14 relative
 *   achievable, and it is why the validation harness can use a tight bar.
 *
 * TAIL SUMMATION
 *   The tail is summed term by term, with successive terms produced by the
 *   exact ratio recurrence
 *
 *     P(x+1)/P(x) = ((r-x)(n-x)) / ((x+1)(b-n+x+1))
 *
 *   and accumulated with logaddexp. The hypergeometric pmf is unimodal, so
 *   once the ratio drops below 1 the remainder is bounded by the geometric
 *   series term/(1-ratio); we stop when that bound falls 50 nats below the
 *   running total, which is well past the point where it could perturb the
 *   result.
 *
 *   We always sum the tail that does *not* contain the mode, complementing
 *   with log1mexp when the caller asked for the other one. Summing the larger
 *   tail directly would walk across the entire support — hundreds of
 *   thousands of terms at whole-genome scale — and the accumulated rounding
 *   is enough to push the total above 1, which surfaces as a positive log
 *   p-value. Anchoring at the mode keeps every summation short, so the cost
 *   is bounded by the spread of the distribution rather than by its support.
 *
 * All external entry points return log10, matching knowYourCG's convention of
 * carrying log10.p.value as the primary quantity and recovering p.value as
 * 10^log10.p.value (which is allowed to underflow to 0 — it is a display
 * column, never a sort key or an input to FDR).
 */

#include "hypergeo.h"

#include <math.h>
#include <float.h>
#include <stdlib.h>

#ifndef M_LN10
#define M_LN10 2.30258509299404568402
#endif

#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

/* log(sqrt(2*pi)) and log(2*pi) */
#define KYCG_LN_SQRT_2PI 0.918938533204672741780329736406
#define KYCG_LN_2PI      1.837877066409345483560659472811

/* ---------------------------------------------------------------- helpers */

static inline double logaddexp(double a, double b) {
  if (a == -INFINITY) return b;
  if (b == -INFINITY) return a;
  double hi = a > b ? a : b;
  double lo = a > b ? b : a;
  return hi + log1p(exp(lo - hi));
}

/**
 * stirlerr(n) = log(n!) - log(sqrt(2*pi*n) * (n/e)^n)
 *
 * The error term of Stirling's approximation. Exact tabulated values for
 * half-integers up to 15 (where the asymptotic series is not yet accurate),
 * then a truncated series whose depth is chosen by magnitude. Transcribed
 * from R's src/nmath/stirlerr.c.
 */
static double stirlerr(double n) {
  static const double S0 = 0.083333333333333333333;        /* 1/12    */
  static const double S1 = 0.00277777777777777777778;      /* 1/360   */
  static const double S2 = 0.00079365079365079365079365;   /* 1/1260  */
  static const double S3 = 0.000595238095238095238095238;  /* 1/1680  */
  static const double S4 = 0.0008417508417508417508417508; /* 1/1188  */

  static const double sferr_halves[31] = {
    0.0,                            /* n=0 - placeholder, never indexed */
    0.1534264097200273452913848,    /* 0.5  */
    0.0810614667953272582196702,    /* 1.0  */
    0.0548141210519176538961390,    /* 1.5  */
    0.0413406959554092940938221,    /* 2.0  */
    0.03316287351993628748511048,   /* 2.5  */
    0.02767792568499833914878929,   /* 3.0  */
    0.02374616365629749597132920,   /* 3.5  */
    0.02079067210376509311152277,   /* 4.0  */
    0.01848845053267318523077934,   /* 4.5  */
    0.01664469118982119216319487,   /* 5.0  */
    0.01513497322191737887351255,   /* 5.5  */
    0.01387612882307074799874573,   /* 6.0  */
    0.01281046524292022692424986,   /* 6.5  */
    0.01189670994589177009505572,   /* 7.0  */
    0.01110455975820691732662991,   /* 7.5  */
    0.010411265261972096497478567,  /* 8.0  */
    0.009799416126158803298389475,  /* 8.5  */
    0.009255462182712732917728637,  /* 9.0  */
    0.008768700134139385462952823,  /* 9.5  */
    0.008330563433362871256469318,  /* 10.0 */
    0.007934114564314020547248100,  /* 10.5 */
    0.007573675487951840794972024,  /* 11.0 */
    0.007244554301320383179543912,  /* 11.5 */
    0.006942840107209529865664152,  /* 12.0 */
    0.006665247032707682442354394,  /* 12.5 */
    0.006408994188004207068439631,  /* 13.0 */
    0.006171712263039457647532867,  /* 13.5 */
    0.005951370112758847735624416,  /* 14.0 */
    0.005746216513010115682023589,  /* 14.5 */
    0.005554733551962801371038690   /* 15.0 */
  };

  if (n <= 15.0) {
    double nn = n + n;
    /* Bound the table index explicitly. Every caller today passes a positive
     * n, so a negative index is unreachable -- but this is the only unchecked
     * array index in the file, and the guard above (n <= 15.0) does not stop
     * one on its own. R's nmath has the same shape and the same latent hole. */
    if (nn >= 0.0 && nn <= 30.0 && nn == (double)(int)nn)
      return sferr_halves[(int)nn];
    return lgamma(n + 1.0) - (n + 0.5) * log(n) + n - KYCG_LN_SQRT_2PI;
  }

  double nn = n * n;
  if (n > 500.0) return (S0 - S1 / nn) / n;
  if (n > 80.0)  return (S0 - (S1 - S2 / nn) / nn) / n;
  if (n > 35.0)  return (S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n;
  return (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n;
}

/**
 * bd0(x, np) = x*log(x/np) + np - x
 *
 * The "deviance part": x times the Kullback-Leibler divergence between x and
 * np. Written directly it cancels catastrophically when x is close to np,
 * which is precisely the common case, so that regime uses the odd-power
 * Taylor series in v = (x-np)/(x+np) instead. From R's src/nmath/bd0.c.
 */
static double bd0(double x, double np) {
  if (!isfinite(x) || !isfinite(np) || np == 0.0) return NAN;

  if (fabs(x - np) < 0.1 * (x + np)) {
    double v = (x - np) / (x + np);
    double s = (x - np) * v;
    if (fabs(s) < DBL_MIN) return s;
    double ej = 2.0 * x * v;
    v = v * v;
    for (int j = 1; j < 1000; ++j) {
      ej *= v;                        /* v^(2j+1) */
      double s1 = s + ej / ((j << 1) + 1);
      if (s1 == s) return s1;         /* term vanished */
      s = s1;
    }
  }
  return x * log(x / np) + np - x;
}

/** Natural log of the binomial point mass, Loader-style. */
static double ldbinom_raw(double x, double n, double p, double q) {
  if (p == 0.0) return (x == 0.0) ? 0.0 : -INFINITY;
  if (q == 0.0) return (x == n)   ? 0.0 : -INFINITY;

  if (x == 0.0) {
    if (n == 0.0) return 0.0;
    return (p < 0.1) ? -bd0(n, n * q) - n * p : n * log(q);
  }
  if (x == n) {
    return (q < 0.1) ? -bd0(n, n * p) - n * q : n * log(p);
  }
  if (x < 0.0 || x > n) return -INFINITY;

  double lc = stirlerr(n) - stirlerr(x) - stirlerr(n - x)
            - bd0(x, n * p) - bd0(n - x, n * q);
  double lf = KYCG_LN_2PI + log(x) + log1p(-x / n);
  return lc - 0.5 * lf;
}

/** Natural log of the hypergeometric point mass. Mirrors R's dhyper(). */
static double ldhyper_nat(double x, double r, double b, double n) {
  if (x < 0.0 || x > n || x > r || n - x > b) return -INFINITY;
  if (n == 0.0) return (x == 0.0) ? 0.0 : -INFINITY;

  double N = r + b;
  double p = n / N;
  double q = (N - n) / N;

  return ldbinom_raw(x, r, p, q)
       + ldbinom_raw(n - x, b, p, q)
       - ldbinom_raw(n, N, p, q);
}

/* Stop summing once the bounded remainder is this many nats below the sum. */
#define KYCG_TAIL_CUTOFF 50.0

/**
 * log(1 - exp(l)) for l <= 0, without cancellation at either end.
 * The crossover at -ln(2) is where expm1 and log1p swap which one is stable.
 */
static double log1mexp(double l) {
  if (l == -INFINITY) return 0.0;      /* 1 - 0 */
  if (l >= 0.0) return -INFINITY;      /* 1 - 1, clamped */
  if (l > -M_LN2) return log(-expm1(l));
  return log1p(-exp(l));
}

/** Most probable value of X; the summation always starts from this side. */
static double hyper_mode(double r, double b, double n) {
  return floor((n + 1.0) * (r + 1.0) / (r + b + 2.0));
}

/** Raw upward sum: natural log of sum_{x >= x0} P(X = x). */
static double lsum_up(double x0, double r, double b, double n) {
  double hi = fmin(n, r);
  if (x0 > hi) return -INFINITY;

  double lterm = ldhyper_nat(x0, r, b, n);
  double lsum = lterm;

  for (double x = x0; x < hi; x += 1.0) {
    double num = (r - x) * (n - x);
    double den = (x + 1.0) * (b - n + x + 1.0);
    if (num <= 0.0 || den <= 0.0) break;
    double lratio = log(num) - log(den);
    lterm += lratio;
    lsum = logaddexp(lsum, lterm);
    if (lratio < 0.0) {
      /* Past the mode: remaining mass < lterm/(1-ratio). */
      double lrem = lterm - log1p(-exp(lratio));
      if (lrem < lsum - KYCG_TAIL_CUTOFF) break;
    }
  }
  return lsum;
}

/** Raw downward sum: natural log of sum_{x <= x0} P(X = x). */
static double lsum_down(double x0, double r, double b, double n) {
  double lo = fmax(0.0, n - b);
  if (x0 < lo) return -INFINITY;

  double lterm = ldhyper_nat(x0, r, b, n);
  double lsum = lterm;

  for (double x = x0; x > lo; x -= 1.0) {
    double num = x * (b - n + x);
    double den = (r - x + 1.0) * (n - x + 1.0);
    if (num <= 0.0 || den <= 0.0) break;
    double lratio = log(num) - log(den);
    lterm += lratio;
    lsum = logaddexp(lsum, lterm);
    if (lratio < 0.0) {
      double lrem = lterm - log1p(-exp(lratio));
      if (lrem < lsum - KYCG_TAIL_CUTOFF) break;
    }
  }
  return lsum;
}

/*
 * Both tails are evaluated by summing the *smaller* one and complementing if
 * necessary. Summing the larger tail directly would walk from one end of the
 * support to the other — hundreds of thousands of terms at whole-genome
 * scale — and the accumulated rounding is enough to push the total above 1,
 * yielding a positive log p-value. Starting from the mode and moving outward
 * keeps every summation short and monotone, and the complement is taken with
 * log1mexp so no precision is lost converting back.
 */

/** Natural log of P(X >= x0). */
static double lphyper_upper_nat(double x0, double r, double b, double n) {
  double lo = fmax(0.0, n - b);
  double hi = fmin(n, r);
  if (x0 <= lo) return 0.0;            /* the whole support: log(1) */
  if (x0 > hi)  return -INFINITY;      /* beyond the support */

  if (x0 <= hyper_mode(r, b, n))       /* upper tail is the big one */
    return log1mexp(lsum_down(x0 - 1.0, r, b, n));

  return fmin(lsum_up(x0, r, b, n), 0.0);
}

/** Natural log of P(X <= x0). */
static double lphyper_lower_nat(double x0, double r, double b, double n) {
  double lo = fmax(0.0, n - b);
  double hi = fmin(n, r);
  if (x0 >= hi) return 0.0;
  if (x0 < lo)  return -INFINITY;

  if (x0 >= hyper_mode(r, b, n))       /* lower tail is the big one */
    return log1mexp(lsum_up(x0 + 1.0, r, b, n));

  return fmin(lsum_down(x0, r, b, n), 0.0);
}

/* ------------------------------------------------------------------- API */

double kycg_ldhyper(double x, double r, double b, double n) {
  return ldhyper_nat(x, r, b, n) / M_LN10;
}

double kycg_lphyper_upper(double x, double r, double b, double n) {
  return lphyper_upper_nat(x, r, b, n) / M_LN10;
}

double kycg_lphyper_lower(double x, double r, double b, double n) {
  return lphyper_lower_nat(x, r, b, n) / M_LN10;
}

double kycg_hypergeo_log10p(uint64_t nDQ, uint64_t nQ, uint64_t nD,
                            uint64_t nU, kycg_alt_t alt) {
  /* knowYourCG parameterizes phyper with m = nQ white, n = nU - nQ black,
   * k = nD draws; see calculate_fisher_pvalue() in testEnrichment.R. */
  if (nU < nQ || nU < nD || nDQ > nQ || nDQ > nD) return NAN;

  double x = (double)nDQ;
  double r = (double)nQ;
  double b = (double)nU - (double)nQ;
  double n = (double)nD;

  switch (alt) {
  case KYCG_ALT_GREATER:
    return lphyper_upper_nat(x, r, b, n) / M_LN10;
  case KYCG_ALT_LESS:
    return lphyper_lower_nat(x, r, b, n) / M_LN10;
  case KYCG_ALT_TWO_SIDED: {
    double up = lphyper_upper_nat(x, r, b, n) / M_LN10;
    double lw = lphyper_lower_nat(x, r, b, n) / M_LN10;
    double v = fmin(up, lw) + log10(2.0);
    return fmin(v, 0.0);
  }
  default:
    return NAN;
  }
}

/* --------------------------------------------------------- BH correction */

typedef struct {
  double v;
  size_t i;
} bh_pair_t;

static int bh_cmp(const void *a, const void *b) {
  double x = ((const bh_pair_t *)a)->v;
  double y = ((const bh_pair_t *)b)->v;
  /* NaN sorts last so it never participates in the step-up. */
  if (isnan(x)) return isnan(y) ? 0 : 1;
  if (isnan(y)) return -1;
  if (x < y) return -1;
  if (x > y) return 1;
  return 0;
}

void kycg_bh_log10(const double *log10p, const size_t *idx, size_t n_idx,
                   double *log10_fdr) {
  if (!n_idx) return;

  bh_pair_t *ord = malloc(n_idx * sizeof(bh_pair_t));
  if (!ord) {
    /* Fill what we promised rather than leaving the caller's buffer as we
     * found it. kycg_apply_fdr pre-fills NaN so today this changes nothing,
     * but the header documents that this function writes every index, and a
     * caller with an uninitialized buffer would otherwise read garbage and
     * report it as an FDR. */
    for (size_t j = 0; j < n_idx; ++j) log10_fdr[idx[j]] = NAN;
    return;
  }

  /* Only finite p-values take part; NaN rows get NaN FDR. */
  size_t m = 0;
  for (size_t j = 0; j < n_idx; ++j) {
    double v = log10p[idx[j]];
    if (isnan(v)) {
      log10_fdr[idx[j]] = NAN;
    } else {
      ord[m].v = v;
      ord[m].i = idx[j];
      ++m;
    }
  }
  if (!m) { free(ord); return; }

  qsort(ord, m, sizeof(bh_pair_t), bh_cmp);

  /* Step-up from the largest p-value down, in log10:
   *   fdr_(j) = min( 1, min_{l >= j} p_(l) * m / l )
   * Seeding the running minimum at log10(1) = 0 is the pmin(1, ...) clamp. */
  double log10_m = log10((double)m);
  double cummin = 0.0;
  for (size_t j = m; j-- > 0; ) {
    double raw = ord[j].v + log10_m - log10((double)(j + 1));
    if (raw < cummin) cummin = raw;
    log10_fdr[ord[j].i] = cummin;
  }

  free(ord);
}
