// SPDX-License-Identifier: LicenseRef-CHOP-Academic-BSD-2-Clause
/**
 * This file is part of kycg.
 *
 * Copyright (C) 2026-present The Children's Hospital of Philadelphia.
 *
 * Use of this software is available to academic and non-profit institutions
 * for research purposes subject to the terms of the 2-Clause BSD License.
 * For use or transfers of the software to commercial entities, please inquire
 * with Dr. Wanding Zhou at zhouw3@chop.edu.
 *
 * See the LICENSE file at the root of this repository for the full terms and
 * the warranty disclaimer.
 */

/**
 * Set enrichment score of a continuous query: `kycg test -c`.
 *
 * GOAL
 *   knowYourCG's testEnrichmentSEA(): does a set crowd the small or the large
 *   values of a continuous query (a beta value, a score)? It is a
 *   Kolmogorov-Smirnov-style running sum, as in GSEA.
 *
 * THE WALK
 *   Sort the rows that have a value, smallest first. Walk them, adding 1/k at
 *   each of the set's k members and subtracting 1/(n-k) at each of the other
 *   n-k rows. The walk ends at zero. Its maximum is es_small, minus its
 *   minimum is es_large, and the estimate is the larger of the two, signed
 *   the GSEA way: positive when the set crowds the large values.
 *   knowYourCG 1.8.x had the opposite sign; 1.9.3 flipped it to this one.
 *
 *   The walk falls between members and rises at them, so its maximum sits
 *   just after a member and its minimum just before one (or at an end, where
 *   it is zero). With h members passed after i rows the walk is
 *   h/k - (i-h)/(n-k), which can be evaluated at those points directly. One
 *   pass over the sorted query therefore serves any number of sets.
 *
 * THE P-VALUE
 *   KYCG_SEA_PERMS random member sets of the same size, each walked the same
 *   way. Under the null es_small and es_large have the same distribution, so
 *   both go into one pool and both directions are judged against it.
 *   knowYourCG once kept them apart, and its p-value then depended on the
 *   sign of the query (zhou-lab/knowYourCG#6). When nothing in the pool
 *   reaches the observed score, a normal fitted to the pool gives the tail.
 *
 *   The null depends only on (n, k). It is seeded from that pair, so it is
 *   reproducible and order-independent, and cached on it.
 */

#include "sea.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#include "wzmisc.h"   /* YAME: wzfatal */

/* ------------------------------------------------------------ the query */

static int ent_cmp(const void *a, const void *b) {
  const kycg_sea_ent_t *x = a, *y = b;
  if (x->val < y->val) return -1;
  if (x->val > y->val) return 1;
  return (x->row > y->row) - (x->row < y->row);   /* ties in row order */
}

int kycg_sea_query_init(kycg_sea_query_t *q, cdata_t *c) {
  memset(q, 0, sizeof(*q));
  if (c->fmt != '3' && c->fmt != '4') return -1;

  q->nrows = c->n;
  q->ent = malloc((c->n ? c->n : 1) * sizeof(kycg_sea_ent_t));
  if (!q->ent) wzfatal("[%s:%d] Cannot allocate the query.\n", __func__, __LINE__);

  if (c->fmt == '3') {
    for (uint64_t i = 0; i < c->n; ++i) {
      uint64_t mu = f3_get_mu(c, i);
      if (!MU2cov(mu)) continue;               /* no coverage: missing */
      q->ent[q->n].val = MU2beta(mu);
      q->ent[q->n++].row = i;
    }
  } else {
    const float *v = (const float *)c->s;
    for (uint64_t i = 0; i < c->n; ++i) {
      if (!(v[i] >= 0.0f)) continue;           /* negative is missing */
      q->ent[q->n].val = v[i];
      q->ent[q->n++].row = i;
    }
  }

  qsort(q->ent, q->n, sizeof(kycg_sea_ent_t), ent_cmp);
  return 0;
}

void kycg_sea_query_free(kycg_sea_query_t *q) {
  free(q->ent);
  memset(q, 0, sizeof(*q));
}

/* ------------------------------------------------------------- the walk */

kycg_sea_es_t kycg_sea_es(const uint64_t *pos, uint64_t k, uint64_t n) {
  kycg_sea_es_t es = {0, 0};
  if (!k || k >= n) return es;
  double up = 1.0 / (double)k, down = 1.0 / (double)(n - k);
  double mx = 0, mn = 0;
  for (uint64_t j = 0; j < k; ++j) {
    double miss = (double)(pos[j] - j);        /* non-members before it */
    double before = (double)j * up - miss * down;
    double after  = (double)(j + 1) * up - miss * down;
    if (before < mn) mn = before;
    if (after > mx) mx = after;
  }
  es.es_small = mx;
  es.es_large = -mn;
  return es;
}

double kycg_sea_estimate(kycg_sea_es_t es) {
  return es.es_large > es.es_small ? es.es_large : -es.es_small;
}

/* -------------------------------------------------------------- the null */

/* splitmix64 seeds xoshiro256**; both public domain (Vigna). */
static uint64_t splitmix64(uint64_t *x) {
  uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

typedef struct { uint64_t s[4]; } rng_t;

static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

static uint64_t rng_next(rng_t *r) {
  uint64_t *s = r->s;
  uint64_t out = rotl(s[1] * 5, 7) * 9;
  uint64_t t = s[1] << 17;
  s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
  s[2] ^= t; s[3] = rotl(s[3], 45);
  return out;
}

/* Uniform on [0, bound), by multiply-shift; the bias is under 2^-32 for any
 * bound a row space can have. */
static inline uint64_t rng_below(rng_t *r, uint64_t bound) {
  return (uint64_t)(((unsigned __int128)rng_next(r) * bound) >> 64);
}

/**
 * One random member set, walked. Floyd's algorithm draws min(k, n-k) distinct
 * rows into a bitmap -- the members, or the non-members when those are fewer
 * -- and a scan of the bitmap walks them in order. The extremes then sit at
 * the drawn rows either way: around a member as in kycg_sea_es(), or, for a
 * drawn non-member, the maximum just before it and the minimum just after.
 */
static kycg_sea_es_t null_draw(rng_t *r, uint64_t *bits, uint64_t n, uint64_t k) {
  int members = k <= n - k;
  uint64_t m = members ? k : n - k;

  for (uint64_t j = n - m; j < n; ++j) {
    uint64_t t = rng_below(r, j + 1);
    if (bits[t >> 6] >> (t & 63) & 1) t = j;
    bits[t >> 6] |= 1ULL << (t & 63);
  }

  double up = 1.0 / (double)k, down = 1.0 / (double)(n - k);
  double mx = 0, mn = 0;
  uint64_t j = 0;                       /* drawn rows passed */
  for (uint64_t w = 0; w < (n + 63) / 64; ++w) {
    uint64_t word = bits[w];
    bits[w] = 0;                        /* clear as we go, for the next draw */
    while (word) {
      uint64_t p = (w << 6) + (uint64_t)__builtin_ctzll(word);
      word &= word - 1;
      double before, after;
      if (members) {
        double miss = (double)(p - j);
        before = (double)j * up - miss * down;
        after  = (double)(j + 1) * up - miss * down;
        if (before < mn) mn = before;
        if (after > mx) mx = after;
      } else {
        double hit = (double)(p - j);
        before = hit * up - (double)j * down;
        after  = hit * up - (double)(j + 1) * down;
        if (before > mx) mx = before;
        if (after < mn) mn = after;
      }
      ++j;
    }
  }

  kycg_sea_es_t es = { mx, -mn };
  return es;
}

static int dbl_cmp(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}

void kycg_sea_null_build(kycg_sea_null_t *nl, uint64_t n, uint64_t k) {
  memset(nl, 0, sizeof(*nl));
  nl->n = n; nl->k = k;
  nl->m = 2 * KYCG_SEA_PERMS;
  nl->draw = malloc(nl->m * sizeof(double));
  uint64_t *bits = calloc((n + 63) / 64 + 1, sizeof(uint64_t));
  if (!nl->draw || !bits)
    wzfatal("[%s:%d] Cannot allocate the null.\n", __func__, __LINE__);

  /* Seeded from the pair alone: a fixed constant, then n, then k. */
  uint64_t x = 0x6b79636753454121ULL ^ (n * 0x9e3779b97f4a7c15ULL) ^ k;
  rng_t r;
  for (int i = 0; i < 4; ++i) r.s[i] = splitmix64(&x);

  for (size_t i = 0; i < KYCG_SEA_PERMS; ++i) {
    kycg_sea_es_t es = null_draw(&r, bits, n, k);
    nl->draw[2 * i]     = es.es_small;
    nl->draw[2 * i + 1] = es.es_large;
  }
  free(bits);

  double s = 0, ss = 0;
  for (size_t i = 0; i < nl->m; ++i) s += nl->draw[i];
  nl->mean = s / (double)nl->m;
  for (size_t i = 0; i < nl->m; ++i) {
    double d = nl->draw[i] - nl->mean;
    ss += d * d;
  }
  nl->sd = sqrt(ss / (double)(nl->m - 1));   /* R's sd() */
  qsort(nl->draw, nl->m, sizeof(double), dbl_cmp);
}

void kycg_sea_null_free(kycg_sea_null_t *nl) {
  free(nl->draw);
  memset(nl, 0, sizeof(*nl));
}

/* log10 of P(Z > z) for a standard normal, without underflow. Past z = 35
 * erfc() leaves the double range, and the asymptotic series takes over. */
static double log10_norm_upper(double z) {
  if (z < 35) return log10(0.5 * erfc(z / sqrt(2.0)));
  double z2 = z * z;
  double ln = -0.5 * z2 - log(z) - 0.5 * log(2 * M_PI)
              + log1p(-1 / z2 + 3 / (z2 * z2));
  return ln / log(10.0);
}

double kycg_sea_log10_p(const kycg_sea_null_t *nl, double es) {
  /* the first draw >= es; everything from there on counts */
  size_t lo = 0, hi = nl->m;
  while (lo < hi) {
    size_t mid = lo + (hi - lo) / 2;
    if (nl->draw[mid] < es) lo = mid + 1; else hi = mid;
  }
  size_t at_or_above = nl->m - lo;
  if (at_or_above) return log10((double)at_or_above / (double)nl->m);

  if (!(nl->sd > 0)) return -INFINITY;     /* a null with no spread */
  return log10_norm_upper((es - nl->mean) / nl->sd);
}

const kycg_sea_null_t *kycg_sea_cache_get(kycg_sea_cache_t *c,
                                          uint64_t n, uint64_t k) {
  for (size_t i = 0; i < c->n; ++i)
    if (c->a[i].n == n && c->a[i].k == k) return &c->a[i];
  if (c->n == c->m) {
    c->m = c->m ? c->m * 2 : 16;
    c->a = realloc(c->a, c->m * sizeof(kycg_sea_null_t));
    if (!c->a) wzfatal("[%s:%d] Cannot allocate.\n", __func__, __LINE__);
  }
  kycg_sea_null_build(&c->a[c->n], n, k);
  return &c->a[c->n++];
}

void kycg_sea_cache_free(kycg_sea_cache_t *c) {
  for (size_t i = 0; i < c->n; ++i) kycg_sea_null_free(&c->a[i]);
  free(c->a);
  memset(c, 0, sizeof(*c));
}

/* ------------------------------------------------------ sets of a mask */

/* Running state of one set during the pass. */
typedef struct {
  uint64_t k, h;
  double   mx, mn, sum;
} walk_t;

static inline void walk_hit(walk_t *w, uint64_t i, uint64_t n, double val) {
  double miss = (double)(i - w->h);
  double up = 1.0 / (double)w->k, down = 1.0 / (double)(n - w->k);
  double before = (double)w->h * up - miss * down;
  double after  = (double)(w->h + 1) * up - miss * down;
  if (before < w->mn) w->mn = before;
  if (after > w->mx) w->mx = after;
  ++w->h;
  w->sum += val;
}

static void walk_finish(const walk_t *w, uint64_t n, kycg_sea_set_t *s) {
  s->k = w->k;
  s->sum = w->sum;
  if (w->k && w->k < n) {
    s->es.es_small = w->mx;
    s->es.es_large = -w->mn;
  }
}

size_t kycg_sea_walk(const kycg_sea_query_t *q, cdata_t *mask,
                     kycg_sea_set_t **out) {
  uint64_t n = q->n;
  *out = NULL;

  if (mask->fmt == '0' || mask->fmt == '6') {
    int f6 = mask->fmt == '6';
    walk_t w = {0};
    kycg_sea_set_t *s = calloc(1, sizeof(kycg_sea_set_t));
    if (!s) wzfatal("[%s:%d] Cannot allocate.\n", __func__, __LINE__);

    for (uint64_t r = 0; r < q->nrows; ++r)
      if (f6 ? (FMT6_IN_UNI(*mask, r) && FMT6_IN_SET(*mask, r))
             : FMT0_IN_SET(*mask, r)) ++s->nD;
    for (uint64_t i = 0; i < n; ++i) {
      uint64_t r = q->ent[i].row;
      if (f6 ? (FMT6_IN_UNI(*mask, r) && FMT6_IN_SET(*mask, r))
             : FMT0_IN_SET(*mask, r)) ++w.k;
    }
    if (w.k && w.k < n) {
      for (uint64_t i = 0; i < n; ++i) {
        uint64_t r = q->ent[i].row;
        if (f6 ? (FMT6_IN_UNI(*mask, r) && FMT6_IN_SET(*mask, r))
               : FMT0_IN_SET(*mask, r))
          walk_hit(&w, i, n, q->ent[i].val);
      }
    } else {
      for (uint64_t i = 0; i < n; ++i) {
        uint64_t r = q->ent[i].row;
        if (f6 ? (FMT6_IN_UNI(*mask, r) && FMT6_IN_SET(*mask, r))
               : FMT0_IN_SET(*mask, r)) w.sum += q->ent[i].val;
      }
    }
    walk_finish(&w, n, s);
    *out = s;
    return 1;
  }

  if (mask->fmt == '2') {
    if (!mask->aux) fmt2_set_aux(mask);
    f2_aux_t *aux = (f2_aux_t *)mask->aux;
    uint64_t nk = aux->nk;
    walk_t *w = calloc(nk ? nk : 1, sizeof(walk_t));
    kycg_sea_set_t *s = calloc(nk ? nk : 1, sizeof(kycg_sea_set_t));
    if (!w || !s) wzfatal("[%s:%d] Cannot allocate.\n", __func__, __LINE__);

    for (uint64_t r = 0; r < q->nrows; ++r) {
      uint64_t st = f2_get_uint64(mask, r);
      if (st >= nk) wzfatal("[%s:%d] State data is corrupted.\n", __func__, __LINE__);
      ++s[st].nD;
    }
    for (uint64_t i = 0; i < n; ++i) ++w[f2_get_uint64(mask, q->ent[i].row)].k;
    for (uint64_t i = 0; i < n; ++i) {
      walk_t *ws = &w[f2_get_uint64(mask, q->ent[i].row)];
      if (ws->k < n) walk_hit(ws, i, n, q->ent[i].val);
      else ws->sum += q->ent[i].val;
    }
    for (uint64_t st = 0; st < nk; ++st) {
      walk_finish(&w[st], n, &s[st]);
      s[st].name = strdup(aux->keys[st]);
    }
    free(w);
    *out = s;
    return nk;
  }

  return (size_t)-1;
}
