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
 * Tests for the set enrichment score behind `kycg test -c`.
 *
 * The walk is checked against the obvious implementation -- a cumulative sum
 * over every row, which is what knowYourCG does -- because the fast form only
 * evaluates it at the members and is easy to get wrong by one. The null is
 * checked for what the p-value relies on: it is reproducible, sorted, and the
 * same in both directions, which is the property zhou-lab/knowYourCG#6 found
 * missing.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sea.h"

static int n_fail = 0;

static void check(const char *what, double got, double want, double tol) {
  int ok = (isnan(got) && isnan(want)) || fabs(got - want) <= tol;
  if (!ok) {
    printf("  FAIL %-46s got %.17g  want %.17g\n", what, got, want);
    ++n_fail;
  }
}

static void check_true(const char *what, int cond) {
  if (!cond) { printf("  FAIL %s\n", what); ++n_fail; }
}

/* knowYourCG's form: s = +1/k or -1/(n-k) per row, cumsum, max and -min. */
static kycg_sea_es_t es_brute(const int *member, uint64_t n) {
  uint64_t k = 0;
  for (uint64_t i = 0; i < n; ++i) k += member[i] != 0;
  kycg_sea_es_t es = {0, 0};
  double cs = 0, mx = -1e300, mn = 1e300;
  for (uint64_t i = 0; i < n; ++i) {
    cs += member[i] ? 1.0 / k : -1.0 / (n - k);
    if (cs > mx) mx = cs;
    if (cs < mn) mn = cs;
  }
  es.es_small = mx > 0 ? mx : 0;   /* the walk also touches 0 at its ends */
  es.es_large = mn < 0 ? -mn : 0;
  return es;
}

static kycg_sea_es_t es_fast(const int *member, uint64_t n) {
  uint64_t *pos = malloc(n * sizeof(uint64_t)), k = 0;
  for (uint64_t i = 0; i < n; ++i) if (member[i]) pos[k++] = i;
  kycg_sea_es_t es = kycg_sea_es(pos, k, n);
  free(pos);
  return es;
}

static void test_walk(void) {
  /* random sets of every density, against the full cumulative sum */
  srand(3);
  int member[500];
  for (int trial = 0; trial < 200; ++trial) {
    uint64_t n = 2 + rand() % 499;
    double dens = (double)(rand() % 100) / 100.0;
    uint64_t k = 0;
    for (uint64_t i = 0; i < n; ++i) {
      member[i] = (double)rand() / RAND_MAX < dens;
      k += member[i];
    }
    if (!k || k == n) continue;
    kycg_sea_es_t a = es_fast(member, n), b = es_brute(member, n);
    check("walk max matches the cumulative sum", a.es_small, b.es_small, 1e-12);
    check("walk min matches the cumulative sum", a.es_large, b.es_large, 1e-12);
  }

  /* the extremes, by hand: members all at the small end, then the large */
  int lo[10] = {1, 1, 1, 0, 0, 0, 0, 0, 0, 0};
  int hi[10] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1};
  kycg_sea_es_t es = es_fast(lo, 10);
  check("all members first: es_small is 1", es.es_small, 1.0, 1e-12);
  check("all members first: estimate is -1", kycg_sea_estimate(es), -1.0, 1e-12);
  es = es_fast(hi, 10);
  check("all members last: es_large is 1", es.es_large, 1.0, 1e-12);
  check("all members last: estimate is +1", kycg_sea_estimate(es), 1.0, 1e-12);

  /* reversing the order swaps the two extremes, so the estimate flips */
  int a[12] = {1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0}, r[12];
  for (int i = 0; i < 12; ++i) r[i] = a[11 - i];
  kycg_sea_es_t f = es_fast(a, 12), g = es_fast(r, 12);
  check("reversal swaps es_small into es_large", f.es_small, g.es_large, 1e-12);
  check("reversal swaps es_large into es_small", f.es_large, g.es_small, 1e-12);
  check("reversal negates the estimate",
        kycg_sea_estimate(f), -kycg_sea_estimate(g), 1e-12);

  /* a tie between the extremes goes to the small side, as in knowYourCG */
  int t[4] = {1, 0, 0, 1};
  es = es_fast(t, 4);
  check("tied extremes are equal", es.es_small, es.es_large, 1e-12);
  check_true("a tie is reported as negative", kycg_sea_estimate(es) < 0);

  /* no walk without both members and non-members */
  uint64_t p[3] = {0, 1, 2};
  es = kycg_sea_es(p, 0, 3);
  check("k = 0 gives no score", es.es_small + es.es_large, 0, 0);
  es = kycg_sea_es(p, 3, 3);
  check("k = n gives no score", es.es_small + es.es_large, 0, 0);
}

static void test_null(void) {
  kycg_sea_null_t a, b;
  kycg_sea_null_build(&a, 1000, 50);
  kycg_sea_null_build(&b, 1000, 50);
  check_true("the null has 2 draws per permutation", a.m == 2 * KYCG_SEA_PERMS);
  check_true("the same (n, k) gives the same null",
             memcmp(a.draw, b.draw, a.m * sizeof(double)) == 0);
  int sorted = 1, nonneg = 1;
  for (size_t i = 0; i < a.m; ++i) {
    if (i && a.draw[i] < a.draw[i - 1]) sorted = 0;
    if (a.draw[i] < 0) nonneg = 0;
  }
  check_true("the null is sorted ascending", sorted);
  check_true("every draw is an extreme, so >= 0", nonneg);
  kycg_sea_null_free(&b);

  /* p-values: the share of the pool at or above the score */
  check("a score of 0 has p = 1", kycg_sea_log10_p(&a, 0), 0, 1e-12);
  double median = a.draw[a.m / 2];
  double lp = kycg_sea_log10_p(&a, median);
  check_true("the null median has p near 1/2", lp > log10(0.4) && lp < log10(0.6));
  double far = kycg_sea_log10_p(&a, 0.9);
  check_true("a score past every draw falls back to the normal tail",
             far < log10(1.0 / a.m) && isfinite(far));
  double farther = kycg_sea_log10_p(&a, 50);
  check_true("the normal tail does not underflow", isfinite(farther) && farther < -300);
  check_true("the tail is monotone", farther < far);

  /* a null with no spread cannot place a score past it */
  kycg_sea_null_t flat = a;
  flat.draw = calloc(a.m, sizeof(double));
  flat.mean = 0; flat.sd = 0;
  check_true("a flat null gives -inf past it", isinf(kycg_sea_log10_p(&flat, 1)));
  free(flat.draw);

  /* Members or non-members drawn, the walk has the same law: a set and its
   * complement walk mirror images. The two sides of the k = n/2 switch
   * must agree in distribution, so their means agree within noise. */
  kycg_sea_null_t s, c;
  kycg_sea_null_build(&s, 400, 100);
  kycg_sea_null_build(&c, 400, 300);
  double se = sqrt(s.sd * s.sd / s.m + c.sd * c.sd / c.m);
  check_true("a set and its complement share a null", fabs(s.mean - c.mean) < 4 * se);
  kycg_sea_null_free(&s);
  kycg_sea_null_free(&c);

  kycg_sea_cache_t cache = {0};
  const kycg_sea_null_t *x = kycg_sea_cache_get(&cache, 1000, 50);
  const kycg_sea_null_t *y = kycg_sea_cache_get(&cache, 1000, 50);
  check_true("the cache returns the stored null", x == y && cache.n == 1);
  check_true("the cached null is the built one",
             memcmp(x->draw, a.draw, a.m * sizeof(double)) == 0);
  kycg_sea_cache_get(&cache, 1000, 51);
  check_true("a new pair adds an entry", cache.n == 2);
  for (int i = 0; i < 20; ++i) kycg_sea_cache_get(&cache, 100, 2 + i);
  check_true("the cache grows", cache.n == 22);
  kycg_sea_cache_free(&cache);
  kycg_sea_null_free(&a);
}

static void test_query_and_walk(void) {
  /* format 4: negative is missing; ties keep row order */
  float v[8] = {0.5f, -1.0f, 0.1f, 0.9f, 0.5f, -1.0f, 0.3f, 0.7f};
  cdata_t c = {0};
  c.fmt = '4'; c.n = 8; c.s = (uint8_t *)v;
  kycg_sea_query_t q;
  check_true("format 4 is a continuous query", kycg_sea_query_init(&q, &c) == 0);
  check_true("missing rows leave the universe", q.n == 6);
  uint64_t want[6] = {2, 6, 0, 4, 7, 3};
  int ord = 1;
  for (int i = 0; i < 6; ++i) if (q.ent[i].row != want[i]) ord = 0;
  check_true("ascending by value, ties in row order", ord);

  /* format 0 mask over rows 3 and 7: the two largest values */
  uint8_t bits[1] = {(1u << 3) | (1u << 7)};
  cdata_t m = {0};
  m.fmt = '0'; m.n = 8; m.s = bits;
  kycg_sea_set_t *s = NULL;
  check_true("a format 0 mask is one set", kycg_sea_walk(&q, &m, &s) == 1);
  check_true("nD counts the set's rows", s[0].nD == 2);
  check_true("k counts members with a value", s[0].k == 2);
  check("sum is over the members", s[0].sum, 0.9f + 0.7f, 1e-6);
  check("members at the top: estimate +1", kycg_sea_estimate(s[0].es), 1.0, 1e-12);
  free(s);

  /* format 6: the set is set AND universe */
  uint8_t b6[2] = {0, 0};
  int rows6[3] = {3, 7, 1};                  /* 1 is missing in the query */
  for (int i = 0; i < 3; ++i) {
    int r = rows6[i];
    b6[r >> 2] |= (uint8_t)(3u << ((r & 3) * 2));
  }
  b6[0] |= (uint8_t)(1u << (2 * 2));         /* row 2: set, not universe */
  cdata_t m6 = {0};
  m6.fmt = '6'; m6.n = 8; m6.s = b6;
  check_true("a format 6 mask is one set", kycg_sea_walk(&q, &m6, &s) == 1);
  check_true("format 6 nD needs both bits", s[0].nD == 3);
  check_true("format 6 k skips the missing row", s[0].k == 2);
  check("format 6 scores like format 0", kycg_sea_estimate(s[0].es), 1.0, 1e-12);
  free(s);

  /* a set with no member that has a value: no walk */
  uint8_t none[1] = {(1u << 1) | (1u << 5)};
  m.s = none;
  kycg_sea_walk(&q, &m, &s);
  check_true("no members with a value: k = 0", s[0].k == 0);
  check("no members: no score", s[0].es.es_small + s[0].es.es_large, 0, 0);
  free(s);

  /* a set that holds every row with a value: no walk either */
  uint8_t all[1] = {0xff};
  m.s = all;
  kycg_sea_walk(&q, &m, &s);
  check_true("every row a member: k = n", s[0].k == q.n);
  check("full overlap: no score", s[0].es.es_small + s[0].es.es_large, 0, 0);
  check("full overlap still sums", s[0].sum, 0.5f + 0.1f + 0.9f + 0.5f + 0.3f + 0.7f, 1e-6);
  free(s);

  /* an unreadable mask format */
  cdata_t bad = {0};
  bad.fmt = '3'; bad.n = 8;
  check_true("a format 3 mask is refused", kycg_sea_walk(&q, &bad, &s) == (size_t)-1);
  kycg_sea_query_free(&q);

  /* format 3: beta from M/U, no coverage is missing */
  uint64_t mu[4] = {(3ULL << 32) | 1, 0, (1ULL << 32) | 3, (2ULL << 32) | 2};
  cdata_t c3 = {0};
  c3.fmt = '3'; c3.n = 4; c3.unit = 8; c3.s = (uint8_t *)mu;
  check_true("format 3 is a continuous query", kycg_sea_query_init(&q, &c3) == 0);
  check_true("format 3 drops uncovered rows", q.n == 3);
  check("format 3 smallest beta", q.ent[0].val, 0.25, 1e-12);
  check("format 3 largest beta", q.ent[2].val, 0.75, 1e-12);
  kycg_sea_query_free(&q);

  /* a count format is not a continuous query */
  cdata_t c0 = {0};
  c0.fmt = '0'; c0.n = 8; c0.s = bits;
  check_true("format 0 is refused as a query", kycg_sea_query_init(&q, &c0) == -1);
}

int main(void) {
  test_walk();
  test_null();
  test_query_and_walk();
  if (n_fail) { printf("%d check(s) failed\n", n_fail); return 1; }
  printf("all checks passed\n");
  return 0;
}
