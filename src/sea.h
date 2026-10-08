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

#ifndef _KYCG_SEA_H
#define _KYCG_SEA_H

#include <stdint.h>
#include <stddef.h>

#include "cdata.h"   /* YAME */

/** Permutations behind every p-value; each one gives two null draws. */
#define KYCG_SEA_PERMS 1000

/**
 * A continuous query, ready to walk: its non-missing rows in ascending order
 * of value, ties kept in row order. That is the order R's sort() produces, and
 * the enrichment score depends on it whenever values tie -- which beta values
 * from counts do all the time.
 */
typedef struct {
  double   val;
  uint64_t row;
} kycg_sea_ent_t;

typedef struct {
  uint64_t        nrows;  /* rows in the record, missing included          */
  uint64_t        n;      /* rows with a value: the universe of the walk   */
  kycg_sea_ent_t *ent;    /* those rows, ascending by (value, row)         */
} kycg_sea_query_t;

/**
 * Build the walk order from a format 3 (M/U counts; beta, missing when there
 * is no coverage) or format 4 (float; negative is missing) record. Returns 0,
 * or -1 for any other format, leaving q empty.
 */
int  kycg_sea_query_init(kycg_sea_query_t *q, cdata_t *c);
void kycg_sea_query_free(kycg_sea_query_t *q);

/**
 * The two extremes of the running sum: +1/k at a member, -1/(n-k) at a
 * non-member, walked from the smallest value up. es_small is its maximum
 * (members crowd the small values), es_large is minus its minimum (members
 * crowd the large values). Both are >= 0.
 */
typedef struct {
  double es_small;
  double es_large;
} kycg_sea_es_t;

/**
 * Extremes of the walk from the sorted 0-based positions of the k members in
 * a universe of n. Evaluated only where they can occur -- just after a member
 * for the maximum, just before one for the minimum -- so the cost is O(k).
 */
kycg_sea_es_t kycg_sea_es(const uint64_t *pos, uint64_t k, uint64_t n);

/**
 * The signed estimate, GSEA's convention: +es_large when the members crowd
 * the large values, -es_small when they crowd the small ones. A tie goes to
 * the small side, as in knowYourCG.
 */
double kycg_sea_estimate(kycg_sea_es_t es);

/**
 * The permutation null for k members among n: KYCG_SEA_PERMS random member
 * sets, each giving its es_small and its es_large. The two have the same
 * distribution (reverse and negate the walk), so both go into one pool of
 * 2 * KYCG_SEA_PERMS draws, sorted ascending.
 *
 * The draws are seeded from (n, k) alone, so the same pair gets the same null
 * in any run and in any order, and a cache keyed on the pair is exact.
 */
typedef struct {
  uint64_t n, k;
  double  *draw;      /* 2 * KYCG_SEA_PERMS, ascending                     */
  size_t   m;
  double   mean, sd;
} kycg_sea_null_t;

void kycg_sea_null_build(kycg_sea_null_t *nl, uint64_t n, uint64_t k);
void kycg_sea_null_free(kycg_sea_null_t *nl);

/**
 * log10 of the one-sided p-value of an observed |estimate|: the share of the
 * pooled null at or above it. When none of the draws reaches it, the tail of
 * a normal fitted to the null takes over, in log space so it never
 * underflows. That fallback is knowYourCG's.
 */
double kycg_sea_log10_p(const kycg_sea_null_t *nl, double es);

/** A cache of nulls, one per (n, k) pair seen. */
typedef struct {
  kycg_sea_null_t *a;
  size_t n, m;
} kycg_sea_cache_t;

const kycg_sea_null_t *kycg_sea_cache_get(kycg_sea_cache_t *c,
                                          uint64_t n, uint64_t k);
void kycg_sea_cache_free(kycg_sea_cache_t *c);

/** One set's walk against a query. */
typedef struct {
  char         *name;     /* state key for a format 2 mask, else NULL      */
  uint64_t      nD;       /* rows in the set, missing query rows included  */
  uint64_t      k;        /* members with a query value                    */
  double        sum;      /* sum of the members' query values              */
  kycg_sea_es_t es;
} kycg_sea_set_t;

/**
 * Walk every set a mask record carries -- one for format 0/1/6, one per state
 * for format 2 -- in a single pass over the query. Returns the number of sets
 * and a malloc'd array (free each name, then the array), or (size_t)-1 for a
 * mask format it cannot read.
 */
size_t kycg_sea_walk(const kycg_sea_query_t *q, cdata_t *mask,
                     kycg_sea_set_t **out);

#endif /* _KYCG_SEA_H */
