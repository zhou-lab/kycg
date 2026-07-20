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

#ifndef _KYCG_ENRICH_H
#define _KYCG_ENRICH_H

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>

#include "hypergeo.h"

/**
 * One enrichment test: one query sample against one knowledgebase record.
 *
 * The four counts come straight from YAME's stats_t (n_u/n_q/n_m/n_o);
 * everything else on this struct is computed here. Strings are owned by the
 * struct and released by kycg_results_free().
 */
typedef struct {
  char *query_file; /* the .cg the sample came from, as displayed           */
  char *query;      /* sample name within the .cg                          */
  char *db_file;    /* knowledgebase file, as displayed                     */
  char *db;         /* record name within the knowledgebase                 */
  char *group;      /* knowledgebase grouping; part of the FDR stratum      */

  uint64_t nU;      /* universe                                             */
  uint64_t nQ;      /* query                                                */
  uint64_t nD;      /* database / mask                                      */
  uint64_t nDQ;     /* overlap                                              */

  double estimate;  /* log2 odds ratio                                      */
  double log10_p;
  double log10_fdr;

  double cf_jaccard;
  double cf_mcc;
  double cf_overlap;
  double cf_npmi;
  double cf_dice;

  double beta;      /* mean beta over the record; NAN when unavailable      */
  double depth;     /* mean depth over the record; NAN when unavailable     */
} kycg_result_t;

/**
 * Fill estimate and the five cf_* coefficients from the four counts.
 * Mirrors knowYourCG's calculate_odds_ratio() and calculate_effect_sizes().
 */
void kycg_effect_sizes(kycg_result_t *r);

/**
 * Compute log10_p for a row from its counts under the given alternative.
 */
void kycg_result_pvalue(kycg_result_t *r, kycg_alt_t alt);

/**
 * Benjamini-Hochberg across `n` results, writing log10_fdr.
 *
 * The correction is always applied within one query sample, because that is
 * the unit knowYourCG corrects over: testEnrichment() takes a single query
 * and calls p.adjust() on the result frame for that query alone. Pooling
 * several samples into one correction would change every reported FDR.
 *
 * When by_group is nonzero the stratum is narrowed further to one
 * knowledgebase, which is knowYourCG's default (mtc_by_group=TRUE). So the
 * stratum is (query_file, query, group), or (query_file, query) when
 * by_group is zero.
 */
void kycg_apply_fdr(kycg_result_t *res, size_t n, int by_group);

/**
 * Sort in knowYourCG's reporting order — ascending log10 p-value, ties broken
 * by descending |estimate|, matching order(log10.p.value, -abs(estimate)) —
 * applied within each query sample, whose rows are kept contiguous.
 */
void kycg_sort_results(kycg_result_t *res, size_t n);

/**
 * The display group for a knowledgebase file name.
 *
 * knowYourCG sets a DB's group to the knowledgebase file name, then strips it
 * for display: "KYCG.MSA.CGI.20220904" renders as "CGI" (drop the KYCG prefix,
 * the platform, and the trailing date). Names not matching that convention
 * are returned unchanged. Caller frees.
 */
char *kycg_display_group(const char *db_file);

/** Header line for the TSV emitted by kycg_write_results(). */
void kycg_write_header(FILE *out);

/** One tab-separated row per result; NaN renders as "NA" (YAME convention). */
void kycg_write_results(FILE *out, const kycg_result_t *res, size_t n);

void kycg_results_free(kycg_result_t *res, size_t n);

#endif /* _KYCG_ENRICH_H */
