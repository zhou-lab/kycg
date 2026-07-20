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

#ifndef _KYCG_HYPERGEO_H
#define _KYCG_HYPERGEO_H

#include <stdint.h>
#include <stddef.h>

/* Alternative hypothesis for the hypergeometric tail test. */
typedef enum {
  KYCG_ALT_GREATER = 0,   /* enrichment:  P(X >= nDQ) */
  KYCG_ALT_LESS,          /* depletion:   P(X <= nDQ) */
  KYCG_ALT_TWO_SIDED,     /* min(2*min(upper,lower), 1) */
} kycg_alt_t;

/**
 * log10 of the hypergeometric point mass P(X = x), where X is the number of
 * white balls drawn in n draws without replacement from an urn holding r
 * white and b black balls.
 */
double kycg_ldhyper(double x, double r, double b, double n);

/**
 * log10 of the upper tail P(X >= x). Equivalent to R's
 *   phyper(x-1, r, b, n, lower.tail=FALSE, log.p=TRUE) / log(10)
 */
double kycg_lphyper_upper(double x, double r, double b, double n);

/**
 * log10 of the lower tail P(X <= x). Equivalent to R's
 *   phyper(x, r, b, n, lower.tail=TRUE, log.p=TRUE) / log(10)
 */
double kycg_lphyper_lower(double x, double r, double b, double n);

/**
 * The knowYourCG enrichment p-value, in log10, from the four counts.
 *
 *   nDQ - overlap (query AND database)
 *   nQ  - query size
 *   nD  - database size
 *   nU  - universe size
 *
 * Mirrors knowYourCG's calculate_fisher_pvalue(), which parameterizes phyper
 * as (m = nQ, n = nU - nQ, k = nD).
 */
double kycg_hypergeo_log10p(uint64_t nDQ, uint64_t nQ, uint64_t nD,
                            uint64_t nU, kycg_alt_t alt);

/**
 * Benjamini-Hochberg step-up correction carried out entirely in log10 space.
 *
 * `log10p[i]` is read for every i in idx[0..n_idx), and the corresponding
 * corrected value is written to `log10_fdr[idx[j]]`. Passing a subset of
 * indices is how FDR gets stratified by group without copying rows.
 *
 * Working in log space is not an optimization: p-values recovered as
 * 10^log10p routinely underflow to exactly 0 for real enrichment results, and
 * a linear-space BH would then report FDR = 0 for every one of them.
 */
void kycg_bh_log10(const double *log10p, const size_t *idx, size_t n_idx,
                   double *log10_fdr);

#endif /* _KYCG_HYPERGEO_H */
