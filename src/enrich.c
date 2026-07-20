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
 * Effect sizes, FDR stratification, ordering, and TSV emission for the
 * enrichment table.
 *
 * GOAL
 *   Turn YAME's four counts per (query, knowledgebase record) pair into the
 *   row knowYourCG's testEnrichment() returns, and write it as tidy TSV that
 *   cinderplot can plot without further arithmetic.
 *
 * COUNT ARITHMETIC
 *   The three derived cells are
 *     nDmQ  = nD - nDQ              in the database, not the query
 *     nQmD  = nQ - nDQ              in the query, not the database
 *     nUmDQ = nU - nQ - nD + nDQ    in neither
 *   All three are computed in double, never in uint64_t: the nUmDQ expression
 *   subtracts before it adds, so unsigned arithmetic would wrap for any row
 *   where nQ + nD exceeds nU before the +nDQ correction lands.
 *
 * FDR IS STRATIFIED
 *   knowYourCG corrects within group, not globally (mtc_by_group=TRUE is the
 *   default). The group is the knowledgebase *file*: R sets
 *   attr(db, "group") <- nm in getDBs.R, and determine_group() falls back to
 *   the MFile column when no group attribute survived. Those are the same
 *   thing, which is why kycg can derive the stratum from the .cm file name
 *   alone and needs no external metadata table to reproduce R's numbers.
 *
 * MISSING VALUES
 *   NaN is the in-memory missing marker and renders as "NA" on output, per
 *   the YAME convention that makes the TSV directly read.table-able.
 */

#include "enrich.h"

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

/* ------------------------------------------------------------ effect sizes */

void kycg_effect_sizes(kycg_result_t *r) {
  double nU  = (double)r->nU;
  double nQ  = (double)r->nQ;
  double nD  = (double)r->nD;
  double nDQ = (double)r->nDQ;

  double nDmQ  = nD - nDQ;
  double nQmD  = nQ - nDQ;
  double nUmDQ = nU - nQ - nD + nDQ;

  /* Odds ratio, with knowYourCG's clamping: an infinite ratio saturates at
   * DBL_MAX, an exactly-zero ratio at DBL_MIN, and 0/0 is genuinely missing. */
  double odds = (nDQ * nUmDQ) / (nQmD * nDmQ);
  if (isinf(odds))      odds = DBL_MAX;
  else if (odds == 0.0) odds = DBL_MIN;

  r->estimate = isnan(odds) ? NAN : log2(odds);

  r->cf_jaccard = nDQ / (nD + nQmD);

  double mcc_num = nDQ * nUmDQ - nQmD * nDmQ;
  double mcc_den = sqrt(nD * (nU - nD) * nQ * (nU - nQ));
  r->cf_mcc = mcc_num / mcc_den;

  r->cf_overlap = nDQ / (nD < nQ ? nD : nQ);

  r->cf_npmi = (log2(nD) + log2(nQ) - 2.0 * log2(nU))
             / (log2(nDQ) - log2(nU)) - 1.0;

  r->cf_dice = (2.0 * nDQ) / (nD + nQ);
}

void kycg_result_pvalue(kycg_result_t *r, kycg_alt_t alt) {
  r->log10_p = kycg_hypergeo_log10p(r->nDQ, r->nQ, r->nD, r->nU, alt);
}

/* -------------------------------------------------------------------- FDR */

void kycg_apply_fdr(kycg_result_t *res, size_t n, int by_group) {
  if (!n) return;

  double *log10p = malloc(n * sizeof(double));
  double *fdr    = malloc(n * sizeof(double));
  size_t *idx    = malloc(n * sizeof(size_t));
  if (!log10p || !fdr || !idx) { free(log10p); free(fdr); free(idx); return; }

  for (size_t i = 0; i < n; ++i) {
    log10p[i] = res[i].log10_p;
    fdr[i] = NAN;
  }

  if (!by_group) {
    for (size_t i = 0; i < n; ++i) idx[i] = i;
    kycg_bh_log10(log10p, idx, n, fdr);
  } else {
    /* One pass per distinct group. Groups are few (one per knowledgebase
     * file), so the quadratic scan is cheaper than building a hash. */
    char **seen = malloc(n * sizeof(char *));
    size_t n_seen = 0;
    for (size_t i = 0; i < n; ++i) {
      const char *g = res[i].group ? res[i].group : "";
      int already = 0;
      for (size_t s = 0; s < n_seen; ++s) {
        if (strcmp(seen[s], g) == 0) { already = 1; break; }
      }
      if (already) continue;
      seen[n_seen++] = (char *)g;

      size_t n_idx = 0;
      for (size_t j = 0; j < n; ++j) {
        const char *gj = res[j].group ? res[j].group : "";
        if (strcmp(gj, g) == 0) idx[n_idx++] = j;
      }
      kycg_bh_log10(log10p, idx, n_idx, fdr);
    }
    free(seen);
  }

  for (size_t i = 0; i < n; ++i) res[i].log10_fdr = fdr[i];

  free(log10p); free(fdr); free(idx);
}

/* ---------------------------------------------------------------- ordering */

static int result_cmp(const void *a, const void *b) {
  const kycg_result_t *x = a, *y = b;

  /* Ascending log10 p-value; NaN last. */
  int xn = isnan(x->log10_p), yn = isnan(y->log10_p);
  if (xn || yn) {
    if (xn && yn) return 0;
    return xn ? 1 : -1;
  }
  if (x->log10_p < y->log10_p) return -1;
  if (x->log10_p > y->log10_p) return 1;

  /* Ties: descending |estimate| (R sorts ascending on -abs(estimate)). */
  double ax = isnan(x->estimate) ? -INFINITY : fabs(x->estimate);
  double ay = isnan(y->estimate) ? -INFINITY : fabs(y->estimate);
  if (ax > ay) return -1;
  if (ax < ay) return 1;
  return 0;
}

void kycg_sort_results(kycg_result_t *res, size_t n) {
  qsort(res, n, sizeof(kycg_result_t), result_cmp);
}

/* ------------------------------------------------------------ group naming */

char *kycg_display_group(const char *db_file) {
  if (!db_file) return NULL;

  /* Only the "KYCG.<platform>.<group...>.<date>" convention is rewritten. */
  if (strncmp(db_file, "KYCG.", 5) != 0) return strdup(db_file);

  const char *rest = db_file + 5;

  /* Split on '.', then drop the first field (platform) and the last (date),
   * rejoining whatever is between them. */
  size_t n_dots = 0;
  for (const char *p = rest; *p; ++p) if (*p == '.') ++n_dots;
  if (n_dots < 2) return strdup(db_file);   /* not enough fields to trim */

  const char *beg = strchr(rest, '.');
  if (!beg) return strdup(db_file);
  ++beg;
  const char *end = strrchr(rest, '.');
  if (!end || end <= beg) return strdup(db_file);

  size_t len = (size_t)(end - beg);
  char *out = malloc(len + 1);
  if (!out) return NULL;
  memcpy(out, beg, len);
  out[len] = '\0';
  return out;
}

/* ------------------------------------------------------------------ output */

void kycg_write_header(FILE *out) {
  fputs("query\tdb_file\tdb\tgroup\t"
        "nU\tnQ\tnD\toverlap\t"
        "estimate\tlog10_p\tp_value\tfdr\tneglog10_fdr\t"
        "cf_jaccard\tcf_mcc\tcf_overlap\tcf_npmi\tcf_dice\t"
        "beta\tdepth\n", out);
}

/* Print a double, or "NA" if it is not a number.
 *
 * Negative zero is normalized to zero: a saturated tail legitimately yields
 * -0.0 here, and "-0" in a results table reads as a typo or a sign error to
 * everyone who sees it. */
static void put_d(FILE *out, double v, const char *fmt) {
  if (isnan(v)) { fputs("\tNA", out); return; }
  if (v == 0.0) v = 0.0;
  fputc('\t', out);
  fprintf(out, fmt, v);
}

void kycg_write_results(FILE *out, const kycg_result_t *res, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    const kycg_result_t *r = &res[i];

    fprintf(out, "%s\t%s\t%s\t%s\t%" PRIu64 "\t%" PRIu64
                 "\t%" PRIu64 "\t%" PRIu64,
            r->query   ? r->query   : "NA",
            r->db_file ? r->db_file : "NA",
            r->db      ? r->db      : "NA",
            r->group   ? r->group   : "NA",
            r->nU, r->nQ, r->nD, r->nDQ);

    put_d(out, r->estimate, "%.10g");
    /* log10_p and neglog10_fdr carry full double precision: they are the
     * primary sort keys, they are what downstream validation compares, and at
     * deep-tail magnitudes (1e5 and beyond) fewer digits would quantize the
     * value far more coarsely than the computation's own error. */
    put_d(out, r->log10_p, "%.15g");

    /* p_value is a convenience column only. It is recovered as 10^log10_p and
     * is *expected* to underflow to 0 for strong results; log10_p is the
     * quantity to sort or threshold on. */
    if (isnan(r->log10_p)) fputs("\tNA", out);
    else fprintf(out, "\t%g", pow(10.0, r->log10_p));

    if (isnan(r->log10_fdr)) {
      fputs("\tNA\tNA", out);
    } else {
      /* neglog10_fdr is precomputed because cinderplot has no -log10()
       * transform in aesthetics, and it is the y-axis of three plots. */
      fprintf(out, "\t%g", pow(10.0, r->log10_fdr));
      put_d(out, -r->log10_fdr, "%.15g");
    }

    put_d(out, r->cf_jaccard, "%.10g");
    put_d(out, r->cf_mcc, "%.10g");
    put_d(out, r->cf_overlap, "%.10g");
    put_d(out, r->cf_npmi, "%.10g");
    put_d(out, r->cf_dice, "%.10g");

    put_d(out, r->beta, "%.3f");
    put_d(out, r->depth, "%.3f");

    fputc('\n', out);
  }
}

void kycg_results_free(kycg_result_t *res, size_t n) {
  if (!res) return;
  for (size_t i = 0; i < n; ++i) {
    free(res[i].query);
    free(res[i].db_file);
    free(res[i].db);
    free(res[i].group);
  }
  free(res);
}
