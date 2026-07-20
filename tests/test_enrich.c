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
 * Tests for FDR stratification and effect-size edge cases.
 *
 * The stratification tests exist because the failure they guard against is
 * invisible in the common case. kycg accumulates every (query, knowledgebase
 * record) pair into one table before correcting, so it is easy to correct
 * across the whole table by accident. With a single-sample query that is
 * indistinguishable from correct behavior, and with an exactly duplicated
 * second sample it is *still* indistinguishable — pooling doubles both the
 * test count and every rank, and the two cancel. Only a second query with
 * genuinely different p-values separates the two implementations, which is
 * what these cases construct.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "enrich.h"

static int n_fail = 0;

static void check(const char *what, double got, double want) {
  int ok = (isnan(got) && isnan(want)) || fabs(got - want) <= 1e-12;
  if (!ok) {
    printf("  FAIL %-46s got %.17g  want %.17g\n", what, got, want);
    ++n_fail;
  }
}

static void check_true(const char *what, int cond) {
  if (!cond) { printf("  FAIL %s\n", what); ++n_fail; }
}

/* Release the strings in a stack-allocated row array.
 * kycg_results_free() also free()s the array itself, so it cannot be used on
 * the fixed-size arrays these tests build. */
static void free_rows(kycg_result_t *res, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    free(res[i].query_file); free(res[i].query);
    free(res[i].group); free(res[i].db_file); free(res[i].db);
  }
}

/* Build a row carrying only what the FDR code reads. */
static kycg_result_t row(const char *qf, const char *q, const char *g,
                         double log10_p) {
  kycg_result_t r = {0};
  r.query_file = strdup(qf);
  r.query      = strdup(q);
  r.group      = strdup(g);
  r.db_file    = strdup(g);
  r.db         = strdup("x");
  r.log10_p    = log10_p;
  r.log10_fdr  = NAN;
  r.estimate   = 0.0;
  return r;
}

/* BH over one stratum of 4, computed by hand for the reference:
 *   sorted p:      -8, -4, -2, -1   (log10)
 *   raw = p + log10(4) - log10(rank)
 *   rank 1: -8 + 0.60206 - 0        = -7.39794
 *   rank 2: -4 + 0.60206 - 0.30103  = -3.69897
 *   rank 3: -2 + 0.60206 - 0.47712  = -1.87506
 *   rank 4: -1 + 0.60206 - 0.60206  = -1
 * step-up from the bottom leaves each unchanged (already increasing). */
static const double L4[4] = {-8.0, -4.0, -2.0, -1.0};
static const double E4[4] = {-7.397940008672037, -3.6989700043360187,
                             -1.8750612633917 , -1.0};

int main(void) {
  /* ---- stratification: two queries, different p-values ----------------
   * Each query gets its own family of 4 tests. If the correction were
   * pooled, m would be 8 rather than 4 and every value would shift. */
  {
    kycg_result_t res[8];
    for (int i = 0; i < 4; ++i) res[i]     = row("a.cg", "1", "kb", L4[i]);
    /* Second query, deliberately NOT a copy: a duplicate would let a pooled
     * implementation pass by coincidence. */
    for (int i = 0; i < 4; ++i) res[4 + i] = row("b.cg", "1", "kb",
                                                 L4[i] * 0.5 - 0.25);

    kycg_apply_fdr(res, 8, 1);

    for (int i = 0; i < 4; ++i) {
      char lbl[64];
      snprintf(lbl, sizeof(lbl), "stratified by query_file[%d]", i);
      check(lbl, res[i].log10_fdr, E4[i]);
    }
    free_rows(res, 8);
  }

  /* ---- same file, two sample names --------------------------------- */
  {
    kycg_result_t res[8];
    for (int i = 0; i < 4; ++i) res[i]     = row("m.cg", "s1", "kb", L4[i]);
    for (int i = 0; i < 4; ++i) res[4 + i] = row("m.cg", "s2", "kb",
                                                 L4[i] * 0.5 - 0.25);
    kycg_apply_fdr(res, 8, 1);
    for (int i = 0; i < 4; ++i) {
      char lbl[64];
      snprintf(lbl, sizeof(lbl), "stratified by sample[%d]", i);
      check(lbl, res[i].log10_fdr, E4[i]);
    }
    free_rows(res, 8);
  }

  /* ---- knowledgebase strata, and -G collapsing them ----------------- */
  {
    kycg_result_t res[8];
    for (int i = 0; i < 4; ++i) res[i]     = row("q.cg", "1", "kbA", L4[i]);
    for (int i = 0; i < 4; ++i) res[4 + i] = row("q.cg", "1", "kbB",
                                                 L4[i] * 0.5 - 0.25);
    kycg_apply_fdr(res, 8, 1);
    for (int i = 0; i < 4; ++i) {
      char lbl[64];
      snprintf(lbl, sizeof(lbl), "stratified by group[%d]", i);
      check(lbl, res[i].log10_fdr, E4[i]);
    }

    /* With by_group off the two knowledgebases merge into one family of 8,
     * so the m factor doubles and the smallest p-value must move. */
    kycg_apply_fdr(res, 8, 0);
    check_true("global FDR differs from per-group",
               fabs(res[0].log10_fdr - E4[0]) > 1e-6);
    check("global FDR uses m=8", res[0].log10_fdr,
          -8.0 + log10(8.0) - log10(1.0));
    free_rows(res, 8);
  }

  /* ---- effect-size clamping, per knowYourCG's calculate_odds_ratio() -- */
  {
    /* nDQ = 0 -> odds ratio 0 -> clamped to DBL_MIN -> log2 = -1022 */
    kycg_result_t r = {0};
    r.nU = 1000; r.nQ = 100; r.nD = 100; r.nDQ = 0;
    kycg_effect_sizes(&r);
    check("estimate clamps at log2(DBL_MIN)", r.estimate, -1022.0);
    check("cf_jaccard with no overlap", r.cf_jaccard, 0.0);
    check("cf_dice with no overlap", r.cf_dice, 0.0);

    /* Perfect containment: nQmD = 0 makes the ratio infinite. */
    kycg_result_t s = {0};
    s.nU = 1000; s.nQ = 100; s.nD = 200; s.nDQ = 100;
    kycg_effect_sizes(&s);
    check_true("estimate saturates near log2(DBL_MAX)",
               s.estimate > 1023.0 && s.estimate <= 1024.0);
    check("cf_overlap is 1 under containment", s.cf_overlap, 1.0);

    /* An empty knowledgebase leaves the ratio at 0/0, which is missing, not
     * infinite -- R maps NaN to NA here rather than clamping. */
    kycg_result_t t = {0};
    t.nU = 1000; t.nQ = 100; t.nD = 0; t.nDQ = 0;
    kycg_effect_sizes(&t);
    check_true("estimate is NA for an empty knowledgebase", isnan(t.estimate));
  }

  if (n_fail) {
    printf("  %d check(s) FAILED\n", n_fail);
    return 1;
  }
  printf("  all checks passed\n");
  return 0;
}
