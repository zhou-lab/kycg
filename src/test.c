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
 * `kycg test` — set enrichment of a CpG query against a knowledgebase.
 *
 * GOAL
 *   Replace the knowYourCG sequencing workflow, which is a shell pipeline
 *   around `yame summary` plus three things R contributed: a hypergeometric
 *   test, FDR correction, and a plot. This command does the first two and
 *   emits TSV for the third.
 *
 * INPUTS
 *   query .cg   one or more samples, any YAME query format (fmt 0/2/3/4/6)
 *   mask  .cm   one or more knowledgebase records (-m), fmt0 sets or fmt2
 *               categorical states
 *
 * SEMANTICS
 *   For every (query sample, knowledgebase record) pair, YAME's summarize1()
 *   yields the four counts of the 2x2 contingency table — universe, query,
 *   database, overlap. kycg turns each into a one-sided hypergeometric tail
 *   probability in log10, six effect-size coefficients, and a
 *   group-stratified BH false discovery rate.
 *
 * ROW SPACES
 *   A .cg and a .cm are only comparable if they index the same reference row
 *   list; row i must mean the same CpG in both. kycg asserts that the two
 *   record lengths agree and fails loudly otherwise. It deliberately does no
 *   more than that: it does not infer platforms and does not try to detect
 *   two files that happen to share a row count but come from different row
 *   spaces. Per DESIGN.md section 2, keeping the row spaces matched is the
 *   user's responsibility, and the machinery to police it costs more than it
 *   is worth.
 *
 * STRUCTURE
 *   The iteration over query records and mask records — including the
 *   seekable vs in-memory split for the mask file — deliberately mirrors
 *   YAME's main_summary() in src/summary.c, so that the two stay comparable
 *   when either side changes. What differs is that results are accumulated
 *   rather than printed per pair: FDR is a property of the whole table, so
 *   nothing can be emitted until every test has been run.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <inttypes.h>

#include "kycg.h"
#include "enrich.h"

/* YAME (submodule) */
#include "cfile.h"
#include "cdata.h"
#include "summary.h"
#include "wzmisc.h"

typedef struct {
  char *fname_out;
  char *fname_snames;
  kycg_alt_t alt;
  int by_group;      /* stratify BH by knowledgebase; on by default */
  int in_memory;
  int no_header;
  int full_name;
} test_conf_t;

/**
 * One opened knowledgebase file.
 *
 * -m is repeatable so a query can be tested against an entire store in a
 * single pass, one -m per knowledgebase file. That matters because
 * reading the query is the expensive part: testing against 30 knowledgebases
 * one process at a time decompresses the query 30 times, and pooling them here
 * decompresses it once. It also puts every result in one table, which is what
 * FDR wants, since the correction runs per (query, knowledgebase) stratum and
 * the ordering is global.
 */
typedef struct {
  char       *fname;      /* path as given                                  */
  const char *disp;       /* what appears in the db_file column             */
  cfile_t     cf;
  snames_t    snames;
  cdata_t    *mem;        /* records held in memory, or NULL if seekable    */
  uint64_t    n_mem;
} mask_src_t;

static int usage(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg test [options] -m <knowledgebase.cm> <query.cg> [...]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Test each query sample for enrichment against each record of a\n");
  fprintf(stderr, "knowledgebase, using a one-sided hypergeometric tail test.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -m FILE   knowledgebase (.cm) to test against [required];\n");
  fprintf(stderr, "              repeatable -- pass a whole store to test them all\n");
  fprintf(stderr, "    -a STR    alternative: greater|less|two.sided [greater]\n");
  fprintf(stderr, "    -G        correct FDR globally instead of within knowledgebase\n");
  fprintf(stderr, "    -s FILE   sample names for the query\n");
  fprintf(stderr, "    -M        load the knowledgebase into memory\n");
  fprintf(stderr, "    -F        report full file paths instead of basenames\n");
  fprintf(stderr, "    -H        suppress the header line\n");
  fprintf(stderr, "    -o FILE   write to FILE instead of stdout\n");
  fprintf(stderr, "    -h        this help\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "The query and the knowledgebase must index the same reference row\n");
  fprintf(stderr, "list; kycg checks that their lengths agree but cannot verify that\n");
  fprintf(stderr, "they describe the same row space.\n");
  fprintf(stderr, "\n");
  return 1;
}

/* Growable result table. */
typedef struct {
  kycg_result_t *a;
  size_t n, m;
} results_t;

static void results_push(results_t *v, kycg_result_t r) {
  if (v->n == v->m) {
    v->m = v->m ? v->m * 2 : 256;
    v->a = realloc(v->a, v->m * sizeof(kycg_result_t));
    if (!v->a) wzfatal("[%s:%d] Cannot allocate result table.\n", __func__, __LINE__);
  }
  v->a[v->n++] = r;
}

/**
 * Run one (query record, mask record) pair and append its rows.
 *
 * summarize1() may return more than one stat per pair: an fmt2 mask carries
 * several categorical states, each of which is its own knowledgebase record
 * and so its own test.
 */
static void run_pair(results_t *v, cdata_t *c_qry, cdata_t *c_mask,
                     const char *sq, const char *sm,
                     const char *fname_qry, const char *fname_mask,
                     test_conf_t *conf) {

  /* The only row-space guard kycg offers. YAME's per-format handlers make the
   * same comparison, but they report it from deep in the library; surfacing it
   * here lets us say which files were involved and what to do about it. */
  if (c_mask->n && c_qry->n != c_mask->n) {
    wzfatal("[%s:%d] Row count mismatch: query '%s' record '%s' has %" PRIu64
            " rows but knowledgebase '%s' record '%s' has %" PRIu64 " rows.\n"
            "These files index different reference row lists and cannot be "
            "compared.\n",
            __func__, __LINE__, fname_qry, sq, c_qry->n, fname_mask, sm,
            c_mask->n);
  }

  /* summarize1 takes non-const char* but does not modify the names. */
  config_t ycfg = {0};
  ycfg.fname_mask = (char *)fname_mask;   /* non-NULL selects the masked path */

  uint64_t n_st = 0;
  stats_t *st = summarize1(c_qry, c_mask, &n_st, (char *)sm, (char *)sq, &ycfg);

  for (uint64_t i = 0; i < n_st; ++i) {
    kycg_result_t r = {0};

    r.query_file = strdup(fname_qry);
    r.query      = strdup(st[i].sq ? st[i].sq : sq);
    r.db_file    = strdup(fname_mask);
    r.db         = strdup(st[i].sm ? st[i].sm : sm);

    /* The FDR stratum is the knowledgebase file. knowYourCG reaches the same
     * grouping from either direction: attr(db,"group") is set to the
     * knowledgebase name, and determine_group() otherwise falls back to
     * MFile. See enrich.c. */
    r.group = strdup(fname_mask);

    r.nU  = st[i].n_u;
    r.nQ  = st[i].n_q;
    r.nD  = st[i].n_m;
    r.nDQ = st[i].n_o;

    /* Beta and depth are pass-through from YAME, matching the conventions in
     * format_stats_and_clean(): a negative beta means "not applicable to this
     * format", and a zero depth sum means no depth was carried. */
    r.beta = (st[i].beta >= 0) ? st[i].beta : NAN;
    if (st[i].sum_depth) {
      uint64_t denom = st[i].n_m ? st[i].n_m : st[i].n_u;
      r.depth = denom ? (double)st[i].sum_depth / (double)denom : NAN;
    } else {
      r.depth = NAN;
    }

    kycg_effect_sizes(&r);
    kycg_result_pvalue(&r, conf->alt);
    r.log10_fdr = NAN;    /* filled once the whole table is known */

    results_push(v, r);
  }

  for (uint64_t i = 0; i < n_st; ++i) { free(st[i].sm); free(st[i].sq); }
  free(st);
}

int main_test(int argc, char *argv[]) {
  test_conf_t conf = {0};
  conf.alt = KYCG_ALT_GREATER;
  conf.by_group = 1;

  /* -m is repeatable; collect the paths as they arrive. */
  char **mask_names = NULL;
  size_t n_masks = 0;

  int c;
  while ((c = getopt(argc, argv, "m:a:Gs:MFHo:h")) >= 0) {
    switch (c) {
    case 'm':
      mask_names = realloc(mask_names, (n_masks + 1) * sizeof(char *));
      if (!mask_names) wzfatal("[%s:%d] Cannot allocate.\n", __func__, __LINE__);
      mask_names[n_masks++] = strdup(optarg);
      break;
    case 'o': conf.fname_out = strdup(optarg); break;
    case 's': conf.fname_snames = strdup(optarg); break;
    case 'G': conf.by_group = 0; break;
    case 'M': conf.in_memory = 1; break;
    case 'F': conf.full_name = 1; break;
    case 'H': conf.no_header = 1; break;
    case 'a':
      if (strcmp(optarg, "greater") == 0)        conf.alt = KYCG_ALT_GREATER;
      else if (strcmp(optarg, "less") == 0)      conf.alt = KYCG_ALT_LESS;
      else if (strcmp(optarg, "two.sided") == 0) conf.alt = KYCG_ALT_TWO_SIDED;
      else { usage(); wzfatal("Unknown alternative: %s.\n", optarg); }
      break;
    case 'h': return usage();
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind >= argc) {
    usage();
    wzfatal("Please supply a query file.\n");
  }
  if (!n_masks) {
    usage();
    wzfatal("Please supply at least one knowledgebase with -m.\n");
  }

  mask_src_t *masks = calloc(n_masks, sizeof(mask_src_t));
  if (!masks) wzfatal("[%s:%d] Cannot allocate.\n", __func__, __LINE__);

  for (size_t i = 0; i < n_masks; ++i) {
    masks[i].fname  = mask_names[i];
    masks[i].disp   = conf.full_name ? mask_names[i]
                                     : get_basename(mask_names[i]);
    masks[i].cf     = open_cfile(mask_names[i]);
    masks[i].snames = loadSampleNamesFromIndex(mask_names[i]);

    /* An unseekable knowledgebase (a pipe) must be slurped, since it is
     * re-read once per query record. Mirrors main_summary()'s handling. */
    int unseekable = bgzf_seek(masks[i].cf.fh, 0, SEEK_SET);
    if (conf.in_memory || unseekable) {
      for (;;) {
        cdata_t m = read_cdata1(&masks[i].cf);
        if (m.n == 0) break;
        prepare_mask(&m);
        masks[i].mem = realloc(masks[i].mem,
                               (masks[i].n_mem + 1) * sizeof(cdata_t));
        if (!masks[i].mem)
          wzfatal("[%s:%d] Cannot allocate mask table.\n", __func__, __LINE__);
        masks[i].mem[masks[i].n_mem++] = m;
      }
    }
  }
  free(mask_names);

  results_t v = {0};

  for (int j = optind; j < argc; ++j) {
    char *fname_qry = argv[j];
    cfile_t cf_qry = open_cfile(fname_qry);

    snames_t snames_qry;
    if (conf.fname_snames) snames_qry = loadSampleNames(conf.fname_snames, 1);
    else snames_qry = loadSampleNamesFromIndex(fname_qry);

    const char *qry_disp = conf.full_name ? fname_qry : get_basename(fname_qry);

    for (uint64_t kq = 0;; ++kq) {
      cdata_t c_qry = read_cdata1(&cf_qry);
      if (c_qry.n == 0) break;

      if (snames_qry.n && kq >= (unsigned)snames_qry.n) {
        wzfatal("[%s:%d] More records (N=%" PRIu64 ") in '%s' than names in "
                "its index (N=%d).\n",
                __func__, __LINE__, kq + 1, fname_qry, snames_qry.n);
      }

      kstring_t sq = {0};
      if (snames_qry.n) kputs(snames_qry.s[kq], &sq);
      else ksprintf(&sq, "%" PRIu64, kq + 1);

      prepare_mask(&c_qry);

      /* Every knowledgebase, against this one decompressed query record. */
      for (size_t mi = 0; mi < n_masks; ++mi) {
        mask_src_t *ms = &masks[mi];

        if (ms->n_mem) {                  /* records already in memory */
          for (uint64_t km = 0; km < ms->n_mem; ++km) {
            kstring_t sm = {0};
            if (ms->snames.n && km < (unsigned)ms->snames.n)
              kputs(ms->snames.s[km], &sm);
            else ksprintf(&sm, "%" PRIu64, km + 1);

            run_pair(&v, &c_qry, &ms->mem[km], sq.s, sm.s,
                     qry_disp, ms->disp, &conf);
            free(sm.s);
          }
        } else {                          /* re-scan the seekable file */
          if (bgzf_seek(ms->cf.fh, 0, SEEK_SET) != 0) {
            wzfatal("[%s:%d] Cannot seek knowledgebase '%s'.\n",
                    __func__, __LINE__, ms->fname);
          }
          for (uint64_t km = 0;; ++km) {
            cdata_t c_mask = read_cdata1(&ms->cf);
            if (c_mask.n == 0) break;
            prepare_mask(&c_mask);

            kstring_t sm = {0};
            if (ms->snames.n && km < (unsigned)ms->snames.n)
              kputs(ms->snames.s[km], &sm);
            else ksprintf(&sm, "%" PRIu64, km + 1);

            run_pair(&v, &c_qry, &c_mask, sq.s, sm.s,
                     qry_disp, ms->disp, &conf);
            free(sm.s);
            free_cdata(&c_mask);
          }
        }
      }

      free(sq.s);
      free_cdata(&c_qry);
      c_qry.s = NULL;
    }

    bgzf_close(cf_qry.fh);
    cleanSampleNames2(snames_qry);
  }

  /* FDR is a property of the whole table, so it waits until every pair has
   * been tested; then the table is ordered for reporting. */
  kycg_apply_fdr(v.a, v.n, conf.by_group);
  kycg_sort_results(v.a, v.n);

  FILE *out = stdout;
  if (conf.fname_out) {
    out = fopen(conf.fname_out, "w");
    if (!out) wzfatal("[%s:%d] Cannot open '%s' for writing.\n",
                      __func__, __LINE__, conf.fname_out);
  }

  if (!conf.no_header) kycg_write_header(out);
  kycg_write_results(out, v.a, v.n);

  if (conf.fname_out) fclose(out);

  kycg_results_free(v.a, v.n);
  for (size_t i = 0; i < n_masks; ++i) {
    for (uint64_t k = 0; k < masks[i].n_mem; ++k) free_cdata(&masks[i].mem[k]);
    free(masks[i].mem);
    bgzf_close(masks[i].cf.fh);
    cleanSampleNames2(masks[i].snames);
    free(masks[i].fname);
  }
  free(masks);
  free(conf.fname_out);
  free(conf.fname_snames);

  return 0;
}
