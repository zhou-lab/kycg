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
 *   spaces. By design, keeping the row spaces matched is the user's
 *   responsibility, and the machinery to police it costs more than it is
 *   worth.
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
#include <stdarg.h>

#include "kycg.h"
#include "args.h"
#include "enrich.h"
#include "sea.h"
#include "store.h"
#include "registry.h"
#include "ordering.h"
#include "yame_ui.h"

/* YAME (submodule) */
#include "cfile.h"
#include "cdata.h"
#include "summary.h"
#include "wzmisc.h"

typedef struct {
  char *fname_out;
  char *fname_snames;
  char *platform;    /* -p: queries are probe ID lists on this array */
  kycg_alt_t alt;
  int by_group;      /* stratify BH by knowledgebase; on by default */
  int in_memory;
  int no_header;
  int full_name;
  int continuous;    /* -c: set enrichment score of a continuous query */
  int alt_given;
  kycg_sea_cache_t nulls;
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
  FILE *o = stderr;
  fprintf(o, "\n");
  fprintf(o, "  %skycg test%s  %s— set enrichment of a CpG query against a "
             "knowledgebase%s\n\n",
          KYCG_H_TITLE, KYCG_H_OFF, KYCG_H_NOTE, KYCG_H_OFF);

  fprintf(o, "%sUsage%s\n", KYCG_H_TITLE, KYCG_H_OFF);
  fprintf(o, "    kycg test [options] -m %s<knowledgebase>%s %s<query.cg>%s [...]\n\n",
          KYCG_H_KEY, KYCG_H_OFF, KYCG_H_KEY, KYCG_H_OFF);
  fprintf(o, "    %sA knowledgebase is a path, or a set named the way fetch names it:%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "      %s-m mm10:CGI%s          one set from the store\n",
          KYCG_H_KEY, KYCG_H_OFF);
  fprintf(o, "      %s-m mm10:CGI,ChromHMM%s two of them\n",
          KYCG_H_KEY, KYCG_H_OFF);
  fprintf(o, "      %s-m mm10%s              everything cached for that target\n\n",
          KYCG_H_KEY, KYCG_H_OFF);
  fprintf(o, "    %sOne row per (query sample, knowledgebase record): the 2x2 counts,%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %sa one-sided hypergeometric tail in log10, six effect sizes, and a%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %sfalse discovery rate corrected within knowledgebase.%s\n\n",
          KYCG_H_NOTE, KYCG_H_OFF);

  fprintf(o, "%sOptions%s\n", KYCG_H_TITLE, KYCG_H_OFF);
  struct { const char *f, *d; } opt[] = {
    {"-m SPEC", "path or target[:sets]; repeatable [required]"},
    {"-p PLAT", "queries are probe ID lists for this array (- reads stdin)"},
    {"-a STR",  "alternative: greater|less|two.sided [greater]"},
    {"-c",      "continuous query: set enrichment score, as GSEA"},
    {"-G",      "correct FDR globally instead of within knowledgebase"},
    {"-s FILE", "sample names for the query"},
    {"-M",      "load the knowledgebase into memory"},
    {"-F",      "report full file paths instead of basenames"},
    {"-H",      "suppress the header line"},
    {"-o FILE", "write to FILE instead of stdout"},
    {"-h",      "this help"},
  };
  for (size_t i = 0; i < sizeof(opt)/sizeof(opt[0]); ++i)
    fprintf(o, "    %s%-8s%s %s\n", KYCG_H_KEY, opt[i].f, KYCG_H_OFF, opt[i].d);
  fprintf(o, "\n");

  fprintf(o, "    %sWith -p a query is a text file of probe IDs, one per line (the%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %sfirst column if there are several). The universe is every probe%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %son the array, as when the list is packed into a .cg:%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "      %skycg test -p EPIC -m EPIC:CGI,ChromHMM hits.txt%s\n\n",
          KYCG_H_KEY, KYCG_H_OFF);

  fprintf(o, "    %sWith -c the query is a value per row (format 3 beta, or format 4),%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %sand each set is scored by where its rows fall in that order:%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %sestimate > 0 at the large values, < 0 at the small ones.%s\n\n",
          KYCG_H_NOTE, KYCG_H_OFF);

  fprintf(o, "    %sOmit -m on a terminal to pick from the store, filtered to the%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %ssets whose row count matches the query.%s\n\n",
          KYCG_H_NOTE, KYCG_H_OFF);

  fprintf(o, "    %s%s p_value and fdr underflow to 0 for strong results. Sort and%s\n",
          KYCG_H_WARN, yame_ui_unicode() ? "!" : "!", KYCG_H_OFF);
  fprintf(o, "      %sthreshold on log10_p and neglog10_fdr, which stay in log space.%s\n\n",
          KYCG_H_NOTE, KYCG_H_OFF);

  fprintf(o, "    %sQuery and knowledgebase must index the same reference row list;%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %skycg checks their lengths agree but cannot verify they describe%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %sthe same row space.%s\n\n", KYCG_H_NOTE, KYCG_H_OFF);
  return 1;
}

/**
 * Row count of a file's first record, or 0 if it cannot be read.
 *
 * Only the first record is touched. A knowledgebase's records all index the
 * same row list by construction, so one is enough to identify the row space,
 * and the cost is proportional to that record rather than to the file -- a
 * 54 MB multi-record .cm is probed as cheaply as a 3 KB one.
 *
 * The count has to come from prepare_mask() rather than from the raw header:
 * cdata_t.n means bytes for some formats and units for others, and reconciling
 * that by hand is exactly the trap around cdata_nbytes().
 */
static uint64_t first_record_rows(const char *path) {
  /* open_cfile() exits the process on a file it cannot open, which would be a
   * hostile way to react to one stray file while scanning a whole store. */
  FILE *probe = fopen(path, "rb");
  if (!probe) return 0;
  fclose(probe);

  cfile_t cf = open_cfile((char *)path);
  cdata_t c = read_cdata1(&cf);

  uint64_t n = 0;
  if (c.n) {
    prepare_mask(&c);
    n = c.n;
  }

  free_cdata(&c);
  bgzf_close(cf.fh);
  return n;
}

/* Group digits for readability: 21867837 -> "21,867,837". */
static const char *commify(uint64_t v, char *buf, size_t n) {
  char raw[32];
  snprintf(raw, sizeof(raw), "%" PRIu64, v);
  size_t len = strlen(raw), out = 0;
  for (size_t i = 0; i < len && out + 2 < n; ++i) {
    if (i && (len - i) % 3 == 0) buf[out++] = ',';
    buf[out++] = raw[i];
  }
  buf[out] = '\0';
  return buf;
}

/* Append one child row to an expansion. */

/* ------------------------------------------------ the no -m picker */


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
/**
 * The -c form of run_pair: walk every set the mask record carries against the
 * sorted query and score each against its permutation null. The row is the
 * same shape as a count test's, so the table reads the same either way: nU
 * and nQ are the rows with a value, nD the set, overlap the set's rows with a
 * value, beta their mean value. The effect-size columns do not apply.
 */
static void run_pair_sea(results_t *v, const kycg_sea_query_t *q,
                         cdata_t *c_mask, const char *sq, const char *sm,
                         const char *fname_qry, const char *fname_mask,
                         test_conf_t *conf) {
  kycg_sea_set_t *s = NULL;
  size_t ns = kycg_sea_walk(q, c_mask, &s);
  if (ns == (size_t)-1) {
    wzfatal("kycg test -c: knowledgebase '%s' record '%s' is format %c. A "
            "set must be format 0, 1, 2 or 6.\n", fname_mask, sm, c_mask->fmt);
  }

  for (size_t i = 0; i < ns; ++i) {
    kycg_result_t r = {0};
    r.query_file = strdup(fname_qry);
    r.query      = strdup(sq);
    r.db_file    = strdup(fname_mask);
    r.db         = strdup(s[i].name ? s[i].name : sm);
    r.group      = strdup(conf->by_group ? fname_mask : "(global)");

    r.nU  = q->n;
    r.nQ  = q->n;
    r.nD  = s[i].nD;
    r.nDQ = s[i].k;

    /* No member has a value, or every row with a value is a member: there is
     * no walk. knowYourCG reports such a set as estimate 0, p-value 1. */
    if (!s[i].k || s[i].k >= q->n) {
      r.estimate = 0;
      r.log10_p  = 0;
    } else {
      r.estimate = kycg_sea_estimate(s[i].es);
      const kycg_sea_null_t *nl = kycg_sea_cache_get(&conf->nulls, q->n, s[i].k);
      r.log10_p = kycg_sea_log10_p(nl, fabs(r.estimate));
    }

    r.cf_jaccard = r.cf_mcc = r.cf_overlap = r.cf_npmi = r.cf_dice = NAN;
    r.beta  = s[i].k ? s[i].sum / (double)s[i].k : NAN;
    r.depth = NAN;
    r.log10_fdr = NAN;
    results_push(v, r);
    free(s[i].name);
  }
  free(s);
}

static void run_pair(results_t *v, cdata_t *c_qry, cdata_t *c_mask,
                     const char *sq, const char *sm,
                     const char *fname_qry, const char *fname_mask,
                     test_conf_t *conf, const kycg_sea_query_t *q) {

  /* The only row-space guard kycg offers. YAME's per-format handlers make the
   * same comparison, but they report it from deep in the library; surfacing it
   * here lets us say which files were involved and what to do about it. */
  if (c_mask->n && c_qry->n != c_mask->n) {
    wzfatal("kycg test: row count mismatch. Query '%s' record '%s' has %"
            PRIu64 " rows but knowledgebase '%s' record '%s' has %" PRIu64
            " rows.\nThese files index different reference row lists and "
            "cannot be compared.\n",
            fname_qry, sq, c_qry->n, fname_mask, sm, c_mask->n);
  }

  if (conf->continuous) {
    run_pair_sea(v, q, c_mask, sq, sm, fname_qry, fname_mask, conf);
    return;
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
    r.group = strdup(conf->by_group ? fname_mask : "(global)");

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

/**
 * Read a probe ID list into a format-0 query over the platform's rows.
 *
 * An array result is usually a list of probe IDs, and packing it into a .cg
 * only to have kycg unpack it again is a detour. The record built here is the
 * one `yame pack -f b` would write for the same list, so the statistics
 * cannot differ between the two routes.
 *
 * One ID per line, first field (tab, comma or space separated). Blank lines,
 * # comments and a Probe_ID header are skipped; duplicates count once. IDs
 * the platform does not have are left out with one note on stderr, since a
 * list exported from a larger array routinely carries a few.
 */
static void read_id_list(const char *fname, const kycg_ordering_t *ord,
                         const char *platform, cdata_t *out) {
  FILE *in = strcmp(fname, "-") == 0 ? stdin : fopen(fname, "r");
  if (!in) wzfatal("kycg test: cannot open '%s'.\n", fname);

  cdata_t c = {0};
  c.fmt = '0';
  c.compressed = 1;
  c.n = ord->n;
  c.s = wzcalloc((c.n >> 3) + 1, 1);

  uint64_t n_ids = 0, n_found = 0, n_missing = 0;
  char missing_eg[128] = {0};
  char *line = NULL;
  size_t cap = 0;
  while (getline(&line, &cap, in) >= 0) {
    char *id = line + strspn(line, " \t");
    id[strcspn(id, "\t, \r\n")] = '\0';
    if (!*id || *id == '#' || strcmp(id, "Probe_ID") == 0) continue;
    ++n_ids;
    long row = kycg_ordering_find(ord, id);
    if (row < 0) {
      if (!n_missing++) snprintf(missing_eg, sizeof(missing_eg), "%s", id);
      continue;
    }
    if (!FMT0_IN_SET(c, (uint64_t)row)) ++n_found;
    FMT0_SET(c, (uint64_t)row);
  }
  free(line);
  if (in != stdin) fclose(in);

  if (!n_found) {
    wzfatal("kycg test: none of the %" PRIu64 " IDs in '%s' are %s probes.\n"
            "Check the platform given with -p.\n", n_ids, fname, platform);
  }
  if (n_missing) {
    fprintf(stderr, "kycg test: %" PRIu64 " of %" PRIu64 " IDs in '%s' are "
            "not %s probes and were left out (e.g. %s).\n",
            n_missing, n_ids, fname, platform, missing_eg);
  }
  *out = c;
}

/** Test one query record against every knowledgebase. */
static void run_record(results_t *v, cdata_t *c_qry, const char *sq,
                       const char *fname_qry, const char *qry_disp,
                       mask_src_t *masks, size_t n_masks, test_conf_t *conf) {
  prepare_mask(c_qry);

  kycg_sea_query_t q = {0};
  if (conf->continuous && kycg_sea_query_init(&q, c_qry) != 0) {
    wzfatal("kycg test -c: query '%s' record '%s' is format %c. A "
            "continuous query must be format 3 (M/U counts, scored as "
            "beta) or format 4 (one value per row).\n",
            fname_qry, sq, c_qry->fmt);
  }

  /* Every knowledgebase, against this one decompressed query record. */
  for (size_t mi = 0; mi < n_masks; ++mi) {
    mask_src_t *ms = &masks[mi];

    if (ms->n_mem) {                  /* records already in memory */
      for (uint64_t km = 0; km < ms->n_mem; ++km) {
        kstring_t sm = {0};
        if (ms->snames.n && km < (unsigned)ms->snames.n)
          kputs(ms->snames.s[km], &sm);
        else ksprintf(&sm, "%" PRIu64, km + 1);

        run_pair(v, c_qry, &ms->mem[km], sq, sm.s,
                 qry_disp, ms->disp, conf, &q);
        free(sm.s);
      }
    } else {                          /* re-scan the seekable file */
      if (bgzf_seek(ms->cf.fh, 0, SEEK_SET) != 0) {
        wzfatal("[%s:%d] Cannot seek knowledgebase '%s'.\n",
                __func__, __LINE__, ms->fname);
      }
      uint64_t n_mask_rec = 0;
      for (uint64_t km = 0;; ++km) {
        cdata_t c_mask = read_cdata1(&ms->cf);
        if (c_mask.n == 0) break;
        ++n_mask_rec;
        prepare_mask(&c_mask);

        kstring_t sm = {0};
        if (ms->snames.n && km < (unsigned)ms->snames.n)
          kputs(ms->snames.s[km], &sm);
        else ksprintf(&sm, "%" PRIu64, km + 1);

        run_pair(v, c_qry, &c_mask, sq, sm.s,
                 qry_disp, ms->disp, conf, &q);
        free(sm.s);
        free_cdata(&c_mask);
      }

      /* Zero records is not an empty result, it is a file that could not
       * be read as a knowledgebase -- truncated, corrupt, or simply not a
       * .cm. Reporting success with an empty table would let a pipeline
       * treat "the knowledgebase failed to load" as "nothing was
       * enriched", which are opposite conclusions. */
      if (!n_mask_rec) {
        wzfatal("[%s:%d] No records in knowledgebase '%s'. The file is "
                "empty or not a readable .cm.\n",
                __func__, __LINE__, ms->fname);
      }
    }
  }

  kycg_sea_query_free(&q);
}

int main_test(int argc, char *argv[]) {
  test_conf_t conf = {0};
  conf.alt = KYCG_ALT_GREATER;
  conf.by_group = 1;

  /* -m is repeatable; collect the paths as they arrive. */
  char **mask_names = NULL;
  size_t n_masks = 0;

  int c;
  /* Options may follow the target; BSD getopt would stop at it. */
  kycg_permute_args(argc, argv, "m:p:a:cGs:MFHo:h");
  while ((c = getopt(argc, argv, "m:p:a:cGs:MFHo:h")) >= 0) {
    switch (c) {
    case 'm': {
      /* A spec, not just a path: "mm10:CGI" names a set the same way fetch
       * does, and resolves against the store. An existing file is taken
       * literally, so plain paths are unaffected. */
      char **paths = NULL;
      size_t n = kycg_resolve_or_offer(optarg, "test", &paths);
      /* resolve_or_offer has already said what went wrong and, on a terminal,
       * already offered to fix it. Repeating that with a usage dump would
       * bury the one line that matters. */
      if (!n) return 1;
      mask_names = realloc(mask_names, (n_masks + n) * sizeof(char *));
      if (!mask_names) wzfatal("[%s:%d] Cannot allocate.\n", __func__, __LINE__);
      for (size_t i = 0; i < n; ++i) mask_names[n_masks++] = paths[i];
      free(paths);   /* the strings are now owned by mask_names */
      break;
    }
    case 'o': conf.fname_out = strdup(optarg); break;
    case 'p': conf.platform = strdup(optarg); break;
    case 's': conf.fname_snames = strdup(optarg); break;
    case 'G': conf.by_group = 0; break;
    case 'M': conf.in_memory = 1; break;
    case 'F': conf.full_name = 1; break;
    case 'H': conf.no_header = 1; break;
    case 'c': conf.continuous = 1; break;
    case 'a':
      conf.alt_given = 1;
      if (strcmp(optarg, "greater") == 0)        conf.alt = KYCG_ALT_GREATER;
      else if (strcmp(optarg, "less") == 0)      conf.alt = KYCG_ALT_LESS;
      else if (strcmp(optarg, "two.sided") == 0) conf.alt = KYCG_ALT_TWO_SIDED;
      else { usage(); wzfatal("Unknown alternative: %s.\n", optarg); }
      break;
    case 'h': usage(); return 0;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind >= argc) {
    usage();
    wzfatal("Please supply a query file.\n");
  }
  /* A set enrichment score is already one-sided, toward whichever end the set
   * crowds; there is no alternative to choose. */
  if (conf.continuous && conf.alt_given)
    wzfatal("kycg test: -a does not apply with -c. The estimate's sign says "
            "which end of the query the set crowds.\n");
  if (conf.platform && conf.continuous)
    wzfatal("kycg test: -p takes a list of probe IDs, which has no values "
            "to score; -c needs a format 3 or 4 query.\n");
  if (conf.platform && conf.fname_snames)
    wzfatal("kycg test: -s does not apply with -p; each ID list is one "
            "query.\n");

  /* -p: the ordering turns probe IDs into rows. Loaded once, checked against
   * the row count this build pins, as annotate does. */
  kycg_ordering_t ord = {0};
  if (conf.platform) {
    const kycg_array_reg_t *ar = NULL;
    for (const kycg_array_reg_t *r = KYCG_ARRAY_REGISTRY; r->platform; ++r)
      if (strcmp(r->platform, conf.platform) == 0) { ar = r; break; }
    if (!ar)
      wzfatal("kycg test: '%s' is not an array platform this build knows. "
              "Run `kycg fetch` to see what is available.\n", conf.platform);
    if (kycg_ordering_load(kycg_store_root(NULL), conf.platform,
                           conf.platform, "test", &ord) != 0) return 1;
    if (ord.n != ar->rows)
      wzfatal("kycg test: the ordering for %s has %zu probes but this build "
              "pins %" PRIu64 ".\n", conf.platform, ord.n, ar->rows);
  }
  /* No -m: offer the store in the same tree `kycg fetch` uses. Only on a
   * terminal -- a pipeline that forgot -m must fail loudly rather than wait
   * for an answer nobody is there to give. */
  if (!n_masks && yame_ui_interactive()) {
    if (optind >= argc) {
      usage();
      wzfatal("Please supply a query file.\n");
    }

    uint64_t qrows = conf.platform ? ord.n : first_record_rows(argv[optind]);
    if (!qrows) wzfatal("Cannot read a record from query '%s'.\n", argv[optind]);

    /* Targets are filtered by the row counts pinned in the registry, so this
     * is a comparison against a table rather than a scan of the store. The
     * previous picker opened every .cm on disk to ask the same question and
     * took about a second to do it. */
    kycg_pick_target_t tg[32];
    size_t n_tg = 0;
    for (const kycg_seq_reg_t *r = KYCG_SEQ_REGISTRY; r->genome && n_tg < 32; ++r)
      if (r->rows == qrows) {
        tg[n_tg].name = r->genome; tg[n_tg].kind = "whole genome";
        tg[n_tg].rows = r->rows; ++n_tg;
      }
    for (const kycg_array_reg_t *r = KYCG_ARRAY_REGISTRY; r->platform && n_tg < 32; ++r)
      if (r->rows == qrows) {
        tg[n_tg].name = r->platform; tg[n_tg].kind = "array";
        tg[n_tg].rows = r->rows; ++n_tg;
      }

    if (!n_tg) {
      char qb[32];
      fprintf(stderr,
              "No knowledgebase collection indexes the same row list as '%s'\n"
              "(%s rows). Sequencing queries need sets for their genome and array\n"
              "queries need sets for their platform; a .cm from one row space is\n"
              "meaningless in the other. Run `kycg fetch` to see what exists.\n",
              argv[optind], commify(qrows, qb, sizeof(qb)));
      return 1;
    }

    char title[256], qb[32];
    snprintf(title, sizeof(title),
             "Knowledgebases for %s rows -- space to choose, t to test",
             commify(qrows, qb, sizeof(qb)));

    char **chosen = NULL;
    size_t n_chosen = kycg_pick_sets(tg, n_tg, title, 't', "test", NULL, NULL, &chosen);
    if (n_chosen == (size_t)-1 || !n_chosen) {
      kycg_free_specs(chosen, n_chosen == (size_t)-1 ? 0 : n_chosen);
      wzfatal("No knowledgebase selected.\n");
    }

    for (size_t i = 0; i < n_chosen; ++i) {
      char **paths = NULL;
      size_t np = kycg_resolve_spec(chosen[i], NULL, &paths);
      if (!np) {
        fprintf(stderr, "  %sskipping %s: not in the store%s\n",
                yame_ui_yellow(), chosen[i], yame_ui_reset());
        continue;
      }
      for (size_t j = 0; j < np; ++j) {
        mask_names = realloc(mask_names, (n_masks + 1) * sizeof(char *));
        if (!mask_names) wzfatal("[%s:%d] Cannot allocate.\n", __func__, __LINE__);
        mask_names[n_masks++] = paths[j];
      }
      free(paths);
    }
    kycg_free_specs(chosen, n_chosen);

    if (!n_masks) wzfatal("Nothing selected is in the store.\n");
  }

  if (!n_masks) {
    /* The one missing flag, said alone: the 35-line help would bury it. */
    wzfatal("kycg test: please supply at least one knowledgebase with -m "
            "(see kycg test -h).\n");
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
    const char *qry_disp = conf.full_name ? fname_qry : get_basename(fname_qry);

    if (conf.platform) {                /* a probe ID list: one record */
      cdata_t c_qry = {0};
      read_id_list(fname_qry, &ord, conf.platform, &c_qry);
      run_record(&v, &c_qry, "1", fname_qry, qry_disp, masks, n_masks, &conf);
      free_cdata(&c_qry);
      continue;
    }

    cfile_t cf_qry = open_cfile(fname_qry);

    snames_t snames_qry;
    if (conf.fname_snames) snames_qry = loadSampleNames(conf.fname_snames, 1);
    else snames_qry = loadSampleNamesFromIndex(fname_qry);

    uint64_t n_qry_rec = 0;
    for (uint64_t kq = 0;; ++kq) {
      cdata_t c_qry = read_cdata1(&cf_qry);
      if (c_qry.n == 0) break;
      ++n_qry_rec;

      if (snames_qry.n && kq >= (unsigned)snames_qry.n) {
        wzfatal("[%s:%d] More records (N=%" PRIu64 ") in '%s' than names in "
                "its index (N=%d).\n",
                __func__, __LINE__, kq + 1, fname_qry, snames_qry.n);
      }

      kstring_t sq = {0};
      if (snames_qry.n) kputs(snames_qry.s[kq], &sq);
      else ksprintf(&sq, "%" PRIu64, kq + 1);

      run_record(&v, &c_qry, sq.s, fname_qry, qry_disp, masks, n_masks, &conf);

      free(sq.s);
      free_cdata(&c_qry);
      c_qry.s = NULL;
    }

    if (!n_qry_rec) {
      wzfatal("[%s:%d] No records in query '%s'. The file is empty or not a "
              "readable .cg.\n", __func__, __LINE__, fname_qry);
    }

    bgzf_close(cf_qry.fh);
    cleanSampleNames2(snames_qry);
  }
  kycg_ordering_free(&ord);

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
  kycg_sea_cache_free(&conf.nulls);
  free(conf.fname_out);
  free(conf.fname_snames);
  free(conf.platform);

  return 0;
}
