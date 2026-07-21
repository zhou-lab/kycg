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
#include "store.h"
#include "registry.h"
#include "ui.h"

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
    {"-a STR",  "alternative: greater|less|two.sided [greater]"},
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

  fprintf(o, "    %sOmit -m on a terminal to pick from the store, filtered to the%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %ssets whose row count matches the query.%s\n\n",
          KYCG_H_NOTE, KYCG_H_OFF);

  fprintf(o, "    %s%s p_value and fdr underflow to 0 for strong results. Sort and%s\n",
          KYCG_H_WARN, kycg_ui_unicode() ? "!" : "!", KYCG_H_OFF);
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
static void kid_push_local(kycg_ui_kids_t *k, unsigned char style,
                           const char *key, const char *fmt, ...) {
  char **v = realloc(k->rows, (k->n + 1) * sizeof(char *));
  if (!v) return;
  k->rows = v;
  char **kv = realloc(k->keys, (k->n + 1) * sizeof(char *));
  if (!kv) return;
  k->keys = kv;
  k->keys[k->n] = key ? strdup(key) : NULL;
  unsigned char *sv = realloc(k->styles, k->n + 1);
  if (!sv) return;
  k->styles = sv;
  k->styles[k->n] = style;

  char buf[1024];
  va_list ap;
  va_start(ap, fmt);
  vsnprintf(buf, sizeof(buf), fmt, ap);
  va_end(ap);
  k->rows[k->n] = strdup(buf);
  if (k->rows[k->n]) ++k->n;
}

/* ------------------------------------------------ the no -m picker */

/**
 * State for offering the store as a tree, the same widget `kycg fetch` uses.
 *
 * Roots are the collections whose pinned row count matches the query; children
 * are the sets actually cached under them, which by construction share that
 * row space. So the list cannot contain a choice that `kycg test` would then
 * refuse.
 */
typedef struct {
  const char *root;        /* store root */
  uint64_t    qrows;

  char        **rows;      /* root row text, tab-separated */
  unsigned char *styles;
  char        **names;     /* the target name for each root */
  size_t        n, m;

  char        **chosen;    /* paths picked, filled on accept */
  size_t        n_chosen, m_chosen;
} pickctx_t;

static void pick_add_target(pickctx_t *p, const char *name, const char *kind,
                            uint64_t rows) {
  char **paths = NULL;
  size_t n_cached = kycg_resolve_spec(name, NULL, &paths);
  kycg_free_specs(paths, n_cached);

  if (p->n == p->m) {
    size_t want = p->m ? p->m * 2 : 16;
    p->rows   = realloc(p->rows,   want * sizeof(char *));
    p->names  = realloc(p->names,  want * sizeof(char *));
    p->styles = realloc(p->styles, want);
    if (!p->rows || !p->names || !p->styles) return;
    p->m = want;
  }

  char buf[512], rb[32];
  snprintf(buf, sizeof(buf), "%s\t%s\t%s\t%zu",
           name, kind, commify(rows, rb, sizeof(rb)), n_cached);
  p->rows[p->n] = strdup(buf);
  p->names[p->n] = strdup(name);
  /* Dim a collection with nothing in it: it is listed so the user learns it
   * exists and could be fetched, not because it can be tested against. */
  p->styles[p->n] = n_cached ? KYCG_ROW_HAVE : KYCG_ROW_MISSING;
  ++p->n;
}

static void pick_expand(void *ctx, const char *row, kycg_ui_kids_t *out) {
  (void)ctx;

  char target[128];
  const char *tab = strchr(row, '\t');
  size_t len = tab ? (size_t)(tab - row) : strlen(row);
  if (len >= sizeof(target)) len = sizeof(target) - 1;
  memcpy(target, row, len);
  target[len] = '\0';

  /* Everything the collection publishes, not just what is here: seeing what
   * is missing is half the point, since f can fetch it and t can then test
   * against it without leaving. */
  size_t n = 0;
  kycg_catalogue_t *cat = kycg_catalogue(target, NULL, &n);
  if (!cat) {
    kid_push_local(out, KYCG_ROW_MISSING, NULL,
                   "catalogue unavailable - try: kycg fetch %s", target);
    return;
  }

  for (size_t i = 0; i < n; ++i) {
    char key[512], setn[256];
    snprintf(key, sizeof(key), "%s:%s", target, cat[i].name);
    const char *dot = strchr(cat[i].name, '.');
    size_t l = dot ? (size_t)(dot - cat[i].name) : strlen(cat[i].name);
    if (l >= sizeof(setn)) l = sizeof(setn) - 1;
    memcpy(setn, cat[i].name, l);
    setn[l] = '\0';
    kid_push_local(out, cat[i].cached ? KYCG_ROW_HAVE : KYCG_ROW_MISSING, key,
                   "%-22.22s %-32.32s %s", setn, cat[i].name,
                   cat[i].cached ? "cached" : "-");
  }
  kycg_catalogue_free(cat, n);
}

static void pick_accept(void *ctx, const char *root, const char *key) {
  (void)root;
  pickctx_t *p = ctx;
  if (p->n_chosen == p->m_chosen) {
    size_t want = p->m_chosen ? p->m_chosen * 2 : 16;
    char **v = realloc(p->chosen, want * sizeof(char *));
    if (!v) return;
    p->chosen = v; p->m_chosen = want;
  }
  p->chosen[p->n_chosen] = strdup(key);
  if (p->chosen[p->n_chosen]) ++p->n_chosen;
}

/** f in the picker: fetch whatever is checked but not yet here, then stay. */
static void pick_commit_fetch(void *ctx) {
  pickctx_t *p = ctx;
  if (p->n_chosen) kycg_fetch_specs(p->chosen, p->n_chosen, NULL);
  for (size_t i = 0; i < p->n_chosen; ++i) free(p->chosen[i]);
  p->n_chosen = 0;
}

/** i in the picker: describe the set under the cursor, same panel as fetch. */
static int pick_key(void *ctx, char key, const char *root,
                    const char *child_key) {
  (void)ctx;
  if (key == 'i') return kycg_kb_show_info(root, child_key);
  return 0;
}

static void pick_free(pickctx_t *p) {
  for (size_t i = 0; i < p->n; ++i) { free(p->rows[i]); free(p->names[i]); }
  free(p->rows); free(p->names); free(p->styles);
  for (size_t i = 0; i < p->n_chosen; ++i) free(p->chosen[i]);
  free(p->chosen);
  memset(p, 0, sizeof(*p));
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
  /* Options may follow the target; BSD getopt would stop at it. */
  kycg_permute_args(argc, argv, "m:a:Gs:MFHo:h");
  while ((c = getopt(argc, argv, "m:a:Gs:MFHo:h")) >= 0) {
    switch (c) {
    case 'm': {
      /* A spec, not just a path: "mm10:CGI" names a set the same way fetch
       * does, and resolves against the store. An existing file is taken
       * literally, so plain paths are unaffected. */
      char **paths = NULL;
      size_t n = kycg_resolve_spec(optarg, NULL, &paths);
      if (!n) {
        usage();
        wzfatal("No knowledgebase matches '%s'.\n"
                "Give a path, or a target like mm10:CGI -- and fetch it first "
                "if it is not in the store.\n", optarg);
      }
      mask_names = realloc(mask_names, (n_masks + n) * sizeof(char *));
      if (!mask_names) wzfatal("[%s:%d] Cannot allocate.\n", __func__, __LINE__);
      for (size_t i = 0; i < n; ++i) mask_names[n_masks++] = paths[i];
      free(paths);   /* the strings are now owned by mask_names */
      break;
    }
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
  /* No -m: offer the store in the same tree `kycg fetch` uses. Only on a
   * terminal -- a pipeline that forgot -m must fail loudly rather than wait
   * for an answer nobody is there to give. */
  if (!n_masks && kycg_ui_interactive()) {
    if (optind >= argc) {
      usage();
      wzfatal("Please supply a query file.\n");
    }

    uint64_t qrows = first_record_rows(argv[optind]);
    if (!qrows) wzfatal("Cannot read a record from query '%s'.\n", argv[optind]);

    /* Targets are filtered by the row counts pinned in the registry, so this
     * is a comparison against a table rather than a scan of the store. The
     * previous picker opened every .cm on disk to ask the same question and
     * took about a second to do it. */
    pickctx_t pc = {0};
    pc.root = kycg_store_root(NULL);
    pc.qrows = qrows;

    for (const kycg_seq_reg_t *r = KYCG_SEQ_REGISTRY; r->genome; ++r)
      if (r->rows == qrows) pick_add_target(&pc, r->genome, "whole genome", r->rows);
    for (const kycg_array_reg_t *r = KYCG_ARRAY_REGISTRY; r->platform; ++r)
      if (r->rows == qrows) pick_add_target(&pc, r->platform, "array", r->rows);

    if (!pc.n) {
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

    kycg_ui_tree_t spec = {0};
    spec.title = title;
    spec.header = "target\tkind\trows\tcached_sets";
    spec.roots = pc.rows;
    spec.root_styles = pc.styles;
    spec.n_roots = pc.n;
    spec.expand = pick_expand;
    /* Two verbs on one screen: fetch what is missing, then test against it.
     * f keeps the browser open (it has a commit); t ends it and the selection
     * is what gets tested. */
    spec.actions[0].key = 'f';
    spec.actions[0].verb = "fetch";
    spec.actions[0].accept = pick_accept;
    spec.actions[0].commit = pick_commit_fetch;
    spec.actions[1].key = 't';
    spec.actions[1].verb = "test";
    spec.actions[1].accept = pick_accept;
    spec.actions[1].commit = NULL;
    spec.n_actions = 2;
    /* Cached sets are the ones worth testing, so they stay checkable. */
    spec.have_selectable = 1;
    /* r and i are the same callbacks the fetch browser uses, so a set is
     * recommended and described identically whichever tree you reached it
     * through. */
    spec.recommend = kycg_kb_recommended;
    spec.on_key = pick_key;
    spec.hint = "i info";
    spec.ctx = &pc;

    int rc = kycg_ui_tree(&spec);
    if (rc != 2 || !pc.n_chosen) {   /* 2 = the test action */
      pick_free(&pc);
      wzfatal("No knowledgebase selected.\n");
    }

    for (size_t i = 0; i < pc.n_chosen; ++i) {
      char **paths = NULL;
      size_t np = kycg_resolve_spec(pc.chosen[i], NULL, &paths);
      if (!np) {
        fprintf(stderr, "  %sskipping %s: not in the store%s\n",
                kycg_ui_yellow(), pc.chosen[i], kycg_ui_reset());
        continue;
      }
      for (size_t j = 0; j < np; ++j) {
        mask_names = realloc(mask_names, (n_masks + 1) * sizeof(char *));
        if (!mask_names) wzfatal("[%s:%d] Cannot allocate.\n", __func__, __LINE__);
        mask_names[n_masks++] = paths[j];
      }
      free(paths);
    }
    pick_free(&pc);
    if (!n_masks) wzfatal("Nothing selected is in the store.\n");
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
