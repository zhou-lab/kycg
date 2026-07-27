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
 * `kycg annotate` — label rows of a TSV by what knowledgebase sets their
 * probes fall in.
 *
 * GOAL
 *   Reproduce knowYourCG's annoProbes(): given probe IDs and a set of
 *   knowledgebases, say which sets each probe belongs to. This is the lookup
 *   half of the same data `kycg test` aggregates over -- testing asks whether
 *   a group of probes is enriched somewhere, annotating asks where one probe
 *   actually is. Analysts want both, usually in that order.
 *
 * WHY A TSV RATHER THAN A LIST OF IDS
 *   A probe list on its own is rarely what anyone has. It is a differential
 *   methylation table, a QC report, a spreadsheet of hits -- and what they
 *   want back is that same table with columns added, not a separate mapping
 *   to join by hand. So this reads a TSV, keeps every column and every row in
 *   order, and appends.
 *
 * HOW A PROBE ID BECOMES A ROW
 *   This is the whole problem. A .cm is positional: it carries no probe IDs,
 *   only a bit or a state per row. The mapping lives in the platform's
 *   ordering file, where a probe's row index *is* its line number after the
 *   header. That file is fetched alongside any array set precisely because a
 *   set without it is a column of anonymous bits.
 *
 *   The ordering is LC_ALL=C-sorted by Probe_ID, which the InfiniumAnnotation
 *   build guarantees, so lookup is a binary search over it rather than a hash
 *   table built at startup. kycg checks that sortedness on load rather than
 *   trusting it: if the guarantee ever broke, silently wrong annotations are
 *   the worst possible failure, and the check costs one pass.
 *
 * PLATFORM IS NAMED, NEVER INFERRED
 *   `-m MSA:CGI` names the platform. A bare path needs `-p`. knowYourCG
 *   guesses the platform from probe ID patterns; kycg does not, for the same
 *   reason it does not guess row spaces anywhere else -- a wrong guess here
 *   produces a full table of confident, wrong answers.
 *
 * OUTPUT
 *   One column per knowledgebase file, holding the labels that file gives the
 *   probe: the state for a categorical set (Island, Shore, ...), or the names
 *   of whichever records contain it for a bitset, comma-joined. A probe in no
 *   set gets NA, as does a probe absent from the ordering -- the two are
 *   different, so the second is also counted and reported on stderr.
 *
 *   -i switches to an indicator matrix: one column per set, 0 or 1. That is
 *   annoProbes(indicator=TRUE), and it is what you want feeding a model
 *   rather than reading.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <inttypes.h>
#include <zlib.h>

#include "kycg.h"
#include "args.h"
#include "store.h"
#include "registry.h"
#include "ui.h"

/* YAME (submodule) */
#include "cfile.h"
#include "cdata.h"
#include "summary.h"
#include "wzmisc.h"

/* ------------------------------------------------------------- the ordering */

typedef struct {
  char  **id;      /* probe IDs, in row order (which is sorted order) */
  size_t  n;
} ordering_t;

static void ordering_free(ordering_t *o) {
  for (size_t i = 0; i < o->n; ++i) free(o->id[i]);
  free(o->id);
  o->id = NULL; o->n = 0;
}

/**
 * Read <store>/<platform>/<platform>.ordering.tsv.gz.
 *
 * YAME v1.33 keys the store on the browser hierarchy -- the address a user
 * navigates to (`<platform>/...`) -- not on the upstream repo that publishes
 * the file, so the old `InfiniumAnnotation/` prefix is gone. The ordering sits
 * at the platform unit root (matching the remote and where sesame-cli lands
 * it). The pre-v1.33 `InfiniumAnnotation/<platform>/` spot is still tried as a
 * fallback for one release, so a store fetched by an older kycg/yame keeps
 * working until it is re-fetched.
 *
 * The row index of a probe is its line number after the header; nothing else
 * in the file matters here. The M/U/col columns describe the bead design and
 * are irrelevant to set membership.
 */
static int ordering_load(const char *store, const char *platform,
                         ordering_t *out) {
  char path[4096], legacy[4096];
  snprintf(path, sizeof(path),
           "%s/%s/%s.ordering.tsv.gz",
           store, platform, platform);
  snprintf(legacy, sizeof(legacy),
           "%s/InfiniumAnnotation/%s/%s.ordering.tsv.gz",
           store, platform, platform);

  gzFile fp = gzopen(path, "rb");
  if (!fp) fp = gzopen(legacy, "rb");   /* one-release fallback (pre-v1.33) */
  if (!fp) {
    fprintf(stderr,
            "kycg annotate: cannot read the probe ordering for %s.\n"
            "  expected: %s\n"
            "  It is fetched with any set for that platform:\n"
            "      yame fetch %s\n", platform, path, platform);
    return -1;
  }

  size_t cap = 1u << 16, n = 0;
  char **v = malloc(cap * sizeof(char *));
  if (!v) { gzclose(fp); return -1; }

  char line[4096];
  int lineno = 0;
  int unsorted_at = -1;
  while (gzgets(fp, line, sizeof(line))) {
    ++lineno;
    if (lineno == 1) continue;              /* header */

    char *tab = strchr(line, '\t');
    if (tab) *tab = '\0';
    size_t len = strcspn(line, "\r\n");
    line[len] = '\0';
    if (!line[0]) continue;

    if (n == cap) {
      size_t want = cap * 2;
      char **nv = realloc(v, want * sizeof(char *));
      if (!nv) { for (size_t i = 0; i < n; ++i) free(v[i]); free(v); gzclose(fp); return -1; }
      v = nv; cap = want;
    }
    v[n] = strdup(line);
    if (!v[n]) break;

    /* Verify the sortedness the binary search depends on, rather than
     * trusting it. A broken ordering would otherwise yield a full table of
     * confident wrong answers, which is worse than any error. */
    if (n && unsorted_at < 0 && strcmp(v[n - 1], v[n]) > 0) unsorted_at = (int)n;
    ++n;
  }
  gzclose(fp);

  if (unsorted_at >= 0) {
    fprintf(stderr,
            "kycg annotate: %s.ordering.tsv.gz is not sorted (at line %d).\n"
            "This build assumes the LC_ALL=C ordering InfiniumAnnotation\n"
            "publishes; refusing to guess row indices from it.\n",
            platform, unsorted_at + 2);
    for (size_t i = 0; i < n; ++i) free(v[i]);
    free(v);
    return -1;
  }

  out->id = v;
  out->n = n;
  return 0;
}

/** Row index of `id`, or -1. Binary search: the ordering is sorted. */
static long ordering_find(const ordering_t *o, const char *id) {
  size_t lo = 0, hi = o->n;
  while (lo < hi) {
    size_t mid = lo + (hi - lo) / 2;
    int c = strcmp(o->id[mid], id);
    if (c == 0) return (long)mid;
    if (c < 0) lo = mid + 1; else hi = mid;
  }
  return -1;
}

/* ----------------------------------------------------------- knowledgebases */

/**
 * The distinct states a categorical record can take.
 *
 * fmt2 stores a key table alongside the per-row values, so the column set for
 * indicator mode comes from the data rather than from a list kycg would have
 * to keep in step with it.
 */
static const char *f2_key(cdata_t *c, uint64_t k) {
  if (!c->aux) fmt2_set_aux(c);
  f2_aux_t *aux = (f2_aux_t *)c->aux;
  return (k < aux->nk) ? aux->keys[k] : "?";
}

/** One .cm, with every record held in memory so rows can be probed at random. */
typedef struct {
  char     *path;
  char      label[256];   /* column name: the set name */
  cdata_t  *rec;
  uint64_t  n_rec;
  snames_t  snames;
} kb_t;

/* "CGI.20220904.cm" -> "CGI"; a path is reduced to its basename first. */
static void kb_label_of(const char *path, char *out, size_t n) {
  const char *base = strrchr(path, '/');
  base = base ? base + 1 : path;
  const char *dot = strchr(base, '.');
  size_t len = dot ? (size_t)(dot - base) : strlen(base);
  if (len >= n) len = n - 1;
  memcpy(out, base, len);
  out[len] = '\0';
}

/** Name of record `k`, or its 1-based index when the file names no records. */
static const char *kb_record_name(const kb_t *kb, uint64_t k, char *buf,
                                  size_t n) {
  if (kb->snames.n && k < (unsigned)kb->snames.n) return kb->snames.s[k];
  snprintf(buf, n, "%" PRIu64, k + 1);
  return buf;
}

static int kb_load(const char *path, uint64_t want_rows, kb_t *out) {
  memset(out, 0, sizeof(*out));
  out->path = strdup(path);
  kb_label_of(path, out->label, sizeof(out->label));
  out->snames = loadSampleNamesFromIndex((char *)path);

  cfile_t cf = open_cfile((char *)path);
  size_t cap = 8;
  out->rec = malloc(cap * sizeof(cdata_t));
  if (!out->rec) { bgzf_close(cf.fh); return -1; }

  for (;;) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    prepare_mask(&c);

    /* Same assertion kycg test makes, for the same reason: a .cm from another
     * row space is not a worse annotation, it is a meaningless one. */
    if (c.n != want_rows) {
      fprintf(stderr,
              "kycg annotate: '%s' indexes %" PRIu64 " rows but the platform "
              "ordering has %" PRIu64 ".\nThey describe different row spaces.\n",
              path, c.n, want_rows);
      free_cdata(&c);
      bgzf_close(cf.fh);
      return -1;
    }

    if (out->n_rec == cap) {
      size_t w = cap * 2;
      cdata_t *nr = realloc(out->rec, w * sizeof(cdata_t));
      if (!nr) { free_cdata(&c); break; }
      out->rec = nr; cap = w;
    }
    out->rec[out->n_rec++] = c;
  }
  bgzf_close(cf.fh);

  if (!out->n_rec) {
    fprintf(stderr, "kycg annotate: no records in '%s'.\n", path);
    return -1;
  }
  return 0;
}

static void kb_free(kb_t *kb) {
  for (uint64_t k = 0; k < kb->n_rec; ++k) free_cdata(&kb->rec[k]);
  free(kb->rec);
  free(kb->path);
  cleanSampleNames2(kb->snames);
}

/**
 * Labels this knowledgebase gives row `row`, joined with `sep` into `out`.
 *
 * Categorical records contribute their state at that row; bitsets contribute
 * the record's name when the bit is set. A row in nothing yields "".
 */
static void kb_labels_at(const kb_t *kb, long row, const char *sep,
                         char *out, size_t n) {
  out[0] = '\0';
  size_t used = 0;

  for (uint64_t k = 0; k < kb->n_rec; ++k) {
    cdata_t *c = &kb->rec[k];
    const char *lab = NULL;
    char idx[32];

    if (c->fmt == '0') {
      if (FMT0_IN_SET(*c, (uint64_t)row)) lab = kb_record_name(kb, k, idx, sizeof(idx));
    } else if (c->fmt == '2') {
      char *s = f2_get_string(c, (uint64_t)row);
      /* A categorical set covers every row, so an empty state is a genuine
       * "unlabelled here" rather than a missing value. */
      if (s && *s) lab = s;
    }
    if (!lab) continue;

    size_t want = strlen(lab) + (used ? strlen(sep) : 0);
    if (used + want + 1 >= n) break;
    if (used) { strcpy(out + used, sep); used += strlen(sep); }
    strcpy(out + used, lab);
    used += strlen(lab);
  }
}

/* --------------------------------------------------------------- the input */

static int usage(void) {
  FILE *o = stderr;
  fprintf(o, "\n");
  fprintf(o, "  %skycg annotate%s  %s— label a TSV's probes by knowledgebase "
             "membership%s\n\n",
          KYCG_H_TITLE, KYCG_H_OFF, KYCG_H_NOTE, KYCG_H_OFF);

  fprintf(o, "%sUsage%s\n", KYCG_H_TITLE, KYCG_H_OFF);
  fprintf(o, "    kycg annotate -m %s<knowledgebase>%s %s<input.tsv>%s\n\n",
          KYCG_H_KEY, KYCG_H_OFF, KYCG_H_KEY, KYCG_H_OFF);
  fprintf(o, "    %sReads a TSV with a column of probe IDs and writes it back with%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %sone column added per knowledgebase, naming the sets each probe%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %sfalls in. Every input column and row is preserved, in order.%s\n\n",
          KYCG_H_NOTE, KYCG_H_OFF);

  fprintf(o, "%sOptions%s\n", KYCG_H_TITLE, KYCG_H_OFF);
  struct { const char *f, *d; } opt[] = {
    {"-m SPEC", "path or platform[:sets]; repeatable [required]"},
    {"-p PLAT", "platform, when -m is a plain path"},
    {"-c COL",  "probe ID column: a name or 1-based index [Probe_ID, else 1]"},
    {"-i",      "indicator columns (0/1) per set instead of joined labels"},
    {"-s SEP",  "separator when a probe is in several sets [,]"},
    {"-H",      "the input has no header line"},
    {"-o FILE", "write to FILE instead of stdout"},
    {"-h",      "this help"},
  };
  for (size_t i = 0; i < sizeof(opt)/sizeof(opt[0]); ++i)
    fprintf(o, "    %s%-8s%s %s\n", KYCG_H_KEY, opt[i].f, KYCG_H_OFF, opt[i].d);
  fprintf(o, "\n");

  fprintf(o, "    %sArrays only: a probe ID means nothing without the platform%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %sordering that gives it a row. The platform is named, never%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %sinferred from the IDs.%s\n\n", KYCG_H_NOTE, KYCG_H_OFF);
  return 1;
}

/* Split a line on tabs in place; returns the field count. */
static size_t split_tabs(char *line, char **fld, size_t max) {
  size_t n = 0;
  char *p = line;
  while (n < max) {
    fld[n++] = p;
    char *t = strchr(p, '\t');
    if (!t) break;
    *t = '\0';
    p = t + 1;
  }
  return n;
}

int main_annotate(int argc, char *argv[]) {
  const char *specs[64];
  size_t n_specs = 0;
  const char *platform = NULL, *colspec = NULL, *sep = ",", *out_path = NULL;
  int indicator = 0, no_header = 0;

  int c;
  kycg_permute_args(argc, argv, "m:p:c:s:o:iHh");
  while ((c = getopt(argc, argv, "m:p:c:s:o:iHh")) >= 0) {
    switch (c) {
    case 'm':
      if (n_specs < 64) specs[n_specs++] = optarg;
      break;
    case 'p': platform = optarg; break;
    case 'c': colspec = optarg; break;
    case 's': sep = optarg; break;
    case 'o': out_path = optarg; break;
    case 'i': indicator = 1; break;
    case 'H': no_header = 1; break;
    case 'h': return usage();
    default: return usage();
    }
  }

  if (optind >= argc) { usage(); wzfatal("Please supply an input TSV.\n"); }

  /* No -m on a terminal: offer the store, exactly as `kycg test` does. Only
   * array platforms are listed -- a probe ID has no meaning without the
   * ordering that gives it a row, so a whole genome is not a choice that
   * could work. Off a terminal this stays an error rather than a wait. */
  char **picked = NULL;
  size_t n_picked = 0;
  if (!n_specs) {
    if (!kycg_ui_interactive()) {
      usage();
      wzfatal("Please supply -m.\n");
    }

    kycg_pick_target_t tg[32];
    size_t n_tg = 0;
    for (const kycg_array_reg_t *r = KYCG_ARRAY_REGISTRY;
         r->platform && n_tg < 32; ++r) {
      tg[n_tg].name = r->platform;
      tg[n_tg].kind = "array";
      tg[n_tg].rows = r->rows;
      ++n_tg;
    }

    n_picked = kycg_pick_sets(tg, n_tg,
                              "Knowledgebases to annotate with -- "
                              "space to choose, a to annotate",
                              'a', "annotate", NULL, NULL, &picked);
    if (n_picked == (size_t)-1) {
      wzfatal("This terminal cannot host the browser; name a set with -m.\n");
    }
    if (!n_picked) wzfatal("No knowledgebase selected.\n");

    for (size_t i = 0; i < n_picked && n_specs < 64; ++i)
      specs[n_specs++] = picked[i];
  }

  /* The platform comes from the spec when it names one, so `-m MSA:CGI` needs
   * nothing further. A plain path carries no platform and cannot be guessed
   * from, so it must be given. */
  char plat_buf[128];
  if (!platform) {
    const char *colon = strchr(specs[0], ':');
    if (colon && !kycg_store_is_file(specs[0])) {
      size_t len = (size_t)(colon - specs[0]);
      if (len >= sizeof(plat_buf)) len = sizeof(plat_buf) - 1;
      memcpy(plat_buf, specs[0], len);
      plat_buf[len] = '\0';
      platform = plat_buf;
    }
  }
  if (!platform) {
    usage();
    wzfatal("Cannot tell which platform '%s' belongs to; name it with -p.\n",
            specs[0]);
  }

  const kycg_array_reg_t *ar = NULL;
  for (const kycg_array_reg_t *r = KYCG_ARRAY_REGISTRY; r->platform; ++r)
    if (strcmp(r->platform, platform) == 0) { ar = r; break; }
  if (!ar) {
    wzfatal("'%s' is not an array platform this build knows. "
            "Run `kycg fetch` to see what is available.\n", platform);
  }

  const char *store = kycg_store_root(NULL);
  ordering_t ord = {0};
  if (ordering_load(store, platform, &ord) != 0) return 1;
  if (ord.n != ar->rows) {
    fprintf(stderr,
            "kycg annotate: the ordering for %s has %zu probes but this build "
            "pins %" PRIu64 ".\n", platform, ord.n, ar->rows);
    ordering_free(&ord);
    return 1;
  }

  /* Resolve every -m to concrete files, then hold them all in memory: the
   * input is streamed once and each row probes every knowledgebase. */
  kb_t *kbs = NULL;
  size_t n_kb = 0, kb_cap = 0;
  for (size_t i = 0; i < n_specs; ++i) {
    char **paths = NULL;
    size_t np = kycg_resolve_or_offer(specs[i], "annotate", &paths);
    if (!np) {
      for (size_t k = 0; k < n_kb; ++k) kb_free(&kbs[k]);
      free(kbs);
      ordering_free(&ord);
      kycg_free_specs(picked, n_picked);
      return 1;
    }
    for (size_t j = 0; j < np; ++j) {
      if (n_kb == kb_cap) {
        size_t w = kb_cap ? kb_cap * 2 : 8;
        kb_t *nk = realloc(kbs, w * sizeof(kb_t));
        if (!nk) break;
        kbs = nk; kb_cap = w;
      }
      if (kb_load(paths[j], ord.n, &kbs[n_kb]) != 0) {
        kycg_free_specs(paths, np);
        for (size_t k = 0; k < n_kb; ++k) kb_free(&kbs[k]);
        free(kbs);
        ordering_free(&ord);
        return 1;
      }
      ++n_kb;
    }
    kycg_free_specs(paths, np);
  }

  if (!n_kb) {
    fprintf(stderr, "kycg annotate: none of the requested sets are in the "
                    "store; nothing to annotate with.\n");
    ordering_free(&ord);
    kycg_free_specs(picked, n_picked);
    return 1;
  }

  FILE *in = strcmp(argv[optind], "-") == 0 ? stdin : fopen(argv[optind], "r");
  if (!in) wzfatal("Cannot open '%s'.\n", argv[optind]);
  FILE *out = out_path ? fopen(out_path, "w") : stdout;
  if (!out) wzfatal("Cannot write '%s'.\n", out_path);

  /* ------------------------------------------------------------- streaming */

  enum { MAXF = 512 };
  char *fld[MAXF];
  char *line = NULL;
  size_t linecap = 0;
  ssize_t len;

  long  col = -1;            /* 0-based index of the probe ID column */
  int   first = 1;
  uint64_t n_rows = 0, n_missing = 0;

  while ((len = getline(&line, &linecap, in)) > 0) {
    while (len && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = '\0';

    /* A copy, because splitting writes NULs into the line and the output has
     * to reproduce the original text. */
    char *keep = strdup(line);
    if (!keep) break;

    size_t nf = split_tabs(line, fld, MAXF);

    if (first && !no_header) {
      first = 0;

      /* Locate the probe column: an explicit -c, else a column literally
       * named Probe_ID, else the first. Named beats positional so a table
       * whose probe column moved still annotates correctly. */
      if (colspec) {
        char *end = NULL;
        long v = strtol(colspec, &end, 10);
        if (end && !*end && v >= 1) col = v - 1;
        else {
          for (size_t i = 0; i < nf; ++i)
            if (strcmp(fld[i], colspec) == 0) { col = (long)i; break; }
          if (col < 0) {
            fprintf(stderr, "kycg annotate: no column named '%s' in the header.\n",
                    colspec);
            free(keep); free(line); fclose(in);
            goto fail;
          }
        }
      } else {
        for (size_t i = 0; i < nf; ++i)
          if (strcmp(fld[i], "Probe_ID") == 0) { col = (long)i; break; }
        if (col < 0) col = 0;
      }

      fputs(keep, out);
      if (indicator) {
        for (size_t k = 0; k < n_kb; ++k) {
          char idx[32];
          for (uint64_t r = 0; r < kbs[k].n_rec; ++r) {
            cdata_t *cd = &kbs[k].rec[r];
            if (cd->fmt == '2') {
              /* One column per state, discovered from the data. */
              uint64_t nkeys = fmt2_get_keys_n(cd);
              for (uint64_t s = 0; s < nkeys; ++s) {
                cdata_t tmp = *cd;
                (void)tmp;
                fprintf(out, "\t%s:%s", kbs[k].label, f2_key(cd, s));
              }
            } else {
              fprintf(out, "\t%s:%s", kbs[k].label,
                      kb_record_name(&kbs[k], r, idx, sizeof(idx)));
            }
          }
        }
      } else {
        for (size_t k = 0; k < n_kb; ++k) fprintf(out, "\t%s", kbs[k].label);
      }
      fputc('\n', out);
      free(keep);
      continue;
    }
    first = 0;

    if (col < 0) col = colspec ? strtol(colspec, NULL, 10) - 1 : 0;

    const char *id = ((size_t)col < nf) ? fld[col] : "";
    long row = *id ? ordering_find(&ord, id) : -1;
    if (row < 0) ++n_missing;
    ++n_rows;

    fputs(keep, out);
    for (size_t k = 0; k < n_kb; ++k) {
      if (indicator) {
        char idx[32];
        for (uint64_t r = 0; r < kbs[k].n_rec; ++r) {
          cdata_t *cd = &kbs[k].rec[r];
          if (cd->fmt == '2') {
            uint64_t nkeys = fmt2_get_keys_n(cd);
            const char *here = (row >= 0) ? f2_get_string(cd, (uint64_t)row) : NULL;
            for (uint64_t s = 0; s < nkeys; ++s) {
              const char *st = f2_key(cd, s);
              fprintf(out, "\t%s", (row < 0) ? "NA"
                      : ((here && st && strcmp(here, st) == 0) ? "1" : "0"));
            }
          } else {
            (void)idx;
            fprintf(out, "\t%s", (row < 0) ? "NA"
                    : (FMT0_IN_SET(*cd, (uint64_t)row) ? "1" : "0"));
          }
        }
      } else {
        if (row < 0) { fputs("\tNA", out); continue; }
        char labels[4096];
        kb_labels_at(&kbs[k], row, sep, labels, sizeof(labels));
        fprintf(out, "\t%s", labels[0] ? labels : "NA");
      }
    }
    fputc('\n', out);
    free(keep);
  }

  free(line);
  if (in != stdin) fclose(in);
  if (out != stdout) fclose(out);

  /* A probe absent from the ordering and a probe in no set both read as NA,
   * so the count that distinguishes them goes to stderr where it cannot
   * corrupt the table. */
  if (n_missing) {
    fprintf(stderr,
            "kycg annotate: %" PRIu64 " of %" PRIu64 " rows had a probe ID not "
            "in the %s ordering; those are NA.\n", n_missing, n_rows, platform);
  }

  for (size_t k = 0; k < n_kb; ++k) kb_free(&kbs[k]);
  free(kbs);
  ordering_free(&ord);
  kycg_free_specs(picked, n_picked);
  return 0;

fail:
  for (size_t k = 0; k < n_kb; ++k) kb_free(&kbs[k]);
  free(kbs);
  ordering_free(&ord);
  kycg_free_specs(picked, n_picked);
  return 1;
}
