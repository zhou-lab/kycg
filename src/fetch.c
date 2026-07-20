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
 * `kycg fetch` — assemble the local knowledgebase store.
 *
 * GOAL
 *   Enrichment testing is only as good as the sets you test against, and until
 *   now those sets have been scattered across three channels with three access
 *   patterns: array knowledgebases in per-platform directories of the
 *   InfiniumAnnotation git repo, sequencing knowledgebases in per-genome Zenodo
 *   records, and a whole-genome CpG reference alongside them. Nothing tied them
 *   together, and nothing verified them. This command is the tie.
 *
 * THE STORE
 *   Everything lands under one directory ($KYCG_DATA_DIR, else ~/.cache/kycg):
 *
 *     <store>/<PLATFORM>/KYCG/<Set>.<date>.cm[.idx]   arrays
 *     <store>/<PLATFORM>/KYCG/SHA256SUMS              as published
 *     <store>/<genome>/<Set>.<date>.cm[.idx]          sequencing
 *     <store>/<genome>/cpg_nocontig.cr                the reference row list
 *     <store>/<genome>/SHA256SUMS                     written by us
 *
 *   The layout is deliberately flat and boring: a `.cm` in the store is an
 *   ordinary file that `kycg test -m` takes by path. There is no database, no
 *   lockfile, and no local bookkeeping — the store is re-verifiable at any time
 *   with `shasum -a 256 -c SHA256SUMS`, using no kycg code at all.
 *
 * TRUST
 *   Both channels verify against a digest compiled into this binary by
 *   tools/make_registry.sh; see src/registry.h and src/digest.c for why the two
 *   channels use different hashes. Downloads land on a temporary path and are
 *   renamed only after their digest matches, so an interrupted or corrupted
 *   fetch can never leave a file in the store that later reads as valid. A file
 *   already present with the right digest is not re-downloaded, which is what
 *   makes re-running fetch after a tag bump cheap.
 *
 * NEVER IMPLICITLY, NEVER INTERACTIVELY
 *   Network access happens here and only here. No other kycg subcommand
 *   downloads anything, and this one never prompts — a prompt would hang a
 *   Nextflow job or a Docker build with no indication of why. Policy lifted
 *   verbatim from sesame-cli, which reached it the same way.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <inttypes.h>

#include "kycg.h"
#include "digest.h"
#include "registry.h"

#ifdef KYCG_HAVE_CURL
#include <curl/curl.h>
#endif

/* ------------------------------------------------------------ store layout */

static const char *store_root(const char *override) {
  static char buf[4096];
  if (override && *override) return override;

  const char *env = getenv("KYCG_DATA_DIR");
  if (env && *env) return env;

  const char *home = getenv("HOME");
  if (!home || !*home) home = ".";
  snprintf(buf, sizeof(buf), "%s/.cache/kycg", home);
  return buf;
}

/* mkdir -p. Returns 0 on success. */
static int mkdir_p(const char *path) {
  char tmp[4096];
  snprintf(tmp, sizeof(tmp), "%s", path);
  size_t n = strlen(tmp);
  if (n && tmp[n-1] == '/') tmp[n-1] = '\0';

  for (char *p = tmp + 1; *p; ++p) {
    if (*p != '/') continue;
    *p = '\0';
    if (mkdir(tmp, 0755) != 0 && errno != EEXIST) return -1;
    *p = '/';
  }
  if (mkdir(tmp, 0755) != 0 && errno != EEXIST) return -1;
  return 0;
}

static int file_exists(const char *path) {
  struct stat st;
  return stat(path, &st) == 0 && S_ISREG(st.st_mode);
}

/* Human-readable byte count, into a caller-supplied buffer. */
static const char *human(uint64_t bytes, char *buf, size_t n) {
  const char *unit[] = {"B", "KB", "MB", "GB"};
  double v = (double)bytes;
  int u = 0;
  while (v >= 1024.0 && u < 3) { v /= 1024.0; ++u; }
  snprintf(buf, n, "%.*f %s", (u == 0 || v >= 100) ? 0 : 1, v, unit[u]);
  return buf;
}

/* ------------------------------------------------------------- set naming */

/**
 * The "set name" of a knowledgebase file: everything before the first dot.
 * "ChromHMM.20220303.cm" -> "ChromHMM". This is what users type for --only,
 * because the dates are an implementation detail of how the sets are versioned
 * and nobody remembers them.
 */
static void set_name_of(const char *fname, char *out, size_t n) {
  const char *dot = strchr(fname, '.');
  size_t len = dot ? (size_t)(dot - fname) : strlen(fname);
  if (len >= n) len = n - 1;
  memcpy(out, fname, len);
  out[len] = '\0';
}

/**
 * Does this file pass the --only filter? `only` is a comma-separated list;
 * a token matches either the whole file name or its set name.
 */
static int passes_filter(const char *fname, const char *only) {
  if (!only || !*only) return 1;

  char setn[256];
  set_name_of(fname, setn, sizeof(setn));

  const char *p = only;
  while (*p) {
    const char *comma = strchr(p, ',');
    size_t len = comma ? (size_t)(comma - p) : strlen(p);
    if (len) {
      if ((strlen(fname) == len && strncasecmp(fname, p, len) == 0) ||
          (strlen(setn)  == len && strncasecmp(setn,  p, len) == 0))
        return 1;
    }
    if (!comma) break;
    p = comma + 1;
  }
  return 0;
}

/* ------------------------------------------------------------- networking */

#ifdef KYCG_HAVE_CURL

typedef struct { char *s; size_t n, m; } membuf_t;

static size_t mem_write(void *data, size_t sz, size_t nm, void *ud) {
  membuf_t *b = ud;
  size_t add = sz * nm;
  if (b->n + add + 1 > b->m) {
    size_t want = (b->n + add + 1) * 2;
    char *p = realloc(b->s, want);
    if (!p) return 0;
    b->s = p; b->m = want;
  }
  memcpy(b->s + b->n, data, add);
  b->n += add;
  b->s[b->n] = '\0';
  return add;
}

static CURL *new_handle(const char *url) {
  CURL *h = curl_easy_init();
  if (!h) return NULL;
  curl_easy_setopt(h, CURLOPT_URL, url);
  curl_easy_setopt(h, CURLOPT_FOLLOWLOCATION, 1L);
  curl_easy_setopt(h, CURLOPT_FAILONERROR, 1L);
  curl_easy_setopt(h, CURLOPT_NOSIGNAL, 1L);
  curl_easy_setopt(h, CURLOPT_CONNECTTIMEOUT, 30L);
  curl_easy_setopt(h, CURLOPT_USERAGENT, "kycg/" KYCG_VERSION);
  return h;
}

/** Fetch a small file into memory. Returns NULL on failure; caller frees. */
static char *http_get_mem(const char *url, size_t *len) {
  CURL *h = new_handle(url);
  if (!h) return NULL;

  membuf_t b = {0};
  curl_easy_setopt(h, CURLOPT_WRITEFUNCTION, mem_write);
  curl_easy_setopt(h, CURLOPT_WRITEDATA, &b);

  CURLcode rc = curl_easy_perform(h);
  curl_easy_cleanup(h);

  if (rc != CURLE_OK) { free(b.s); return NULL; }
  if (len) *len = b.n;
  return b.s;
}

/* Progress rendering, TTY only -- a redirected log should not collect
 * thousands of carriage returns. */
typedef struct { const char *label; int tty; } prog_t;

static int on_xfer(void *ud, curl_off_t dltotal, curl_off_t dlnow,
                   curl_off_t ultotal, curl_off_t ulnow) {
  (void)ultotal; (void)ulnow;
  prog_t *p = ud;
  if (!p->tty || dltotal <= 0) return 0;
  int pct = (int)((100.0 * (double)dlnow) / (double)dltotal);
  char a[32], b[32];
  fprintf(stderr, "\r  %-38s %3d%%  %s / %s   ", p->label, pct,
          human((uint64_t)dlnow, a, sizeof(a)),
          human((uint64_t)dltotal, b, sizeof(b)));
  fflush(stderr);
  return 0;
}

/** Download to `path`. Returns 0 on success. */
static int http_get_file(const char *url, const char *path, const char *label) {
  FILE *fp = fopen(path, "wb");
  if (!fp) return -1;

  CURL *h = new_handle(url);
  if (!h) { fclose(fp); return -1; }

  prog_t pr = { label, isatty(fileno(stderr)) };
  curl_easy_setopt(h, CURLOPT_WRITEDATA, fp);
  curl_easy_setopt(h, CURLOPT_XFERINFOFUNCTION, on_xfer);
  curl_easy_setopt(h, CURLOPT_XFERINFODATA, &pr);
  curl_easy_setopt(h, CURLOPT_NOPROGRESS, 0L);

  CURLcode rc = curl_easy_perform(h);
  curl_easy_cleanup(h);
  fclose(fp);

  if (pr.tty) { fprintf(stderr, "\r%*s\r", 78, ""); fflush(stderr); }

  if (rc != CURLE_OK) { unlink(path); return -1; }
  return 0;
}

#endif /* KYCG_HAVE_CURL */

/* ---------------------------------------------------------------- fetching */

typedef struct {
  const char *store;
  const char *only;
  const char *tag;
  int dry_run;
  int force;
} fetch_conf_t;

/* Tallies for the closing summary. */
typedef struct {
  uint64_t n_got, n_skip, n_fail;
  uint64_t bytes_got;
} tally_t;

#ifdef KYCG_HAVE_CURL

/**
 * Fetch one file into `dir` and verify it.
 *
 * `want_sha` or `want_md5` (exactly one non-NULL) is the expected digest. The
 * download goes to a ".part" sibling and is renamed in only after the digest
 * matches, so the store never contains an unverified file even briefly.
 */
static int fetch_one(const char *url, const char *dir, const char *fname,
                     const char *want_sha, const char *want_md5,
                     const fetch_conf_t *conf, tally_t *t) {
  char path[4096], part[4200], got[65];
  snprintf(path, sizeof(path), "%s/%s", dir, fname);

  /* Already present and correct? Then there is nothing to do -- this is what
   * makes re-running fetch after a tag bump cost only the changed files. */
  if (!conf->force && file_exists(path)) {
    int ok = 0;
    if (want_sha && kycg_sha256_file(path, got) == 0)
      ok = kycg_digest_equal(got, want_sha);
    else if (want_md5) {
      char m[33];
      if (kycg_md5_file(path, m) == 0) ok = kycg_digest_equal(m, want_md5);
    }
    if (ok) { ++t->n_skip; return 0; }
  }

  snprintf(part, sizeof(part), "%s.part", path);
  if (http_get_file(url, part, fname) != 0) {
    fprintf(stderr, "  ! %-40s download failed\n", fname);
    ++t->n_fail;
    return -1;
  }

  int ok = 0;
  if (want_sha && kycg_sha256_file(part, got) == 0)
    ok = kycg_digest_equal(got, want_sha);
  else if (want_md5) {
    char m[33];
    if (kycg_md5_file(part, m) == 0) ok = kycg_digest_equal(m, want_md5);
  }

  if (!ok) {
    unlink(part);
    fprintf(stderr, "  ! %-40s DIGEST MISMATCH -- discarded\n", fname);
    ++t->n_fail;
    return -1;
  }

  struct stat st;
  uint64_t sz = (stat(part, &st) == 0) ? (uint64_t)st.st_size : 0;

  if (rename(part, path) != 0) {
    unlink(part);
    fprintf(stderr, "  ! %-40s could not be moved into the store\n", fname);
    ++t->n_fail;
    return -1;
  }

  char hb[32];
  fprintf(stderr, "  + %-40s %s\n", fname, human(sz, hb, sizeof(hb)));
  ++t->n_got;
  t->bytes_got += sz;
  return 0;
}

#endif /* KYCG_HAVE_CURL */

/* A parsed SHA256SUMS line. */
typedef struct { char sha[65]; char name[512]; } sums_ent_t;

/** Parse "<64 hex>  <name>" lines. Returns a malloc'd array; sets *n. */
static sums_ent_t *parse_sums(const char *text, size_t *n) {
  size_t cap = 64, cnt = 0;
  sums_ent_t *v = malloc(cap * sizeof(sums_ent_t));
  if (!v) return NULL;

  const char *p = text;
  while (*p) {
    const char *eol = strchr(p, '\n');
    size_t len = eol ? (size_t)(eol - p) : strlen(p);

    if (len > 66) {
      if (cnt == cap) {
        cap *= 2;
        sums_ent_t *nv = realloc(v, cap * sizeof(sums_ent_t));
        if (!nv) { free(v); return NULL; }
        v = nv;
      }
      memcpy(v[cnt].sha, p, 64);
      v[cnt].sha[64] = '\0';

      /* Skip the separator (two spaces, or " *" for binary mode). */
      const char *q = p + 64;
      while ((size_t)(q - p) < len && (*q == ' ' || *q == '*')) ++q;
      size_t nlen = len - (size_t)(q - p);
      if (nlen >= sizeof(v[cnt].name)) nlen = sizeof(v[cnt].name) - 1;
      memcpy(v[cnt].name, q, nlen);
      v[cnt].name[nlen] = '\0';
      if (v[cnt].name[0]) ++cnt;
    }

    if (!eol) break;
    p = eol + 1;
  }
  *n = cnt;
  return v;
}

/* ------------------------------------------------------------ registry ops */

static const kycg_array_reg_t *find_array(const char *name) {
  for (const kycg_array_reg_t *r = KYCG_ARRAY_REGISTRY; r->platform; ++r)
    if (strcasecmp(r->platform, name) == 0) return r;
  return NULL;
}

static const kycg_seq_reg_t *find_seq(const char *name) {
  for (const kycg_seq_reg_t *r = KYCG_SEQ_REGISTRY; r->genome; ++r)
    if (strcasecmp(r->genome, name) == 0) return r;
  return NULL;
}

/* -------------------------------------------------------- array platforms */

#ifdef KYCG_HAVE_CURL

static int fetch_array(const kycg_array_reg_t *reg, const fetch_conf_t *conf,
                       tally_t *t) {
  if (!reg->sums_sha256) {
    fprintf(stderr,
            "kycg fetch: platform '%s' is not published at tag %s.\n",
            reg->platform, conf->tag);
    return 1;
  }

  char url[4096];
  snprintf(url, sizeof(url), "%s/%s/%s/KYCG/%s",
           KYCG_IA_BASE_URL, conf->tag, reg->platform, KYCG_IA_SUMS_FILE);

  size_t len = 0;
  char *sums = http_get_mem(url, &len);
  if (!sums) {
    fprintf(stderr, "kycg fetch: cannot reach %s\n", url);
    return 1;
  }

  /* The anchor. Everything downloaded below is trusted only because this
   * one comparison held. */
  char got[65];
  kycg_sha256_buf(sums, len, got);
  if (!kycg_digest_equal(got, reg->sums_sha256)) {
    fprintf(stderr,
            "kycg fetch: SHA256SUMS for %s does not match the digest pinned in\n"
            "this build (tag %s).\n"
            "  expected %s\n"
            "  got      %s\n"
            "Refusing to fetch. This build pins tag %s; if upstream has moved,\n"
            "regenerate src/registry.h with tools/make_registry.sh and rebuild.\n",
            reg->platform, conf->tag, reg->sums_sha256, got, conf->tag);
    free(sums);
    return 1;
  }

  size_t n_ent = 0;
  sums_ent_t *ent = parse_sums(sums, &n_ent);
  if (!ent) { free(sums); return 1; }

  char dir[4096];
  snprintf(dir, sizeof(dir), "%s/%s/KYCG", store_root(conf->store),
           reg->platform);

  if (conf->dry_run) {
    fprintf(stderr, "Would fetch into %s:\n", dir);
    size_t n = 0;
    for (size_t i = 0; i < n_ent; ++i) {
      if (!passes_filter(ent[i].name, conf->only)) continue;
      fprintf(stderr, "  %s\n", ent[i].name);
      ++n;
    }
    fprintf(stderr, "%zu file(s). Sizes are not published for this channel.\n", n);
    free(ent); free(sums);
    return 0;
  }

  if (mkdir_p(dir) != 0) {
    fprintf(stderr, "kycg fetch: cannot create %s\n", dir);
    free(ent); free(sums);
    return 1;
  }

  fprintf(stderr, "%s KYCG sets (tag %s) -> %s\n", reg->platform, conf->tag, dir);

  for (size_t i = 0; i < n_ent; ++i) {
    if (!passes_filter(ent[i].name, conf->only)) continue;
    snprintf(url, sizeof(url), "%s/%s/%s/KYCG/%s",
             KYCG_IA_BASE_URL, conf->tag, reg->platform, ent[i].name);
    fetch_one(url, dir, ent[i].name, ent[i].sha, NULL, conf, t);
  }

  /* Keep the manifest, so the store can be re-verified without kycg. */
  char sp[4200];
  snprintf(sp, sizeof(sp), "%s/%s", dir, KYCG_IA_SUMS_FILE);
  FILE *fp = fopen(sp, "wb");
  if (fp) { fwrite(sums, 1, len, fp); fclose(fp); }

  free(ent);
  free(sums);
  return 0;
}

/* -------------------------------------------------------- sequencing sets */

static int fetch_seq(const kycg_seq_reg_t *reg, const fetch_conf_t *conf,
                     tally_t *t) {
  char dir[4096];
  snprintf(dir, sizeof(dir), "%s/%s", store_root(conf->store), reg->genome);

  if (conf->dry_run) {
    uint64_t total = 0;
    size_t n = 0;
    fprintf(stderr, "Would fetch into %s:\n", dir);
    for (const kycg_zfile_t *f = reg->files; f->name; ++f) {
      if (!passes_filter(f->name, conf->only)) continue;
      char hb[32];
      fprintf(stderr, "  %-44s %s\n", f->name, human(f->size, hb, sizeof(hb)));
      total += f->size;
      ++n;
    }
    char hb[32];
    fprintf(stderr, "%zu file(s), %s total.\n", n, human(total, hb, sizeof(hb)));
    return 0;
  }

  if (mkdir_p(dir) != 0) {
    fprintf(stderr, "kycg fetch: cannot create %s\n", dir);
    return 1;
  }

  fprintf(stderr, "%s knowledgebases (Zenodo %s) -> %s\n",
          reg->genome, reg->record, dir);

  for (const kycg_zfile_t *f = reg->files; f->name; ++f) {
    if (!passes_filter(f->name, conf->only)) continue;
    char url[4096];
    snprintf(url, sizeof(url), "%s/%s/files/%s",
             KYCG_ZENODO_BASE, reg->record, f->name);
    fetch_one(url, dir, f->name, NULL, f->md5, conf, t);
  }

  /* Zenodo publishes md5, but the store should be checkable with the same
   * `shasum -a 256 -c SHA256SUMS` incantation as the array side, so we compute
   * and write one over whatever is now on disk. */
  char sp[4200];
  snprintf(sp, sizeof(sp), "%s/%s", dir, KYCG_IA_SUMS_FILE);
  FILE *fp = fopen(sp, "wb");
  if (fp) {
    for (const kycg_zfile_t *f = reg->files; f->name; ++f) {
      char path[4200], sha[65];
      snprintf(path, sizeof(path), "%s/%s", dir, f->name);
      if (file_exists(path) && kycg_sha256_file(path, sha) == 0)
        fprintf(fp, "%s  %s\n", sha, f->name);
    }
    fclose(fp);
  }

  return 0;
}

#endif /* KYCG_HAVE_CURL */

/* ------------------------------------------------------------------ usage */

static int usage(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg fetch [options] <platform|genome>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Download and verify knowledgebases into a local store.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Targets:\n");
  fprintf(stderr, "    arrays    ");
  for (const kycg_array_reg_t *r = KYCG_ARRAY_REGISTRY; r->platform; ++r)
    fprintf(stderr, "%s ", r->platform);
  fprintf(stderr, "\n    genomes   ");
  for (const kycg_seq_reg_t *r = KYCG_SEQ_REGISTRY; r->genome; ++r)
    fprintf(stderr, "%s ", r->genome);
  fprintf(stderr, "\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -d DIR    store directory [$KYCG_DATA_DIR, else ~/.cache/kycg]\n");
  fprintf(stderr, "    -o SETS   comma-separated subset, by set name or file name\n");
  fprintf(stderr, "              e.g. -o CGI,ChromHMM,TFBS\n");
  fprintf(stderr, "    -n        dry run: list what would be fetched, download nothing\n");
  fprintf(stderr, "    -f        re-download even if present and verified\n");
  fprintf(stderr, "    -t TAG    InfiniumAnnotation tag, arrays only [%s]\n", KYCG_IA_TAG);
  fprintf(stderr, "    -h        this help\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Fetched sets are ordinary files; pass one to `kycg test -m`.\n");
  fprintf(stderr, "kycg never downloads outside this command and never prompts.\n");
  fprintf(stderr, "\n");
  return 1;
}

int main_fetch(int argc, char *argv[]) {
  fetch_conf_t conf = {0};
  conf.tag = KYCG_IA_TAG;

  int c;
  while ((c = getopt(argc, argv, "d:o:t:nfh")) >= 0) {
    switch (c) {
    case 'd': conf.store = optarg; break;
    case 'o': conf.only = optarg; break;
    case 't': conf.tag = optarg; break;
    case 'n': conf.dry_run = 1; break;
    case 'f': conf.force = 1; break;
    case 'h': return usage();
    default: return usage();
    }
  }

  if (optind >= argc) {
    usage();
    fprintf(stderr, "kycg fetch: please name a platform or genome.\n");
    return 1;
  }

#ifndef KYCG_HAVE_CURL
  fprintf(stderr,
          "kycg fetch: this build has no network support.\n"
          "kycg was compiled without libcurl, so fetch is unavailable. Install\n"
          "libcurl development headers and rebuild, or populate the store by\n"
          "hand -- fetched files are ordinary .cm files and `kycg test -m` takes\n"
          "any path. See `kycg list` for the expected layout.\n");
  return 1;
#else
  curl_global_init(CURL_GLOBAL_DEFAULT);

  tally_t t = {0};
  int rc = 0;

  for (int j = optind; j < argc; ++j) {
    const char *target = argv[j];
    const kycg_array_reg_t *ar = find_array(target);
    const kycg_seq_reg_t   *sr = find_seq(target);

    if (ar)      rc |= fetch_array(ar, &conf, &t);
    else if (sr) rc |= fetch_seq(sr, &conf, &t);
    else {
      fprintf(stderr,
              "kycg fetch: '%s' is not a known platform or genome.\n"
              "Run `kycg list` to see what is available.\n", target);
      rc = 1;
    }
  }

  if (!conf.dry_run && (t.n_got || t.n_skip || t.n_fail)) {
    char hb[32];
    fprintf(stderr, "\n%" PRIu64 " fetched (%s), %" PRIu64 " already current",
            t.n_got, human(t.bytes_got, hb, sizeof(hb)), t.n_skip);
    if (t.n_fail) fprintf(stderr, ", %" PRIu64 " FAILED", t.n_fail);
    fprintf(stderr, ".\n");
  }

  curl_global_cleanup();
  return (rc || t.n_fail) ? 1 : 0;
#endif
}

/* -------------------------------------------------------------- kycg list */

static int list_usage(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg list [options] [target ...]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Show available knowledgebase collections and what is cached.\n");
  fprintf(stderr, "With a target, list the individual sets it carries.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -d DIR    store directory [$KYCG_DATA_DIR, else ~/.cache/kycg]\n");
  fprintf(stderr, "    -h        this help\n");
  fprintf(stderr, "\n");
  return 1;
}

/* Count the .cm files already in a store directory. */
static uint64_t count_cached(const char *dir) {
  char sums[4200];
  snprintf(sums, sizeof(sums), "%s/%s", dir, KYCG_IA_SUMS_FILE);

  FILE *fp = fopen(sums, "rb");
  if (!fp) return 0;

  uint64_t n = 0;
  char line[1024];
  while (fgets(line, sizeof(line), fp)) {
    char *nm = strstr(line, "  ");
    if (!nm) continue;
    nm += 2;
    size_t len = strlen(nm);
    while (len && (nm[len-1] == '\n' || nm[len-1] == '\r')) nm[--len] = '\0';
    if (len > 3 && strcmp(nm + len - 3, ".cm") == 0) {
      char path[4300];
      snprintf(path, sizeof(path), "%s/%s", dir, nm);
      if (file_exists(path)) ++n;
    }
  }
  fclose(fp);
  return n;
}

int main_list(int argc, char *argv[]) {
  const char *store = NULL;
  int c;
  while ((c = getopt(argc, argv, "d:h")) >= 0) {
    switch (c) {
    case 'd': store = optarg; break;
    case 'h': return list_usage();
    default: return list_usage();
    }
  }

  const char *root = store_root(store);

  /* A named target: enumerate its sets. */
  if (optind < argc) {
    for (int j = optind; j < argc; ++j) {
      const char *target = argv[j];
      const kycg_seq_reg_t *sr = find_seq(target);
      const kycg_array_reg_t *ar = find_array(target);

      if (sr) {
        /* The Zenodo file list is compiled in, so this works offline. */
        char dir[4096];
        snprintf(dir, sizeof(dir), "%s/%s", root, sr->genome);
        printf("# %s -- Zenodo %s (doi %s)\n", sr->genome, sr->record, sr->doi);
        printf("set\tfile\tsize\tcached\n");
        for (const kycg_zfile_t *f = sr->files; f->name; ++f) {
          char setn[256], path[4300], hb[32];
          set_name_of(f->name, setn, sizeof(setn));
          snprintf(path, sizeof(path), "%s/%s", dir, f->name);
          printf("%s\t%s\t%s\t%s\n", setn, f->name,
                 human(f->size, hb, sizeof(hb)),
                 file_exists(path) ? "yes" : "no");
        }
      } else if (ar) {
        /* Array file lists are not compiled in -- the whole point of anchoring
         * on SHA256SUMS is that upstream can add a set without a kycg rebuild.
         * So we can only enumerate what is already cached. */
        char dir[4096];
        snprintf(dir, sizeof(dir), "%s/%s/KYCG", root, ar->platform);
        char sums[4300];
        snprintf(sums, sizeof(sums), "%s/%s", dir, KYCG_IA_SUMS_FILE);

        printf("# %s -- InfiniumAnnotation tag %s\n", ar->platform, KYCG_IA_TAG);
        FILE *fp = fopen(sums, "rb");
        if (!fp) {
          printf("# not fetched yet; run: kycg fetch %s\n", ar->platform);
          continue;
        }
        printf("set\tfile\tcached\n");
        char line[1024];
        while (fgets(line, sizeof(line), fp)) {
          char *nm = strstr(line, "  ");
          if (!nm) continue;
          nm += 2;
          size_t len = strlen(nm);
          while (len && (nm[len-1] == '\n' || nm[len-1] == '\r')) nm[--len] = '\0';
          if (len < 4 || strcmp(nm + len - 3, ".cm") != 0) continue;
          char setn[256], path[4300];
          set_name_of(nm, setn, sizeof(setn));
          snprintf(path, sizeof(path), "%s/%s", dir, nm);
          printf("%s\t%s\t%s\n", setn, nm, file_exists(path) ? "yes" : "no");
        }
        fclose(fp);
      } else {
        fprintf(stderr, "kycg list: '%s' is not a known platform or genome.\n",
                target);
        return 1;
      }
    }
    return 0;
  }

  /* No target: the overview. */
  printf("# store: %s\n", root);
  printf("target\tkind\tsource\tcached_sets\n");

  for (const kycg_array_reg_t *r = KYCG_ARRAY_REGISTRY; r->platform; ++r) {
    char dir[4096];
    snprintf(dir, sizeof(dir), "%s/%s/KYCG", root, r->platform);
    printf("%s\tarray\tInfiniumAnnotation@%s\t%" PRIu64 "\n",
           r->platform, KYCG_IA_TAG, count_cached(dir));
  }

  for (const kycg_seq_reg_t *r = KYCG_SEQ_REGISTRY; r->genome; ++r) {
    char dir[4096];
    uint64_t avail = 0;
    for (const kycg_zfile_t *f = r->files; f->name; ++f) {
      size_t len = strlen(f->name);
      if (len > 3 && strcmp(f->name + len - 3, ".cm") == 0) ++avail;
    }
    snprintf(dir, sizeof(dir), "%s/%s", root, r->genome);

    uint64_t have = 0;
    for (const kycg_zfile_t *f = r->files; f->name; ++f) {
      size_t len = strlen(f->name);
      if (len <= 3 || strcmp(f->name + len - 3, ".cm") != 0) continue;
      char path[4300];
      snprintf(path, sizeof(path), "%s/%s", dir, f->name);
      if (file_exists(path)) ++have;
    }
    printf("%s\tsequencing\tzenodo:%s\t%" PRIu64 "/%" PRIu64 "\n",
           r->genome, r->record, have, avail);
  }

  return 0;
}
