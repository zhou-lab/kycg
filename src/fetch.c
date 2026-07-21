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
 *   Enrichment testing is only as good as the sets you test against, and those
 *   sets were scattered across three channels with three access patterns and
 *   no verification. This command is the tie.
 *
 * TWO WAYS IN
 *   kycg fetch hg38                 everything for a target
 *   kycg fetch hg38:CGI,ChromHMM    a named subset, no questions asked
 *   kycg list                       browse, check what you want, fetch it
 *
 *   The colon form exists so a pipeline can name exactly what it wants on one
 *   line. The browser exists because nobody memorizes 33 set names -- and
 *   since browsing the catalogue and choosing from it are the same activity,
 *   `kycg list` is where both happen. `kycg fetch` with no target simply opens
 *   it; there is no second guided flow to drift out of step with the first.
 *   See main_list() and fetch_picked() below.
 *
 * PLAN, CONFIRM, EXECUTE
 *   Nothing downloads until a plan has been built and shown: which files, how
 *   large, where they land, and what is already present. This matters because
 *   the hg38 collection is 363 MB and the difference between wanting all of it
 *   and wanting two sets is easy to express and easy to get wrong.
 *
 *   Building the plan for the Zenodo channel is free -- the file list is
 *   compiled in. Building it for the array channel costs one small request for
 *   SHA256SUMS, since the whole point of anchoring on that manifest is that
 *   kycg does not carry the file list and upstream can add a set without a
 *   rebuild. So the array path reaches the network before the confirmation, by
 *   a few kilobytes, and says so while it does.
 *
 * PROMPTING WITHOUT BREAKING AUTOMATION
 *   DESIGN.md originally said "never prompt", because a prompt hangs a
 *   Nextflow job or a Docker build forever with no indication of why. That
 *   guarantee is preserved exactly, by gating every question on an interactive
 *   terminal (see ui.c): off a TTY, an explicit target proceeds without asking
 *   and a missing target is an error rather than a wait. -y forces the same
 *   behavior on a TTY.
 *
 *   What has not changed: kycg still downloads in this command and nowhere
 *   else. No other subcommand touches the network.
 *
 * TRUST
 *   Both channels verify against a digest compiled into this binary by
 *   tools/make_registry.sh; see src/registry.h and src/digest.c for why the two
 *   channels use different hashes. Downloads land on a ".part" sibling and are
 *   renamed only after their digest matches, so an interrupted or corrupted
 *   fetch can never leave a file in the store that later reads as valid.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <inttypes.h>
#include <sys/stat.h>

#include "kycg.h"
#include "digest.h"
#include "registry.h"
#include "store.h"
#include "ui.h"

#ifdef KYCG_HAVE_CURL
#include <curl/curl.h>
#endif

/* ------------------------------------------------------------- set naming */

/**
 * The "set name" of a knowledgebase file: everything before the first dot.
 * "ChromHMM.20220303.cm" -> "ChromHMM". This is what users select and what
 * -o matches, because the dates are how sets are versioned and nobody
 * remembers them.
 */
static void set_name_of(const char *fname, char *out, size_t n) {
  const char *dot = strchr(fname, '.');
  size_t len = dot ? (size_t)(dot - fname) : strlen(fname);
  if (len >= n) len = n - 1;
  memcpy(out, fname, len);
  out[len] = '\0';
}

/** Does this file pass the subset filter? A token matches the whole file name
 *  or its set name. NULL/empty filter passes everything. */
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

/* Group digits: 21867837 -> "21,867,837". Row counts are the one number here
 * a user compares by eye against `kycg info`, so they get separators. */
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

/* ------------------------------------------------------------- registry ops */

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

/* --------------------------------------------------------------- the plan */

typedef struct {
  char     name[512];
  char     url[4096];
  char     sha[65];       /* "" when this channel does not use sha256 */
  char     md5[33];       /* "" when this channel does not use md5    */
  uint64_t size;          /* 0 = not published by this channel        */
  int      have;          /* already present and digest-verified      */
} plan_item_t;

typedef struct {
  plan_item_t *a;
  size_t       n, m;
  char         dir[4096];
  char         target[128];
  char         source[256];
  int          sizes_known;
  char        *sums_text;  /* array channel: the manifest, to write out */
  size_t       sums_len;
} plan_t;

static void plan_free(plan_t *p) {
  free(p->a);
  free(p->sums_text);
  memset(p, 0, sizeof(*p));
}

static plan_item_t *plan_add(plan_t *p) {
  if (p->n == p->m) {
    size_t want = p->m ? p->m * 2 : 64;
    plan_item_t *v = realloc(p->a, want * sizeof(plan_item_t));
    if (!v) return NULL;
    p->a = v; p->m = want;
  }
  plan_item_t *it = &p->a[p->n++];
  memset(it, 0, sizeof(*it));
  return it;
}

/**
 * Mark items already present with the right digest.
 *
 * This hashes what is on disk rather than trusting size or mtime, so a
 * truncated or edited file is re-fetched rather than silently kept. For a
 * fully-populated 363 MB store that costs about a second, which is a fair
 * price for the summary being accurate.
 */
static void plan_check_present(plan_t *p, int force) {
  if (force) return;

  for (size_t i = 0; i < p->n; ++i) {
    plan_item_t *it = &p->a[i];
    char path[4600];
    snprintf(path, sizeof(path), "%s/%s", p->dir, it->name);
    if (!kycg_store_is_file(path)) continue;

    char got[65];
    if (it->sha[0]) {
      if (kycg_sha256_file(path, got) == 0 && kycg_digest_equal(got, it->sha))
        it->have = 1;
    } else if (it->md5[0]) {
      char m[33];
      if (kycg_md5_file(path, m) == 0 && kycg_digest_equal(m, it->md5))
        it->have = 1;
    }
  }
}

/** The brew-style summary shown before anything is downloaded. */
static void plan_show(const plan_t *p) {
  size_t n_todo = 0, n_have = 0;
  uint64_t bytes_todo = 0;

  for (size_t i = 0; i < p->n; ++i) {
    if (p->a[i].have) { ++n_have; continue; }
    ++n_todo;
    bytes_todo += p->a[i].size;
  }

  fprintf(stderr, "\n%s==>%s %sKnowledgebase sets to download%s\n",
          kycg_ui_cyan(), kycg_ui_reset(), kycg_ui_bold(), kycg_ui_reset());
  fprintf(stderr, "    %s%s  %s  %s  %s  %s%s\n\n",
          kycg_ui_dim(), p->target, kycg_ui_bullet(), p->source,
          kycg_ui_bullet(), p->dir, kycg_ui_reset());

  if (!n_todo) {
    fprintf(stderr, "    %severything selected is already present and "
                    "verified.%s\n\n", kycg_ui_dim(), kycg_ui_reset());
    return;
  }

  for (size_t i = 0; i < p->n; ++i) {
    const plan_item_t *it = &p->a[i];
    if (it->have) continue;
    if (it->size) {
      char hb[24];
      fprintf(stderr, "    %-46s %s%8s%s\n", it->name,
              kycg_ui_dim(), kycg_ui_human(it->size, hb, sizeof(hb)),
              kycg_ui_reset());
    } else {
      fprintf(stderr, "    %s\n", it->name);
    }
  }

  fputc('\n', stderr);
  if (p->sizes_known) {
    char hb[24];
    fprintf(stderr, "    %s%zu file(s), %s%s", kycg_ui_bold(), n_todo,
            kycg_ui_human(bytes_todo, hb, sizeof(hb)), kycg_ui_reset());
  } else {
    /* The array channel publishes no sizes; saying so beats inventing them. */
    fprintf(stderr, "    %s%zu file(s)%s%s (sizes not published by this "
                    "channel)%s", kycg_ui_bold(), n_todo, kycg_ui_reset(),
            kycg_ui_dim(), kycg_ui_reset());
  }
  if (n_have)
    fprintf(stderr, "%s   %zu already present, skipped%s",
            kycg_ui_dim(), n_have, kycg_ui_reset());
  fprintf(stderr, "\n\n");
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

/* Drives the spinner and bar from libcurl's transfer callback. No thread is
 * needed: the callback fires often enough to animate on its own. */
static int on_xfer(void *ud, curl_off_t dltotal, curl_off_t dlnow,
                   curl_off_t ultotal, curl_off_t ulnow) {
  (void)ultotal; (void)ulnow;
  kycg_prog_update((kycg_prog_t *)ud, (uint64_t)dlnow, (uint64_t)dltotal);
  return 0;
}

static int http_get_file(const char *url, const char *path, kycg_prog_t *pr) {
  FILE *fp = fopen(path, "wb");
  if (!fp) return -1;

  CURL *h = new_handle(url);
  if (!h) { fclose(fp); return -1; }

  curl_easy_setopt(h, CURLOPT_WRITEDATA, fp);
  curl_easy_setopt(h, CURLOPT_XFERINFOFUNCTION, on_xfer);
  curl_easy_setopt(h, CURLOPT_XFERINFODATA, pr);
  curl_easy_setopt(h, CURLOPT_NOPROGRESS, 0L);

  CURLcode rc = curl_easy_perform(h);
  curl_easy_cleanup(h);
  fclose(fp);

  if (rc != CURLE_OK) { unlink(path); return -1; }
  return 0;
}

#endif /* KYCG_HAVE_CURL */

/* ---------------------------------------------------------- sums parsing */

typedef struct { char sha[65]; char name[512]; } sums_ent_t;

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

/* ------------------------------------------------------------ plan builds */

typedef struct {
  const char *store;
  const char *only;
  const char *tag;
  int dry_run;
  int force;
  int assume_yes;
} fetch_conf_t;

#ifdef KYCG_HAVE_CURL

/** Build the plan for an array platform. Costs one small manifest request. */
static int build_plan_array(const kycg_array_reg_t *reg,
                            const fetch_conf_t *conf, plan_t *plan) {
  if (!reg->sums_sha256) {
    fprintf(stderr, "kycg fetch: platform '%s' is not published at tag %s.\n",
            reg->platform, conf->tag);
    return -1;
  }

  char url[4096];
  snprintf(url, sizeof(url), "%s/%s/%s/KYCG/%s",
           KYCG_IA_BASE_URL, conf->tag, reg->platform, KYCG_IA_SUMS_FILE);

  if (kycg_ui_fancy())
    fprintf(stderr, "%s  reading %s manifest...%s\r",
            kycg_ui_dim(), reg->platform, kycg_ui_reset());

  size_t len = 0;
  char *sums = http_get_mem(url, &len);

  if (kycg_ui_fancy()) fputs("\r\033[2K", stderr);

  if (!sums) {
    fprintf(stderr, "kycg fetch: cannot reach %s\n", url);
    return -1;
  }

  /* The anchor. Everything below is trusted only because this held. */
  char got[65];
  kycg_sha256_buf(sums, len, got);
  if (!kycg_digest_equal(got, reg->sums_sha256)) {
    fprintf(stderr,
            "%s%s%s SHA256SUMS for %s does not match the digest pinned in this "
            "build (tag %s).\n"
            "  expected %s\n  got      %s\n"
            "Refusing to fetch. If upstream has moved, regenerate\n"
            "src/registry.h with tools/make_registry.sh and rebuild.\n",
            kycg_ui_red(), kycg_ui_cross(), kycg_ui_reset(),
            reg->platform, conf->tag, reg->sums_sha256, got);
    free(sums);
    return -1;
  }

  size_t n_ent = 0;
  sums_ent_t *ent = parse_sums(sums, &n_ent);
  if (!ent) { free(sums); return -1; }

  snprintf(plan->dir, sizeof(plan->dir), "%s/%s/KYCG",
           kycg_store_root(conf->store), reg->platform);
  snprintf(plan->target, sizeof(plan->target), "%s", reg->platform);
  snprintf(plan->source, sizeof(plan->source), "InfiniumAnnotation %s",
           conf->tag);
  plan->sizes_known = 0;      /* this channel publishes no sizes */
  plan->sums_text = sums;
  plan->sums_len = len;

  for (size_t i = 0; i < n_ent; ++i) {
    if (!passes_filter(ent[i].name, conf->only)) continue;
    plan_item_t *it = plan_add(plan);
    if (!it) break;
    snprintf(it->name, sizeof(it->name), "%s", ent[i].name);
    snprintf(it->url, sizeof(it->url), "%s/%s/%s/KYCG/%s",
             KYCG_IA_BASE_URL, conf->tag, reg->platform, ent[i].name);
    snprintf(it->sha, sizeof(it->sha), "%s", ent[i].sha);
  }

  free(ent);
  return 0;
}

#endif /* KYCG_HAVE_CURL */

/** Build the plan for a genome. Entirely offline: the list is compiled in. */
static int build_plan_seq(const kycg_seq_reg_t *reg, const fetch_conf_t *conf,
                          plan_t *plan) {
  snprintf(plan->dir, sizeof(plan->dir), "%s/%s",
           kycg_store_root(conf->store), reg->genome);
  snprintf(plan->target, sizeof(plan->target), "%s", reg->genome);
  snprintf(plan->source, sizeof(plan->source), "Zenodo %s", reg->record);
  plan->sizes_known = 1;

  for (const kycg_zfile_t *f = reg->files; f->name; ++f) {
    if (!passes_filter(f->name, conf->only)) continue;
    plan_item_t *it = plan_add(plan);
    if (!it) break;
    snprintf(it->name, sizeof(it->name), "%s", f->name);
    snprintf(it->url, sizeof(it->url), "%s/%s/files/%s",
             KYCG_ZENODO_BASE, reg->record, f->name);
    snprintf(it->md5, sizeof(it->md5), "%s", f->md5);
    it->size = f->size;
  }
  return 0;
}

/* -------------------------------------------------------------- execution */

typedef struct {
  uint64_t n_got, n_skip, n_fail;
  uint64_t bytes_got;
} tally_t;

#ifdef KYCG_HAVE_CURL

static int execute_plan(const plan_t *plan, tally_t *t) {
  if (kycg_store_mkdir_p(plan->dir) != 0) {
    fprintf(stderr, "kycg fetch: cannot create %s\n", plan->dir);
    return -1;
  }

  for (size_t i = 0; i < plan->n; ++i) {
    const plan_item_t *it = &plan->a[i];

    if (it->have) {
      ++t->n_skip;
      continue;
    }

    char path[4600], part[4700];
    snprintf(path, sizeof(path), "%s/%s", plan->dir, it->name);
    snprintf(part, sizeof(part), "%s.part", path);

    kycg_prog_t pr;
    kycg_prog_begin(&pr, it->name, it->size);

    if (http_get_file(it->url, part, &pr) != 0) {
      kycg_prog_done(&pr, "download failed", 0);
      ++t->n_fail;
      continue;
    }

    int ok = 0;
    char got[65];
    if (it->sha[0]) {
      if (kycg_sha256_file(part, got) == 0) ok = kycg_digest_equal(got, it->sha);
    } else if (it->md5[0]) {
      char m[33];
      if (kycg_md5_file(part, m) == 0) ok = kycg_digest_equal(m, it->md5);
    }

    if (!ok) {
      unlink(part);
      kycg_prog_done(&pr, "digest mismatch - discarded", 0);
      ++t->n_fail;
      continue;
    }

    struct stat st;
    uint64_t sz = (stat(part, &st) == 0) ? (uint64_t)st.st_size : 0;

    if (rename(part, path) != 0) {
      unlink(part);
      kycg_prog_done(&pr, "could not move into the store", 0);
      ++t->n_fail;
      continue;
    }

    char hb[24];
    kycg_prog_done(&pr, kycg_ui_human(sz, hb, sizeof(hb)), 1);
    ++t->n_got;
    t->bytes_got += sz;
  }

  /* Keep a manifest so the store re-verifies with shasum and no kycg code. */
  char sp[4700];
  snprintf(sp, sizeof(sp), "%s/%s", plan->dir, KYCG_IA_SUMS_FILE);

  if (plan->sums_text) {
    FILE *fp = fopen(sp, "wb");
    if (fp) { fwrite(plan->sums_text, 1, plan->sums_len, fp); fclose(fp); }
  } else {
    /* Zenodo publishes md5 only, so we compute the sha256 side ourselves. */
    FILE *fp = fopen(sp, "wb");
    if (fp) {
      for (size_t i = 0; i < plan->n; ++i) {
        char path[4600], sha[65];
        snprintf(path, sizeof(path), "%s/%s", plan->dir, plan->a[i].name);
        if (kycg_store_is_file(path) && kycg_sha256_file(path, sha) == 0)
          fprintf(fp, "%s  %s\n", sha, plan->a[i].name);
      }
      fclose(fp);
    }
  }

  return 0;
}

/* ----------------------------------------------------------- the guided run */

#endif /* KYCG_HAVE_CURL */

/* ------------------------------------------------------------------ usage */

static int usage(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg fetch [options] [<target>[:<sets>] ...]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Download and verify knowledgebases into a local store.\n");
  fprintf(stderr, "With no target on a terminal, opens the `kycg list` browser,\n");
  fprintf(stderr, "where sets are checked with space and fetched with f.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Examples:\n");
  fprintf(stderr, "    kycg fetch                     browse and pick\n");
  fprintf(stderr, "    kycg fetch hg38                every set for a target\n");
  fprintf(stderr, "    kycg fetch hg38:CGI,ChromHMM   just those sets\n");
  fprintf(stderr, "    kycg fetch -n hg38             show the plan, download nothing\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Targets:\n");
  fprintf(stderr, "    genomes   ");
  for (const kycg_seq_reg_t *r = KYCG_SEQ_REGISTRY; r->genome; ++r)
    fprintf(stderr, "%s ", r->genome);
  fprintf(stderr, "\n    arrays    ");
  for (const kycg_array_reg_t *r = KYCG_ARRAY_REGISTRY; r->platform; ++r)
    fprintf(stderr, "%s ", r->platform);
  fprintf(stderr, "\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -d DIR    store directory [$KYCG_DATA_DIR, else ~/.cache/kycg]\n");
  fprintf(stderr, "    -o SETS   subset by set name; same as the :SETS suffix\n");
  fprintf(stderr, "    -y        assume yes; do not ask to confirm\n");
  fprintf(stderr, "    -n        dry run: show the plan, download nothing\n");
  fprintf(stderr, "    -f        re-download even if present and verified\n");
  fprintf(stderr, "    -t TAG    InfiniumAnnotation tag, arrays only [%s]\n", KYCG_IA_TAG);
  fprintf(stderr, "    -h        this help\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Fetched sets are ordinary files; pass one to `kycg test -m`.\n");
  fprintf(stderr, "kycg downloads here and nowhere else, and never asks anything\n");
  fprintf(stderr, "when stdin is not a terminal.\n");
  fprintf(stderr, "\n");
  return 1;
}

/* ------------------------------------------------------------------ driver */

int main_fetch(int argc, char *argv[]) {
  fetch_conf_t conf = {0};
  conf.tag = KYCG_IA_TAG;

  int c;
  while ((c = getopt(argc, argv, "d:o:t:nfyh")) >= 0) {
    switch (c) {
    case 'd': conf.store = optarg; break;
    case 'o': conf.only = optarg; break;
    case 't': conf.tag = optarg; break;
    case 'n': conf.dry_run = 1; break;
    case 'f': conf.force = 1; break;
    case 'y': conf.assume_yes = 1; break;
    case 'h': return usage();
    default: return usage();
    }
  }

#ifndef KYCG_HAVE_CURL
  (void)argc;
  fprintf(stderr,
          "kycg fetch: this build has no network support.\n"
          "kycg was compiled without libcurl, so fetch is unavailable. Install\n"
          "libcurl development headers and rebuild, or populate the store by\n"
          "hand -- fetched files are ordinary .cm files and `kycg test -m` takes\n"
          "any path. See `kycg list` for the expected layout.\n");
  return 1;
#else
  curl_global_init(CURL_GLOBAL_DEFAULT);

  /* Collect targets, either from argv or from the guided run. */
  const char *argv_targets[64];
  int n_targets = 0;

  if (optind < argc) {
    for (int j = optind; j < argc && n_targets < 64; ++j)
      argv_targets[n_targets++] = argv[j];
  } else if (kycg_ui_interactive()) {
    /* No target on a terminal: hand over to the browser in `kycg list`, which
     * is the interactive fetch surface -- browse, check, fetch, one screen.
     * Keeping a second guided flow here would be a worse copy of it. */
    char *lav[4];
    int lac = 0;
    lav[lac++] = "list";
    if (conf.store) { lav[lac++] = "-d"; lav[lac++] = (char *)conf.store; }
    lav[lac] = NULL;
    optind = 1;
    curl_global_cleanup();
    return main_list(lac, lav);
  } else {
    usage();
    fprintf(stderr,
            "kycg fetch: no target given, and stdin is not a terminal so there\n"
            "is nobody to ask. Name a target explicitly, e.g. `kycg fetch hg38`.\n");
    curl_global_cleanup();
    return 1;
  }

  tally_t t = {0};
  int rc = 0;

  for (int j = 0; j < n_targets; ++j) {
    /* Split "hg38:CGI,ChromHMM" into target and subset. */
    char target[128];
    const char *spec = argv_targets[j];
    const char *colon = strchr(spec, ':');
    const char *only = conf.only;

    if (colon) {
      size_t len = (size_t)(colon - spec);
      if (len >= sizeof(target)) len = sizeof(target) - 1;
      memcpy(target, spec, len);
      target[len] = '\0';
      only = colon + 1;
      if (!*only) only = NULL;      /* "hg38:" means everything */
    } else {
      snprintf(target, sizeof(target), "%s", spec);
    }

    fetch_conf_t tc = conf;
    tc.only = only;

    const kycg_array_reg_t *ar = find_array(target);
    const kycg_seq_reg_t   *sr = find_seq(target);

    if (!ar && !sr) {
      fprintf(stderr,
              "kycg fetch: '%s' is not a known platform or genome.\n"
              "Run `kycg list` to see what is available.\n", target);
      rc = 1;
      continue;
    }

    plan_t plan = {0};
    int prc = sr ? build_plan_seq(sr, &tc, &plan)
                 : build_plan_array(ar, &tc, &plan);
    if (prc != 0) { plan_free(&plan); rc = 1; continue; }

    if (!plan.n) {
      fprintf(stderr, "kycg fetch: nothing in '%s' matches that selection.\n",
              target);
      plan_free(&plan);
      rc = 1;
      continue;
    }

    plan_check_present(&plan, tc.force);
    plan_show(&plan);

    if (tc.dry_run) { plan_free(&plan); continue; }

    size_t n_todo = 0;
    for (size_t i = 0; i < plan.n; ++i) if (!plan.a[i].have) ++n_todo;
    if (!n_todo) { t.n_skip += plan.n; plan_free(&plan); continue; }

    /* The confirmation. Skipped when nobody can answer, which is what keeps
     * this safe to run inside a pipeline. */
    if (!tc.assume_yes && kycg_ui_interactive()) {
      if (!kycg_ui_confirm("Proceed?", 1)) {
        fprintf(stderr, "Cancelled.\n");
        plan_free(&plan);
        continue;
      }
      fputc('\n', stderr);
    }

    if (execute_plan(&plan, &t) != 0) rc = 1;
    plan_free(&plan);
  }

  if (!conf.dry_run && (t.n_got || t.n_skip || t.n_fail)) {
    char hb[24];
    fprintf(stderr, "\n%s%s%s %" PRIu64 " fetched (%s)",
            kycg_ui_green(), kycg_ui_check(), kycg_ui_reset(),
            t.n_got, kycg_ui_human(t.bytes_got, hb, sizeof(hb)));
    if (t.n_skip) fprintf(stderr, ", %" PRIu64 " already current", t.n_skip);
    if (t.n_fail)
      fprintf(stderr, ", %s%" PRIu64 " FAILED%s",
              kycg_ui_red(), t.n_fail, kycg_ui_reset());
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
  fprintf(stderr, "Browse knowledgebase collections and fetch from them.\n");
  fprintf(stderr, "On a terminal this is an interactive tree: arrows move,\n");
  fprintf(stderr, "right unfolds a target, space checks a set, f fetches the\n");
  fprintf(stderr, "checked ones, q quits. Redirect stdout for plain TSV.\n");
  fprintf(stderr, "With a target named, lists the individual sets it carries.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -d DIR    store directory [$KYCG_DATA_DIR, else ~/.cache/kycg]\n");
  fprintf(stderr, "    -h        this help\n");
  fprintf(stderr, "\n");
  return 1;
}

/*
 * `kycg list` output is buffered rather than printed directly, so it can be
 * handed to the in-place browser when someone is watching and written as plain
 * TSV when it is not. The distinction is stdout: a redirected stdout means the
 * caller wants data, and turning that into a full-screen widget would break
 * every script that pipes this into cut or awk.
 */
typedef struct { char **a; unsigned char *st; size_t n, m; } rows_t;

static void rows_push(rows_t *r, unsigned char style, const char *fmt, ...) {
  if (r->n == r->m) {
    size_t want = r->m ? r->m * 2 : 64;
    char **v = realloc(r->a, want * sizeof(char *));
    if (!v) return;
    r->a = v;
    unsigned char *sv = realloc(r->st, want);
    if (!sv) return;
    r->st = sv;
    r->m = want;
  }
  r->st[r->n] = style;
  char buf[2048];
  va_list ap;
  va_start(ap, fmt);
  vsnprintf(buf, sizeof(buf), fmt, ap);
  va_end(ap);
  r->a[r->n++] = strdup(buf);
}

static void rows_free(rows_t *r) {
  for (size_t i = 0; i < r->n; ++i) free(r->a[i]);
  free(r->a);
  free(r->st);
  memset(r, 0, sizeof(*r));
}

/** Browse if stdout is a terminal, else write TSV. `comment` is a leading
 *  '#' line preserved in plain mode and folded into the title when browsing. */
static void rows_emit(rows_t *r, const char *comment, const char *header) {
  if (isatty(STDOUT_FILENO) &&
      kycg_ui_browse(comment ? comment : "kycg list", header,
                     (const char **)r->a, r->st, r->n) == 0) {
    rows_free(r);
    return;
  }
  if (comment) printf("# %s\n", comment);
  if (header) printf("%s\n", header);
  for (size_t i = 0; i < r->n; ++i) printf("%s\n", r->a[i]);
  rows_free(r);
}

/*
 * `kycg list` is also the fetch picker, so its expand and accept callbacks
 * share one context: the store root they read, and the selection they build.
 */
typedef struct {
  char **target;   /* parallel arrays, one entry per checked set */
  char **file;
  size_t n, m;
} picks_t;

typedef struct {
  char root[4096];   /* the store being browsed; 'd' can change it */
  picks_t picks;
  rows_t *rows;      /* the overview, rewritten in place after any change */
  char  *title;      /* the widget's title buffer, rewritten with the root */
  size_t title_sz;
} listctx_t;

static void picks_add(picks_t *p, const char *target, const char *file) {
  if (p->n == p->m) {
    size_t want = p->m ? p->m * 2 : 32;
    char **t = realloc(p->target, want * sizeof(char *));
    char **f = realloc(p->file, want * sizeof(char *));
    if (!t || !f) { free(t); free(f); return; }
    p->target = t; p->file = f; p->m = want;
  }
  p->target[p->n] = strdup(target);
  p->file[p->n] = strdup(file);
  ++p->n;
}

static void picks_free(picks_t *p) {
  for (size_t i = 0; i < p->n; ++i) { free(p->target[i]); free(p->file[i]); }
  free(p->target); free(p->file);
  memset(p, 0, sizeof(*p));
}

/** Records one checked row; the target is the first field of its parent. */
static void on_pick(void *ctx, const char *root, const char *key) {
  listctx_t *lc = ctx;
  char target[128];
  const char *tab = strchr(root, '\t');
  size_t len = tab ? (size_t)(tab - root) : strlen(root);
  if (len >= sizeof(target)) len = sizeof(target) - 1;
  memcpy(target, root, len);
  target[len] = '\0';
  picks_add(&lc->picks, target, key);
}

/**
 * A platform's SHA256SUMS, fetched and verified, memoized for this run.
 *
 * Returns a malloc'd copy the caller frees, or NULL when there is no network
 * (or no libcurl). The digest is checked against the compiled anchor exactly
 * as a fetch would: a catalogue is worth no more trust than the files in it.
 */
static char *array_manifest(const kycg_array_reg_t *ar, size_t *len_out) {
#ifndef KYCG_HAVE_CURL
  (void)ar; (void)len_out;
  return NULL;
#else
  /* One slot per platform; the registry is small and fixed. */
  static struct { const char *plat; char *text; size_t len; } memo[16];
  static size_t n_memo = 0;

  for (size_t i = 0; i < n_memo; ++i) {
    if (memo[i].plat != ar->platform) continue;
    if (!memo[i].text) return NULL;
    char *dup = malloc(memo[i].len + 1);
    if (!dup) return NULL;
    memcpy(dup, memo[i].text, memo[i].len + 1);
    if (len_out) *len_out = memo[i].len;
    return dup;
  }

  char *text = NULL;
  size_t len = 0;

  if (ar->sums_sha256) {
    char url[4096];
    snprintf(url, sizeof(url), "%s/%s/%s/KYCG/%s",
             KYCG_IA_BASE_URL, KYCG_IA_TAG, ar->platform, KYCG_IA_SUMS_FILE);
    text = http_get_mem(url, &len);

    if (text) {
      char got[65];
      kycg_sha256_buf(text, len, got);
      if (!kycg_digest_equal(got, ar->sums_sha256)) { free(text); text = NULL; }
    }
  }

  if (n_memo < 16) {
    memo[n_memo].plat = ar->platform;
    memo[n_memo].len = len;
    memo[n_memo].text = NULL;
    if (text) {
      memo[n_memo].text = malloc(len + 1);
      if (memo[n_memo].text) memcpy(memo[n_memo].text, text, len + 1);
    }
    ++n_memo;
  }

  if (len_out) *len_out = len;
  return text;
#endif
}

/* Append one preformatted child line to an expansion. */
static void kid_push(kycg_ui_kids_t *k, unsigned char style, const char *key,
                     const char *fmt, ...) {
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

/**
 * Children of one row of the `kycg list` overview: the sets that target holds.
 *
 * Widths are fixed rather than measured, because the tree renders children one
 * parent at a time — columns that re-measured on each expansion would shift
 * under rows already on screen.
 *
 * Whole-genome targets can be listed offline, since the file list is compiled
 * in. Array platforms cannot: anchoring on SHA256SUMS is what lets upstream
 * add a set without a kycg rebuild, so the set list only exists once fetched.
 */
static void expand_target(void *ctx, const char *row, kycg_ui_kids_t *out) {
  const char *root = ((listctx_t *)ctx)->root;

  char target[128];
  const char *tab = strchr(row, '\t');
  size_t len = tab ? (size_t)(tab - row) : strlen(row);
  if (len >= sizeof(target)) len = sizeof(target) - 1;
  memcpy(target, row, len);
  target[len] = '\0';

  const kycg_seq_reg_t *sr = find_seq(target);
  const kycg_array_reg_t *ar = find_array(target);

  if (sr) {
    char dir[4096];
    snprintf(dir, sizeof(dir), "%s/%s", root, sr->genome);
    for (const kycg_zfile_t *fl = sr->files; fl->name; ++fl) {
      size_t l = strlen(fl->name);
      if (l <= 3 || strcmp(fl->name + l - 3, ".cm") != 0) continue;
      char setn[256], path[4400], hb[24];
      set_name_of(fl->name, setn, sizeof(setn));
      snprintf(path, sizeof(path), "%s/%s", dir, fl->name);
      int have = kycg_store_is_file(path);
      kid_push(out, have ? KYCG_ROW_HAVE : KYCG_ROW_MISSING, fl->name,
               "%-22.22s %-32.32s %9s  %s", setn, fl->name,
               kycg_ui_human(fl->size, hb, sizeof(hb)),
               have ? "cached" : "-");
    }
    return;
  }

  if (ar) {
    char dir[4096], sums[4400];
    snprintf(dir, sizeof(dir), "%s/%s/KYCG", root, ar->platform);
    snprintf(sums, sizeof(sums), "%s/%s", dir, KYCG_IA_SUMS_FILE);

    /* The catalogue lives in the manifest, and an unfetched platform has no
     * local copy of it -- which made every array look as though it had no
     * knowledgebases at all, when it simply had none *here*. So pull the
     * manifest (a couple of kilobytes, verified against the compiled anchor)
     * to list what exists. It is cached in memory for the session and never
     * written to the store, because writing it would claim files are present
     * that are not. */
    char *text = NULL;
    size_t len = 0;

    FILE *fp = fopen(sums, "rb");
    if (fp) {
      fseek(fp, 0, SEEK_END);
      long sz = ftell(fp);
      fseek(fp, 0, SEEK_SET);
      if (sz > 0) {
        text = malloc((size_t)sz + 1);
        if (text) { len = fread(text, 1, (size_t)sz, fp); text[len] = '\0'; }
      }
      fclose(fp);
    } else {
      text = array_manifest(ar, &len);
    }

    if (!text) {
      kid_push(out, KYCG_ROW_MISSING, NULL,
               "catalogue unavailable - run: kycg fetch %s", ar->platform);
      return;
    }

    size_t n_ent = 0;
    sums_ent_t *ent = parse_sums(text, &n_ent);
    for (size_t i = 0; ent && i < n_ent; ++i) {
      const char *nm = ent[i].name;
      size_t l = strlen(nm);
      if (l < 4 || strcmp(nm + l - 3, ".cm") != 0) continue;
      char setn[256], path[4400];
      set_name_of(nm, setn, sizeof(setn));
      snprintf(path, sizeof(path), "%s/%s", dir, nm);
      int have = kycg_store_is_file(path);
      kid_push(out, have ? KYCG_ROW_HAVE : KYCG_ROW_MISSING, nm,
               "%-22.22s %-32.32s %9s  %s", setn, nm, "",
               have ? "cached" : "-");
    }
    free(ent);
    free(text);
  }
}

/**
 * Fetch what the picker checked: one plan per target, confirmed once.
 *
 * The selection arrives as (target, file name) pairs, which is turned back
 * into the same comma-separated subset string the -o flag takes, so this path
 * and `kycg fetch hg38:CGI` build their plans through identical code. Index
 * sidecars are added alongside their .cm, since a set's index is part of the
 * set and picking one row should not leave half of it behind.
 */
static int fetch_picked(const picks_t *picks, const char *store) {
#ifndef KYCG_HAVE_CURL
  (void)picks; (void)store;
  fprintf(stderr,
          "kycg: this build has no network support (compiled without "
          "libcurl), so the selection cannot be fetched.\n");
  return 1;
#else
  curl_global_init(CURL_GLOBAL_DEFAULT);

  /* Four lines at the foot of the browser: what is happening, the last file
   * finished, the one in flight, and the prompt or tally. The catalogue stays
   * visible above, so fetching reads as part of browsing rather than an
   * errand it sent you on. */
  kycg_ui_panel_open(4);

  tally_t t = {0};
  int rc = 0;
  char **done = calloc(picks->n, sizeof(char *));
  size_t n_done = 0;

  for (size_t i = 0; i < picks->n; ++i) {
    /* One plan per distinct target. */
    int seen = 0;
    for (size_t d = 0; d < n_done; ++d)
      if (strcmp(done[d], picks->target[i]) == 0) { seen = 1; break; }
    if (seen) continue;
    done[n_done++] = picks->target[i];

    size_t cap = 1;
    for (size_t j = 0; j < picks->n; ++j)
      if (strcmp(picks->target[j], picks->target[i]) == 0)
        cap += strlen(picks->file[j]) * 2 + 12;

    char *only = malloc(cap);
    if (!only) { rc = 1; break; }
    only[0] = '\0';
    for (size_t j = 0; j < picks->n; ++j) {
      if (strcmp(picks->target[j], picks->target[i]) != 0) continue;
      if (only[0]) strcat(only, ",");
      strcat(only, picks->file[j]);
      strcat(only, ",");
      strcat(only, picks->file[j]);
      strcat(only, ".idx");
    }

    fetch_conf_t conf = {0};
    conf.tag = KYCG_IA_TAG;
    conf.store = store;
    conf.only = only;

    const kycg_seq_reg_t   *sr = find_seq(picks->target[i]);
    const kycg_array_reg_t *ar = find_array(picks->target[i]);

    plan_t plan = {0};
    int prc = sr ? build_plan_seq(sr, &conf, &plan)
                 : (ar ? build_plan_array(ar, &conf, &plan) : -1);
    free(only);
    if (prc != 0 || !plan.n) { plan_free(&plan); rc = 1; continue; }

    plan_check_present(&plan, 0);

    size_t n_todo = 0;
    uint64_t bytes = 0;
    for (size_t k = 0; k < plan.n; ++k)
      if (!plan.a[k].have) { ++n_todo; bytes += plan.a[k].size; }
    if (!n_todo) { t.n_skip += plan.n; plan_free(&plan); continue; }

    char hb[24];
    kycg_ui_panel_line(0, "  %s%s%s  %s  %zu file(s)%s%s%s  %s  %s%s%s",
                       kycg_ui_bold(), plan.target, kycg_ui_reset(),
                       kycg_ui_bullet(), n_todo,
                       plan.sizes_known ? ", " : "",
                       plan.sizes_known ? kycg_ui_human(bytes, hb, sizeof(hb)) : "",
                       "", kycg_ui_bullet(),
                       kycg_ui_dim(), plan.dir, kycg_ui_reset());

    if (!kycg_ui_panel_confirm(3, "Proceed?", 1)) {
      kycg_ui_panel_line(3, "  %scancelled%s", kycg_ui_dim(), kycg_ui_reset());
      plan_free(&plan);
      continue;
    }
    kycg_ui_panel_line(3, " ");

    if (execute_plan(&plan, &t) != 0) rc = 1;
    plan_free(&plan);
  }

  free(done);

  if (t.n_got || t.n_skip || t.n_fail) {
    char hb[24], line[512];
    int o = snprintf(line, sizeof(line), "%" PRIu64 " fetched (%s)",
                     t.n_got, kycg_ui_human(t.bytes_got, hb, sizeof(hb)));
    if (t.n_skip)
      o += snprintf(line + o, sizeof(line) - (size_t)o,
                    ", %" PRIu64 " already current", t.n_skip);
    if (t.n_fail)
      snprintf(line + o, sizeof(line) - (size_t)o,
               ", %" PRIu64 " FAILED", t.n_fail);

    kycg_ui_panel_line(0, "  %s%s%s %s%s",
                       t.n_fail ? kycg_ui_red() : kycg_ui_green(),
                       t.n_fail ? kycg_ui_cross() : kycg_ui_check(),
                       kycg_ui_reset(), line, kycg_ui_reset());
    kycg_ui_panel_pause(3, "press any key to return to the browser");
  }

  kycg_ui_panel_close();
  curl_global_cleanup();
  return (rc || t.n_fail) ? 1 : 0;
#endif
}

static uint64_t count_cached(const char *dir);

/**
 * Fill `rows` with one line per target, counting what the store holds.
 *
 * Called for the initial view and again after anything changes it -- a fetch,
 * or a different store directory. Rebuilding is cheap: it stats files, and
 * the target list is fixed and short.
 */
static void build_overview(const char *root, rows_t *rows) {
  for (const kycg_seq_reg_t *r = KYCG_SEQ_REGISTRY; r->genome; ++r) {
    char dir[4096];
    snprintf(dir, sizeof(dir), "%s/%s", root, r->genome);

    uint64_t avail = 0, have = 0;
    for (const kycg_zfile_t *f = r->files; f->name; ++f) {
      size_t len = strlen(f->name);
      if (len <= 3 || strcmp(f->name + len - 3, ".cm") != 0) continue;
      ++avail;
      char path[4400];
      snprintf(path, sizeof(path), "%s/%s", dir, f->name);
      if (kycg_store_is_file(path)) ++have;
    }
    /* Green means "something here is usable", dim means "nothing yet". How
     * much is in the count. Distinguishing complete from partial by colour
     * would read as a different meaning on the array rows below, where the
     * total is unknowable by design -- anchoring on SHA256SUMS is what keeps
     * the file list out of the binary -- so the rule stays the same for both. */
    char rb[32];
    rows_push(rows, have ? KYCG_ROW_HAVE : KYCG_ROW_MISSING,
              "%s\twhole genome\t%s\tzenodo:%s\t%" PRIu64 "/%" PRIu64,
              r->genome, commify(r->rows, rb, sizeof(rb)),
              r->record, have, avail);
  }

  for (const kycg_array_reg_t *r = KYCG_ARRAY_REGISTRY; r->platform; ++r) {
    char dir[4096];
    snprintf(dir, sizeof(dir), "%s/%s/KYCG", root, r->platform);
    uint64_t nc = count_cached(dir);
    /* Arrays report only what is cached, not a total: the catalogue lives in
     * the SHA256SUMS manifest, which is the point of anchoring on it. Unfold
     * a platform and the browser fetches that manifest to show the rest. */
    char rb[32];
    rows_push(rows, nc ? KYCG_ROW_HAVE : KYCG_ROW_MISSING,
              "%s\tarray\t%s\tInfiniumAnnotation@%s\t%" PRIu64,
              r->platform, commify(r->rows, rb, sizeof(rb)), KYCG_IA_TAG, nc);
  }
}

/**
 * Rewrite the overview in place.
 *
 * The tree holds the addresses of the row and style arrays, so this refills
 * them rather than replacing them. The target list is fixed, so the counts
 * change but the row count never does.
 */
static void refresh_overview(listctx_t *lc) {
  for (size_t i = 0; i < lc->rows->n; ++i) free(lc->rows->a[i]);
  lc->rows->n = 0;
  build_overview(lc->root, lc->rows);
}

/**
 * `f` in the browser: fetch everything checked, drawing into the panel so the
 * catalogue stays on screen, then refresh the counts.
 */
static void on_commit(void *ctx) {
  listctx_t *lc = ctx;
  if (lc->picks.n) fetch_picked(&lc->picks, lc->root);
  picks_free(&lc->picks);
  refresh_overview(lc);
}

/**
 * `d` in the browser: point it at a different store.
 *
 * The store location is the one thing that cannot be changed from inside
 * otherwise -- every other question the browser answers is about its contents.
 */
static int on_list_key(void *ctx, char key) {
  listctx_t *lc = ctx;
  if (key != 'd') return 0;

  char buf[4096];
  snprintf(buf, sizeof(buf), "%s", lc->root);

  kycg_ui_panel_open(3);
  kycg_ui_panel_line(0, "  %s%s%s", kycg_ui_dim(),
                     "store directory (enter to accept, esc to cancel)",
                     kycg_ui_reset());
  int ok = kycg_ui_panel_ask(1, "store:", buf, sizeof(buf));
  kycg_ui_panel_close();

  if (!ok || !buf[0] || strcmp(buf, lc->root) == 0) return 0;

  snprintf(lc->root, sizeof(lc->root), "%s", buf);
  /* The title names the store, so it goes stale with it. */
  if (lc->title) snprintf(lc->title, lc->title_sz, "store: %s", lc->root);
  refresh_overview(lc);
  return 1;
}

static uint64_t count_cached(const char *dir) {
  char sums[4600];
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
      char path[4700];
      snprintf(path, sizeof(path), "%s/%s", dir, nm);
      if (kycg_store_is_file(path)) ++n;
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

  const char *root = kycg_store_root(store);

  if (optind < argc) {
    for (int j = optind; j < argc; ++j) {
      const char *target = argv[j];
      const kycg_seq_reg_t *sr = find_seq(target);
      const kycg_array_reg_t *ar = find_array(target);

      if (sr) {
        char dir[4096], title[256];
        snprintf(dir, sizeof(dir), "%s/%s", root, sr->genome);
        char rb[32];
        snprintf(title, sizeof(title), "%s -- %s rows -- Zenodo %s (doi %s)",
                 sr->genome, commify(sr->rows, rb, sizeof(rb)),
                 sr->record, sr->doi);

        rows_t rows = {0};
        for (const kycg_zfile_t *f = sr->files; f->name; ++f) {
          char setn[256], path[4400], hb[24];
          set_name_of(f->name, setn, sizeof(setn));
          snprintf(path, sizeof(path), "%s/%s", dir, f->name);
          int have = kycg_store_is_file(path);
          rows_push(&rows, have ? KYCG_ROW_HAVE : KYCG_ROW_MISSING,
                    "%s\t%s\t%s\t%s", setn, f->name,
                    kycg_ui_human(f->size, hb, sizeof(hb)),
                    have ? "yes" : "no");
        }
        rows_emit(&rows, title, "set\tfile\tsize\tcached");
      } else if (ar) {
        char dir[4096], title[256];
        snprintf(dir, sizeof(dir), "%s/%s/KYCG", root, ar->platform);
        char sums[4400];
        snprintf(sums, sizeof(sums), "%s/%s", dir, KYCG_IA_SUMS_FILE);
        char rb[32];
        snprintf(title, sizeof(title), "%s -- %s rows -- InfiniumAnnotation %s",
                 ar->platform, commify(ar->rows, rb, sizeof(rb)), KYCG_IA_TAG);

        /* Read the local manifest if the platform is fetched. If it is not,
         * pull the catalogue only when someone is watching -- a redirected
         * stdout means a script is reading, and a script must not trigger a
         * download it did not ask for. */
        char *text = NULL;
        size_t tlen = 0;
        FILE *fp = fopen(sums, "rb");
        if (fp) {
          fseek(fp, 0, SEEK_END);
          long sz = ftell(fp);
          fseek(fp, 0, SEEK_SET);
          if (sz > 0) {
            text = malloc((size_t)sz + 1);
            if (text) { tlen = fread(text, 1, (size_t)sz, fp); text[tlen] = '\0'; }
          }
          fclose(fp);
        } else if (isatty(STDOUT_FILENO)) {
          text = array_manifest(ar, &tlen);
        }

        if (!text) {
          printf("# %s\n# not fetched yet; run: kycg fetch %s\n",
                 title, ar->platform);
          continue;
        }

        rows_t rows = {0};
        size_t n_ent = 0;
        sums_ent_t *ent = parse_sums(text, &n_ent);
        for (size_t i = 0; ent && i < n_ent; ++i) {
          const char *nm = ent[i].name;
          size_t len = strlen(nm);
          if (len < 4 || strcmp(nm + len - 3, ".cm") != 0) continue;
          char setn[256], path[4400];
          set_name_of(nm, setn, sizeof(setn));
          snprintf(path, sizeof(path), "%s/%s", dir, nm);
          int have = kycg_store_is_file(path);
          rows_push(&rows, have ? KYCG_ROW_HAVE : KYCG_ROW_MISSING,
                    "%s\t%s\t%s", setn, nm, have ? "yes" : "no");
        }
        free(ent);
        free(text);
        rows_emit(&rows, title, "set\tfile\tcached");
      } else {
        fprintf(stderr, "kycg list: '%s' is not a known platform or genome.\n",
                target);
        return 1;
      }
    }
    return 0;
  }

  rows_t rows = {0};
  build_overview(root, &rows);

  char title[4200];
  snprintf(title, sizeof(title), "store: %s", root);

  /* On a terminal the overview unfolds and doubles as the fetch picker: the
   * cached_sets column says how many a target holds, the obvious next question
   * is which ones, and the one after that is "get me those". Answering all
   * three in one screen beats re-running with a target named and a -o list
   * typed from memory. */
  if (isatty(STDOUT_FILENO)) {
    listctx_t lc = {0};
    snprintf(lc.root, sizeof(lc.root), "%s", root);
    lc.rows = &rows;
    lc.title = title;
    lc.title_sz = sizeof(title);

    kycg_ui_tree_t spec = {0};
    spec.title = title;
    spec.header = "target\tkind\trows\tsource\tcached_sets";
    spec.roots = rows.a;
    spec.root_styles = rows.st;
    spec.n_roots = rows.n;
    spec.expand = expand_target;
    spec.accept = on_pick;
    spec.commit = on_commit;
    spec.on_key = on_list_key;
    spec.hint = "d store";
    spec.ctx = &lc;

    int rc = kycg_ui_tree(&spec);

    if (rc >= 0) {          /* the widget ran; plain output is not wanted */
      picks_free(&lc.picks);
      rows_free(&rows);
      return 0;
    }
    picks_free(&lc.picks);  /* -1: terminal cannot host it, print plainly */
  }

  rows_emit(&rows, title, "target\tkind\trows\tsource\tcached_sets");

  return 0;
}
