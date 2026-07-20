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
 * THREE WAYS IN
 *   kycg fetch hg38                 everything for a target
 *   kycg fetch hg38:CGI,ChromHMM    a named subset, no questions asked
 *   kycg fetch                      guided: location, target, sets, confirm
 *
 *   The colon form exists so a pipeline can name exactly what it wants on one
 *   line. The guided form exists because nobody memorizes 33 set names.
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

/** Distinct set names in a plan, with a size note per set. Caller frees. */
static size_t collect_sets(const plan_t *plan, char ***names_out,
                           char ***notes_out) {
  char **names = calloc(plan->n, sizeof(char *));
  uint64_t *sz = calloc(plan->n, sizeof(uint64_t));
  size_t *cnt = calloc(plan->n, sizeof(size_t));
  size_t n = 0;
  if (!names || !sz || !cnt) { free(names); free(sz); free(cnt); return 0; }

  for (size_t i = 0; i < plan->n; ++i) {
    char s[256];
    set_name_of(plan->a[i].name, s, sizeof(s));

    size_t k;
    for (k = 0; k < n; ++k) if (strcmp(names[k], s) == 0) break;
    if (k == n) { names[n] = strdup(s); ++n; }
    sz[k] += plan->a[i].size;
    cnt[k] += 1;
  }

  char **notes = calloc(n, sizeof(char *));
  for (size_t k = 0; k < n; ++k) {
    char buf[128], hb[24];
    if (sz[k]) kycg_ui_human(sz[k], hb, sizeof(hb));
    else snprintf(hb, sizeof(hb), "size n/a");
    snprintf(buf, sizeof(buf), "%s", hb);
    notes[k] = strdup(buf);
  }

  free(sz); free(cnt);
  *names_out = names;
  *notes_out = notes;
  return n;
}

/**
 * The guided run: where, what, which sets.
 *
 * Only reachable on an interactive terminal; main_fetch() checks before
 * calling. Writes the chosen target into `target_out` and the chosen set list
 * into a malloc'd string returned through `only_out` (NULL meaning all).
 */
static int wizard(fetch_conf_t *conf, char *target_out, size_t target_sz,
                  char **only_out) {
  fprintf(stderr, "\n%s%skycg fetch%s %s- guided setup%s\n",
          kycg_ui_bold(), kycg_ui_cyan(), kycg_ui_reset(),
          kycg_ui_dim(), kycg_ui_reset());

  /* 1. Where. */
  char *dir = kycg_ui_ask("\nWhere should knowledgebases be stored?",
                          kycg_store_root(conf->store));
  if (!dir) return -1;
  conf->store = dir;   /* leaked deliberately: lives until process exit */

  /* 2. What. */
  const char *items[32];
  const char *notes[32];
  char notebuf[32][64];
  size_t n = 0;

  for (const kycg_seq_reg_t *r = KYCG_SEQ_REGISTRY; r->genome && n < 32; ++r) {
    size_t sets = 0;
    uint64_t bytes = 0;
    for (const kycg_zfile_t *f = r->files; f->name; ++f) {
      size_t l = strlen(f->name);
      if (l > 3 && strcmp(f->name + l - 3, ".cm") == 0) ++sets;
      bytes += f->size;
    }
    char hb[24];
    snprintf(notebuf[n], sizeof(notebuf[n]), "sequencing  %zu sets, %s",
             sets, kycg_ui_human(bytes, hb, sizeof(hb)));
    items[n] = r->genome;
    notes[n] = notebuf[n];
    ++n;
  }
  for (const kycg_array_reg_t *r = KYCG_ARRAY_REGISTRY; r->platform && n < 32; ++r) {
    snprintf(notebuf[n], sizeof(notebuf[n]), "array");
    items[n] = r->platform;
    notes[n] = notebuf[n];
    ++n;
  }

  long pick = kycg_ui_choose("Which collection?", items, notes, n);
  if (pick < 0) return -1;
  snprintf(target_out, target_sz, "%s", items[pick]);

  /* 3. Which sets. Requires a plan, which for arrays means the manifest. */
  plan_t probe = {0};
  const kycg_seq_reg_t *sr = find_seq(target_out);
  const kycg_array_reg_t *ar = find_array(target_out);

  fetch_conf_t probe_conf = *conf;
  probe_conf.only = NULL;

  int rc = sr ? build_plan_seq(sr, &probe_conf, &probe)
              : build_plan_array(ar, &probe_conf, &probe);
  if (rc != 0) { plan_free(&probe); return -1; }

  char **snames = NULL, **snotes = NULL;
  size_t n_sets = collect_sets(&probe, &snames, &snotes);
  plan_free(&probe);

  if (!n_sets) return -1;

  int *flags = kycg_ui_multiselect("Which sets?",
                                   (const char **)snames,
                                   (const char **)snotes, n_sets, 1);
  if (!flags) {
    for (size_t i = 0; i < n_sets; ++i) { free(snames[i]); free(snotes[i]); }
    free(snames); free(snotes);
    return -1;
  }

  /* All selected -> no filter at all, which keeps the plan honest about
   * files (like .idx sidecars) whose set name never appears in the menu. */
  size_t chosen = 0;
  for (size_t i = 0; i < n_sets; ++i) chosen += (size_t)(flags[i] != 0);

  if (chosen == n_sets) {
    *only_out = NULL;
  } else {
    size_t cap = 1;
    for (size_t i = 0; i < n_sets; ++i)
      if (flags[i]) cap += strlen(snames[i]) + 1;
    char *only = malloc(cap);
    only[0] = '\0';
    for (size_t i = 0; i < n_sets; ++i) {
      if (!flags[i]) continue;
      if (only[0]) strcat(only, ",");
      strcat(only, snames[i]);
    }
    *only_out = only;
  }

  free(flags);
  for (size_t i = 0; i < n_sets; ++i) { free(snames[i]); free(snotes[i]); }
  free(snames); free(snotes);
  return 0;
}

#endif /* KYCG_HAVE_CURL */

/* ------------------------------------------------------------------ usage */

static int usage(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg fetch [options] [<target>[:<sets>] ...]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Download and verify knowledgebases into a local store.\n");
  fprintf(stderr, "With no target on a terminal, asks where, what, and which sets.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Examples:\n");
  fprintf(stderr, "    kycg fetch                     guided setup\n");
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
  char wiz_target[128];
  char *wiz_only = NULL;
  const char *argv_targets[64];
  int n_targets = 0;

  if (optind < argc) {
    for (int j = optind; j < argc && n_targets < 64; ++j)
      argv_targets[n_targets++] = argv[j];
  } else if (kycg_ui_interactive()) {
    if (wizard(&conf, wiz_target, sizeof(wiz_target), &wiz_only) != 0) {
      fprintf(stderr, "\nNothing to do.\n");
      curl_global_cleanup();
      return 1;
    }
    if (wiz_only) conf.only = wiz_only;
    argv_targets[n_targets++] = wiz_target;
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

  free(wiz_only);
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
        char dir[4096];
        snprintf(dir, sizeof(dir), "%s/%s", root, sr->genome);
        printf("# %s -- Zenodo %s (doi %s)\n", sr->genome, sr->record, sr->doi);
        printf("set\tfile\tsize\tcached\n");
        for (const kycg_zfile_t *f = sr->files; f->name; ++f) {
          char setn[256], path[4400], hb[24];
          set_name_of(f->name, setn, sizeof(setn));
          snprintf(path, sizeof(path), "%s/%s", dir, f->name);
          printf("%s\t%s\t%s\t%s\n", setn, f->name,
                 kycg_ui_human(f->size, hb, sizeof(hb)),
                 kycg_store_is_file(path) ? "yes" : "no");
        }
      } else if (ar) {
        char dir[4096];
        snprintf(dir, sizeof(dir), "%s/%s/KYCG", root, ar->platform);
        char sums[4400];
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
          char setn[256], path[4400];
          set_name_of(nm, setn, sizeof(setn));
          snprintf(path, sizeof(path), "%s/%s", dir, nm);
          printf("%s\t%s\t%s\n", setn, nm,
                 kycg_store_is_file(path) ? "yes" : "no");
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
    printf("%s\tsequencing\tzenodo:%s\t%" PRIu64 "/%" PRIu64 "\n",
           r->genome, r->record, have, avail);
  }

  return 0;
}
