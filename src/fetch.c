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
 *   kycg fetch hg38:CGI,ChromHMM    a named subset, no questions asked
 *   kycg fetch                      browse, check what you want, fetch it
 *
 *   The colon form exists so a pipeline can name exactly what it wants on one
 *   line. The browser exists because nobody memorizes 32 set names -- and
 *   since browsing the catalogue and choosing from it are the same activity,
 *   one command does both. There was a separate `kycg list` for the browsing
 *   half; it was removed once that stopped being a separate act, so there is
 *   no second entry point to drift out of step with this one.
 *   See browse_catalogue() and fetch_picked() below.
 *
 * PLAN, CONFIRM, EXECUTE
 *   Nothing downloads until a plan has been built and shown: which files, how
 *   large, where they land, and what is already present. This matters because
 *   the hg38 collection is 363 MB and the difference between wanting all of it
 *   and wanting two sets is easy to express and easy to get wrong.
 *
 *   Building a plan costs one small request for that collection's SHA256SUMS,
 *   since the whole point of anchoring on that manifest is that
 *   kycg does not carry the file list and upstream can add a set without a
 *   rebuild. So the array path reaches the network before the confirmation, by
 *   a few kilobytes, and says so while it does.
 *
 * PROMPTING WITHOUT BREAKING AUTOMATION
 *   The original rule here was "never prompt", because a prompt hangs a
 *   Nextflow job or a Docker build forever with no indication of why. That
 *   guarantee is preserved exactly, by gating every question on an interactive
 *   terminal (see ui.c): off a TTY, an explicit target proceeds without asking
 *   and a missing target is an error rather than a wait. -f forces the same
 *   behavior on a TTY.
 *
 *   What has not changed: kycg still downloads in this command and nowhere
 *   else. No other subcommand touches the network.
 *
 * TRUST
 *   Both channels verify against a sha256 compiled into this binary by
 *   tools/make_registry.sh; see src/registry.h. Downloads land on a ".part"
 *   sibling and are
 *   renamed only after their digest matches, so an interrupted or corrupted
 *   fetch can never leave a file in the store that later reads as valid.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <fcntl.h>
#include <inttypes.h>
#include <sys/stat.h>
#include <dirent.h>
#include <time.h>

#include "kycg.h"
#include "args.h"
#include "digest.h"
#include "registry.h"
#include "store.h"
#include "ui.h"
#include "kbinfo.h"

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

/**
 * Is this a file a user would select -- a knowledgebase set, or the reference
 * row list that gives the whole-genome sets their coordinates?
 *
 * Excludes .idx sidecars, which are not independently useful: they are fetched
 * alongside the .cm they index. Includes .cr, which is not a set but is the
 * thing every sequencing analysis is positioned against, and which was
 * previously fetchable only by typing its name.
 */
static int is_selectable(const char *name) {
  size_t l = strlen(name);
  if (l > 3 && strcmp(name + l - 3, ".cm") == 0) return 1;
  if (l > 3 && strcmp(name + l - 3, ".cr") == 0) return 1;
  return 0;
}

/**
 * Does this file pass the subset filter?
 *
 * Returns 0 for no, KYCG_MATCH_SET when a token matched the set name, and
 * KYCG_MATCH_EXACT when it named the file outright. The distinction matters
 * because a set name can cover several published versions -- mm10 carries
 * ChromHMM.20220318 and ChromHMM.20220414 -- and those should collapse to one,
 * while a file named in full must be taken at its word. NULL/empty passes
 * everything, as a set-name match.
 */
#define KYCG_MATCH_SET   1
#define KYCG_MATCH_EXACT 2

static int passes_filter(const char *fname, const char *only) {
  if (!only || !*only) return KYCG_MATCH_SET;

  char setn[256];
  set_name_of(fname, setn, sizeof(setn));

  const char *p = only;
  while (*p) {
    const char *comma = strchr(p, ',');
    size_t len = comma ? (size_t)(comma - p) : strlen(p);
    if (len) {
      if (strlen(fname) == len && strncasecmp(fname, p, len) == 0)
        return KYCG_MATCH_EXACT;
      if (strlen(setn) == len && strncasecmp(setn, p, len) == 0)
        return KYCG_MATCH_SET;
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


/* ------------------------------------------------------- one collection */

/**
 * A fetchable collection, whichever channel it comes from.
 *
 * The two channels differ only in where their files sit and what they are
 * called. Both publish a SHA256SUMS at a pinned tag, both verify against an
 * anchor compiled in here, and both name their sets the same way -- so past
 * this struct there is one code path. That was not true while the whole-genome
 * sets came from Zenodo, which had no manifest and forced a parallel
 * implementation pinned file-by-file against md5.
 */
typedef struct coll_s {
  const char         *target;      /* "hg38" or "MSA", as the user types it */
  char                base[1024];  /* URL prefix, no trailing slash         */
  char                dir[4096];   /* where it lands in the store           */
  char                source[256]; /* provenance, for display               */
  const char         *anchor;      /* sha256 of SHA256SUMS; NULL = unpinned */
  const kycg_fsize_t *sizes;       /* display only, may be NULL/incomplete  */
  uint64_t            n_sets;      /* sets published at the pinned tag       */
  int                 unpinned_tag;/* -t named a tag this build cannot verify */

  /* The companion: the file that gives a row index its identity. For a genome
   * that is cpg_nocontig.cr, the coordinate list every sequencing analysis is
   * positioned against; for an array it is the probe ordering, without which a
   * set is a column of anonymous bits. Neither is a knowledgebase and neither
   * is optional, so both ride along with any fetch rather than being offered
   * as a choice. */
  char                comp_name[256];  /* "" when there is none              */
  char                comp_base[1024]; /* may differ from base               */
  const char         *comp_anchor;     /* NULL: listed in the main manifest  */
} coll_t;

struct coll_s;
#define KYCG_MF_LOCAL 0
#define KYCG_MF_FETCH 1

static char *coll_manifest(const struct coll_s *c, size_t *len_out, int mode);
static int coll_companion_sha(const struct coll_s *c, char sha_out[65], int mode);

/**
 * Resolve a target name to a collection. Returns 0 on success.
 *
 * `tag` overrides the array channel's pinned tag (the -t flag); NULL uses the
 * one compiled in. It does not apply to whole genomes, whose tags are per
 * repository and come from the registry.
 */
static int coll_for(const char *target, const char *store, const char *tag,
                    coll_t *c) {
  memset(c, 0, sizeof(*c));
  if (!tag || !*tag) tag = KYCG_IA_TAG;

  const kycg_seq_reg_t *sr = find_seq(target);
  if (sr) {
    c->target = sr->genome;
    snprintf(c->base, sizeof(c->base), "%s/%s/raw/%s",
             KYCG_KB_BASE_URL, sr->repo, sr->tag);
    snprintf(c->dir, sizeof(c->dir), "%s/%s", store, sr->genome);
    snprintf(c->source, sizeof(c->source), "%s %s", sr->repo, sr->tag);
    c->anchor = sr->sums_sha256;
    c->sizes = sr->sizes;
    c->n_sets = sr->n_sets;
    /* Listed in the same manifest as the sets, so no second anchor. */
    snprintf(c->comp_name, sizeof(c->comp_name), "cpg_nocontig.cr");
    snprintf(c->comp_base, sizeof(c->comp_base), "%s", c->base);
    c->comp_anchor = NULL;
    return 0;
  }

  const kycg_array_reg_t *ar = find_array(target);
  if (ar) {
    c->target = ar->platform;
    snprintf(c->base, sizeof(c->base), "%s/%s/%s/KYCG",
             KYCG_IA_BASE_URL, tag, ar->platform);
    snprintf(c->dir, sizeof(c->dir), "%s/%s/KYCG", store, ar->platform);
    snprintf(c->source, sizeof(c->source), "InfiniumAnnotation %s", tag);
    /* The anchor is only valid for the tag it was generated against. Asking
     * for a different one must fail the manifest check rather than quietly
     * fetch something this build cannot verify. */
    c->unpinned_tag = (strcmp(tag, KYCG_IA_TAG) != 0);
    c->anchor = ar->sums_sha256;
    c->sizes = NULL;      /* this channel publishes no sizes */
    c->n_sets = ar->n_sets;
    /* One level up, under the platform's own manifest. */
    if (ar->plat_sums_sha256) {
      snprintf(c->comp_name, sizeof(c->comp_name), "%s.ordering.tsv.gz",
               ar->platform);
      snprintf(c->comp_base, sizeof(c->comp_base), "%s/%s/%s",
               KYCG_IA_BASE_URL, tag, ar->platform);
      c->comp_anchor = ar->plat_sums_sha256;
    }
    return 0;
  }

  return -1;
}

/** Published size of a file, or 0 when unknown. Display only. */
static uint64_t coll_size_of(const coll_t *c, const char *name) {
  if (!c->sizes) return 0;
  for (const kycg_fsize_t *f = c->sizes; f->name; ++f)
    if (strcmp(f->name, name) == 0) return f->size;
  return 0;
}

/* --------------------------------------------------- spec -> store paths */

void kycg_free_specs(char **v, size_t n) {
  if (!v) return;
  for (size_t i = 0; i < n; ++i) free(v[i]);
  free(v);
}

size_t kycg_resolve_spec(const char *spec, const char *store, char ***out) {
  *out = NULL;
  if (!spec || !*spec) return 0;

  /* An existing path wins outright. Nothing else can then be ambiguous, and a
   * file whose name happens to contain a colon still works. */
  if (kycg_store_is_file(spec)) {
    char **v = malloc(sizeof(char *));
    if (!v) return 0;
    v[0] = strdup(spec);
    *out = v;
    return 1;
  }

  char target[128];
  const char *only = NULL;
  const char *colon = strchr(spec, ':');
  size_t len = colon ? (size_t)(colon - spec) : strlen(spec);
  if (len >= sizeof(target)) len = sizeof(target) - 1;
  memcpy(target, spec, len);
  target[len] = '\0';
  if (colon && colon[1]) only = colon + 1;

  coll_t c;
  if (coll_for(target, kycg_store_root(store), NULL, &c) != 0) return 0;

  /* Read the directory rather than the manifest: what is testable is what is
   * actually here, and a manifest lists what upstream publishes. */
  DIR *d = opendir(c.dir);
  if (!d) return 0;

  char **v = NULL;
  unsigned char *how = NULL;
  size_t n = 0, m = 0;
  struct dirent *e;
  while ((e = readdir(d))) {
    size_t l = strlen(e->d_name);
    if (e->d_name[0] == '.') continue;
    if (l <= 3 || strcmp(e->d_name + l - 3, ".cm") != 0) continue;
    int kind = passes_filter(e->d_name, only);
    if (!kind) continue;

    if (n == m) {
      size_t want = m ? m * 2 : 16;
      /* Adopt each block as it succeeds. Freeing them on partial failure
       * released the block realloc had just returned while v/how still held
       * the invalidated originals -- the caller then got a freed pointer to
       * walk and free again. */
      char **nv = realloc(v, want * sizeof(char *));
      if (!nv) break;
      v = nv;
      unsigned char *nh = realloc(how, want);
      if (!nh) break;
      how = nh;
      m = want;
    }
    char path[4600];
    snprintf(path, sizeof(path), "%s/%s", c.dir, e->d_name);
    v[n] = strdup(path);
    if (!v[n]) break;
    how[n] = (unsigned char)kind;
    ++n;
  }
  closedir(d);

  /* Stable order, so a run is reproducible and its output diffable. Sorting
   * by path also puts versions of a set next to each other, newest last,
   * because the date is a fixed-width field in the name. */
  if (n > 1) {
    for (size_t i = 0; i + 1 < n; ++i)          /* keep `how` with its path */
      for (size_t j = i + 1; j < n; ++j)
        if (strcmp(v[i], v[j]) > 0) {
          char *tp = v[i]; v[i] = v[j]; v[j] = tp;
          unsigned char th = how[i]; how[i] = how[j]; how[j] = th;
        }
  }

  /* One file per set name: a set that has been republished would otherwise be
   * tested twice over, which reads as two independent knowledgebases and is
   * corrected as two families. The newest wins; a file named in full is never
   * dropped, so a version can always be pinned by writing it out. */
  size_t k = 0;
  for (size_t i = 0; i < n; ++i) {
    if (how[i] == KYCG_MATCH_SET && k > 0) {
      char a[256], b[256];
      const char *pa = strrchr(v[k-1], '/'), *pb = strrchr(v[i], '/');
      set_name_of(pa ? pa + 1 : v[k-1], a, sizeof(a));
      set_name_of(pb ? pb + 1 : v[i],   b, sizeof(b));
      if (strcmp(a, b) == 0 && how[k-1] == KYCG_MATCH_SET) {
        fprintf(stderr,
                "%skycg: %s has more than one version; using %s%s\n",
                kycg_ui_dim(), b, pb ? pb + 1 : v[i], kycg_ui_reset());
        free(v[k-1]);
        v[k-1] = v[i];        /* the later one sorts newer */
        how[k-1] = how[i];
        continue;
      }
    }
    v[k] = v[i];
    how[k] = how[i];
    ++k;
  }
  n = k;

  free(how);
  *out = v;
  return n;
}

/* --------------------------------------------------------------- the plan */

typedef struct {
  char     name[512];
  char     url[4096];
  char     sha[65];       /* always set: both channels publish sha256 */
  uint64_t size;          /* 0 = not published by this channel        */
  int      have;          /* already present and digest-verified      */
} plan_item_t;

typedef struct {
  plan_item_t *a;
  size_t       n, m;
  size_t       n_sets;     /* items that are knowledgebases, not the companion */
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
  curl_easy_setopt(h, CURLOPT_FAILONERROR, 1L);
  curl_easy_setopt(h, CURLOPT_NOSIGNAL, 1L);
  curl_easy_setopt(h, CURLOPT_CONNECTTIMEOUT, 30L);
  curl_easy_setopt(h, CURLOPT_USERAGENT, "kycg/" KYCG_VERSION);

  /* Redirects have to be followed: github.com/<repo>/raw/... legitimately
   * lands on raw.githubusercontent.com. But they are confined to https and
   * bounded in number.
   *
   * Every byte is digest-checked against a compiled-in anchor, so a downgrade
   * could not substitute content -- what it could do is move the transfer to
   * cleartext, exposing which knowledgebases are being fetched, and an
   * unbounded chain is a way to waste a client's time. Every runtime URL is
   * built from compiled-in https constants (KYCG_IA_BASE_URL,
   * KYCG_KB_BASE_URL), so nothing legitimate ever needs another scheme. */
  curl_easy_setopt(h, CURLOPT_FOLLOWLOCATION, 1L);
  curl_easy_setopt(h, CURLOPT_MAXREDIRS, 10L);
#if defined(CURLOPT_REDIR_PROTOCOLS_STR) && LIBCURL_VERSION_NUM >= 0x075500
  curl_easy_setopt(h, CURLOPT_REDIR_PROTOCOLS_STR, "https");
#elif defined(CURLOPT_REDIR_PROTOCOLS)
  /* Pre-7.85 spelling; the bitmask form is deprecated but is what older
   * libcurl understands, and bioconda still builds against those. */
  curl_easy_setopt(h, CURLOPT_REDIR_PROTOCOLS, (long)CURLPROTO_HTTPS);
#endif
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

/*
 * Download to `path`, which is a temporary name the caller renames on success.
 *
 * O_EXCL, and the caller makes the name unique per process: the temp name used
 * to be a fixed "<file>.part", so two concurrent fetches into one store wrote
 * the same file. The likely outcome was a digest mismatch on both sides, but
 * there was a window where one process verified the digest, the other
 * truncated and rewrote the file, and the first renamed the other's partial
 * content into place as verified. Both are fetching identical bytes from the
 * same origin so the practical risk was nil, but a temp file that cannot
 * collide costs nothing.
 */
static int http_get_file(const char *url, const char *path, kycg_prog_t *pr) {
  int fd = open(path, O_WRONLY | O_CREAT | O_EXCL, 0644);
  if (fd < 0) return -1;
  FILE *fp = fdopen(fd, "wb");
  if (!fp) { close(fd); unlink(path); return -1; }

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
      if (v[cnt].name[0] && kycg_store_safe_name(v[cnt].name)) ++cnt;
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
  int direct;      /* -f: fetch now, no browser, no questions */
  int redownload;  /* -r: re-fetch even what is present and verified */
} fetch_conf_t;

/**
 * Build the plan for a collection: one manifest request, then one item per
 * file that survives the subset filter.
 *
 * This is the whole of what used to be two functions. Both channels publish a
 * manifest of sha256 digests at a pinned tag, so "which files, and what should
 * each one hash to" has one answer and one implementation.
 */
static int build_plan(const coll_t *c, const fetch_conf_t *conf, plan_t *plan) {
  if (c->unpinned_tag) {
    fprintf(stderr,
            "kycg fetch: this build pins InfiniumAnnotation %s, so it holds no\n"
            "digest for the tag you asked for and cannot verify anything fetched\n"
            "from it. Regenerate src/registry.h with tools/make_registry.sh and\n"
            "rebuild to move tags.\n", KYCG_IA_TAG);
    return -1;
  }
  if (!c->anchor) {
    fprintf(stderr, "kycg fetch: '%s' has no published manifest in this build.\n",
            c->target);
    return -1;
  }

  if (kycg_ui_fancy() && !kycg_ui_panel_active())
    fprintf(stderr, "%s  reading %s manifest...%s\r",
            kycg_ui_dim(), c->target, kycg_ui_reset());

  size_t len = 0;
  char *sums = coll_manifest(c, &len, KYCG_MF_FETCH);

  if (kycg_ui_fancy() && !kycg_ui_panel_active()) fputs("\r\033[2K", stderr);

  if (!sums) {
    fprintf(stderr,
            "%s%s%s cannot read the manifest for '%s'.\n"
            "Either the network is unavailable, or its SHA256SUMS does not match\n"
            "the digest pinned in this build. Refusing to fetch. If upstream has\n"
            "moved, regenerate src/registry.h with tools/make_registry.sh.\n",
            kycg_ui_red(), kycg_ui_cross(), kycg_ui_reset(), c->target);
    return -1;
  }

  size_t n_ent = 0;
  sums_ent_t *ent = parse_sums(sums, &n_ent);
  if (!ent) { free(sums); return -1; }

  snprintf(plan->dir, sizeof(plan->dir), "%s", c->dir);
  snprintf(plan->target, sizeof(plan->target), "%s", c->target);
  snprintf(plan->source, sizeof(plan->source), "%s", c->source);
  plan->sums_text = sums;
  plan->sums_len = len;
  plan->sizes_known = (c->sizes != NULL);

  for (size_t i = 0; i < n_ent; ++i) {
    /* The companion ignores the subset filter: asking for one set still means
     * asking for the thing that makes that set interpretable. */
    int is_comp = (c->comp_name[0] && strcmp(ent[i].name, c->comp_name) == 0);
    if (!is_comp && !passes_filter(ent[i].name, conf->only)) continue;
    if (!is_comp) ++plan->n_sets;
    plan_item_t *it = plan_add(plan);
    if (!it) break;
    snprintf(it->name, sizeof(it->name), "%s", ent[i].name);
    snprintf(it->url, sizeof(it->url), "%s/%s", c->base, ent[i].name);
    snprintf(it->sha, sizeof(it->sha), "%s", ent[i].sha);
    it->size = coll_size_of(c, ent[i].name);
  }

  /* A companion under a manifest of its own is not in the list above. */
  if (c->comp_name[0] && c->comp_anchor) {
    char csha[65];
    if (coll_companion_sha(c, csha, KYCG_MF_FETCH)) {
      plan_item_t *it = plan_add(plan);
      if (it) {
        snprintf(it->name, sizeof(it->name), "%s", c->comp_name);
        snprintf(it->url, sizeof(it->url), "%s/%s", c->comp_base, c->comp_name);
        snprintf(it->sha, sizeof(it->sha), "%s", csha);
      }
    }
  }

  free(ent);
  return 0;
}

/* -------------------------------------------------------------- execution */

typedef struct {
  uint64_t n_got, n_skip, n_fail;
  uint64_t bytes_got;
} tally_t;

#ifdef KYCG_HAVE_CURL

/*
 * Remove temp files a previous run left behind.
 *
 * A failed download unlinks its own temp file, so these only appear when a
 * process was killed outright. The name carries a pid to keep concurrent
 * fetches apart, which means they would otherwise accumulate one per killed
 * run instead of being overwritten. Only files older than a day are swept, so
 * a fetch running right now in another process is never touched.
 */
static void sweep_stale_parts(const char *dir) {
  DIR *d = opendir(dir);
  if (!d) return;

  time_t now = time(NULL);
  struct dirent *e;
  while ((e = readdir(d))) {
    size_t l = strlen(e->d_name);
    if (l < 6 || strcmp(e->d_name + l - 5, ".part") != 0) continue;

    char p[4800];
    snprintf(p, sizeof(p), "%s/%s", dir, e->d_name);
    struct stat st;
    if (stat(p, &st) == 0 && S_ISREG(st.st_mode) &&
        now - st.st_mtime > 24 * 60 * 60)
      unlink(p);
  }
  closedir(d);
}

static int execute_plan(const plan_t *plan, tally_t *t) {
  if (kycg_store_mkdir_p(plan->dir) != 0) {
    fprintf(stderr, "kycg fetch: cannot create %s\n", plan->dir);
    return -1;
  }
  sweep_stale_parts(plan->dir);

  for (size_t i = 0; i < plan->n; ++i) {
    const plan_item_t *it = &plan->a[i];

    if (it->have) {
      ++t->n_skip;
      continue;
    }

    char path[4600], part[4750];
    snprintf(path, sizeof(path), "%s/%s", plan->dir, it->name);
    /* Per-process temp name; see http_get_file. Still ends in ".part", which
     * is what keeps it invisible to the store's ".cm" enumeration. */
    snprintf(part, sizeof(part), "%s.%ld.part", path, (long)getpid());

    kycg_prog_t pr;
    kycg_prog_begin(&pr, it->name, it->size);

    if (http_get_file(it->url, part, &pr) != 0) {
      kycg_prog_done(&pr, "download failed", 0);
      ++t->n_fail;
      continue;
    }

    /* No digest, no rename. Both channels always publish one, so an item
     * without a sha is a bug rather than a permissive case. */
    int ok = 0;
    char got[65];
    if (it->sha[0]) {
      if (kycg_sha256_file(part, got) == 0) ok = kycg_digest_equal(got, it->sha);
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
    /* build_plan always sets sums_text, so this is belt and braces: rebuild a
     * manifest from what actually landed rather than leave the store without
     * one. */
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
  FILE *o = stderr;
  fprintf(o, "\n");
  fprintf(o, "  %skycg fetch%s  %s— browse, choose and download knowledgebases%s\n\n",
          KYCG_H_TITLE, KYCG_H_OFF, KYCG_H_NOTE, KYCG_H_OFF);

  fprintf(o, "%sUsage%s\n", KYCG_H_TITLE, KYCG_H_OFF);
  fprintf(o, "    kycg fetch [options] [%s<target>%s[:%s<sets>%s] ...]\n\n",
          KYCG_H_KEY, KYCG_H_OFF, KYCG_H_KEY, KYCG_H_OFF);

  fprintf(o, "    %sOn a terminal this opens the catalogue -- a tree of every%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %scollection. Naming a target opens it already checked, so you see%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %swhat will be downloaded, can narrow it, and press f to start.%s\n\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %sWith stdout redirected it prints the catalogue as TSV; with a%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %starget it downloads. -f forces that on a terminal too.%s\n\n",
          KYCG_H_NOTE, KYCG_H_OFF);

  fprintf(o, "%sExamples%s\n", KYCG_H_TITLE, KYCG_H_OFF);
  fprintf(o, "    kycg fetch %-22s browse everything\n", "");
  fprintf(o, "    kycg fetch %s%-22s%s open with hg38 checked, f to start\n",
          KYCG_H_KEY, "hg38", KYCG_H_OFF);
  fprintf(o, "    kycg fetch %s%-22s%s open with just those checked\n",
          KYCG_H_KEY, "hg38:CGI,ChromHMM", KYCG_H_OFF);
  fprintf(o, "    kycg fetch %s%-22s%s download it, no questions\n",
          KYCG_H_KEY, "-f hg38", KYCG_H_OFF);
  fprintf(o, "    kycg fetch %s%-22s%s TSV catalogue, for scripts\n\n",
          KYCG_H_KEY, "| cut -f1", KYCG_H_OFF);

  fprintf(o, "%sIn the tree%s\n", KYCG_H_TITLE, KYCG_H_OFF);
  {
    /* Two fixed columns; padding the keys and the descriptions separately is
     * what keeps the right column from drifting with the text on the left. */
    struct { const char *k, *d; } key[] = {
      {"arrows, j k", "move"},          {"space", "check a set"},
      {"right, l",    "unfold"},        {"f",     "fetch what is checked"},
      {"left, h",     "fold"},          {"d",     "change the store"},
      {"home, end",   "jump"},          {"q",     "quit"},
    };
    for (size_t i = 0; i < sizeof(key)/sizeof(key[0]); i += 2)
      fprintf(o, "    %s%-12s%s %-10s  %s%-6s%s %s\n",
              KYCG_H_KEY, key[i].k,   KYCG_H_OFF, key[i].d,
              KYCG_H_KEY, key[i+1].k, KYCG_H_OFF, key[i+1].d);
  }
  fprintf(o, "\n");

  fprintf(o, "%sTargets%s\n", KYCG_H_TITLE, KYCG_H_OFF);
  fprintf(o, "    %-9s ", "genomes");
  for (const kycg_seq_reg_t *r = KYCG_SEQ_REGISTRY; r->genome; ++r)
    fprintf(o, "%s%s%s ", KYCG_H_KEY, r->genome, KYCG_H_OFF);
  fprintf(o, "\n    %-9s ", "arrays");
  for (const kycg_array_reg_t *r = KYCG_ARRAY_REGISTRY; r->platform; ++r)
    fprintf(o, "%s%s%s ", KYCG_H_KEY, r->platform, KYCG_H_OFF);
  fprintf(o, "\n\n");

  fprintf(o, "%sOptions%s\n", KYCG_H_TITLE, KYCG_H_OFF);
  struct { const char *f, *d; } opt[] = {
    {"-d DIR", "store directory [$KYCG_DATA_DIR, else ~/.cache/kycg]"},
    {"-o SETS", "subset by set name; same as the :SETS suffix"},
    {"-f", "download now: no browser, no questions"},
    {"-r", "re-download even what is present and verified"},
    {"-t TAG", "InfiniumAnnotation tag, arrays only"},
    {"-h", "this help"},
  };
  for (size_t i = 0; i < sizeof(opt)/sizeof(opt[0]); ++i)
    fprintf(o, "    %s%-8s%s %s\n", KYCG_H_KEY, opt[i].f, KYCG_H_OFF, opt[i].d);
  fprintf(o, "\n");

  fprintf(o, "    %sThe reference row list -- a genome's .cr, a platform's probe%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %sordering -- always comes too; it is what gives a row its identity.%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %sNothing is downloaded outside this command, and nothing is asked%s\n",
          KYCG_H_NOTE, KYCG_H_OFF);
  fprintf(o, "    %swhen stdin is not a terminal.%s\n\n", KYCG_H_NOTE, KYCG_H_OFF);
  return 1;
}

/* The browser lives further down this file; `fetch` with no target, and
 * `fetch <target>` on a terminal, both hand off to it. */
static int browse_catalogue(int argc, char *argv[]);

int main_fetch(int argc, char *argv[]) {
  fetch_conf_t conf = {0};
  conf.tag = KYCG_IA_TAG;

  int c;
  /* Options may follow the target; BSD getopt would stop at it. */
  kycg_permute_args(argc, argv, "d:o:t:frh");
  while ((c = getopt(argc, argv, "d:o:t:frh")) >= 0) {
    switch (c) {
    case 'd': conf.store = optarg; break;
    case 'o': conf.only = optarg; break;
    case 't': conf.tag = optarg; break;
    case 'f': conf.direct = 1; break;
    case 'r': conf.redownload = 1; break;
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
          "any path. See `kycg fetch` for the expected layout.\n");
  return 1;
#else
  curl_global_init(CURL_GLOBAL_DEFAULT);

  /* Collect targets, either from argv or from the guided run. */
  const char *argv_targets[64];
  int n_targets = 0;

  /* No target, or -l: this is a request to see the catalogue rather than to
   * download one. On a terminal that is the browser -- where checking a set
   * and fetching it happen anyway -- and off one it is plain TSV, so
   * `kycg fetch | cut -f1` still works. */
  /* No target: the catalogue. A tree on a terminal, TSV when redirected, so
   * `kycg fetch | cut -f1` still works. */
  if (optind >= argc) {
    char *lav[8];
    int lac = 0;
    lav[lac++] = "fetch";
    if (conf.store) { lav[lac++] = "-d"; lav[lac++] = (char *)conf.store; }
    for (int j = optind; j < argc && lac < 7; ++j) lav[lac++] = argv[j];
    lav[lac] = NULL;
    optind = 1;
    curl_global_cleanup();
    return browse_catalogue(lac, lav);
  }

  /* A named target on a terminal opens the catalogue with that target chosen,
   * so the download is seen before it starts and can be narrowed in the same
   * screen. -f skips straight to fetching, and so does a non-terminal, which
   * is what keeps scripts working. */
  if (!conf.direct && kycg_ui_interactive()) {
    char *lav[8];
    int lac = 0;
    lav[lac++] = "fetch";
    if (conf.store) { lav[lac++] = "-d"; lav[lac++] = (char *)conf.store; }
    for (int j = optind; j < argc && lac < 7; ++j) lav[lac++] = argv[j];
    lav[lac] = NULL;
    optind = 1;
    curl_global_cleanup();
    return browse_catalogue(lac, lav);
  }

  for (int j = optind; j < argc && n_targets < 64; ++j)
    argv_targets[n_targets++] = argv[j];

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

    coll_t coll;
    if (coll_for(target, kycg_store_root(tc.store), tc.tag, &coll) != 0) {
      fprintf(stderr,
              "kycg fetch: '%s' is not a known platform or genome.\n"
              "Run `kycg fetch` to see what is available.\n", target);
      rc = 1;
      continue;
    }

    plan_t plan = {0};
    int prc = build_plan(&coll, &tc, &plan);
    if (prc != 0) { plan_free(&plan); rc = 1; continue; }

    if (!plan.n) {
      fprintf(stderr, "kycg fetch: nothing in '%s' matches that selection.\n",
              target);
      plan_free(&plan);
      rc = 1;
      continue;
    }

    /* A named subset that matched no set is a failed request. The plan is not
     * empty in that case -- the companion reference is always in it -- so
     * without this a typo like "mm10:CGl" fetched only the .cr and reported
     * "1 fetched", leaving the user believing CGI was downloaded. */
    if (only && *only && !plan.n_sets) {
      fprintf(stderr,
              "kycg fetch: no set in '%s' matches '%s'.\n"
              "Run `kycg fetch %s` to see what it publishes.\n",
              target, only, target);
      plan_free(&plan);
      rc = 1;
      continue;
    }

    plan_check_present(&plan, tc.redownload);
    plan_show(&plan);

    size_t n_todo = 0;
    for (size_t i = 0; i < plan.n; ++i) if (!plan.a[i].have) ++n_todo;
    if (!n_todo) { t.n_skip += plan.n; plan_free(&plan); continue; }

    /* No confirmation here: reaching this point means either -f or no
     * terminal, and both say "do it". Anyone who wants to look first gets the
     * catalogue, which is the same screen with the same numbers. */
    if (execute_plan(&plan, &t) != 0) rc = 1;
    plan_free(&plan);
  }

  if (t.n_got || t.n_skip || t.n_fail) {
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

/* ------------------------------------------------ the catalogue browser */

static int browse_usage(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg fetch [options] [target ...]\n");
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
 * `kycg fetch` output is buffered rather than printed directly, so it can be
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
      kycg_ui_browse(comment ? comment : "kycg fetch", header,
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
 * the browser is also the fetch picker, so its expand and accept callbacks
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
  const char *only;  /* subset named on the command line, or NULL */
} listctx_t;

static void picks_add(picks_t *p, const char *target, const char *file) {
  if (p->n == p->m) {
    size_t want = p->m ? p->m * 2 : 32;
    char **t = realloc(p->target, want * sizeof(char *));
    if (!t) return;
    p->target = t;
    char **f = realloc(p->file, want * sizeof(char *));
    if (!f) return;
    p->file = f;
    p->m = want;
  }
  /* Both strdups must land before n advances, or picks_free walks a NULL. */
  char *t = strdup(target), *f = strdup(file);
  if (!t || !f) { free(t); free(f); return; }
  p->target[p->n] = t;
  p->file[p->n] = f;
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
 * A platform's SHA256SUMS: the catalogue of what that platform publishes.
 *
 * Looked for in three places, cheapest first -- the local store, this run's
 * memo, then the network. `mode` is KYCG_MF_LOCAL to stop before the network
 * (the overview redraws too often to reach it) or KYCG_MF_FETCH to go all the
 * way.
 *
 * The store's copy is accepted only if it hashes to the pinned anchor. That
 * single check makes staleness self-healing: a manifest written by an older
 * kycg -- which recorded only the files it fetched, not the whole catalogue --
 * fails it and is replaced, with no refresh button to remember to press. It is
 * the same test a freshly downloaded manifest has to pass, applied to the copy
 * already on disk.
 *
 * Returns a malloc'd copy the caller frees, or NULL if unavailable. A fetched
 * manifest is verified against the compiled anchor exactly as a download would
 * be: a catalogue deserves no more trust than the files it lists.
 */
static struct { const char *plat; char *text; size_t len; } g_manifest[16];
static size_t g_manifest_n = 0;

static char *coll_manifest(const struct coll_s *c, size_t *len_out, int mode) {
  if (len_out) *len_out = 0;

  /* 1. The store's own copy, written when the collection was fetched. It is
   *    authoritative and always current, so it is never memoized. */
  char sums[4400];
  snprintf(sums, sizeof(sums), "%s/%s", c->dir, KYCG_IA_SUMS_FILE);
  FILE *fp = fopen(sums, "rb");
  if (fp) {
    fseek(fp, 0, SEEK_END);
    long sz = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char *text = NULL;
    size_t len = 0;
    if (sz > 0 && (text = malloc((size_t)sz + 1))) {
      len = fread(text, 1, (size_t)sz, fp);
      text[len] = '\0';
    }
    fclose(fp);
    if (text) {
      char got[65];
      kycg_sha256_buf(text, len, got);
      if (c->anchor && kycg_digest_equal(got, c->anchor)) {
        if (len_out) *len_out = len;
        return text;
      }
      /* Does not match what this build pins -- written by an older kycg, or
       * against a different tag. Not trustworthy as a catalogue. */
      free(text);
    }
  }

  /* 2. Something this run already pulled. */
  for (size_t i = 0; i < g_manifest_n; ++i) {
    if (g_manifest[i].plat != c->anchor) continue;
    if (!g_manifest[i].text) return NULL;
    char *dup = malloc(g_manifest[i].len + 1);
    if (!dup) return NULL;
    memcpy(dup, g_manifest[i].text, g_manifest[i].len + 1);
    if (len_out) *len_out = g_manifest[i].len;
    return dup;
  }

  if (mode == KYCG_MF_LOCAL) return NULL;

#ifndef KYCG_HAVE_CURL
  return NULL;
#else
  char *text = NULL;
  size_t len = 0;

  if (c->anchor) {
    char url[5200];
    snprintf(url, sizeof(url), "%s/%s", c->base, KYCG_IA_SUMS_FILE);
    text = http_get_mem(url, &len);
    if (text) {
      char got[65];
      kycg_sha256_buf(text, len, got);
      if (!kycg_digest_equal(got, c->anchor)) { free(text); text = NULL; }
    }
  }

  /* Remember the outcome either way: a failed lookup memoized as NULL stops
   * a dead platform being retried on every redraw. */
  if (g_manifest_n < 16) {
    g_manifest[g_manifest_n].plat = c->anchor;
    g_manifest[g_manifest_n].len = len;
    g_manifest[g_manifest_n].text = NULL;
    if (text) {
      g_manifest[g_manifest_n].text = malloc(len + 1);
      if (g_manifest[g_manifest_n].text)
        memcpy(g_manifest[g_manifest_n].text, text, len + 1);
    }
    ++g_manifest_n;
  }

  /* Replace a rejected local copy, so the store stops carrying a manifest that
   * disagrees with the tag this build pins. */
  if (text && kycg_store_is_file(sums)) {
    FILE *out = fopen(sums, "wb");
    if (out) { fwrite(text, 1, len, out); fclose(out); }
  }

  if (len_out) *len_out = len;
  return text;
#endif
}

/**
 * Digest of the companion file, or 0 if it cannot be established.
 *
 * A genome's companion is listed in the same manifest as its sets; an array's
 * lives under the platform manifest, so that one is fetched and verified in
 * its own right. Either way the file is trusted exactly as a set is.
 */
static int coll_companion_sha(const struct coll_s *c, char sha_out[65], int mode) {
  if (!c->comp_name[0]) return 0;

  coll_t src = *c;
  if (c->comp_anchor) {          /* a manifest of its own */
    snprintf(src.base, sizeof(src.base), "%s", c->comp_base);
    src.anchor = c->comp_anchor;
    /* dir has to move with base. It did not, so coll_manifest's write-back of
     * a re-verified copy put the *platform* manifest at the *KYCG* directory's
     * SHA256SUMS -- self-healing, since the next run rejects it against the
     * anchor and re-downloads, but in the meantime `shasum -c SHA256SUMS` in
     * that directory failed for no visible reason and count_cached lost its
     * catalogue. The companion lives one level up from KYCG/. */
    size_t dl = strlen(src.dir);
    if (dl > 5 && strcmp(src.dir + dl - 5, "/KYCG") == 0)
      src.dir[dl - 5] = '\0';
  }

  size_t len = 0;
  char *text = coll_manifest(&src, &len, mode);
  if (!text) return 0;

  size_t n_ent = 0;
  sums_ent_t *ent = parse_sums(text, &n_ent);
  int found = 0;
  for (size_t i = 0; ent && i < n_ent; ++i) {
    if (strcmp(ent[i].name, c->comp_name) != 0) continue;
    snprintf(sha_out, 65, "%s", ent[i].sha);
    found = 1;
    break;
  }
  free(ent);
  free(text);
  return found;
}

/** How many sets a collection publishes, or 0 if the catalogue is not to hand. */
static uint64_t coll_set_total(const coll_t *c) {
  size_t len = 0;
  char *text = coll_manifest(c, &len, KYCG_MF_LOCAL);
  if (!text) return 0;

  size_t n_ent = 0;
  sums_ent_t *ent = parse_sums(text, &n_ent);
  uint64_t n = 0;
  for (size_t i = 0; ent && i < n_ent; ++i) {
    size_t l = strlen(ent[i].name);
    if (l > 3 && strcmp(ent[i].name + l - 3, ".cm") == 0) ++n;
  }
  free(ent);
  free(text);
  return n;
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
 * Children of one row of the `kycg fetch` overview: the sets that target holds.
 *
 * Widths are fixed rather than measured, because the tree renders children one
 * parent at a time — columns that re-measured on each expansion would shift
 * under rows already on screen.
 *
 * Whole-genome targets can be listed offline, since the file list is compiled
 * in. Array platforms cannot: anchoring on SHA256SUMS is what lets upstream
 * add a set without a kycg rebuild, so the set list only exists once fetched.
 */
static void refresh_overview(listctx_t *lc);

static void expand_target(void *ctx, const char *row, kycg_ui_kids_t *out) {
  listctx_t *lc = ctx;

  char target[128];
  const char *tab = strchr(row, '\t');
  size_t len = tab ? (size_t)(tab - row) : strlen(row);
  if (len >= sizeof(target)) len = sizeof(target) - 1;
  memcpy(target, row, len);
  target[len] = '\0';

  coll_t c;
  if (coll_for(target, lc->root, NULL, &c) != 0) return;

  /* The catalogue lives in the manifest. A collection with nothing fetched has
   * no local copy, which used to make it look as though it published nothing
   * at all -- so pull it (a couple of kilobytes, verified against the compiled
   * anchor) to list what exists. Memoized for the run, never written to the
   * store: writing it would claim files are present that are not. */
  int had = (coll_set_total(&c) != 0);
  size_t mlen = 0;
  char *text = coll_manifest(&c, &mlen, KYCG_MF_FETCH);

  /* The catalogue is what the overview needs to show a denominator, so the
   * row above updates the moment it becomes knowable. */
  if (text && !had) refresh_overview(lc);

  if (!text) {
    kid_push(out, KYCG_ROW_MISSING, NULL,
             "catalogue unavailable - run: kycg fetch %s", c.target);
    return;
  }

  /* The companion goes first and in red: it is not a choice, it comes with
   * whatever you pick, and burying it alphabetically among things that are
   * choices would misrepresent it. */
  if (c.comp_name[0]) {
    char cpath[4400], hb[24];
    snprintf(cpath, sizeof(cpath), "%s/%s", c.dir, c.comp_name);
    int chave = kycg_store_is_file(cpath);
    uint64_t csz = coll_size_of(&c, c.comp_name);
    char setn[256];
    set_name_of(c.comp_name, setn, sizeof(setn));
    kid_push(out, chave ? KYCG_ROW_HAVE : KYCG_ROW_REQUIRED, NULL,
             "%-22.22s %-32.32s %9s  %s", setn, c.comp_name,
             csz ? kycg_ui_human(csz, hb, sizeof(hb)) : "",
             chave ? "cached" : "always fetched");
  }

  size_t n_ent = 0;
  sums_ent_t *ent = parse_sums(text, &n_ent);
  for (size_t i = 0; ent && i < n_ent; ++i) {
    const char *nm = ent[i].name;
    if (!is_selectable(nm)) continue;
    /* Already shown above, out of alphabetical order and on purpose. */
    if (c.comp_name[0] && strcmp(nm, c.comp_name) == 0) continue;

    char setn[256], path[4400], hb[24];
    set_name_of(nm, setn, sizeof(setn));
    snprintf(path, sizeof(path), "%s/%s", c.dir, nm);
    int have = kycg_store_is_file(path);

    uint64_t sz = coll_size_of(&c, nm);
    kid_push(out, have ? KYCG_ROW_HAVE : KYCG_ROW_MISSING, nm,
             "%-22.22s %-32.32s %9s  %s", setn, nm,
             sz ? kycg_ui_human(sz, hb, sizeof(hb)) : "",
             have ? "cached" : "-");
  }
  free(ent);
  free(text);
}

static uint64_t count_cached(const coll_t *c);

/**
 * Fill `rows` with one line per target, counting what the store holds.
 *
 * Called for the initial view and again after anything changes it -- a fetch,
 * or a different store directory. Rebuilding is cheap: it stats files, and
 * the target list is fixed and short.
 */
static void build_overview(const char *root, rows_t *rows) {
  for (const kycg_seq_reg_t *r = KYCG_SEQ_REGISTRY; r->genome; ++r) {
    coll_t c;
    if (coll_for(r->genome, root, NULL, &c) != 0) continue;

    uint64_t have = count_cached(&c);
    /* Pinned, not discovered: the tag is immutable, so its set count is a
     * fixed fact about this build. The overview therefore shows have/total
     * immediately, without a manifest and without touching the network. */
    uint64_t nt = c.n_sets;

    char rb[32], cnt[64];
    if (nt) snprintf(cnt, sizeof(cnt), "%" PRIu64 "/%" PRIu64, have, nt);
    else    snprintf(cnt, sizeof(cnt), "%" PRIu64, have);

    /* Green means "something here is usable", dim means "nothing yet". How
     * much is in the count. */
    rows_push(rows, have ? KYCG_ROW_HAVE : KYCG_ROW_MISSING,
              "%s\twhole genome\t%s\t%s\t%s",
              r->genome, commify(r->rows, rb, sizeof(rb)), c.source, cnt);
  }

  for (const kycg_array_reg_t *r = KYCG_ARRAY_REGISTRY; r->platform; ++r) {
    coll_t c;
    if (coll_for(r->platform, root, NULL, &c) != 0) continue;

    uint64_t nc = count_cached(&c);
    uint64_t nt = c.n_sets;

    char rb[32], cnt[64];
    if (nt) snprintf(cnt, sizeof(cnt), "%" PRIu64 "/%" PRIu64, nc, nt);
    else    snprintf(cnt, sizeof(cnt), "%" PRIu64, nc);
    rows_push(rows, nc ? KYCG_ROW_HAVE : KYCG_ROW_MISSING,
              "%s\tarray\t%s\t%s\t%s",
              r->platform, commify(r->rows, rb, sizeof(rb)), c.source, cnt);
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
 * Fetch what the picker checked: one plan per target, confirmed once.
 *
 * The selection arrives as (target, file name) pairs and is turned back into
 * the same comma-separated subset string the -o flag takes, so this path and
 * `kycg fetch hg38:CGI` build their plans through identical code. Index
 * sidecars are added alongside their .cm, since a set's index is part of the
 * set and picking one row should not leave half of it behind.
 *
 * Renders into the browser's panel rather than the normal screen, so the
 * catalogue stays visible above the download.
 */
static int fetch_picked(const picks_t *picks, const char *store) {
#ifndef KYCG_HAVE_CURL
  (void)picks; (void)store;
  fprintf(stderr, "kycg: this build has no network support.\n");
  return 1;
#else
  curl_global_init(CURL_GLOBAL_DEFAULT);
  kycg_ui_panel_open(4);

  tally_t t = {0};
  int rc = 0;
  char **done = calloc(picks->n, sizeof(char *));
  if (!done) return 1;
  size_t n_done = 0;

  for (size_t i = 0; i < picks->n; ++i) {
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

    coll_t coll;
    plan_t plan = {0};
    int prc = (coll_for(picks->target[i], kycg_store_root(store), NULL, &coll) == 0)
              ? build_plan(&coll, &conf, &plan) : -1;
    free(only);
    if (prc != 0 || !plan.n) { plan_free(&plan); rc = 1; continue; }

    plan_check_present(&plan, conf.redownload);

    size_t n_todo = 0;
    uint64_t bytes = 0;
    for (size_t k = 0; k < plan.n; ++k)
      if (!plan.a[k].have) { ++n_todo; bytes += plan.a[k].size; }
    if (!n_todo) { t.n_skip += plan.n; plan_free(&plan); continue; }

    char hb[24];
    kycg_ui_panel_line(0, "  %s%s%s  %s  %zu file(s)%s%s  %s  %s%s%s",
                       kycg_ui_bold(), plan.target, kycg_ui_reset(),
                       kycg_ui_bullet(), n_todo,
                       bytes ? ", " : "",
                       bytes ? kycg_ui_human(bytes, hb, sizeof(hb)) : "",
                       kycg_ui_bullet(),
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
 * Which collection a target belongs to, for looking up a recommendation.
 * Genomes recommend by name; every array platform shares one list, since the
 * sets are the same annotation projected onto different probe orderings.
 */
static const char *reco_collection(const char *target) {
  if (find_seq(target)) return target;
  if (find_array(target)) return "array";
  return NULL;
}

/** Is this set part of its collection's recommended selection? */
int kycg_kb_recommended(void *ctx, const char *root, const char *key) {
  (void)ctx;
  char target[128], setn[256];
  const char *tab = strchr(root, '\t');
  size_t len = tab ? (size_t)(tab - root) : strlen(root);
  if (len >= sizeof(target)) len = sizeof(target) - 1;
  memcpy(target, root, len);
  target[len] = '\0';

  const char *coll = reco_collection(target);
  if (!coll) return 0;

  const char *file = strchr(key, ':');
  file = file ? file + 1 : key;
  set_name_of(file, setn, sizeof(setn));
  return kycg_kbinfo_recommended(setn, coll);
}

/* ------------------------------------------------- provenance panel layout */

/*
 * A rendered panel: fully formatted lines, ANSI already embedded.
 *
 * The panel API needs its height at open time, but the fields being shown are
 * sentences of no fixed length that have to be wrapped to the terminal. So the
 * whole thing is laid out into this buffer first and the height read off it.
 * Wrapping happens on word boundaries and the labelled fields hang-indent, so
 * a three-line citation still reads as one field.
 */
#define INFO_MAX_LINES 64

typedef struct {
  char *line[INFO_MAX_LINES];
  int   n;
  int   cols;      /* usable width; the pane is narrower than the terminal */
} info_lay_t;

static void lay_push(info_lay_t *L, const char *s) {
  if (L->n >= INFO_MAX_LINES) return;
  L->line[L->n] = strdup(s ? s : "");
  if (L->line[L->n]) ++L->n;
}

static void lay_free(info_lay_t *L) {
  for (int i = 0; i < L->n; ++i) free(L->line[i]);
  L->n = 0;
}

static void lay_head(info_lay_t *L, const char *setn, const char *title) {
  char buf[1024];
  snprintf(buf, sizeof(buf), "  %s%s%s  %s%s%s", kycg_ui_bold(), setn,
           kycg_ui_reset(), kycg_ui_cyan(), title ? title : "",
           kycg_ui_reset());
  lay_push(L, buf);
  lay_push(L, "");
}

/*
 * Wrap `text` into the panel, under an optional dim label.
 *
 * `label` NULL means running prose at the left margin; otherwise the label is
 * printed once in a fixed-width gutter and continuation lines align under the
 * text rather than under the label.
 */
static void lay_wrap(info_lay_t *L, const char *label, const char *text) {
  if (!text || !*text) return;

  const int gutter = label ? 14 : 2;   /* "  processing  " is the widest */
  int avail = (L->cols > 0 ? L->cols : kycg_ui_cols()) - gutter - 2;
  if (avail < 20) avail = 20;

  const char *p = text;
  int first = 1;
  while (*p) {
    while (*p == ' ') ++p;
    if (!*p) break;

    /* Longest prefix that fits, broken at the last space; a single word
     * longer than the line is emitted whole and allowed to be truncated,
     * which beats hyphenating a DOI. */
    size_t rest = strlen(p), take = rest;
    if (rest > (size_t)avail) {
      size_t brk = 0;
      for (size_t i = 0; i < (size_t)avail; ++i) if (p[i] == ' ') brk = i;
      take = brk ? brk : (size_t)avail;
    }

    char buf[1024], head[64];
    if (label && first)
      snprintf(head, sizeof(head), "  %s%-*s%s", kycg_ui_dim(), gutter - 4,
               label, kycg_ui_reset());
    else
      snprintf(head, sizeof(head), "%*s", gutter, "");

    snprintf(buf, sizeof(buf), "%s%s%.*s", head, label && first ? "  " : "",
             (int)take, p);
    lay_push(L, buf);

    p += take;
    first = 0;
  }
  lay_push(L, "");
}

/** `i` in the browser: what is this set, and where did it come from. */
/**
 * The `i` pane: what the set under the cursor is, and where it came from.
 *
 * A detail callback rather than a modal panel, so the arrow keys keep working
 * while it is open and the text follows the cursor. It runs on every redraw,
 * which is affordable because kycg_kbinfo_find is a scan of a compiled-in
 * table and the wrapping is a few hundred bytes of formatting.
 */
void kycg_kb_detail(void *ctx, const char *root, const char *child_key,
                    int cols, kycg_ui_detail_t *out) {
  (void)ctx;
  (void)root;                 /* the set name is enough to describe a set */
  out->rows = NULL;
  out->n = 0;

  info_lay_t L = {0};
  L.cols = cols > 20 ? cols : 20;

  if (!child_key) {
    /* A collection row: the columns above already describe it, so the pane
     * collapses rather than showing a placeholder. With the pane on by
     * default that placeholder would have been the first thing anyone saw,
     * occupying half the screen to say nothing. */
    return;
  } else {
    char setn[256];
    const char *file = strchr(child_key, ':');
    file = file ? file + 1 : child_key;
    set_name_of(file, setn, sizeof(setn));

    const kycg_kbinfo_t *k = kycg_kbinfo_find(setn);
    if (!k) {
      char buf[512];
      snprintf(buf, sizeof(buf), "  %s%s%s  %s(nothing recorded about this set)%s",
               kycg_ui_bold(), setn, kycg_ui_reset(),
               kycg_ui_dim(), kycg_ui_reset());
      lay_push(&L, buf);
    } else {
      lay_head(&L, setn, k->title);
      lay_wrap(&L, NULL, k->biology);
      lay_wrap(&L, "source", k->source);
      lay_wrap(&L, "citation", k->citation);
      lay_wrap(&L, "processing", k->processing);
    }
  }

  out->rows = malloc((size_t)(L.n ? L.n : 1) * sizeof(char *));
  if (!out->rows) { lay_free(&L); return; }
  for (int i = 0; i < L.n; ++i) out->rows[i] = L.line[i];
  out->n = (size_t)L.n;
  /* Ownership moves to the caller; do not free the strings here. */
  L.n = 0;
  lay_free(&L);
}

/* --------------------------------------------- choosing sets from the store */

/*
 * The store offered as a tree, shared by `kycg test` and `kycg annotate`.
 *
 * Both need the same thing -- let the user see what a collection publishes,
 * fetch what is missing, and hand back what was chosen -- and they differ only
 * in which collections are worth offering. So the caller supplies the target
 * list and the verb, and everything below is common. Keeping one copy is the
 * point: two pickers would drift, and the one that drifted would be the one
 * used less.
 */

/**
 * State for offering the store as a tree, the same widget `kycg fetch` uses.
 *
 * Roots are the collections the caller offered; children are the sets actually
 * cached under them, which by construction share that collection's row space.
 * So the list cannot contain a choice the caller would then have to refuse.
 */
typedef struct {
  const char *root;        /* store root */

  char        **rows;      /* root row text, tab-separated */
  unsigned char *styles;
  char        **names;     /* the target name for each root */
  size_t        n, m;

  char        **chosen;    /* paths picked, filled on accept */
  size_t        n_chosen, m_chosen;
} pickctx_t;

static void pk_add_target(pickctx_t *p, const char *name, const char *kind,
                            uint64_t rows) {
  char **paths = NULL;
  size_t n_cached = kycg_resolve_spec(name, NULL, &paths);
  kycg_free_specs(paths, n_cached);

  if (p->n == p->m) {
    /* `p->x = realloc(p->x, ...)` loses the original on failure, and returning
     * without updating p->m left the arrays grown but the capacity stale --
     * the next call re-entered this branch, realloc(NULL) handed back a fresh
     * block, every prior name leaked, and pick_free() then freed garbage.
     * Adopt each block only once it has succeeded, and update m last. */
    size_t want = p->m ? p->m * 2 : 16;
    char **nr = realloc(p->rows, want * sizeof(char *));
    if (!nr) return;
    p->rows = nr;
    char **nn = realloc(p->names, want * sizeof(char *));
    if (!nn) return;
    p->names = nn;
    unsigned char *ns = realloc(p->styles, want);
    if (!ns) return;
    p->styles = ns;
    p->m = want;
  }

  char buf[512], rb[32];
  snprintf(buf, sizeof(buf), "%s\t%s\t%s\t%zu",
           name, kind, commify(rows, rb, sizeof(rb)), n_cached);
  /* Both strings must land before n advances: a NULL row reaches
   * kycg_ui_tree(), whose strchr() on it would crash. */
  char *row_s = strdup(buf), *name_s = strdup(name);
  if (!row_s || !name_s) { free(row_s); free(name_s); return; }
  p->rows[p->n] = row_s;
  p->names[p->n] = name_s;
  /* Dim a collection with nothing in it: it is listed so the user learns it
   * exists and could be fetched, not because it can be tested against. */
  p->styles[p->n] = n_cached ? KYCG_ROW_HAVE : KYCG_ROW_MISSING;
  ++p->n;
}

static void pk_expand(void *ctx, const char *row, kycg_ui_kids_t *out) {
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
    kid_push(out, KYCG_ROW_MISSING, NULL,
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
    kid_push(out, cat[i].cached ? KYCG_ROW_HAVE : KYCG_ROW_MISSING, key,
                   "%-22.22s %-32.32s %s", setn, cat[i].name,
                   cat[i].cached ? "cached" : "-");
  }
  kycg_catalogue_free(cat, n);
}

static void pk_accept(void *ctx, const char *root, const char *key) {
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
static void pk_commit_fetch(void *ctx) {
  pickctx_t *p = ctx;
  if (p->n_chosen) kycg_fetch_specs(p->chosen, p->n_chosen, NULL);
  for (size_t i = 0; i < p->n_chosen; ++i) free(p->chosen[i]);
  p->n_chosen = 0;
}

static void pk_free(pickctx_t *p) {
  for (size_t i = 0; i < p->n; ++i) { free(p->rows[i]); free(p->names[i]); }
  free(p->rows); free(p->names); free(p->styles);
  for (size_t i = 0; i < p->n_chosen; ++i) free(p->chosen[i]);
  free(p->chosen);
  memset(p, 0, sizeof(*p));
}
/**
 * Offer `targets` in the browser and return what the user chose.
 *
 * `verb_key`/`verb` label the action that ends the session (t for test, a for
 * annotate); `f` fetches without leaving, so a set that is missing can be
 * downloaded and then used in one sitting.
 *
 * Returns the number of chosen specs and fills *out with malloc'd
 * "target:file" strings, freeable with kycg_free_specs. Returns 0 when the
 * user quit without choosing, and (size_t)-1 when the terminal cannot host
 * the browser at all.
 */
size_t kycg_pick_sets(const kycg_pick_target_t *targets, size_t n_targets,
                      const char *title, char verb_key, const char *verb,
                      char ***out) {
  *out = NULL;
  if (!n_targets) return 0;

  pickctx_t pc = {0};
  pc.root = kycg_store_root(NULL);
  for (size_t i = 0; i < n_targets; ++i)
    pk_add_target(&pc, targets[i].name, targets[i].kind, targets[i].rows);

  kycg_ui_tree_t spec = {0};
  spec.title = title;
  spec.header = "target\tkind\trows\tcached_sets";
  spec.roots = pc.rows;
  spec.root_styles = pc.styles;
  spec.n_roots = pc.n;
  spec.expand = pk_expand;
  /* Two verbs on one screen: fetch what is missing, then act on it. f keeps
   * the browser open (it has a commit); the caller's verb ends it and the
   * selection is what gets used. */
  spec.actions[0].key = 'f';
  spec.actions[0].verb = "fetch";
  spec.actions[0].accept = pk_accept;
  spec.actions[0].commit = pk_commit_fetch;
  spec.actions[1].key = verb_key;
  spec.actions[1].verb = verb;
  spec.actions[1].accept = pk_accept;
  spec.actions[1].commit = NULL;
  spec.n_actions = 2;
  /* Cached sets are the ones worth acting on, so they stay checkable. */
  spec.have_selectable = 1;
  /* The same callbacks the fetch browser uses, so a set is recommended and
   * described identically whichever tree you reached it through. */
  spec.recommend = kycg_kb_recommended;
  spec.detail_key = 'i';
  spec.detail_verb = "info";
  spec.detail = kycg_kb_detail;
  spec.ctx = &pc;

  int rc = kycg_ui_tree(&spec);
  if (rc < 0) { pk_free(&pc); return (size_t)-1; }
  if (rc != 2 || !pc.n_chosen) { pk_free(&pc); return 0; }

  /* Hand the chosen specs over; pk_free must not release them. */
  *out = pc.chosen;
  size_t n = pc.n_chosen;
  pc.chosen = NULL;
  pc.n_chosen = 0;
  pk_free(&pc);
  return n;
}

/** Which rows a named target arrives with already checked. */
static int on_preselect(void *ctx, const char *root, const char *key) {
  (void)root;
  listctx_t *lc = ctx;
  return passes_filter(key, lc->only);
}

/**
 * `d` in the browser: point it at a different store.
 *
 * The store location is the one thing that cannot be changed from inside
 * otherwise -- every other question the browser answers is about its contents.
 */
static int on_list_key(void *ctx, char key, const char *root,
                       const char *child_key) {
  listctx_t *lc = ctx;
  (void)root; (void)child_key;      /* d is about the store, not a row */

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

/**
 * How many .cm sets are actually on disk in a directory.
 *
 * Counted by looking, not by consulting the manifest. What is present is a
 * property of the filesystem, and reading it from a manifest made the number
 * wrong whenever that manifest was missing, partial, or stale -- a store with
 * eight sets in it reported zero.
 */
void kycg_catalogue_free(kycg_catalogue_t *v, size_t n) {
  if (!v) return;
  for (size_t i = 0; i < n; ++i) free(v[i].name);
  free(v);
}

kycg_catalogue_t *kycg_catalogue(const char *target, const char *store,
                                 size_t *n) {
  *n = 0;
  coll_t c;
  if (coll_for(target, kycg_store_root(store), NULL, &c) != 0) return NULL;

  size_t len = 0;
  char *text = coll_manifest(&c, &len, KYCG_MF_FETCH);
  if (!text) return NULL;

  size_t n_ent = 0;
  sums_ent_t *ent = parse_sums(text, &n_ent);
  kycg_catalogue_t *v = calloc(n_ent ? n_ent : 1, sizeof(kycg_catalogue_t));
  if (!v) { free(ent); free(text); *n = 0; return NULL; }
  size_t k = 0;

  for (size_t i = 0; ent && v && i < n_ent; ++i) {
    size_t l = strlen(ent[i].name);
    if (l <= 3 || strcmp(ent[i].name + l - 3, ".cm") != 0) continue;
    char path[4600];
    snprintf(path, sizeof(path), "%s/%s", c.dir, ent[i].name);
    v[k].name = strdup(ent[i].name);
    v[k].cached = kycg_store_is_file(path);
    if (v[k].name) ++k;
  }
  free(ent);
  free(text);

  *n = k;
  return v;
}

int kycg_fetch_specs(char *const *specs, size_t n, const char *store) {
  picks_t p = {0};
  for (size_t i = 0; i < n; ++i) {
    /* Split "mm10:CGI.20220904.cm" the way the picker hands it over. */
    const char *colon = strchr(specs[i], ':');
    if (!colon || !colon[1]) continue;
    char target[128];
    size_t len = (size_t)(colon - specs[i]);
    if (len >= sizeof(target)) len = sizeof(target) - 1;
    memcpy(target, specs[i], len);
    target[len] = '\0';
    picks_add(&p, target, colon + 1);
  }
  int rc = p.n ? fetch_picked(&p, store) : 0;
  picks_free(&p);
  return rc;
}

/**
 * How many of this collection's sets are actually present.
 *
 * Counts the intersection of what is on disk with what the collection
 * publishes -- not every .cm in the directory. A store can legitimately hold a
 * file the pinned tag no longer publishes: fetch mm10 at v1, pin v2, and the
 * ChromHMM build v2 dropped is still sitting there. Counting it produced
 * "29/28", a numerator larger than its own denominator, which reads as a
 * broken counter rather than as a true fact about the store. It also let any
 * unrelated .cm dropped in the directory inflate the count.
 *
 * The denominator is the pinned set count, so the numerator has to be measured
 * against the same catalogue or the two are not comparable.
 *
 * With no local manifest there is nothing to intersect against and every .cm
 * counts. That is the hand-populated store -- unpacked from a tarball, mounted
 * off a shared drive -- which store.c supports deliberately, and reporting
 * zero for it would be a worse answer than counting a stray file.
 */
static uint64_t count_cached(const coll_t *c) {
  size_t len = 0;
  char *text = coll_manifest(c, &len, KYCG_MF_LOCAL);
  size_t n_ent = 0;
  sums_ent_t *ent = text ? parse_sums(text, &n_ent) : NULL;

  DIR *d = opendir(c->dir);
  if (!d) { free(ent); free(text); return 0; }

  uint64_t n = 0;
  struct dirent *e;
  while ((e = readdir(d))) {
    size_t l = strlen(e->d_name);
    if (e->d_name[0] == '.') continue;
    if (l <= 3 || strcmp(e->d_name + l - 3, ".cm") != 0) continue;

    if (!ent) { ++n; continue; }        /* no catalogue to check against */

    for (size_t i = 0; i < n_ent; ++i) {
      if (strcmp(ent[i].name, e->d_name) == 0) { ++n; break; }
    }
  }
  closedir(d);

  free(ent);
  free(text);
  return n;
}

static int browse_catalogue(int argc, char *argv[]) {
  const char *store = NULL;
  int c;
  /* Options may follow the target; BSD getopt would stop at it. */
  kycg_permute_args(argc, argv, "d:h");
  while ((c = getopt(argc, argv, "d:h")) >= 0) {
    switch (c) {
    case 'd': store = optarg; break;
    case 'h': return browse_usage();
    default: return browse_usage();
    }
  }

  const char *root = kycg_store_root(store);

  /* A target named on the command line opens the catalogue on it, already
   * chosen. Same syntax as a fetch, so `hg38:CGI` narrows what is checked. */
  char open_target[128] = {0};
  const char *open_only = NULL;
  if (optind < argc && isatty(STDOUT_FILENO)) {
    const char *spec = argv[optind];
    const char *colon = strchr(spec, ':');
    size_t len = colon ? (size_t)(colon - spec) : strlen(spec);
    if (len >= sizeof(open_target)) len = sizeof(open_target) - 1;
    memcpy(open_target, spec, len);
    open_target[len] = '\0';
    if (colon && colon[1]) open_only = colon + 1;
    if (!find_seq(open_target) && !find_array(open_target)) open_target[0] = '\0';
  }

  if (optind < argc && !open_target[0]) {
    for (int j = optind; j < argc; ++j) {
      /* Same target[:sets] syntax as a fetch, so `-l hg38:CGI` narrows the
       * listing rather than failing on a name it does not recognize. */
      char target[128];
      const char *spec = argv[j];
      const char *only = NULL;
      const char *colon = strchr(spec, ':');
      if (colon) {
        size_t len = (size_t)(colon - spec);
        if (len >= sizeof(target)) len = sizeof(target) - 1;
        memcpy(target, spec, len);
        target[len] = '\0';
        if (colon[1]) only = colon + 1;
      } else {
        snprintf(target, sizeof(target), "%s", spec);
      }

      coll_t c;
      if (coll_for(target, root, NULL, &c) == 0) {
        const kycg_seq_reg_t *sr = find_seq(target);
        const kycg_array_reg_t *ar = find_array(target);
        char title[512], rb[32];
        (void)ar;

        if (sr)
          snprintf(title, sizeof(title), "%s -- %s rows -- %s (archive doi %s)",
                   sr->genome, commify(sr->rows, rb, sizeof(rb)),
                   c.source, sr->doi);
        else
          snprintf(title, sizeof(title), "%s -- %s rows -- %s",
                   ar->platform, commify(ar->rows, rb, sizeof(rb)), c.source);

        /* Read the local manifest if fetched. If not, pull the catalogue only
         * when someone is watching -- a redirected stdout means a script is
         * reading, and a script must not trigger a download it did not ask
         * for. */
        size_t tlen = 0;
        char *text = coll_manifest(&c, &tlen,
                                   isatty(STDOUT_FILENO) ? KYCG_MF_FETCH
                                                         : KYCG_MF_LOCAL);
        if (!text) {
          printf("# %s\n# not fetched yet; run: kycg fetch %s\n",
                 title, c.target);
          continue;
        }

        rows_t rows = {0};
        size_t n_ent = 0;
        sums_ent_t *ent = parse_sums(text, &n_ent);
        for (size_t i = 0; ent && i < n_ent; ++i) {
          const char *nm = ent[i].name;
          if (!is_selectable(nm)) continue;
          if (!passes_filter(nm, only)) continue;
          char setn[256], path[4400], hb[24];
          set_name_of(nm, setn, sizeof(setn));
          snprintf(path, sizeof(path), "%s/%s", c.dir, nm);
          int have = kycg_store_is_file(path);
          uint64_t sz = coll_size_of(&c, nm);
          rows_push(&rows, have ? KYCG_ROW_HAVE : KYCG_ROW_MISSING,
                    "%s\t%s\t%s\t%s", setn, nm,
                    sz ? kycg_ui_human(sz, hb, sizeof(hb)) : "-",
                    have ? "yes" : "no");
        }
        free(ent);
        free(text);
        rows_emit(&rows, title, "set\tfile\tsize\tcached");
      } else {
        fprintf(stderr, "kycg fetch: '%s' is not a known platform or genome.\n",
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
    lc.only = open_only;

    kycg_ui_tree_t spec = {0};
    spec.title = title;
    spec.header = "target\tkind\trows\tsource\tcached_sets";
    spec.roots = rows.a;
    spec.root_styles = rows.st;
    spec.n_roots = rows.n;
    spec.expand = expand_target;
    spec.actions[0].key = 'f';
    spec.actions[0].verb = "fetch";
    spec.actions[0].accept = on_pick;
    spec.actions[0].commit = on_commit;
    spec.n_actions = 1;
    spec.on_key = on_list_key;
    spec.recommend = kycg_kb_recommended;
    spec.detail_key = 'i';
    spec.detail_verb = "info";
    spec.detail = kycg_kb_detail;
    spec.hint = "d store";
    spec.ctx = &lc;
    if (open_target[0]) {
      spec.open_root = open_target;
      spec.preselect = on_preselect;
    }

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
