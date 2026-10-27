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
 * The probe ordering of an array platform: probe ID -> row index.
 *
 * A .cm is positional and carries no probe IDs, so every command that starts
 * from probe IDs (annotate, and test -p) goes through this file. Shared here
 * so the two read it, verify it and search it the same way.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "kycg.h"
#include "ordering.h"

void kycg_ordering_free(kycg_ordering_t *o) {
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
 *
 * The file is asked for by its store path, through the resolver every tool in
 * the suite shares, so the path, the state and the sentence that repairs a
 * missing file are libyame's answers rather than four kycg guesses. The policy
 * on top of them is the suite's: a file this build does not pin is a warning
 * and the run continues (an ordering rarely changes, and stopping an analysis
 * over a digest helps nobody), a file that is not there at all is fatal.
 */
int kycg_ordering_load(const char *store, const char *platform,
                       const char *spec_hint, const char *cmd,
                       kycg_ordering_t *out) {
  char spec[4096], path[4096], advice[1024] = {0};
  snprintf(spec, sizeof(spec), "%s/%s.ordering.tsv.gz", platform, platform);

  yame_store_state_t st = yame_store_resolve(kycg_fetch_cfg(), spec, store,
                                             path, sizeof(path), NULL,
                                             advice, sizeof(advice));
  /* The repair is kycg's own fetch of the sets named: it brings the ordering
   * THIS build pins. libyame's advice names a `yame fetch`, which brings the
   * one the installed yame pins, and that need not be the same release. */
  char fix[1200];
  if (spec_hint && strchr(spec_hint, ':'))
    snprintf(fix, sizeof(fix), "kycg fetch -f %s", spec_hint);
  else
    snprintf(fix, sizeof(fix), "kycg fetch -f %s:<set>", platform);
  (void)advice;

  if (st == YAME_STORE_STALE)
    fprintf(stderr,
            "kycg %s: the %s ordering on disk is from a different "
            "release than this build pins;\n  using it anyway. To replace "
            "it: %s -r\n", cmd, platform, fix);

  gzFile fp = (st == YAME_STORE_CURRENT || st == YAME_STORE_STALE)
            ? gzopen(path, "rb") : NULL;
  if (!fp) {
    /* One-release fallback: a store fetched before YAME v1.33 keyed the
     * ordering under InfiniumAnnotation/, which the registry no longer lists,
     * so the resolver cannot see it. */
    char legacy[4096];
    snprintf(legacy, sizeof(legacy),
             "%s/InfiniumAnnotation/%s/%s.ordering.tsv.gz",
             store, platform, platform);
    fp = gzopen(legacy, "rb");
  }
  if (!fp) {
    fprintf(stderr,
            "kycg %s: cannot read the probe ordering for %s (%s).\n"
            "  It comes with any set fetched for the platform:\n      %s\n",
            cmd, platform, spec, fix);
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
            "kycg %s: %s.ordering.tsv.gz is not sorted (at line %d).\n"
            "This build assumes the LC_ALL=C ordering InfiniumAnnotation\n"
            "publishes; refusing to guess row indices from it.\n",
            cmd, platform, unsorted_at + 2);
    for (size_t i = 0; i < n; ++i) free(v[i]);
    free(v);
    return -1;
  }

  out->id = v;
  out->n = n;
  return 0;
}

/** Row index of `id`, or -1. Binary search: the ordering is sorted. */
long kycg_ordering_find(const kycg_ordering_t *o, const char *id) {
  size_t lo = 0, hi = o->n;
  while (lo < hi) {
    size_t mid = lo + (hi - lo) / 2;
    int c = strcmp(o->id[mid], id);
    if (c == 0) return (long)mid;
    if (c < 0) lo = mid + 1; else hi = mid;
  }
  return -1;
}
