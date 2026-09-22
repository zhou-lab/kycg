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
 * Tests for how the store is enumerated: the walk behind every "offer me what
 * I have" list, and the root resolution behind every path.
 *
 * The walk's rules are not arbitrary. A duplicate .cm is not a cosmetic bug:
 * testing one knowledgebase twice inside a stratum inflates BH's m and shifts
 * every FDR in that stratum. So a symlinked *directory* is not descended into
 * (a link to an ancestor would yield the same file twice) while a symlinked
 * *file* is kept (that is what makes a store assembled from a shared mount
 * work). Both halves are checked here by building a store that contains each.
 */

#include "store.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

static int failures = 0;

static void fail(const char *fmt, const char *a, const char *b) {
  printf("  FAIL ");
  printf(fmt, a, b);
  printf("\n");
  ++failures;
}

static void touch(const char *path) {
  FILE *f = fopen(path, "wb");
  if (f) { fputs("x", f); fclose(f); }
}

/* Is `needle` (a basename) among the n found paths? */
static int found(char **v, size_t n, const char *needle) {
  for (size_t i = 0; i < n; ++i) {
    const char *b = strrchr(v[i], '/');
    if (b && strcmp(b + 1, needle) == 0) return 1;
  }
  return 0;
}

int main(void) {
  char root[] = "/tmp/kycg_storewalk_XXXXXX";
  if (!mkdtemp(root)) { printf("  FAIL cannot make a temp store\n"); return 1; }

  char p[4096];

  /* ---- mkdir -p makes a whole chain, and is idempotent ------------------ */
  snprintf(p, sizeof(p), "%s/mm10/KYCG", root);
  if (kycg_store_mkdir_p(p) != 0) { printf("  FAIL mkdir_p on a new chain\n"); ++failures; }
  if (kycg_store_mkdir_p(p) != 0) { printf("  FAIL mkdir_p on an existing dir\n"); ++failures; }

  snprintf(p, sizeof(p), "%s/EPIC/KYCG", root);
  kycg_store_mkdir_p(p);

  /* ---- a store with everything the walk has to judge -------------------- */
  snprintf(p, sizeof(p), "%s/mm10/KYCG/CGI.20220904.cm", root);      touch(p);
  snprintf(p, sizeof(p), "%s/mm10/KYCG/CGI.20220904.cm.idx", root);  touch(p);
  snprintf(p, sizeof(p), "%s/mm10/cpg_nocontig.cr", root);           touch(p);
  snprintf(p, sizeof(p), "%s/EPIC/KYCG/ChromHMM.20220303.cm", root); touch(p);
  snprintf(p, sizeof(p), "%s/SHA256SUMS", root);                     touch(p);
  /* A dotfile: invisible to the walk by design. */
  snprintf(p, sizeof(p), "%s/mm10/.hidden.cm", root);                touch(p);

  /* A symlinked FILE is legitimate -- a store assembled from a shared mount
   * is mostly these -- and must be found. */
  char target[4096], link[4096];
  snprintf(target, sizeof(target), "%s/EPIC/KYCG/ChromHMM.20220303.cm", root);
  snprintf(link, sizeof(link), "%s/EPIC/KYCG/Linked.20220303.cm", root);
  int have_symlinks = (symlink(target, link) == 0);

  /* A symlinked DIRECTORY pointing at an ancestor: descending it would report
   * every .cm above a second time. */
  char loop[4096];
  snprintf(loop, sizeof(loop), "%s/EPIC/loop", root);
  if (have_symlinks) symlink(root, loop);

  size_t n = 0;
  char **v = kycg_store_find_cm(root, &n);

  if (!found(v, n, "CGI.20220904.cm"))      fail("%s%s", "a .cm in a KYCG dir was missed", "");
  if (!found(v, n, "ChromHMM.20220303.cm")) fail("%s%s", "a .cm under another platform was missed", "");
  if (found(v, n, "CGI.20220904.cm.idx"))   fail("%s%s", "a .cm.idx sidecar was listed as a set", "");
  if (found(v, n, "cpg_nocontig.cr"))       fail("%s%s", "a .cr reference was listed as a set", "");
  if (found(v, n, "SHA256SUMS"))            fail("%s%s", "the manifest was listed as a set", "");
  if (found(v, n, ".hidden.cm"))            fail("%s%s", "a dotfile was listed", "");
  if (have_symlinks && !found(v, n, "Linked.20220303.cm"))
    fail("%s%s", "a symlinked .cm was skipped", "");

  /* No path twice -- the FDR-inflating case. */
  for (size_t i = 0; i + 1 < n; ++i)
    if (strcmp(v[i], v[i + 1]) == 0)
      fail("the walk returned %s twice%s", v[i], "");

  /* Sorted, because the output of `kycg test` over a whole store has to be
   * diffable between runs. */
  for (size_t i = 0; i + 1 < n; ++i)
    if (strcmp(v[i], v[i + 1]) > 0)
      fail("the walk is unsorted: %s before %s", v[i], v[i + 1]);

  /* ---- display paths ---------------------------------------------------- */
  snprintf(p, sizeof(p), "%s/mm10/KYCG/CGI.20220904.cm", root);
  if (strcmp(kycg_store_relative(root, p), "mm10/KYCG/CGI.20220904.cm") != 0)
    fail("relative() gave %s%s", kycg_store_relative(root, p), "");
  /* A path outside the root is returned whole rather than mangled. */
  if (strcmp(kycg_store_relative(root, "/elsewhere/x.cm"), "/elsewhere/x.cm") != 0)
    fail("relative() mangled a path outside the root%s%s", "", "");
  /* The root itself has nothing after it; returning "" would print a blank. */
  if (strcmp(kycg_store_relative(root, root), root) != 0)
    fail("relative() of the root itself was not the root%s%s", "", "");

  kycg_store_free_list(v, n);
  kycg_store_free_list(NULL, 0);          /* must not crash */

  /* ---- is_file distinguishes the three kinds of answer ------------------ */
  snprintf(p, sizeof(p), "%s/SHA256SUMS", root);
  if (!kycg_store_is_file(p))        { printf("  FAIL is_file said no to a file\n"); ++failures; }
  snprintf(p, sizeof(p), "%s/mm10", root);
  if (kycg_store_is_file(p))         { printf("  FAIL is_file said yes to a directory\n"); ++failures; }
  if (kycg_store_is_file("/nonexistent/kycg/x")) {
    printf("  FAIL is_file said yes to a missing path\n"); ++failures; }

  /* ---- an empty store is empty, not an error ---------------------------- */
  {
    char empty[] = "/tmp/kycg_storewalk_e_XXXXXX";
    if (mkdtemp(empty)) {
      size_t en = 1;
      char **ev = kycg_store_find_cm(empty, &en);
      if (en != 0) { printf("  FAIL an empty store reported %zu sets\n", en); ++failures; }
      kycg_store_free_list(ev, en);
      rmdir(empty);
    }
  }
  /* A root that does not exist is the same answer, not a crash. */
  {
    size_t en = 1;
    char **ev = kycg_store_find_cm("/nonexistent/kycg/store", &en);
    if (en != 0) { printf("  FAIL a missing store reported %zu sets\n", en); ++failures; }
    kycg_store_free_list(ev, en);
  }

  /* ---- root resolution: -d beats the env, which beats the default ------- */
  setenv("YAME_DATA_HOME", "/tmp/kycg_env_store", 1);
  if (strcmp(kycg_store_root(NULL), "/tmp/kycg_env_store") != 0)
    fail("root() ignored $YAME_DATA_HOME: %s%s", kycg_store_root(NULL), "");
  if (strcmp(kycg_store_root("/tmp/kycg_flag_store"), "/tmp/kycg_flag_store") != 0)
    fail("root() ignored the -d override: %s%s", kycg_store_root("/tmp/kycg_flag_store"), "");
  unsetenv("YAME_DATA_HOME");
  {
    const char *dflt = kycg_store_root(NULL);
    if (!dflt || !*dflt || dflt[0] != '/')
      fail("root() gave no absolute default: %s%s", dflt ? dflt : "(null)", "");
  }

  if (failures) { printf("  %d check(s) failed\n", failures); return 1; }
  printf("  all checks passed\n");
  return 0;
}
