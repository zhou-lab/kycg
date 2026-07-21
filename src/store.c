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
 * Local knowledgebase store: path resolution and enumeration.
 *
 * The store is a plain directory tree, and everything here treats it as one.
 * Enumeration walks it rather than consulting an index, so a store that was
 * populated by hand -- unpacked from a tarball, copied off a shared drive,
 * symlinked from a lab NFS mount -- works exactly like one built by
 * `kycg fetch`. That property is the reason there is no manifest: any manifest
 * would have to be maintained, and would then be capable of disagreeing with
 * what is actually on disk.
 */

#include "store.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>

int kycg_store_safe_name(const char *s) {
  if (!s || !*s) return 0;
  if (s[0] == '.') return 0;              /* ".", "..", and hidden files */
  for (const char *p = s; *p; ++p) {
    if (*p == '/' || *p == '\\') return 0;
    if ((unsigned char)*p < 0x20) return 0;   /* control chars, incl. newline */
  }
  return 1;
}

const char *kycg_store_root(const char *override) {
  static char buf[4096];
  if (override && *override) return override;

  const char *env = getenv("KYCG_DATA_DIR");
  if (env && *env) return env;

  const char *home = getenv("HOME");
  if (!home || !*home) home = ".";
  snprintf(buf, sizeof(buf), "%s/.cache/kycg", home);
  return buf;
}

int kycg_store_mkdir_p(const char *path) {
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

int kycg_store_is_file(const char *path) {
  struct stat st;
  return stat(path, &st) == 0 && S_ISREG(st.st_mode);
}

static int ends_with(const char *s, const char *suf) {
  size_t ls = strlen(s), lf = strlen(suf);
  return ls >= lf && strcmp(s + ls - lf, suf) == 0;
}

/* Push a malloc'd copy of `path` onto a growable vector. */
static int push(char ***v, size_t *n, size_t *m, const char *path) {
  if (*n == *m) {
    size_t want = *m ? *m * 2 : 64;
    char **p = realloc(*v, want * sizeof(char *));
    if (!p) return -1;
    *v = p; *m = want;
  }
  (*v)[*n] = strdup(path);
  if (!(*v)[*n]) return -1;
  ++*n;
  return 0;
}

/* Recursive walk, bounded in depth because the store is only ever two levels
 * (<target>/ or <target>/KYCG/) and a symlink loop should not be fatal. */
static void walk(const char *dir, int depth, char ***v, size_t *n, size_t *m) {
  if (depth > 4) return;

  DIR *d = opendir(dir);
  if (!d) return;

  struct dirent *e;
  while ((e = readdir(d))) {
    if (e->d_name[0] == '.') continue;

    char path[4096];
    snprintf(path, sizeof(path), "%s/%s", dir, e->d_name);

    /* lstat, not stat: a symlinked directory inside the store would otherwise
     * be descended into, and a link pointing at an ancestor yields the same
     * .cm more than once. A duplicate is not cosmetic here -- testing one
     * knowledgebase twice in a stratum inflates BH's m and shifts every FDR
     * in it. Symlinked *files* are still followed, which is what makes a
     * store assembled from a shared mount work. */
    struct stat st;
    if (lstat(path, &st) != 0) continue;

    if (S_ISDIR(st.st_mode)) {
      walk(path, depth + 1, v, n, m);
    } else if (ends_with(e->d_name, ".cm")) {
      /* .cm.idx is excluded by construction: it does not end in ".cm".
       * S_ISREG is deliberately not required: lstat above means a symlinked
       * .cm reports as a link, and those are legitimate in a store assembled
       * from a shared mount. */
      if (push(v, n, m, path) != 0) return;
    }
  }
  closedir(d);
}

static int cmp_str(const void *a, const void *b) {
  return strcmp(*(const char *const *)a, *(const char *const *)b);
}

char **kycg_store_find_cm(const char *root, size_t *n) {
  char **v = NULL;
  size_t cnt = 0, cap = 0;

  walk(root, 0, &v, &cnt, &cap);
  if (cnt > 1) qsort(v, cnt, sizeof(char *), cmp_str);

  *n = cnt;
  return v;
}

void kycg_store_free_list(char **v, size_t n) {
  if (!v) return;
  for (size_t i = 0; i < n; ++i) free(v[i]);
  free(v);
}

const char *kycg_store_relative(const char *root, const char *path) {
  size_t lr = strlen(root);
  if (strncmp(path, root, lr) == 0) {
    const char *p = path + lr;
    while (*p == '/') ++p;
    if (*p) return p;
  }
  return path;
}
