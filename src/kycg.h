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

#ifndef _KYCG_H
#define _KYCG_H

#include "ui.h"

#define KYCG_VERSION "0.1.0"

/*
 * Shared styling for help text. Each expands to an empty string off a TTY (see
 * ui.c), so a redirected --help is byte-identical to plain text.
 */
#define KYCG_H_TITLE kycg_ui_bold()
#define KYCG_H_KEY   kycg_ui_cyan()
#define KYCG_H_NOTE  kycg_ui_dim()
#define KYCG_H_WARN  kycg_ui_yellow()
#define KYCG_H_OFF   kycg_ui_reset()

/**
 * Resolve a knowledgebase spec to file paths in the store.
 *
 * Accepts what `kycg fetch` accepts -- "mm10", "mm10:CGI", "MSA:CGI,ChromHMM"
 * -- so a set is named the same way whether it is being downloaded or tested
 * against. A spec that is an existing file is returned as-is, which keeps
 * plain paths working and leaves no ambiguity to resolve.
 *
 * Returns the number of paths and fills *out with a malloc'd vector of
 * malloc'd strings (free with kycg_free_specs). Zero means nothing matched;
 * the caller reports that, since only it knows what to say.
 */
size_t kycg_resolve_spec(const char *spec, const char *store, char ***out);
void kycg_free_specs(char **v, size_t n);

/**
 * Fetch these specs into the store. Renders into the widget's panel when one
 * is open, so it can be called from inside a browser without leaving it.
 * Returns 0 on success.
 */
int kycg_fetch_specs(char *const *specs, size_t n, const char *store);

/** One set a collection publishes, and whether it is here. */
typedef struct {
  char *name;
  int   cached;
} kycg_catalogue_t;

/**
 * Every .cm a collection publishes, cached or not. Returns NULL when the
 * catalogue cannot be read; sets *n. Caller frees with kycg_catalogue_free.
 */
kycg_catalogue_t *kycg_catalogue(const char *target, const char *store,
                                 size_t *n);
void kycg_catalogue_free(kycg_catalogue_t *v, size_t n);

/**
 * The two browser callbacks that describe knowledgebases rather than move
 * them. Both are shared with `kycg test`'s picker so that `r` and `i` mean
 * exactly the same thing in both trees -- a second implementation would be
 * free to drift, which is the failure this project exists to avoid.
 *
 * kycg_kb_recommended matches kycg_ui_preselect_fn; kycg_kb_detail matches
 * kycg_ui_detail_fn and backs the `i` pane, which stays open while the cursor
 * moves.
 */
int kycg_kb_recommended(void *ctx, const char *root, const char *key);
void kycg_kb_detail(void *ctx, const char *root, const char *child_key,
                    int cols, kycg_ui_detail_t *out);

/* Subcommand entry points, dispatched from main() on argv[1]. */
/**
 * A collection worth offering in the picker: its name, what kind of row space
 * it indexes, and how many rows that is.
 */
typedef struct {
  const char *name;
  const char *kind;   /* "whole genome" or "array" */
  uint64_t    rows;
} kycg_pick_target_t;

/**
 * Offer these collections in the catalogue browser and return what was chosen.
 *
 * Shared by `kycg test` and `kycg annotate`, which differ only in which
 * collections are worth offering: test filters by row count against its query,
 * annotate offers the array platforms, since a probe ID has no meaning without
 * one. `f` fetches inside the browser so a missing set can be downloaded and
 * used in one sitting; `verb_key` ends it.
 *
 * Returns the count and fills *out with malloc'd "target:file" specs (free
 * with kycg_free_specs). 0 means the user quit; (size_t)-1 means the terminal
 * cannot host the browser.
 */
size_t kycg_pick_sets(const kycg_pick_target_t *targets, size_t n_targets,
                      const char *title, char verb_key, const char *verb,
                      char ***out);

int main_test(int argc, char *argv[]);
int main_info(int argc, char *argv[]);
int main_fetch(int argc, char *argv[]);
int main_annotate(int argc, char *argv[]);

#endif /* _KYCG_H */
