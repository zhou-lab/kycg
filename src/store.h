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

#ifndef _KYCG_STORE_H
#define _KYCG_STORE_H

#include <stddef.h>

/**
 * The local knowledgebase store.
 *
 * Layout is deliberately flat: a .cm in the store is an ordinary file that
 * `kycg test -m` takes by path. There is no database and no lockfile, so the
 * store stays inspectable with ls and verifiable with shasum.
 *
 * The unified asset layout shared across the tool suite:
 *
 *   <root>/InfiniumAnnotation/<PLATFORM>/<PLATFORM>.ordering.tsv.gz   arrays
 *   <root>/InfiniumAnnotation/<PLATFORM>/KYCG/<Set>.<date>.cm        arrays
 *   <root>/KYCGKB/<genome>/<Set>.<date>.cm                          sequencing
 */

/**
 * Resolve the store root. Thin shim over libyame's yame_assets_root:
 * `override`, else $KYCG_DATA_DIR, else $YAME_DATA_HOME, else the shared data
 * tier (${XDG_DATA_HOME:-~/.local/share}/yame) -- so the store kycg reads is
 * the one every tool in the suite fills.
 */
const char *kycg_store_root(const char *override);

/** mkdir -p. Returns 0 on success. */
int kycg_store_mkdir_p(const char *path);

/**
 * Is this manifest entry a name we are willing to write to?
 *
 * The anchor proves a manifest is the one upstream published. It says nothing
 * about whether upstream should get to choose *where* a file lands -- and the
 * name from a manifest line was being joined straight onto the store directory
 * with snprintf("%s/%s"). An entry reading
 *
 *     <64 hex>  ../../../../.ssh/authorized_keys
 *
 * therefore resolved outside the store entirely, and since whoever writes the
 * manifest also serves the bytes, the digest check would pass and the file
 * would be renamed into place. Verified by construction, not theory: joining
 * that name to a real store directory and fopen()ing it writes the file
 * outside the store.
 *
 * A network attacker cannot reach this -- they cannot forge a manifest that
 * satisfies the compiled-in anchor. The realistic route is an upstream
 * repository compromise that a maintainer then pins with make_registry.sh,
 * which hashes whatever it is served without inspecting the names. That is
 * exactly the case where authentication is not authorization.
 *
 * The store layout is flat, so the rule is simply: no separators, no dot-
 * leading names, nothing that could be a traversal component. All 196 entries
 * across the currently published manifests are plain filenames, so this
 * rejects nothing legitimate.
 */
int kycg_store_safe_name(const char *s);

/** Nonzero if `path` is an existing regular file. */
int kycg_store_is_file(const char *path);

/**
 * Every .cm in the store, as full paths, sorted. Returns a malloc'd vector of
 * malloc'd strings; sets *n. Index sidecars (.cm.idx) are excluded — they are
 * not independently testable.
 */
char **kycg_store_find_cm(const char *root, size_t *n);

void kycg_store_free_list(char **v, size_t n);

/**
 * A store path rendered for display: "mm10/ChromHMM.20220414.cm" rather than
 * the absolute path, which is mostly the user's home directory repeated.
 * Returns a pointer into `path`.
 */
const char *kycg_store_relative(const char *root, const char *path);

#endif /* _KYCG_STORE_H */
