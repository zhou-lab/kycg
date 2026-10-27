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


#ifndef _KYCG_ORDERING_H
#define _KYCG_ORDERING_H

#include <stddef.h>

typedef struct {
  char  **id;      /* probe IDs, in row order (which is sorted order) */
  size_t  n;
} kycg_ordering_t;

/**
 * Read <store>/<platform>/<platform>.ordering.tsv.gz into `out`.
 *
 * `spec_hint` is the -m spec the caller was given, used only to word the
 * `kycg fetch` that repairs a missing file; `cmd` names the calling command
 * in messages. Returns 0, or -1 after saying what went wrong.
 */
int  kycg_ordering_load(const char *store, const char *platform,
                        const char *spec_hint, const char *cmd,
                        kycg_ordering_t *out);

/** Row index of `id`, or -1. Binary search: the ordering is sorted. */
long kycg_ordering_find(const kycg_ordering_t *o, const char *id);

void kycg_ordering_free(kycg_ordering_t *o);

#endif /* _KYCG_ORDERING_H */
