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

#ifndef _KYCG_DIGEST_H
#define _KYCG_DIGEST_H

#include <stddef.h>
#include <stdint.h>

/**
 * Hex digest of a file. `out` must hold 65 bytes (sha256) or 33 (md5),
 * including the terminating NUL. Returns 0 on success, -1 if the file could
 * not be read.
 */
int kycg_sha256_file(const char *path, char out[65]);
/* md5 was needed while the whole-genome sets came from Zenodo, whose record
 * API publishes nothing stronger. Both channels now publish SHA256SUMS, so
 * nothing calls this; it is kept because the archival Zenodo deposits still
 * carry md5 and a future verifier against them would want it. */
int kycg_md5_file(const char *path, char out[33]);

/** Hex sha256 of a memory buffer. `out` must hold 65 bytes. */
void kycg_sha256_buf(const void *data, size_t len, char out[65]);

/** Constant-time-ish case-insensitive hex comparison. Nonzero if equal. */
int kycg_digest_equal(const char *a, const char *b);

#endif /* _KYCG_DIGEST_H */
