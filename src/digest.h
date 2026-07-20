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
int kycg_md5_file(const char *path, char out[33]);

/** Hex sha256 of a memory buffer. `out` must hold 65 bytes. */
void kycg_sha256_buf(const void *data, size_t len, char out[65]);

/** Constant-time-ish case-insensitive hex comparison. Nonzero if equal. */
int kycg_digest_equal(const char *a, const char *b);

#endif /* _KYCG_DIGEST_H */
