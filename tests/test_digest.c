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
 * Tests for the digests every download is checked against.
 *
 * The expected values are the published RFC/FIPS test vectors, not values
 * this code produced -- a digest test that checks an implementation against
 * itself would pass just as happily on a broken one. The empty-input and
 * multi-block cases are here because they are where hand-written padding goes
 * wrong: the length is encoded in bits, big-endian for SHA-256 and
 * little-endian for MD5, and a buffer that ends exactly on a block boundary
 * needs a whole extra block of padding.
 *
 * MD5 is dead code today -- nothing fetches from Zenodo any more -- but it is
 * the only verifier for the archival md5 deposits, so it is kept, and kept
 * correct.
 */

#include "digest.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static int failures = 0;

static void want(const char *label, const char *expect, const char *got) {
  if (strcmp(expect, got) != 0) {
    printf("  FAIL %s\n       expected: %s\n       actual:   %s\n",
           label, expect, got);
    ++failures;
  }
}

/* A temp file holding `data`, since the file digests take a path. */
static char *tmp_with(const void *data, size_t len, char *path) {
  snprintf(path, 64, "/tmp/kycg_digest_%d_%p", (int)getpid(), data);
  FILE *f = fopen(path, "wb");
  if (!f) return NULL;
  if (len) fwrite(data, 1, len, f);
  fclose(f);
  return path;
}

int main(void) {
  char out[65], path[64];

  /* ---- SHA-256, FIPS 180-2 vectors ------------------------------------- */
  kycg_sha256_buf("", 0, out);
  want("sha256 of the empty string",
       "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855", out);

  kycg_sha256_buf("abc", 3, out);
  want("sha256(\"abc\")",
       "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad", out);

  /* Two blocks: 56 bytes is one byte past what fits with the length field, so
   * the padding spills into a second block. */
  kycg_sha256_buf("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq", 56, out);
  want("sha256 of a two-block message",
       "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1", out);

  /* Exactly one block (64 bytes): padding needs a whole extra block. */
  {
    char buf[64];
    memset(buf, 'a', sizeof(buf));
    kycg_sha256_buf(buf, sizeof(buf), out);
    want("sha256 of exactly one block",
         "ffe054fe7ae0cb6dc65c3af9b61d5209f439851db43d0ba5997337df154668eb", out);
  }

  /* ---- the file forms read the same bytes ------------------------------ */
  if (tmp_with("abc", 3, path)) {
    if (kycg_sha256_file(path, out) != 0) {
      printf("  FAIL kycg_sha256_file failed on a readable file\n"); ++failures;
    } else {
      want("sha256_file(\"abc\")",
           "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad", out);
    }
    unlink(path);
  }

  /* A file larger than the read buffer, so the streaming loop runs more than
   * once. 200000 'a' bytes: the digest is a published value for that input. */
  {
    size_t n = 200000;
    char *big = malloc(n);
    if (big) {
      memset(big, 'a', n);
      kycg_sha256_buf(big, n, out);
      char expect[65];
      snprintf(expect, sizeof(expect), "%s", out);
      if (tmp_with(big, n, path)) {
        char fout[65];
        kycg_sha256_file(path, fout);
        want("sha256_file matches sha256_buf over several buffers", expect, fout);
        unlink(path);
      }
      free(big);
    }
  }

  /* A file that is not there is an error, not a digest of nothing: silently
   * hashing the empty string would make a missing download verify. */
  if (kycg_sha256_file("/nonexistent/kycg/test/file", out) == 0) {
    printf("  FAIL sha256_file succeeded on a missing file\n"); ++failures;
  }

  /* ---- MD5, RFC 1321 vectors ------------------------------------------- */
  {
    char md5[33];
    if (tmp_with("", 0, path)) {
      kycg_md5_file(path, md5);
      want("md5 of the empty file", "d41d8cd98f00b204e9800998ecf8427e", md5);
      unlink(path);
    }
    if (tmp_with("abc", 3, path)) {
      kycg_md5_file(path, md5);
      want("md5(\"abc\")", "900150983cd24fb0d6963f7d28e17f72", md5);
      unlink(path);
    }
    if (tmp_with("message digest", 14, path)) {
      kycg_md5_file(path, md5);
      want("md5(\"message digest\")", "f96b697d7cb7938d525a2f31aaf161d0", md5);
      unlink(path);
    }
    /* 80 bytes: past one block, and the RFC's longest vector. */
    if (tmp_with("12345678901234567890123456789012345678901234567890"
                 "123456789012345678901234567890", 80, path)) {
      kycg_md5_file(path, md5);
      want("md5 of the 80-byte vector", "57edf4a22be3c955ac49da2e2107b67a", md5);
      unlink(path);
    }
    if (kycg_md5_file("/nonexistent/kycg/test/file", md5) == 0) {
      printf("  FAIL md5_file succeeded on a missing file\n"); ++failures;
    }
  }

  /* ---- the comparison used on every verify ----------------------------- */
  {
    const char *a = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad";
    const char *A = "BA7816BF8F01CFEA414140DE5DAE2223B00361A396177A9CB410FF61F20015AD";
    if (!kycg_digest_equal(a, a)) { printf("  FAIL a digest differs from itself\n"); ++failures; }
    if (!kycg_digest_equal(a, A)) { printf("  FAIL case matters in a hex digest\n"); ++failures; }
    if (kycg_digest_equal(a, "deadbeef")) { printf("  FAIL a short digest matched\n"); ++failures; }
    if (kycg_digest_equal(a, NULL) || kycg_digest_equal(NULL, a) ||
        kycg_digest_equal(NULL, NULL)) {
      printf("  FAIL a NULL digest compared equal\n"); ++failures;
    }
    /* One bit apart: the whole point of the check. */
    if (kycg_digest_equal(a,
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ae")) {
      printf("  FAIL two digests differing in one nibble matched\n"); ++failures;
    }
  }

  if (failures) { printf("  %d check(s) failed\n", failures); return 1; }
  printf("  all checks passed\n");
  return 0;
}
