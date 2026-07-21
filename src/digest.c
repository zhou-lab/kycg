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
 * SHA-256 (FIPS 180-4) and MD5 (RFC 1321).
 *
 * Self-contained on purpose. OpenSSL would be an entire dependency for two
 * hash functions, and CommonCrypto is macOS-only; a few hundred lines is
 * cheaper than either and keeps the bioconda surface at htslib + zlib.
 * The SHA-256 core follows the same structure as sesame-cli's src/sha256.c,
 * which solves this identical problem for the same annotation repository.
 *
 * WHY ONLY SHA-256
 *   Both resource channels publish a SHA256SUMS per directory, so both are
 *   sha256 chains anchored on a digest compiled into the binary by
 *   tools/make_registry.sh. kycg therefore never trusts a digest it fetched at
 *   run time.
 *
 *   This was not always true. The whole-genome sets came from Zenodo, whose
 *   record API publishes an md5 per file and nothing stronger, so that channel
 *   was md5 and this file carried both algorithms. Once KYCGKB_<genome> began
 *   publishing its own SHA256SUMS the md5 half became dead weight and was
 *   removed; Zenodo remains the citable archive but is no longer fetched from.
 */

#include "digest.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define KYCG_DIGEST_BUFSZ 65536

/* ------------------------------------------------------------------ sha256 */

typedef struct {
  uint32_t h[8];
  uint64_t len;
  uint8_t  buf[64];
  size_t   n;
} sha256_ctx;

static const uint32_t SHA256_K[64] = {
  0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
  0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
  0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
  0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
  0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
  0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
  0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
  0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

#define ROR32(x,n) (((x) >> (n)) | ((x) << (32 - (n))))

static void sha256_block(sha256_ctx *c, const uint8_t *p) {
  uint32_t w[64], a, b, cc, d, e, f, g, h;
  int i;

  for (i = 0; i < 16; ++i)
    w[i] = ((uint32_t)p[i*4] << 24) | ((uint32_t)p[i*4+1] << 16) |
           ((uint32_t)p[i*4+2] << 8) | (uint32_t)p[i*4+3];
  for (; i < 64; ++i) {
    uint32_t s0 = ROR32(w[i-15],7) ^ ROR32(w[i-15],18) ^ (w[i-15] >> 3);
    uint32_t s1 = ROR32(w[i-2],17) ^ ROR32(w[i-2],19) ^ (w[i-2] >> 10);
    w[i] = w[i-16] + s0 + w[i-7] + s1;
  }

  a=c->h[0]; b=c->h[1]; cc=c->h[2]; d=c->h[3];
  e=c->h[4]; f=c->h[5]; g=c->h[6];  h=c->h[7];

  for (i = 0; i < 64; ++i) {
    uint32_t S1 = ROR32(e,6) ^ ROR32(e,11) ^ ROR32(e,25);
    uint32_t ch = (e & f) ^ (~e & g);
    uint32_t t1 = h + S1 + ch + SHA256_K[i] + w[i];
    uint32_t S0 = ROR32(a,2) ^ ROR32(a,13) ^ ROR32(a,22);
    uint32_t mj = (a & b) ^ (a & cc) ^ (b & cc);
    uint32_t t2 = S0 + mj;
    h=g; g=f; f=e; e=d+t1; d=cc; cc=b; b=a; a=t1+t2;
  }

  c->h[0]+=a; c->h[1]+=b; c->h[2]+=cc; c->h[3]+=d;
  c->h[4]+=e; c->h[5]+=f; c->h[6]+=g;  c->h[7]+=h;
}

static void sha256_init(sha256_ctx *c) {
  c->h[0]=0x6a09e667; c->h[1]=0xbb67ae85; c->h[2]=0x3c6ef372; c->h[3]=0xa54ff53a;
  c->h[4]=0x510e527f; c->h[5]=0x9b05688c; c->h[6]=0x1f83d9ab; c->h[7]=0x5be0cd19;
  c->len = 0; c->n = 0;
}

static void sha256_update(sha256_ctx *c, const void *data, size_t len) {
  const uint8_t *p = data;
  c->len += len;
  while (len) {
    size_t take = 64 - c->n;
    if (take > len) take = len;
    memcpy(c->buf + c->n, p, take);
    c->n += take; p += take; len -= take;
    if (c->n == 64) { sha256_block(c, c->buf); c->n = 0; }
  }
}

static void sha256_final(sha256_ctx *c, uint8_t out[32]) {
  uint64_t bits = c->len * 8;
  int i;

  c->buf[c->n++] = 0x80;
  if (c->n > 56) {
    memset(c->buf + c->n, 0, 64 - c->n);
    sha256_block(c, c->buf);
    c->n = 0;
  }
  memset(c->buf + c->n, 0, 56 - c->n);
  for (i = 0; i < 8; ++i) c->buf[56+i] = (uint8_t)(bits >> (56 - 8*i));
  sha256_block(c, c->buf);

  for (i = 0; i < 8; ++i) {
    out[i*4]   = (uint8_t)(c->h[i] >> 24);
    out[i*4+1] = (uint8_t)(c->h[i] >> 16);
    out[i*4+2] = (uint8_t)(c->h[i] >> 8);
    out[i*4+3] = (uint8_t)(c->h[i]);
  }
}

static void hexify(const uint8_t *raw, size_t n, char *out) {
  static const char HEX[] = "0123456789abcdef";
  size_t i;
  for (i = 0; i < n; ++i) {
    out[i*2]   = HEX[raw[i] >> 4];
    out[i*2+1] = HEX[raw[i] & 0xf];
  }
  out[n*2] = '\0';
}

void kycg_sha256_buf(const void *data, size_t len, char out[65]) {
  sha256_ctx c;
  uint8_t raw[32];
  sha256_init(&c);
  sha256_update(&c, data, len);
  sha256_final(&c, raw);
  hexify(raw, 32, out);
}

int kycg_sha256_file(const char *path, char out[65]) {
  FILE *fp = fopen(path, "rb");
  if (!fp) return -1;

  sha256_ctx c;
  sha256_init(&c);

  static uint8_t buf[KYCG_DIGEST_BUFSZ];
  size_t got;
  while ((got = fread(buf, 1, sizeof(buf), fp)) > 0) sha256_update(&c, buf, got);

  int err = ferror(fp);
  fclose(fp);
  if (err) return -1;

  uint8_t raw[32];
  sha256_final(&c, raw);
  hexify(raw, 32, out);
  return 0;
}

/* --------------------------------------------------------------------- md5 */

typedef struct {
  uint32_t h[4];
  uint64_t len;
  uint8_t  buf[64];
  size_t   n;
} md5_ctx;

static const uint32_t MD5_K[64] = {
  0xd76aa478,0xe8c7b756,0x242070db,0xc1bdceee,0xf57c0faf,0x4787c62a,0xa8304613,0xfd469501,
  0x698098d8,0x8b44f7af,0xffff5bb1,0x895cd7be,0x6b901122,0xfd987193,0xa679438e,0x49b40821,
  0xf61e2562,0xc040b340,0x265e5a51,0xe9b6c7aa,0xd62f105d,0x02441453,0xd8a1e681,0xe7d3fbc8,
  0x21e1cde6,0xc33707d6,0xf4d50d87,0x455a14ed,0xa9e3e905,0xfcefa3f8,0x676f02d9,0x8d2a4c8a,
  0xfffa3942,0x8771f681,0x6d9d6122,0xfde5380c,0xa4beea44,0x4bdecfa9,0xf6bb4b60,0xbebfbc70,
  0x289b7ec6,0xeaa127fa,0xd4ef3085,0x04881d05,0xd9d4d039,0xe6db99e5,0x1fa27cf8,0xc4ac5665,
  0xf4292244,0x432aff97,0xab9423a7,0xfc93a039,0x655b59c3,0x8f0ccc92,0xffeff47d,0x85845dd1,
  0x6fa87e4f,0xfe2ce6e0,0xa3014314,0x4e0811a1,0xf7537e82,0xbd3af235,0x2ad7d2bb,0xeb86d391
};

static const uint8_t MD5_S[64] = {
  7,12,17,22, 7,12,17,22, 7,12,17,22, 7,12,17,22,
  5, 9,14,20, 5, 9,14,20, 5, 9,14,20, 5, 9,14,20,
  4,11,16,23, 4,11,16,23, 4,11,16,23, 4,11,16,23,
  6,10,15,21, 6,10,15,21, 6,10,15,21, 6,10,15,21
};

#define ROL32(x,n) (((x) << (n)) | ((x) >> (32 - (n))))

static void md5_block(md5_ctx *c, const uint8_t *p) {
  uint32_t m[16], a = c->h[0], b = c->h[1], cc = c->h[2], d = c->h[3];
  int i;

  for (i = 0; i < 16; ++i)
    m[i] = (uint32_t)p[i*4] | ((uint32_t)p[i*4+1] << 8) |
           ((uint32_t)p[i*4+2] << 16) | ((uint32_t)p[i*4+3] << 24);

  for (i = 0; i < 64; ++i) {
    uint32_t f;
    int g;
    if (i < 16)      { f = (b & cc) | (~b & d);      g = i; }
    else if (i < 32) { f = (d & b) | (~d & cc);      g = (5*i + 1) & 15; }
    else if (i < 48) { f = b ^ cc ^ d;               g = (3*i + 5) & 15; }
    else             { f = cc ^ (b | ~d);            g = (7*i) & 15; }

    f = f + a + MD5_K[i] + m[g];
    a = d; d = cc; cc = b;
    b = b + ROL32(f, MD5_S[i]);
  }

  c->h[0]+=a; c->h[1]+=b; c->h[2]+=cc; c->h[3]+=d;
}

static void md5_init(md5_ctx *c) {
  c->h[0]=0x67452301; c->h[1]=0xefcdab89; c->h[2]=0x98badcfe; c->h[3]=0x10325476;
  c->len = 0; c->n = 0;
}

static void md5_update(md5_ctx *c, const void *data, size_t len) {
  const uint8_t *p = data;
  c->len += len;
  while (len) {
    size_t take = 64 - c->n;
    if (take > len) take = len;
    memcpy(c->buf + c->n, p, take);
    c->n += take; p += take; len -= take;
    if (c->n == 64) { md5_block(c, c->buf); c->n = 0; }
  }
}

static void md5_final(md5_ctx *c, uint8_t out[16]) {
  uint64_t bits = c->len * 8;
  int i;

  c->buf[c->n++] = 0x80;
  if (c->n > 56) {
    memset(c->buf + c->n, 0, 64 - c->n);
    md5_block(c, c->buf);
    c->n = 0;
  }
  memset(c->buf + c->n, 0, 56 - c->n);
  /* MD5 length is little-endian, unlike SHA-256's big-endian. */
  for (i = 0; i < 8; ++i) c->buf[56+i] = (uint8_t)(bits >> (8*i));
  md5_block(c, c->buf);

  for (i = 0; i < 4; ++i) {
    out[i*4]   = (uint8_t)(c->h[i]);
    out[i*4+1] = (uint8_t)(c->h[i] >> 8);
    out[i*4+2] = (uint8_t)(c->h[i] >> 16);
    out[i*4+3] = (uint8_t)(c->h[i] >> 24);
  }
}

int kycg_md5_file(const char *path, char out[33]) {
  FILE *fp = fopen(path, "rb");
  if (!fp) return -1;

  md5_ctx c;
  md5_init(&c);

  static uint8_t buf[KYCG_DIGEST_BUFSZ];
  size_t got;
  while ((got = fread(buf, 1, sizeof(buf), fp)) > 0) md5_update(&c, buf, got);

  int err = ferror(fp);
  fclose(fp);
  if (err) return -1;

  uint8_t raw[16];
  md5_final(&c, raw);
  hexify(raw, 16, out);
  return 0;
}

/* ------------------------------------------------------------------ compare */

int kycg_digest_equal(const char *a, const char *b) {
  if (!a || !b) return 0;
  size_t i;
  int diff = 0;
  for (i = 0; a[i] && b[i]; ++i)
    diff |= tolower((unsigned char)a[i]) ^ tolower((unsigned char)b[i]);
  return diff == 0 && a[i] == '\0' && b[i] == '\0';
}
