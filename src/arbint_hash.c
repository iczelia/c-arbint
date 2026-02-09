/*  arbint - portable arbitrary-precision computation library

    Copyright (C) 2026 Kamila Szewczyk (k@iczelia.net)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program. If not, see <https://www.gnu.org/licenses/>.  */

#include "arbint_hash.h"

#include "config.h"

#include "arbint_cpu.h"

#include <string.h>

#define ARBINT_HASH_LIMB_BYTES (ARBINT_LIMB_BITS / 8u)

/* ========== Function Pointer Dispatch ========== */

typedef void (*arbint_sha256_compress_fn_t)(uint32_t state[8],
                                            const uint8_t block[64]);

typedef uint32_t (*arbint_crc32c_fn_t)(uint32_t crc, const uint8_t * data,
                                       size_t len);

static arbint_sha256_compress_fn_t arbint_select_sha256_compress(void) {
#if HAS_SHA_NI_ALWAYS
  return arbint_sha256_compress_shani;
#elif HAS_SHA_NI
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_SHA)
             ? arbint_sha256_compress_shani
             : arbint_sha256_compress_generic;
#else
  return arbint_sha256_compress_generic;
#endif /* HAS_SHA_NI_ALWAYS */
}

static arbint_crc32c_fn_t arbint_select_crc32c(void) {
#if HAS_SSE42_CRC32_ALWAYS
  return arbint_crc32c_sse42;
#elif HAS_SSE42_CRC32
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_CRC32)
             ? arbint_crc32c_sse42
             : arbint_crc32c_generic;
#else
  return arbint_crc32c_generic;
#endif /* HAS_SSE42_CRC32_ALWAYS */
}

static arbint_sha256_compress_fn_t sha256_compress = NULL;
static arbint_crc32c_fn_t crc32c_impl = NULL;

/* ========== SHA-256 Streaming Interface ========== */

/*  SHA-256 initial hash values (FIPS 180-4 section 5.3.3).  */
static void arbint_sha256_init(arbint_sha256_state_t * st) {
  st->h[0] = 0x6a09e667u;
  st->h[1] = 0xbb67ae85u;
  st->h[2] = 0x3c6ef372u;
  st->h[3] = 0xa54ff53au;
  st->h[4] = 0x510e527fu;
  st->h[5] = 0x9b05688cu;
  st->h[6] = 0x1f83d9abu;
  st->h[7] = 0x5be0cd19u;
  st->buf_len = 0u;
  st->total_len = 0u;
}

/*  Feed data into SHA-256 state, compressing complete 64-byte blocks.  */
static void arbint_sha256_update(arbint_sha256_state_t * st,
                                 const uint8_t * data, size_t len) {
  size_t i = 0u;

  st->total_len += (uint64_t) len;

  /*  If we have buffered data, try to fill a block.  */
  if (st->buf_len > 0u) {
    size_t need = 64u - st->buf_len;
    if (len < need) {
      memcpy(st->buf + st->buf_len, data, len);
      st->buf_len += len;
      return;
    }
    memcpy(st->buf + st->buf_len, data, need);
    sha256_compress(st->h, st->buf);
    st->buf_len = 0u;
    i = need;
  }

  /*  Process complete blocks directly from input.  */
  while (i + 64u <= len) {
    sha256_compress(st->h, data + i);
    i += 64u;
  }

  /*  Buffer any remaining bytes.  */
  if (i < len) {
    st->buf_len = len - i;
    memcpy(st->buf, data + i, st->buf_len);
  }
}

/*  Store a 32-bit value as big-endian bytes.  */
static void arbint_sha256_store_be32(uint8_t * dst, uint32_t v) {
  dst[0] = (uint8_t) (v >> 24u);
  dst[1] = (uint8_t) (v >> 16u);
  dst[2] = (uint8_t) (v >> 8u);
  dst[3] = (uint8_t) v;
}

/*  Finalize SHA-256: apply FIPS 180-4 padding and write digest.  */
static void arbint_sha256_final(arbint_sha256_state_t * st,
                                uint8_t digest[32]) {
  uint64_t bit_len = st->total_len * 8u;
  size_t pad_start;
  unsigned i;

  /*  Append 0x80 byte.  */
  st->buf[st->buf_len++] = 0x80u;
  pad_start = st->buf_len;

  /*  If not enough room for 8-byte length, pad and compress, then
      start a fresh block.  */
  if (pad_start > 56u) {
    memset(st->buf + pad_start, 0, 64u - pad_start);
    sha256_compress(st->h, st->buf);
    pad_start = 0u;
  }

  /*  Zero-pad up to byte 56, then append 64-bit big-endian bit count.  */
  memset(st->buf + pad_start, 0, 56u - pad_start);
  st->buf[56] = (uint8_t) (bit_len >> 56u);
  st->buf[57] = (uint8_t) (bit_len >> 48u);
  st->buf[58] = (uint8_t) (bit_len >> 40u);
  st->buf[59] = (uint8_t) (bit_len >> 32u);
  st->buf[60] = (uint8_t) (bit_len >> 24u);
  st->buf[61] = (uint8_t) (bit_len >> 16u);
  st->buf[62] = (uint8_t) (bit_len >> 8u);
  st->buf[63] = (uint8_t) bit_len;
  sha256_compress(st->h, st->buf);

  /*  Write digest as big-endian words.  */
  for (i = 0u; i < 8u; ++i)
    arbint_sha256_store_be32(digest + i * 4u, st->h[i]);
}

/* ========== Canonical Data Feeder ========== */

/*  Callback type: receives (ctx, data, len) for each chunk.  */
typedef void (*arbint_hash_feed_fn_t)(void * ctx, const uint8_t * data,
                                      size_t len);

/*  Count the number of significant bytes in a limb (0 for zero).  */
static size_t arbint_hash_limb_used_bytes(arbint_limb_t x) {
  size_t n = 0u;
  while (x != (arbint_limb_t) 0u) {
    ++n;
    x >>= 8u;
  }
  return n;
}

/*  Feed the canonical representation of op into a callback.
    Canonical form: sign byte (0x00/0x01/0xff) followed by the magnitude
    in little-endian byte order with no leading zero bytes.  */
static void arbint_feed_canonical(const arbint_t op,
                                  arbint_hash_feed_fn_t feed, void * ctx) {
  int sign;
  size_t used;
  const arbint_limb_t * p;
  uint8_t sign_byte;

  sign = (op[0]._sz > 0) - (op[0]._sz < 0);
  used = arbint_abs_sz(op[0]._sz);

  /*  Sign byte.  */
  if (sign == 0)
    sign_byte = 0x00u;
  else if (sign > 0)
    sign_byte = 0x01u;
  else
    sign_byte = 0xffu;
  feed(ctx, &sign_byte, 1u);

  if (used == 0u)
    return;

  p = ARBINT_CLIMBS(op);

  /*  Feed all limbs except the top one as full-width bytes.  */
  {
    size_t i;
    for (i = 0u; i < used - 1u; ++i) {
      arbint_limb_t w = p[i];
      uint8_t tmp[ARBINT_HASH_LIMB_BYTES];
      size_t j;
      for (j = 0u; j < ARBINT_HASH_LIMB_BYTES; ++j) {
        tmp[j] = (uint8_t) (w & (arbint_limb_t) 0xffu);
        w >>= 8u;
      }
      feed(ctx, tmp, ARBINT_HASH_LIMB_BYTES);
    }
  }

  /*  Feed the top limb with only its significant bytes.  */
  {
    arbint_limb_t w = p[used - 1u];
    size_t top_bytes = arbint_hash_limb_used_bytes(w);
    uint8_t tmp[ARBINT_HASH_LIMB_BYTES];
    size_t j;
    for (j = 0u; j < top_bytes; ++j) {
      tmp[j] = (uint8_t) (w & (arbint_limb_t) 0xffu);
      w >>= 8u;
    }
    feed(ctx, tmp, top_bytes);
  }
}

/* ========== SHA-256 Feed Callback ========== */

static void arbint_sha256_feed_cb(void * ctx, const uint8_t * data,
                                  size_t len) {
  arbint_sha256_update((arbint_sha256_state_t *) ctx, data, len);
}

/* ========== CRC32C Feed Callback ========== */

typedef struct arbint_crc32c_ctx {
  uint32_t crc;
  arbint_crc32c_fn_t fn;
} arbint_crc32c_ctx_t;

static void arbint_crc32c_feed_cb(void * ctx, const uint8_t * data,
                                  size_t len) {
  arbint_crc32c_ctx_t * st = (arbint_crc32c_ctx_t *) ctx;
  st->crc = st->fn(st->crc, data, len);
}

/* ========== Public API ========== */

ARBINT_API arbint_err_t arbint_hash_slow(const arbint_t op, uint8_t * out_hash,
                                         size_t hash_len) {
  if (op == NULL || out_hash == NULL)
    return ARBINT_EINVAL;
  if (hash_len == 0u)
    return ARBINT_OK;

  /*  One-time initialization.  */
  if (sha256_compress == NULL)
    sha256_compress = arbint_select_sha256_compress();

  if (hash_len <= 32u) {
    /*  Single SHA-256, possibly truncated.  */
    arbint_sha256_state_t st;
    uint8_t digest[32];
    arbint_sha256_init(&st);
    arbint_feed_canonical(op, arbint_sha256_feed_cb, &st);
    arbint_sha256_final(&st, digest);
    memcpy(out_hash, digest, hash_len);
    return ARBINT_OK;
  }

  /*  Counter mode for hash_len > 32: SHA256(canonical || le32(counter))
      for counter = 0, 1, 2, ... producing 32 bytes each.  */
  {
    size_t offset = 0u;
    uint32_t counter = 0u;

    while (offset < hash_len) {
      arbint_sha256_state_t st;
      uint8_t digest[32];
      uint8_t ctr_bytes[4];
      size_t chunk;

      arbint_sha256_init(&st);
      arbint_feed_canonical(op, arbint_sha256_feed_cb, &st);

      /*  Append 4-byte little-endian counter.  */
      ctr_bytes[0] = (uint8_t) (counter & 0xffu);
      ctr_bytes[1] = (uint8_t) ((counter >> 8u) & 0xffu);
      ctr_bytes[2] = (uint8_t) ((counter >> 16u) & 0xffu);
      ctr_bytes[3] = (uint8_t) ((counter >> 24u) & 0xffu);
      arbint_sha256_update(&st, ctr_bytes, 4u);

      arbint_sha256_final(&st, digest);

      chunk = hash_len - offset;
      if (chunk > 32u)
        chunk = 32u;
      memcpy(out_hash + offset, digest, chunk);
      offset += chunk;
      ++counter;
    }
    return ARBINT_OK;
  }
}

ARBINT_API arbint_err_t arbint_hash_fast(const arbint_t op,
                                         uint32_t * out_hash) {
  arbint_crc32c_ctx_t ctx;

  if (op == NULL || out_hash == NULL)
    return ARBINT_EINVAL;

  /*  One-time initialization.  */
  if (crc32c_impl == NULL)
    crc32c_impl = arbint_select_crc32c();

  ctx.crc = 0xFFFFFFFFu;
  ctx.fn = crc32c_impl;
  arbint_feed_canonical(op, arbint_crc32c_feed_cb, &ctx);

  *out_hash = ctx.crc ^ 0xFFFFFFFFu;
  return ARBINT_OK;
}
