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

#include "test_framework.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

int main(void) {
  ARBINT_TEST_DECLARE_FAILURES();

  arbint_ctx_t ctx;
  arbint_t a, b;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);

  /*  Test 1: NULL argument checks.  */
  {
    uint8_t buf[32];
    uint32_t fast_out;

    CHECK_EQ_I(arbint_hash_slow(NULL, buf, 32u), ARBINT_EINVAL);
    CHECK_EQ_I(arbint_hash_slow(a, NULL, 32u), ARBINT_EINVAL);
    CHECK_EQ_I(arbint_hash_fast(NULL, &fast_out), ARBINT_EINVAL);
    CHECK_EQ_I(arbint_hash_fast(a, NULL), ARBINT_EINVAL);
  }

  /*  Test 2: hash_len == 0 returns OK without writing.  */
  {
    uint8_t buf[1] = {0xAAu};
    CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, buf, 0u), ARBINT_OK);
    CHECK_EQ_I(buf[0], 0xAAu); /*  Untouched.  */
  }

  /*  Test 3: Determinism — hashing the same value twice gives identical
   * output.  */
  {
    uint8_t h1[32];
    uint8_t h2[32];
    uint32_t f1, f2;

    CHECK_EQ_I(arbint_set_i32(a, 12345), ARBINT_OK);

    CHECK_EQ_I(arbint_hash_slow(a, h1, 32u), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, h2, 32u), ARBINT_OK);
    CHECK(memcmp(h1, h2, 32u) == 0);

    CHECK_EQ_I(arbint_hash_fast(a, &f1), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_fast(a, &f2), ARBINT_OK);
    CHECK_EQ_I(f1, f2);
  }

  /*  Test 4: Sign sensitivity — different signs produce different hashes.  */
  {
    uint8_t h_pos[32], h_neg[32], h_zero[32];
    uint32_t f_pos, f_neg, f_zero;

    CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, h_pos, 32u), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_fast(a, &f_pos), ARBINT_OK);

    CHECK_EQ_I(arbint_set_i32(a, -42), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, h_neg, 32u), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_fast(a, &f_neg), ARBINT_OK);

    CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, h_zero, 32u), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_fast(a, &f_zero), ARBINT_OK);

    CHECK(memcmp(h_pos, h_neg, 32u) != 0);
    CHECK(memcmp(h_pos, h_zero, 32u) != 0);
    CHECK(memcmp(h_neg, h_zero, 32u) != 0);
    CHECK_NE_I(f_pos, f_neg);
    CHECK_NE_I(f_pos, f_zero);
    CHECK_NE_I(f_neg, f_zero);
  }

  /*  Test 5: Different values produce different hashes.  */
  {
    uint8_t h1[32], h2[32];
    uint32_t f1, f2;

    CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, 2), ARBINT_OK);

    CHECK_EQ_I(arbint_hash_slow(a, h1, 32u), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(b, h2, 32u), ARBINT_OK);
    CHECK(memcmp(h1, h2, 32u) != 0);

    CHECK_EQ_I(arbint_hash_fast(a, &f1), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_fast(b, &f2), ARBINT_OK);
    CHECK_NE_I(f1, f2);
  }

  /*  Test 6: Truncation — hash_len=16 matches first 16 bytes of hash_len=32.
   */
  {
    uint8_t full[32];
    uint8_t trunc[16];

    CHECK_EQ_I(arbint_set_i32(a, 99999), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, full, 32u), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, trunc, 16u), ARBINT_OK);
    CHECK(memcmp(full, trunc, 16u) == 0);
  }

  /*  Test 7: Counter mode — hash_len=64.
      First 32 bytes should NOT match hash_len=32 (counter mode appends
      counter suffix for len>32, whereas len<=32 has no counter suffix).
      But the two 32-byte halves of the 64-byte output should differ.  */
  {
    uint8_t h64[64];

    CHECK_EQ_I(arbint_set_i32(a, 777), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, h64, 64u), ARBINT_OK);

    /*  The two halves should differ (counter 0 vs counter 1).  */
    CHECK(memcmp(h64, h64 + 32, 32u) != 0);
  }

  /*  Test 8: Counter mode — hash_len=33.
      The first 32 bytes = SHA256(canonical || le32(0)).
      Byte 33 = first byte of SHA256(canonical || le32(1)).
      Verify consistency: 64-byte output's first 33 bytes match.  */
  {
    uint8_t h33[33];
    uint8_t h64[64];

    CHECK_EQ_I(arbint_set_i32(a, 777), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, h33, 33u), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, h64, 64u), ARBINT_OK);
    CHECK(memcmp(h33, h64, 33u) == 0);
  }

  /*  Test 9: Large number — build via squaring, verify deterministic.  */
  {
    uint8_t h1[32], h2[32];
    uint32_t f1, f2;
    int i;

    CHECK_EQ_I(arbint_set_u32(a, 1000000u), ARBINT_OK);
    for (i = 0; i < 10; ++i)
      CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);

    CHECK_EQ_I(arbint_hash_slow(a, h1, 32u), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, h2, 32u), ARBINT_OK);
    CHECK(memcmp(h1, h2, 32u) == 0);

    CHECK_EQ_I(arbint_hash_fast(a, &f1), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_fast(a, &f2), ARBINT_OK);
    CHECK_EQ_I(f1, f2);
  }

  /*  Test 10: Large number with different value gives different hash.  */
  {
    uint8_t h1[32], h2[32];
    uint32_t f1, f2;
    int i;

    CHECK_EQ_I(arbint_set_u32(a, 1000000u), ARBINT_OK);
    CHECK_EQ_I(arbint_set_u32(b, 1000001u), ARBINT_OK);
    for (i = 0; i < 10; ++i) {
      CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
      CHECK_EQ_I(arbint_sqr(b, b), ARBINT_OK);
    }

    CHECK_EQ_I(arbint_hash_slow(a, h1, 32u), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(b, h2, 32u), ARBINT_OK);
    CHECK(memcmp(h1, h2, 32u) != 0);

    CHECK_EQ_I(arbint_hash_fast(a, &f1), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_fast(b, &f2), ARBINT_OK);
    CHECK_NE_I(f1, f2);
  }

  /*  Test 11: hash_fast of zero.  */
  {
    uint32_t h;
    CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_fast(a, &h), ARBINT_OK);
    /*  Just verify it doesn't crash and returns some value.
        CRC32C of single byte 0x00 with init 0xFFFFFFFF is deterministic.  */
    (void) h;
  }

  /*  Test 12: Small values all produce unique fast hashes.  */
  {
    uint32_t hashes[20];
    int i, j;
    for (i = 0; i < 20; ++i) {
      CHECK_EQ_I(arbint_set_i32(a, i), ARBINT_OK);
      CHECK_EQ_I(arbint_hash_fast(a, &hashes[i]), ARBINT_OK);
    }
    for (i = 0; i < 20; ++i) {
      for (j = i + 1; j < 20; ++j)
        CHECK_NE_I(hashes[i], hashes[j]);
    }
  }

  /*  Test 13: Very large hash_len (128 bytes) works and is deterministic.  */
  {
    uint8_t h1[128], h2[128];
    CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, h1, 128u), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, h2, 128u), ARBINT_OK);
    CHECK(memcmp(h1, h2, 128u) == 0);
  }

  /*  Test 14: hash_len=1 works.  */
  {
    uint8_t h1, h2;
    CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, &h1, 1u), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, &h2, 1u), ARBINT_OK);
    CHECK_EQ_I(h1, h2);
  }

  /*  Test 15: Hashing zero with hash_slow.  */
  {
    uint8_t h1[32], h2[32];
    CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, h1, 32u), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, h2, 32u), ARBINT_OK);
    CHECK(memcmp(h1, h2, 32u) == 0);
  }

  /*  Test 16: Known-answer vectors for hash_slow.
      Canonical form is sign byte + LE magnitude (no leading zero bytes).
      Expected digests verified against system sha256sum.  */
  {
    /*  SHA-256("\x01\x2a") = SHA-256 of {sign=+1, mag=42}  */
    static const uint8_t expect_42[32] = {
        0x12u, 0xa0u, 0xf6u, 0x5cu, 0xb2u, 0x57u, 0x38u, 0xc3u,
        0x25u, 0x1fu, 0x2du, 0xdfu, 0xabu, 0x71u, 0x29u, 0xfbu,
        0x80u, 0xdeu, 0x0fu, 0x7fu, 0x05u, 0xe3u, 0xe1u, 0x05u,
        0xccu, 0xacu, 0x2fu, 0x2bu, 0x71u, 0x07u, 0x6eu, 0x9du};
    /*  SHA-256("\x00") = SHA-256 of {sign=0}  */
    static const uint8_t expect_0[32] = {
        0x6eu, 0x34u, 0x0bu, 0x9cu, 0xffu, 0xb3u, 0x7au, 0x98u,
        0x9cu, 0xa5u, 0x44u, 0xe6u, 0xbbu, 0x78u, 0x0au, 0x2cu,
        0x78u, 0x90u, 0x1du, 0x3fu, 0xb3u, 0x37u, 0x38u, 0x76u,
        0x85u, 0x11u, 0xa3u, 0x06u, 0x17u, 0xafu, 0xa0u, 0x1du};
    /*  SHA-256("\xff\x01") = SHA-256 of {sign=-1, mag=1}  */
    static const uint8_t expect_neg1[32] = {
        0x43u, 0x7cu, 0xb4u, 0x3au, 0x30u, 0x22u, 0x6eu, 0x63u,
        0x9bu, 0x33u, 0xd8u, 0x45u, 0x33u, 0xcfu, 0xc3u, 0xddu,
        0xd9u, 0x70u, 0xa4u, 0x00u, 0x55u, 0xb9u, 0xc8u, 0xc0u,
        0xbfu, 0x7eu, 0x37u, 0x95u, 0xcfu, 0x18u, 0x4eu, 0xb6u};
    uint8_t h[32];

    CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, h, 32u), ARBINT_OK);
    CHECK(memcmp(h, expect_42, 32u) == 0);

    CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, h, 32u), ARBINT_OK);
    CHECK(memcmp(h, expect_0, 32u) == 0);

    CHECK_EQ_I(arbint_set_i32(a, -1), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_slow(a, h, 32u), ARBINT_OK);
    CHECK(memcmp(h, expect_neg1, 32u) == 0);
  }

  /*  Test 17: Known-answer vectors for hash_fast (CRC32C).
      CRC32C(init=0xFFFFFFFF, data) ^ 0xFFFFFFFF.  */
  {
    uint32_t h;

    CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_fast(a, &h), ARBINT_OK);
    CHECK_EQ_I(h, 0xa99c4943u);

    CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_fast(a, &h), ARBINT_OK);
    CHECK_EQ_I(h, 0x527d5351u);

    CHECK_EQ_I(arbint_set_i32(a, -1), ARBINT_OK);
    CHECK_EQ_I(arbint_hash_fast(a, &h), ARBINT_OK);
    CHECK_EQ_I(h, 0xa0e9d052u);
  }

  /*  Cleanup  */
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);

  ARBINT_TEST_FINISH("test_hash");
}
