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

ARBINT_TEST_DECLARE_FAILURES();

/* ========== arbint_sizeinbase tests ========== */

static void test_sizeinbase_zero(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  sizeinbase(0, base) == 1 for all valid bases.  */
  CHECK_EQ_I(arbint_sizeinbase(x, 2), 1u);
  CHECK_EQ_I(arbint_sizeinbase(x, 10), 1u);
  CHECK_EQ_I(arbint_sizeinbase(x, 16), 1u);
  CHECK_EQ_I(arbint_sizeinbase(x, 62), 1u);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

static void test_sizeinbase_small(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  1 in any base has 1 digit.  */
  CHECK_EQ_I(arbint_set_i32(x, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_sizeinbase(x, 2), 1u);
  CHECK_EQ_I(arbint_sizeinbase(x, 10), 1u);
  CHECK_EQ_I(arbint_sizeinbase(x, 16), 1u);

  /*  255 = 0xFF = 11111111b.  */
  CHECK_EQ_I(arbint_set_i32(x, 255), ARBINT_OK);
  CHECK_EQ_I(arbint_sizeinbase(x, 2), 8u);
  CHECK_EQ_I(arbint_sizeinbase(x, 16), 2u);
  /*  255 in base 10 is 3 digits; may return 3 or 4 (GMP semantics).  */
  CHECK(arbint_sizeinbase(x, 10) >= 3u);
  CHECK(arbint_sizeinbase(x, 10) <= 4u);

  /*  256 = 0x100 = 100000000b.  */
  CHECK_EQ_I(arbint_set_i32(x, 256), ARBINT_OK);
  CHECK_EQ_I(arbint_sizeinbase(x, 2), 9u);
  CHECK_EQ_I(arbint_sizeinbase(x, 16), 3u);

  /*  1000: 4 digits in base 10.  */
  CHECK_EQ_I(arbint_set_i32(x, 1000), ARBINT_OK);
  CHECK(arbint_sizeinbase(x, 10) >= 4u);
  CHECK(arbint_sizeinbase(x, 10) <= 5u);

  /*  Sign is ignored (magnitude only).  */
  CHECK_EQ_I(arbint_set_i32(x, -255), ARBINT_OK);
  CHECK_EQ_I(arbint_sizeinbase(x, 2), 8u);
  CHECK_EQ_I(arbint_sizeinbase(x, 16), 2u);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

static void test_sizeinbase_powers_of_two(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  unsigned i;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  2^k has exactly k+1 bits.  */
  CHECK_EQ_I(arbint_set_i32(x, 1), ARBINT_OK);
  for (i = 0u; i < 20u; ++i) {
    size_t nb = arbint_nbits(x);
    CHECK_EQ_I(nb, (size_t) (i + 1u));
    CHECK_EQ_I(arbint_sizeinbase(x, 2), nb);
    CHECK_EQ_I(arbint_mul_u32(x, x, 2u), ARBINT_OK);
  }

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

static void test_sizeinbase_base8(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  0o777 = 511 = 9 bits -> 3 octal digits.  */
  CHECK_EQ_I(arbint_set_i32(x, 511), ARBINT_OK);
  CHECK_EQ_I(arbint_sizeinbase(x, 8), 3u);

  /*  0o1000 = 512 = 10 bits -> 4 octal digits (ceil(10/3) = 4).  */
  CHECK_EQ_I(arbint_set_i32(x, 512), ARBINT_OK);
  CHECK_EQ_I(arbint_sizeinbase(x, 8), 4u);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

static void test_sizeinbase_invalid(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(x, 42), ARBINT_OK);

  /*  Invalid base returns 1.  */
  CHECK_EQ_I(arbint_sizeinbase(x, 0), 1u);
  CHECK_EQ_I(arbint_sizeinbase(x, 1), 1u);
  CHECK_EQ_I(arbint_sizeinbase(x, 63), 1u);
  CHECK_EQ_I(arbint_sizeinbase(x, -1), 1u);
  CHECK_EQ_I(arbint_sizeinbase(NULL, 10), 1u);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

static void test_sizeinbase_large(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  size_t nb;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  Build a large number via repeated squaring: 2^1024.  */
  CHECK_EQ_I(arbint_set_i32(x, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(x, x), ARBINT_OK); /* 4 = 2^2 */
  CHECK_EQ_I(arbint_sqr(x, x), ARBINT_OK); /* 2^4 */
  CHECK_EQ_I(arbint_sqr(x, x), ARBINT_OK); /* 2^8 */
  CHECK_EQ_I(arbint_sqr(x, x), ARBINT_OK); /* 2^16 */
  CHECK_EQ_I(arbint_sqr(x, x), ARBINT_OK); /* 2^32 */
  CHECK_EQ_I(arbint_sqr(x, x), ARBINT_OK); /* 2^64 */
  CHECK_EQ_I(arbint_sqr(x, x), ARBINT_OK); /* 2^128 */
  CHECK_EQ_I(arbint_sqr(x, x), ARBINT_OK); /* 2^256 */
  CHECK_EQ_I(arbint_sqr(x, x), ARBINT_OK); /* 2^512 */
  CHECK_EQ_I(arbint_sqr(x, x), ARBINT_OK); /* 2^1024 */

  nb = arbint_nbits(x);
  CHECK_EQ_I(nb, 1025u);
  CHECK_EQ_I(arbint_sizeinbase(x, 2), 1025u);
  CHECK_EQ_I(arbint_sizeinbase(x, 16), 257u);

  /*  2^1024 in base 10 has 309 digits; allow +1.  */
  CHECK(arbint_sizeinbase(x, 10) >= 309u);
  CHECK(arbint_sizeinbase(x, 10) <= 310u);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

/* ========== arbint_shl tests ========== */

static void test_shl_zero(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  0 << k == 0 for any k.  */
  CHECK_EQ_I(arbint_shl(r, a, 0u), ARBINT_OK);
  check_i32_value(r, 0);
  CHECK_EQ_I(arbint_shl(r, a, 100u), ARBINT_OK);
  check_i32_value(r, 0);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_shl_by_zero(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  a << 0 == a.  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(r, a, 0u), ARBINT_OK);
  check_i32_value(r, 42);

  CHECK_EQ_I(arbint_set_i32(a, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(r, a, 0u), ARBINT_OK);
  check_i32_value(r, -7);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_shl_small(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  1 << 1 == 2.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(r, a, 1u), ARBINT_OK);
  check_i32_value(r, 2);

  /*  1 << 10 == 1024.  */
  CHECK_EQ_I(arbint_shl(r, a, 10u), ARBINT_OK);
  check_i32_value(r, 1024);

  /*  3 << 4 == 48.  */
  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(r, a, 4u), ARBINT_OK);
  check_i32_value(r, 48);

  /*  -5 << 3 == -40.  */
  CHECK_EQ_I(arbint_set_i32(a, -5), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(r, a, 3u), ARBINT_OK);
  check_i32_value(r, -40);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_shl_whole_limb(void) {
  arbint_ctx_t ctx;
  arbint_t a, r, expected;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  /*  Shift by exactly LIMB_BITS: should insert one zero limb.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(r, a, ARBINT_LIMB_BITS), ARBINT_OK);

  /*  Verify via multiplication: 1 * 2^LIMB_BITS.  */
  CHECK_EQ_I(arbint_set_i32(expected, 1), ARBINT_OK);
  {
    unsigned i;
    for (i = 0u; i < ARBINT_LIMB_BITS; ++i)
      CHECK_EQ_I(arbint_mul_u32(expected, expected, 2u), ARBINT_OK);
  }
  CHECK_EQ_I(arbint_cmp(r, expected), 0);

  /*  Shift by 2 * LIMB_BITS.  */
  CHECK_EQ_I(arbint_set_i32(a, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(r, a, 2u * ARBINT_LIMB_BITS), ARBINT_OK);
  CHECK_EQ_I(arbint_nbits(r), arbint_nbits(a) + 2u * ARBINT_LIMB_BITS);

  arbint_clear(expected);
  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_shl_cross_reference(void) {
  arbint_ctx_t ctx;
  arbint_t a, shl_result, mul_result;
  unsigned k;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(shl_result, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(mul_result, &ctx), ARBINT_OK);

  /*  For various k, verify shl(a, k) == a * 2^k by repeated doubling.  */
  CHECK_EQ_I(arbint_set_i32(a, 12345), ARBINT_OK);

  for (k = 0u; k < 130u && g_failures == 0; ++k) {
    unsigned i;
    CHECK_EQ_I(arbint_shl(shl_result, a, k), ARBINT_OK);

    /*  Compute a * 2^k via repeated mul_u32(x, 2).  */
    CHECK_EQ_I(arbint_set(mul_result, a), ARBINT_OK);
    for (i = 0u; i < k; ++i)
      CHECK_EQ_I(arbint_mul_u32(mul_result, mul_result, 2u), ARBINT_OK);

    CHECK_EQ_I(arbint_cmp(shl_result, mul_result), 0);
  }

  arbint_clear(mul_result);
  arbint_clear(shl_result);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_shl_aliasing(void) {
  arbint_ctx_t ctx;
  arbint_t a, expected;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  /*  In-place: a <<= 5.  */
  CHECK_EQ_I(arbint_set_i32(a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(expected, 3200), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 5u), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(a, expected), 0);

  /*  In-place with whole-limb shift.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, ARBINT_LIMB_BITS + 3u), ARBINT_OK);
  CHECK_EQ_I(arbint_nbits(a), (size_t) (ARBINT_LIMB_BITS + 4u));

  arbint_clear(expected);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ========== arbint_shr tests ========== */

static void test_shr_zero(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  0 >> k == 0 for any k.  */
  CHECK_EQ_I(arbint_shr(r, a, 0u), ARBINT_OK);
  check_i32_value(r, 0);
  CHECK_EQ_I(arbint_shr(r, a, 100u), ARBINT_OK);
  check_i32_value(r, 0);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_shr_by_zero(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_shr(r, a, 0u), ARBINT_OK);
  check_i32_value(r, 42);

  CHECK_EQ_I(arbint_set_i32(a, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_shr(r, a, 0u), ARBINT_OK);
  check_i32_value(r, -7);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_shr_small(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  1024 >> 1 == 512.  */
  CHECK_EQ_I(arbint_set_i32(a, 1024), ARBINT_OK);
  CHECK_EQ_I(arbint_shr(r, a, 1u), ARBINT_OK);
  check_i32_value(r, 512);

  /*  1024 >> 10 == 1.  */
  CHECK_EQ_I(arbint_shr(r, a, 10u), ARBINT_OK);
  check_i32_value(r, 1);

  /*  1024 >> 11 == 0.  */
  CHECK_EQ_I(arbint_shr(r, a, 11u), ARBINT_OK);
  check_i32_value(r, 0);

  /*  48 >> 4 == 3.  */
  CHECK_EQ_I(arbint_set_i32(a, 48), ARBINT_OK);
  CHECK_EQ_I(arbint_shr(r, a, 4u), ARBINT_OK);
  check_i32_value(r, 3);

  /*  -40 >> 3 == -5 (truncation toward zero).  */
  CHECK_EQ_I(arbint_set_i32(a, -40), ARBINT_OK);
  CHECK_EQ_I(arbint_shr(r, a, 3u), ARBINT_OK);
  check_i32_value(r, -5);

  /*  -7 >> 1 == -3 (truncation toward zero: -7/2 = -3.5 -> -3).  */
  CHECK_EQ_I(arbint_set_i32(a, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_shr(r, a, 1u), ARBINT_OK);
  check_i32_value(r, -3);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_shr_shift_away(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Shifting right by more bits than the number has -> zero.  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_shr(r, a, 100u), ARBINT_OK);
  check_i32_value(r, 0);

  /*  Negative number shifted away -> zero.  */
  CHECK_EQ_I(arbint_set_i32(a, -42), ARBINT_OK);
  CHECK_EQ_I(arbint_shr(r, a, 100u), ARBINT_OK);
  check_i32_value(r, 0);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_shr_whole_limb(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Build a 2-limb number, shift right by exactly LIMB_BITS.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, ARBINT_LIMB_BITS + 5u), ARBINT_OK);
  CHECK_EQ_I(arbint_add_i32(a, a, 7), ARBINT_OK);
  /*  a = 2^(LIMB_BITS+5) + 7.  */

  CHECK_EQ_I(arbint_shr(r, a, ARBINT_LIMB_BITS), ARBINT_OK);
  /*  r = (2^(LIMB_BITS+5) + 7) >> LIMB_BITS = 2^5 = 32.  */
  check_i32_value(r, 32);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_shr_cross_reference(void) {
  arbint_ctx_t ctx;
  arbint_t a, shr_result, div_result;
  unsigned k;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(shr_result, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(div_result, &ctx), ARBINT_OK);

  /*  Verify shr(a, k) == tdiv_q_u32(a, 2^k) for small k.  */
  CHECK_EQ_I(arbint_set_i32(a, 999999), ARBINT_OK);

  for (k = 0u; k < 20u && g_failures == 0; ++k) {
    uint32_t divisor = (uint32_t) 1u << k;
    CHECK_EQ_I(arbint_shr(shr_result, a, k), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_q_u32(div_result, a, divisor), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(shr_result, div_result), 0);
  }

  arbint_clear(div_result);
  arbint_clear(shr_result);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_shr_aliasing(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);

  /*  In-place: a >>= 3.  */
  CHECK_EQ_I(arbint_set_i32(a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_shr(a, a, 3u), ARBINT_OK);
  check_i32_value(a, 12);

  /*  In-place shift to zero.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shr(a, a, 5u), ARBINT_OK);
  check_i32_value(a, 0);

  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ========== shl/shr roundtrip ========== */

static void test_shl_shr_roundtrip(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, c;
  unsigned k;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_set_i32(a, 12345), ARBINT_OK);

  for (k = 0u; k < 200u && g_failures == 0; ++k) {
    /*  (a << k) >> k == a.  */
    CHECK_EQ_I(arbint_shl(b, a, k), ARBINT_OK);
    CHECK_EQ_I(arbint_shr(c, b, k), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(c, a), 0);
  }

  /*  Also test with negative.  */
  CHECK_EQ_I(arbint_set_i32(a, -6789), ARBINT_OK);
  for (k = 0u; k < 200u && g_failures == 0; ++k) {
    CHECK_EQ_I(arbint_shl(b, a, k), ARBINT_OK);
    CHECK_EQ_I(arbint_shr(c, b, k), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(c, a), 0);
  }

  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ========== Large multi-limb shift tests ========== */

static void test_shl_large(void) {
  arbint_ctx_t ctx;
  arbint_t a, r, mul_r;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(mul_r, &ctx), ARBINT_OK);

  /*  Build a multi-limb number via squaring: 3^128.  */
  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  {
    unsigned i;
    for (i = 0u; i < 7u; ++i)
      CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
  }
  /*  a = 3^128, a multi-limb number.  */

  /*  Shift left by 200 bits, verify via separate mul_u32 chain.  */
  CHECK_EQ_I(arbint_shl(r, a, 200u), ARBINT_OK);

  CHECK_EQ_I(arbint_set(mul_r, a), ARBINT_OK);
  {
    unsigned i;
    for (i = 0u; i < 200u; ++i)
      CHECK_EQ_I(arbint_mul_u32(mul_r, mul_r, 2u), ARBINT_OK);
  }
  CHECK_EQ_I(arbint_cmp(r, mul_r), 0);

  arbint_clear(mul_r);
  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_shr_large(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;
  size_t nb_before;
  size_t nb_after;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Build 3^128 then shift it right by 100 bits.  */
  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  {
    unsigned i;
    for (i = 0u; i < 7u; ++i)
      CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
  }

  nb_before = arbint_nbits(a);
  CHECK(nb_before > 100u);

  CHECK_EQ_I(arbint_shr(r, a, 100u), ARBINT_OK);
  nb_after = arbint_nbits(r);
  CHECK_EQ_I(nb_after, nb_before - 100u);

  /*  Roundtrip: shl back by 100 should give a value with same high bits.  */
  {
    arbint_t roundtrip;
    CHECK_EQ_I(arbint_init(roundtrip, &ctx), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(roundtrip, r, 100u), ARBINT_OK);
    /*  roundtrip <= a, and a - roundtrip < 2^100.  */
    CHECK(arbint_cmp(roundtrip, a) <= 0);
    arbint_clear(roundtrip);
  }

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_shl_null(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_shl(NULL, a, 1u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_shl(a, NULL, 1u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_shr(NULL, a, 1u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_shr(a, NULL, 1u), ARBINT_EINVAL);

  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

int main(void) {
  /*  sizeinbase tests.  */
  test_sizeinbase_zero();
  test_sizeinbase_small();
  test_sizeinbase_powers_of_two();
  test_sizeinbase_base8();
  test_sizeinbase_invalid();
  test_sizeinbase_large();

  /*  shl tests.  */
  test_shl_zero();
  test_shl_by_zero();
  test_shl_small();
  test_shl_whole_limb();
  test_shl_cross_reference();
  test_shl_aliasing();
  test_shl_null();
  test_shl_large();

  /*  shr tests.  */
  test_shr_zero();
  test_shr_by_zero();
  test_shr_small();
  test_shr_shift_away();
  test_shr_whole_limb();
  test_shr_cross_reference();
  test_shr_aliasing();
  test_shr_large();

  /*  roundtrip.  */
  test_shl_shr_roundtrip();

  ARBINT_TEST_FINISH("test_shift");
}
