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

ARBINT_TEST_DECLARE_FAILURES();

/* ========== testbit ========== */

static void test_testbit_zero(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  int out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  All bits of zero should be 0.  */
  CHECK_EQ_I(arbint_testbit(x, 0u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);
  CHECK_EQ_I(arbint_testbit(x, 63u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);
  CHECK_EQ_I(arbint_testbit(x, 1000u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

static void test_testbit_positive(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  int out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  x = 0b1010 = 10  */
  CHECK_EQ_I(arbint_set_i32(x, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_testbit(x, 0u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);
  CHECK_EQ_I(arbint_testbit(x, 1u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_testbit(x, 2u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);
  CHECK_EQ_I(arbint_testbit(x, 3u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  /*  Beyond magnitude: should be 0.  */
  CHECK_EQ_I(arbint_testbit(x, 4u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);
  CHECK_EQ_I(arbint_testbit(x, 100u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

static void test_testbit_negative(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  int out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  x = -1: TC is ...1111 (all bits set).  */
  CHECK_EQ_I(arbint_set_i32(x, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_testbit(x, 0u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_testbit(x, 63u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_testbit(x, 1000u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  x = -2: TC is ...11111110.  */
  CHECK_EQ_I(arbint_set_i32(x, -2), ARBINT_OK);
  CHECK_EQ_I(arbint_testbit(x, 0u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);
  CHECK_EQ_I(arbint_testbit(x, 1u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_testbit(x, 100u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  x = -10: TC is ...11110110.  bit0=0, bit1=1, bit2=1, bit3=0.  */
  CHECK_EQ_I(arbint_set_i32(x, -10), ARBINT_OK);
  CHECK_EQ_I(arbint_testbit(x, 0u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);
  CHECK_EQ_I(arbint_testbit(x, 1u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_testbit(x, 2u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_testbit(x, 3u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

static void test_testbit_null(void) {
  CHECK_EQ_I(arbint_testbit(NULL, 0u, NULL), ARBINT_EINVAL);
}

/* ========== setbit ========== */

static void test_setbit_positive(void) {
  arbint_ctx_t ctx;
  arbint_t x;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  Set bit 0 of zero -> 1.  */
  CHECK_EQ_I(arbint_setbit(x, 0u), ARBINT_OK);
  check_i32_value(x, 1);

  /*  Set bit 3 of 1 -> 9 (0b1001).  */
  CHECK_EQ_I(arbint_setbit(x, 3u), ARBINT_OK);
  check_i32_value(x, 9);

  /*  Set already-set bit: no change.  */
  CHECK_EQ_I(arbint_setbit(x, 0u), ARBINT_OK);
  check_i32_value(x, 9);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

static void test_setbit_negative(void) {
  arbint_ctx_t ctx;
  arbint_t x;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  x = -10 (TC: ...11110110).  Set bit 3 -> TC: ...11111110 = -2.  */
  CHECK_EQ_I(arbint_set_i32(x, -10), ARBINT_OK);
  CHECK_EQ_I(arbint_setbit(x, 3u), ARBINT_OK);
  check_i32_value(x, -2);

  /*  x = -1 (TC: all 1s).  Setting any bit is no-op.  */
  CHECK_EQ_I(arbint_set_i32(x, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_setbit(x, 5u), ARBINT_OK);
  check_i32_value(x, -1);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

static void test_setbit_extend(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  size_t nb;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  Set a high bit to force expansion beyond 1 limb.  */
  CHECK_EQ_I(arbint_setbit(x, 128u), ARBINT_OK);
  nb = arbint_nbits(x);
  CHECK_EQ_I(nb, 129u);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

/* ========== clrbit ========== */

static void test_clrbit_positive(void) {
  arbint_ctx_t ctx;
  arbint_t x;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  x = 10 (0b1010).  Clear bit 1 -> 8 (0b1000).  */
  CHECK_EQ_I(arbint_set_i32(x, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_clrbit(x, 1u), ARBINT_OK);
  check_i32_value(x, 8);

  /*  Clear bit 3 -> 0 (0b0000).  */
  CHECK_EQ_I(arbint_clrbit(x, 3u), ARBINT_OK);
  check_i32_value(x, 0);

  /*  Clearing bit of zero: no-op.  */
  CHECK_EQ_I(arbint_clrbit(x, 0u), ARBINT_OK);
  check_i32_value(x, 0);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

static void test_clrbit_negative(void) {
  arbint_ctx_t ctx;
  arbint_t x;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  x = -2 (TC: ...11111110).  Clear bit 1 -> TC: ...11111100 = -4.  */
  CHECK_EQ_I(arbint_set_i32(x, -2), ARBINT_OK);
  CHECK_EQ_I(arbint_clrbit(x, 1u), ARBINT_OK);
  check_i32_value(x, -4);

  /*  x = -1 (TC: all 1s).  Clear bit 0 -> TC: ...11111110 = -2.  */
  CHECK_EQ_I(arbint_set_i32(x, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_clrbit(x, 0u), ARBINT_OK);
  check_i32_value(x, -2);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

/* ========== ctz ========== */

static void test_ctz_basic(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  size_t out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  ctz(0) = EDOM.  */
  CHECK_EQ_I(arbint_ctz(x, &out), ARBINT_EDOM);

  /*  ctz(1) = 0.  */
  CHECK_EQ_I(arbint_set_i32(x, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_ctz(x, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0u);

  /*  ctz(8) = 3.  */
  CHECK_EQ_I(arbint_set_i32(x, 8), ARBINT_OK);
  CHECK_EQ_I(arbint_ctz(x, &out), ARBINT_OK);
  CHECK_EQ_I(out, 3u);

  /*  ctz(-8) = 3.  */
  CHECK_EQ_I(arbint_set_i32(x, -8), ARBINT_OK);
  CHECK_EQ_I(arbint_ctz(x, &out), ARBINT_OK);
  CHECK_EQ_I(out, 3u);

  /*  ctz(2^10 = 1024) = 10.  */
  CHECK_EQ_I(arbint_set_i32(x, 1024), ARBINT_OK);
  CHECK_EQ_I(arbint_ctz(x, &out), ARBINT_OK);
  CHECK_EQ_I(out, 10u);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

static void test_ctz_multilimb(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  size_t out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  Build a value with trailing zeros spanning a limb boundary.
      x = 1 << 128 (bit 128 set, all lower bits zero).  */
  CHECK_EQ_I(arbint_set_i32(x, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(x, x, 128u), ARBINT_OK);
  CHECK_EQ_I(arbint_ctz(x, &out), ARBINT_OK);
  CHECK_EQ_I(out, 128u);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

/* ========== clz ========== */

static void test_clz_basic(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  size_t out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  clz(0) = EDOM.  */
  CHECK_EQ_I(arbint_clz(x, &out), ARBINT_EDOM);

  /*  clz(1): top limb = 1, clz = LIMB_BITS - 1.  */
  CHECK_EQ_I(arbint_set_i32(x, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_clz(x, &out), ARBINT_OK);
  CHECK_EQ_I(out, ARBINT_LIMB_BITS - 1u);

  /*  clz(-1) = 0 (negative TC has MSB set).  */
  CHECK_EQ_I(arbint_set_i32(x, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_clz(x, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0u);

  /*  clz(-100) = 0.  */
  CHECK_EQ_I(arbint_set_i32(x, -100), ARBINT_OK);
  CHECK_EQ_I(arbint_clz(x, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0u);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

/* ========== popcount ========== */

static void test_popcount_basic(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  size_t out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  popcount(0) = 0.  */
  CHECK_EQ_I(arbint_popcount(x, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0u);

  /*  popcount(1) = 1.  */
  CHECK_EQ_I(arbint_set_i32(x, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_popcount(x, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1u);

  /*  popcount(7) = 3.  */
  CHECK_EQ_I(arbint_set_i32(x, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_popcount(x, &out), ARBINT_OK);
  CHECK_EQ_I(out, 3u);

  /*  popcount(255) = 8.  */
  CHECK_EQ_I(arbint_set_i32(x, 255), ARBINT_OK);
  CHECK_EQ_I(arbint_popcount(x, &out), ARBINT_OK);
  CHECK_EQ_I(out, 8u);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

static void test_popcount_negative(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  size_t out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  popcount(-1) within W(-1) = 1 limb.
      TC of -1 in one limb = 0xFFFFFFFFFFFFFFFF -> LIMB_BITS bits set.  */
  CHECK_EQ_I(arbint_set_i32(x, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_popcount(x, &out), ARBINT_OK);
  CHECK_EQ_I(out, ARBINT_LIMB_BITS);

  /*  popcount(-2) within W(-2) = 1 limb.
      TC of -2 = 0xFFFFFFFFFFFFFFFE -> LIMB_BITS - 1 bits.  */
  CHECK_EQ_I(arbint_set_i32(x, -2), ARBINT_OK);
  CHECK_EQ_I(arbint_popcount(x, &out), ARBINT_OK);
  CHECK_EQ_I(out, ARBINT_LIMB_BITS - 1u);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

static void test_popcount_large(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  size_t out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  Build a large number via repeated squaring: 3^256.
      popcount should be consistent between two calls.  */
  CHECK_EQ_I(arbint_set_i32(x, 3), ARBINT_OK);
  {
    int i;
    for (i = 0; i < 8; ++i)
      CHECK_EQ_I(arbint_sqr(x, x), ARBINT_OK);
  }
  CHECK_EQ_I(arbint_popcount(x, &out), ARBINT_OK);
  CHECK(out > 0u);

  size_t nb = arbint_nbits(x);
  CHECK(out <= nb);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

/* ========== hammingdist ========== */

static void test_hammingdist_basic(void) {
  arbint_ctx_t ctx;
  arbint_t a, b;
  size_t out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);

  /*  hamming(0, 0) = 0.  */
  CHECK_EQ_I(arbint_hammingdist(a, b, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0u);

  /*  hamming(0, 7) = 3.  */
  CHECK_EQ_I(arbint_set_i32(b, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_hammingdist(a, b, &out), ARBINT_OK);
  CHECK_EQ_I(out, 3u);

  /*  hamming(x, x) = 0.  */
  CHECK_EQ_I(arbint_set_i32(a, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_hammingdist(a, b, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0u);

  /*  hamming(5, 3) = hamming(0b101, 0b011) = 2 (bits 0 and 2 differ).  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_hammingdist(a, b, &out), ARBINT_OK);
  CHECK_EQ_I(out, 2u);

  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_hammingdist_negative(void) {
  arbint_ctx_t ctx;
  arbint_t a, b;
  size_t out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);

  /*  hamming(-1, -1) = 0.  */
  CHECK_EQ_I(arbint_set_i32(a, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_hammingdist(a, b, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0u);

  /*  hamming(-1, -2) within max(W(-1), W(-2)) = 1 limb.
      TC(-1) = 0xFFFF...FFFF, TC(-2) = 0xFFFF...FFFE.
      XOR = 0x0000...0001 -> popcount = 1.  */
  CHECK_EQ_I(arbint_set_i32(b, -2), ARBINT_OK);
  CHECK_EQ_I(arbint_hammingdist(a, b, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1u);

  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_hammingdist_mixed_sign(void) {
  arbint_ctx_t ctx;
  arbint_t a, b;
  size_t out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);

  /*  hamming(1, -1) within max(W(1), W(-1)) = 1 limb.
      TC(1) = 0x0000...0001, TC(-1) = 0xFFFF...FFFF.
      XOR = 0xFFFF...FFFE -> popcount = LIMB_BITS - 1.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_hammingdist(a, b, &out), ARBINT_OK);
  CHECK_EQ_I(out, ARBINT_LIMB_BITS - 1u);

  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ========== not ========== */

static void test_not_basic(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  not(0) = -1.  */
  CHECK_EQ_I(arbint_not(r, a), ARBINT_OK);
  check_i32_value(r, -1);

  /*  not(-1) = 0.  */
  CHECK_EQ_I(arbint_set_i32(a, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_not(r, a), ARBINT_OK);
  check_i32_value(r, 0);

  /*  not(5) = -6.  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_not(r, a), ARBINT_OK);
  check_i32_value(r, -6);

  /*  not(-6) = 5.  */
  CHECK_EQ_I(arbint_set_i32(a, -6), ARBINT_OK);
  CHECK_EQ_I(arbint_not(r, a), ARBINT_OK);
  check_i32_value(r, 5);

  /*  not(not(x)) = x.  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_not(r, a), ARBINT_OK);
  CHECK_EQ_I(arbint_not(r, r), ARBINT_OK);
  check_i32_value(r, 42);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_not_aliasing(void) {
  arbint_ctx_t ctx;
  arbint_t x;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  In-place not.  */
  CHECK_EQ_I(arbint_set_i32(x, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_not(x, x), ARBINT_OK);
  check_i32_value(x, -101);

  CHECK_EQ_I(arbint_not(x, x), ARBINT_OK);
  check_i32_value(x, 100);

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

/* ========== or ========== */

static void test_or_basic(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  0 | 0 = 0.  */
  CHECK_EQ_I(arbint_or(r, a, b), ARBINT_OK);
  check_i32_value(r, 0);

  /*  5 | 3 = 7.  (0b101 | 0b011 = 0b111)  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_or(r, a, b), ARBINT_OK);
  check_i32_value(r, 7);

  /*  x | 0 = x.  */
  CHECK_EQ_I(arbint_set_i32(b, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_or(r, a, b), ARBINT_OK);
  check_i32_value(r, 5);

  /*  x | x = x.  */
  CHECK_EQ_I(arbint_or(r, a, a), ARBINT_OK);
  check_i32_value(r, 5);

  arbint_clear(r);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_or_negative(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  -1 | x = -1 for any x (TC all 1s).  */
  CHECK_EQ_I(arbint_set_i32(a, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_or(r, a, b), ARBINT_OK);
  check_i32_value(r, -1);

  /*  -4 | -2 = -2.
      TC(-4) = ...11100, TC(-2) = ...11110.
      OR = ...11110 = -2.  */
  CHECK_EQ_I(arbint_set_i32(a, -4), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -2), ARBINT_OK);
  CHECK_EQ_I(arbint_or(r, a, b), ARBINT_OK);
  check_i32_value(r, -2);

  arbint_clear(r);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ========== and ========== */

static void test_and_basic(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  5 & 3 = 1.  (0b101 & 0b011 = 0b001)  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_and(r, a, b), ARBINT_OK);
  check_i32_value(r, 1);

  /*  x & 0 = 0.  */
  CHECK_EQ_I(arbint_set_i32(b, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_and(r, a, b), ARBINT_OK);
  check_i32_value(r, 0);

  /*  x & -1 = x (TC -1 is all 1s).  */
  CHECK_EQ_I(arbint_set_i32(b, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_and(r, a, b), ARBINT_OK);
  check_i32_value(r, 5);

  /*  x & x = x.  */
  CHECK_EQ_I(arbint_and(r, a, a), ARBINT_OK);
  check_i32_value(r, 5);

  arbint_clear(r);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_and_negative(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  -4 & -2 = -4.
      TC(-4) = ...11100, TC(-2) = ...11110.
      AND = ...11100 = -4.  */
  CHECK_EQ_I(arbint_set_i32(a, -4), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -2), ARBINT_OK);
  CHECK_EQ_I(arbint_and(r, a, b), ARBINT_OK);
  check_i32_value(r, -4);

  /*  5 & -4 = 4.
      TC(5) = ...000101, TC(-4) = ...111100.
      AND = ...000100 = 4.  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -4), ARBINT_OK);
  CHECK_EQ_I(arbint_and(r, a, b), ARBINT_OK);
  check_i32_value(r, 4);

  arbint_clear(r);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ========== xor ========== */

static void test_xor_basic(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  5 ^ 3 = 6.  (0b101 ^ 0b011 = 0b110)  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_xor(r, a, b), ARBINT_OK);
  check_i32_value(r, 6);

  /*  x ^ 0 = x.  */
  CHECK_EQ_I(arbint_set_i32(b, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_xor(r, a, b), ARBINT_OK);
  check_i32_value(r, 5);

  /*  x ^ x = 0.  */
  CHECK_EQ_I(arbint_xor(r, a, a), ARBINT_OK);
  check_i32_value(r, 0);

  arbint_clear(r);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_xor_negative(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  x ^ -1 = not(x).  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_xor(r, a, b), ARBINT_OK);
  check_i32_value(r, -6);

  /*  -4 ^ -2.
      TC(-4) = ...11100, TC(-2) = ...11110.
      XOR = ...00010 = 2.  */
  CHECK_EQ_I(arbint_set_i32(a, -4), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -2), ARBINT_OK);
  CHECK_EQ_I(arbint_xor(r, a, b), ARBINT_OK);
  check_i32_value(r, 2);

  arbint_clear(r);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ========== aliasing for binary ops ========== */

static void test_bitwise_aliasing(void) {
  arbint_ctx_t ctx;
  arbint_t a, b;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);

  /*  or: rop == a.  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_or(a, a, b), ARBINT_OK);
  check_i32_value(a, 7);

  /*  and: rop == b.  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_and(b, a, b), ARBINT_OK);
  check_i32_value(b, 1);

  /*  xor: rop == a == b.  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_xor(a, a, a), ARBINT_OK);
  check_i32_value(a, 0);

  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ========== identities ========== */

static void test_bitwise_identities(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r1, r2;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r2, &ctx), ARBINT_OK);

  /*  De Morgan's: not(a & b) = not(a) | not(b).  */
  CHECK_EQ_I(arbint_set_i32(a, 123), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 456), ARBINT_OK);

  CHECK_EQ_I(arbint_and(r1, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_not(r1, r1), ARBINT_OK);

  CHECK_EQ_I(arbint_not(r2, a), ARBINT_OK);
  {
    arbint_t tmp;
    CHECK_EQ_I(arbint_init(tmp, &ctx), ARBINT_OK);
    CHECK_EQ_I(arbint_not(tmp, b), ARBINT_OK);
    CHECK_EQ_I(arbint_or(r2, r2, tmp), ARBINT_OK);
    arbint_clear(tmp);
  }
  CHECK_EQ_I(arbint_cmp(r1, r2), 0);

  /*  De Morgan's: not(a | b) = not(a) & not(b).  */
  CHECK_EQ_I(arbint_or(r1, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_not(r1, r1), ARBINT_OK);

  CHECK_EQ_I(arbint_not(r2, a), ARBINT_OK);
  {
    arbint_t tmp;
    CHECK_EQ_I(arbint_init(tmp, &ctx), ARBINT_OK);
    CHECK_EQ_I(arbint_not(tmp, b), ARBINT_OK);
    CHECK_EQ_I(arbint_and(r2, r2, tmp), ARBINT_OK);
    arbint_clear(tmp);
  }
  CHECK_EQ_I(arbint_cmp(r1, r2), 0);

  /*  Verify with negative operands.  */
  CHECK_EQ_I(arbint_set_i32(a, -123), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 456), ARBINT_OK);

  CHECK_EQ_I(arbint_and(r1, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_not(r1, r1), ARBINT_OK);

  CHECK_EQ_I(arbint_not(r2, a), ARBINT_OK);
  {
    arbint_t tmp;
    CHECK_EQ_I(arbint_init(tmp, &ctx), ARBINT_OK);
    CHECK_EQ_I(arbint_not(tmp, b), ARBINT_OK);
    CHECK_EQ_I(arbint_or(r2, r2, tmp), ARBINT_OK);
    arbint_clear(tmp);
  }
  CHECK_EQ_I(arbint_cmp(r1, r2), 0);

  arbint_clear(r2);
  arbint_clear(r1);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ========== multi-limb binary ops ========== */

static void test_bitwise_multilimb(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r, r2;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r2, &ctx), ARBINT_OK);

  /*  Build multi-limb values via shifting.
      a = 1 << 128, b = 1 << 64.
      or(a, b) should have both bits set.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 128u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(b, b, 64u), ARBINT_OK);

  CHECK_EQ_I(arbint_or(r, a, b), ARBINT_OK);
  {
    size_t nb = arbint_nbits(r);
    CHECK_EQ_I(nb, 129u);
  }
  {
    int bit64;
    int bit128;
    CHECK_EQ_I(arbint_testbit(r, 64u, &bit64), ARBINT_OK);
    CHECK_EQ_I(arbint_testbit(r, 128u, &bit128), ARBINT_OK);
    CHECK_EQ_I(bit64, 1);
    CHECK_EQ_I(bit128, 1);
  }

  /*  and(a, b) = 0 (no overlapping bits).  */
  CHECK_EQ_I(arbint_and(r, a, b), ARBINT_OK);
  check_i32_value(r, 0);

  /*  xor(a | b, a) = b.  */
  CHECK_EQ_I(arbint_or(r, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_xor(r, r, a), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, b), 0);

  /*  Multi-limb negative: not(a) where a has 3 limbs.  */
  CHECK_EQ_I(arbint_not(r, a), ARBINT_OK);
  CHECK_EQ_I(arbint_not(r2, r), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r2, a), 0);

  arbint_clear(r2);
  arbint_clear(r);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ========== stochastic ========== */

static void test_bitwise_stochastic(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r1, r2, r3;
  enum { ITERS = 100 };
  int i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r2, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r3, &ctx), ARBINT_OK);

  /*  For i = 0..ITERS, use a = i*37 - 500, b = i*53 - 300.
      Verify: a ^ b ^ b = a (XOR is self-inverse).  */
  for (i = 0; i < ITERS && g_failures == 0; ++i) {
    int32_t av = (int32_t) (i * 37 - 500);
    int32_t bv = (int32_t) (i * 53 - 300);

    CHECK_EQ_I(arbint_set_i32(a, av), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, bv), ARBINT_OK);

    CHECK_EQ_I(arbint_xor(r1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_xor(r1, r1, b), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(r1, a), 0);

    /*  Verify: (a | b) & a = a & (a | b) = a (absorption-like via OR).
        Actually: a & (a | b) = a.  */
    CHECK_EQ_I(arbint_or(r1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_and(r2, a, r1), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(r2, a), 0);

    /*  Verify: (a & b) | a = a.  */
    CHECK_EQ_I(arbint_and(r1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_or(r2, r1, a), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(r2, a), 0);
  }

  arbint_clear(r3);
  arbint_clear(r2);
  arbint_clear(r1);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ========== null pointer checks ========== */

static void test_null_checks(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  size_t out;
  int iout;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_testbit(NULL, 0u, &iout), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_testbit(a, 0u, NULL), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_setbit(NULL, 0u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_clrbit(NULL, 0u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_ctz(NULL, &out), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_ctz(a, NULL), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_clz(NULL, &out), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_clz(a, NULL), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_popcount(NULL, &out), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_popcount(a, NULL), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_hammingdist(NULL, a, &out), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_hammingdist(a, NULL, &out), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_hammingdist(a, a, NULL), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_or(NULL, a, a), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_or(a, NULL, a), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_or(a, a, NULL), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_and(NULL, a, a), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_xor(NULL, a, a), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_not(NULL, a), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_not(a, NULL), ARBINT_EINVAL);

  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ========== i32/u32 immediate variants ========== */

static void test_and_u32(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  5 & 3u = 1.  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_and_u32(r, a, 3u), ARBINT_OK);
  check_i32_value(r, 1);

  /*  x & 0u = 0.  */
  CHECK_EQ_I(arbint_and_u32(r, a, 0u), ARBINT_OK);
  check_i32_value(r, 0);

  /*  0 & 5u = 0.  */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_and_u32(r, a, 5u), ARBINT_OK);
  check_i32_value(r, 0);

  /*  -4 & 7u.  TC(-4) = ...11100.  ...11100 & 0b111 = 0b100 = 4.  */
  CHECK_EQ_I(arbint_set_i32(a, -4), ARBINT_OK);
  CHECK_EQ_I(arbint_and_u32(r, a, 7u), ARBINT_OK);
  check_i32_value(r, 4);

  /*  -1 & 0xFFu = 255.  TC(-1) = all-ones.  */
  CHECK_EQ_I(arbint_set_i32(a, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_and_u32(r, a, 0xFFu), ARBINT_OK);
  check_u32_value(r, 255u);

  /*  x & UINT32_MAX.  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_and_u32(r, a, 0xFFFFFFFFu), ARBINT_OK);
  check_i32_value(r, 42);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_and_i32(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Positive b delegates to u32.  */
  CHECK_EQ_I(arbint_set_i32(a, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_and_i32(r, a, 3), ARBINT_OK);
  check_i32_value(r, 3);

  /*  (+a) & (-b): TC(-1) = all-ones, so a & -1 = a.  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_and_i32(r, a, -1), ARBINT_OK);
  check_i32_value(r, 42);

  /*  (+a) & (-2): TC(-2) = ...FFFE. Clears bit 0 only.  */
  CHECK_EQ_I(arbint_set_i32(a, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_and_i32(r, a, -2), ARBINT_OK);
  check_i32_value(r, 6);

  /*  (-a) & (-b): (-4) & (-2) = -4.  */
  CHECK_EQ_I(arbint_set_i32(a, -4), ARBINT_OK);
  CHECK_EQ_I(arbint_and_i32(r, a, -2), ARBINT_OK);
  check_i32_value(r, -4);

  /*  0 & (-1) = 0.  */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_and_i32(r, a, -1), ARBINT_OK);
  check_i32_value(r, 0);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_or_u32(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  5 | 3u = 7.  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_or_u32(r, a, 3u), ARBINT_OK);
  check_i32_value(r, 7);

  /*  0 | 5u = 5.  */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_or_u32(r, a, 5u), ARBINT_OK);
  check_i32_value(r, 5);

  /*  x | 0u = x.  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_or_u32(r, a, 0u), ARBINT_OK);
  check_i32_value(r, 42);

  /*  (-4) | 3u.  TC(-4) = ...11100.  ...11100 | 0b011 = ...11111 = -1.  */
  CHECK_EQ_I(arbint_set_i32(a, -4), ARBINT_OK);
  CHECK_EQ_I(arbint_or_u32(r, a, 3u), ARBINT_OK);
  check_i32_value(r, -1);

  /*  (-8) | 1u.  TC(-8) = ...11000.  ...11000 | 1 = ...11001 = -7.  */
  CHECK_EQ_I(arbint_set_i32(a, -8), ARBINT_OK);
  CHECK_EQ_I(arbint_or_u32(r, a, 1u), ARBINT_OK);
  check_i32_value(r, -7);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_or_i32(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Positive b delegates to u32.  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_or_i32(r, a, 3), ARBINT_OK);
  check_i32_value(r, 7);

  /*  (+a) | (-1) = -1.  TC(-1) = all-ones.  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_or_i32(r, a, -1), ARBINT_OK);
  check_i32_value(r, -1);

  /*  0 | (-2) = -2.  */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_or_i32(r, a, -2), ARBINT_OK);
  check_i32_value(r, -2);

  /*  (-4) | (-2) = -2.  */
  CHECK_EQ_I(arbint_set_i32(a, -4), ARBINT_OK);
  CHECK_EQ_I(arbint_or_i32(r, a, -2), ARBINT_OK);
  check_i32_value(r, -2);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_xor_u32(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  5 ^ 3u = 6.  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_xor_u32(r, a, 3u), ARBINT_OK);
  check_i32_value(r, 6);

  /*  0 ^ 5u = 5.  */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_xor_u32(r, a, 5u), ARBINT_OK);
  check_i32_value(r, 5);

  /*  x ^ 0u = x.  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_xor_u32(r, a, 0u), ARBINT_OK);
  check_i32_value(r, 42);

  /*  (-4) ^ 3u.  TC(-4) = ...11100.  ...11100 ^ 011 = ...11111 = -1.  */
  CHECK_EQ_I(arbint_set_i32(a, -4), ARBINT_OK);
  CHECK_EQ_I(arbint_xor_u32(r, a, 3u), ARBINT_OK);
  check_i32_value(r, -1);

  /*  (-1) ^ 1u = -2.  TC(-1) = all-ones, XOR 1 -> ...FFFE.  */
  CHECK_EQ_I(arbint_set_i32(a, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_xor_u32(r, a, 1u), ARBINT_OK);
  check_i32_value(r, -2);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_xor_i32(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  x ^ -1 = not(x).  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_xor_i32(r, a, -1), ARBINT_OK);
  check_i32_value(r, -6);

  /*  0 ^ -1 = -1.  */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_xor_i32(r, a, -1), ARBINT_OK);
  check_i32_value(r, -1);

  /*  (-4) ^ (-2) = 2.  */
  CHECK_EQ_I(arbint_set_i32(a, -4), ARBINT_OK);
  CHECK_EQ_I(arbint_xor_i32(r, a, -2), ARBINT_OK);
  check_i32_value(r, 2);

  /*  (-1) ^ (-1) = 0.  */
  CHECK_EQ_I(arbint_set_i32(a, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_xor_i32(r, a, -1), ARBINT_OK);
  check_i32_value(r, 0);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ========== i32/u32 cross-validation against arbint_t ops ========== */

static void test_imm_cross_validation(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r1, r2;
  enum { ITERS = 100 };
  int i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r2, &ctx), ARBINT_OK);

  /*  For each i, verify and/or/xor_i32 matches and/or/xor with arbint_t.  */
  for (i = 0; i < ITERS && g_failures == 0; ++i) {
    int32_t av = (int32_t) (i * 37 - 500);
    int32_t bv = (int32_t) (i * 53 - 300);

    CHECK_EQ_I(arbint_set_i32(a, av), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, bv), ARBINT_OK);

    CHECK_EQ_I(arbint_and(r1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_and_i32(r2, a, bv), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(r1, r2), 0);

    CHECK_EQ_I(arbint_or(r1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_or_i32(r2, a, bv), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(r1, r2), 0);

    CHECK_EQ_I(arbint_xor(r1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_xor_i32(r2, a, bv), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(r1, r2), 0);
  }

  arbint_clear(r2);
  arbint_clear(r1);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ========== popcount complement identity ========== */

static void test_popcount_complement(void) {
  arbint_ctx_t ctx;
  arbint_t x, nx;
  size_t pc_x;
  size_t pc_nx;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(nx, &ctx), ARBINT_OK);

  /*  Simple case: x = 42, not(42) = -43.  Both single-limb.  */
  CHECK_EQ_I(arbint_set_i32(x, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_not(nx, x), ARBINT_OK);
  CHECK_EQ_I(arbint_popcount(x, &pc_x), ARBINT_OK);
  CHECK_EQ_I(arbint_popcount(nx, &pc_nx), ARBINT_OK);
  CHECK_EQ_I(pc_x + pc_nx, ARBINT_LIMB_BITS);

  /*  Multi-limb: x = 3^128 (large positive).  */
  CHECK_EQ_I(arbint_set_i32(x, 3), ARBINT_OK);
  int i;
  for (i = 0; i < 7; ++i)
    CHECK_EQ_I(arbint_sqr(x, x), ARBINT_OK);
  CHECK_EQ_I(arbint_not(nx, x), ARBINT_OK);
  CHECK_EQ_I(arbint_popcount(x, &pc_x), ARBINT_OK);
  CHECK_EQ_I(arbint_popcount(nx, &pc_nx), ARBINT_OK);
  /*  Both have same number of used limbs (mag_inc from arbint_not
      only grows if all limbs are max, which is extremely unlikely).  */
  size_t nb_x = arbint_nbits(x);
  size_t nb_nx = arbint_nbits(nx);
  size_t max_nb = (nb_x > nb_nx) ? nb_x : nb_nx;
  /*  Round up to limb boundary.  */
  size_t w = ((max_nb + ARBINT_LIMB_BITS - 1u) / ARBINT_LIMB_BITS)
              * ARBINT_LIMB_BITS;
  CHECK_EQ_I(pc_x + pc_nx, w);

  arbint_clear(nx);
  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

/* ========== hamming symmetry ========== */

static void test_hammingdist_symmetry(void) {
  arbint_ctx_t ctx;
  arbint_t a, b;
  size_t h1;
  size_t h2;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);

  /*  hamming(a, b) == hamming(b, a) for all sign combos.  */

  /*  (+, +).  */
  CHECK_EQ_I(arbint_set_i32(a, 123), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 456), ARBINT_OK);
  CHECK_EQ_I(arbint_hammingdist(a, b, &h1), ARBINT_OK);
  CHECK_EQ_I(arbint_hammingdist(b, a, &h2), ARBINT_OK);
  CHECK_EQ_I(h1, h2);

  /*  (+, -).  */
  CHECK_EQ_I(arbint_set_i32(a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -50), ARBINT_OK);
  CHECK_EQ_I(arbint_hammingdist(a, b, &h1), ARBINT_OK);
  CHECK_EQ_I(arbint_hammingdist(b, a, &h2), ARBINT_OK);
  CHECK_EQ_I(h1, h2);

  /*  (-, -).  */
  CHECK_EQ_I(arbint_set_i32(a, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -13), ARBINT_OK);
  CHECK_EQ_I(arbint_hammingdist(a, b, &h1), ARBINT_OK);
  CHECK_EQ_I(arbint_hammingdist(b, a, &h2), ARBINT_OK);
  CHECK_EQ_I(h1, h2);

  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ========== multi-limb signed bitwise ops ========== */

static void test_bitwise_multilimb_signed(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r1, r2;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r2, &ctx), ARBINT_OK);

  /*  Build multi-limb negative values via squaring + negate.
      a = -(3^64), b = -(5^64).  Both multi-limb.  */
  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 5), ARBINT_OK);
  {
    int i;
    for (i = 0; i < 6; ++i) {
      CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
      CHECK_EQ_I(arbint_sqr(b, b), ARBINT_OK);
    }
  }
  CHECK_EQ_I(arbint_neg(a, a), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(b, b), ARBINT_OK);

  /*  (-a) & (-b): verify identity (a & b) | (a ^ b) = (a | b).  */
  {
    arbint_t t1, t2;
    CHECK_EQ_I(arbint_init(t1, &ctx), ARBINT_OK);
    CHECK_EQ_I(arbint_init(t2, &ctx), ARBINT_OK);

    CHECK_EQ_I(arbint_and(r1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_xor(r2, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_or(t1, r1, r2), ARBINT_OK);
    CHECK_EQ_I(arbint_or(t2, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(t1, t2), 0);

    arbint_clear(t2);
    arbint_clear(t1);
  }

  /*  not(not(a)) = a for multi-limb negative.  */
  CHECK_EQ_I(arbint_not(r1, a), ARBINT_OK);
  CHECK_EQ_I(arbint_not(r2, r1), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r2, a), 0);

  /*  hamming symmetry for multi-limb negative.  */
  {
    size_t h1;
    size_t h2;
    CHECK_EQ_I(arbint_hammingdist(a, b, &h1), ARBINT_OK);
    CHECK_EQ_I(arbint_hammingdist(b, a, &h2), ARBINT_OK);
    CHECK_EQ_I(h1, h2);
    CHECK(h1 > 0u);
  }

  arbint_clear(r2);
  arbint_clear(r1);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ========== stochastic with full identity check ========== */

static void test_stochastic_identity(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, ab_and, ab_xor, ab_or, lhs;
  enum { ITERS = 200 };
  int i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(ab_and, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(ab_xor, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(ab_or, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(lhs, &ctx), ARBINT_OK);

  /*  Verify (a & b) | (a ^ b) == (a | b) for random signed pairs.
      Covers all 4 sign combinations.  */
  for (i = 0; i < ITERS && g_failures == 0; ++i) {
    int32_t av = (int32_t) (i * 73 - 7000);
    int32_t bv = (int32_t) (i * 97 - 10000);

    CHECK_EQ_I(arbint_set_i32(a, av), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, bv), ARBINT_OK);

    CHECK_EQ_I(arbint_and(ab_and, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_xor(ab_xor, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_or(lhs, ab_and, ab_xor), ARBINT_OK);
    CHECK_EQ_I(arbint_or(ab_or, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(lhs, ab_or), 0);
  }

  arbint_clear(lhs);
  arbint_clear(ab_or);
  arbint_clear(ab_xor);
  arbint_clear(ab_and);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ========== i32/u32 multi-limb cross-validation ========== */

static void test_imm_multilimb(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r1, r2;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r2, &ctx), ARBINT_OK);

  /*  Build multi-limb a = 3^128.  */
  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  {
    int i;
    for (i = 0; i < 7; ++i)
      CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
  }

  /*  and_u32 with multi-limb positive.  */
  CHECK_EQ_I(arbint_set_u32(b, 0xFFu), ARBINT_OK);
  CHECK_EQ_I(arbint_and(r1, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_and_u32(r2, a, 0xFFu), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r1, r2), 0);

  /*  or_u32 with multi-limb positive.  */
  CHECK_EQ_I(arbint_set_u32(b, 0xFFu), ARBINT_OK);
  CHECK_EQ_I(arbint_or(r1, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_or_u32(r2, a, 0xFFu), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r1, r2), 0);

  /*  xor_u32 with multi-limb positive.  */
  CHECK_EQ_I(arbint_set_u32(b, 0xFFu), ARBINT_OK);
  CHECK_EQ_I(arbint_xor(r1, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_xor_u32(r2, a, 0xFFu), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r1, r2), 0);

  /*  Negate a for multi-limb negative tests.  */
  CHECK_EQ_I(arbint_neg(a, a), ARBINT_OK);

  CHECK_EQ_I(arbint_set_u32(b, 0xFFu), ARBINT_OK);
  CHECK_EQ_I(arbint_and(r1, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_and_u32(r2, a, 0xFFu), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r1, r2), 0);

  CHECK_EQ_I(arbint_or(r1, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_or_u32(r2, a, 0xFFu), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r1, r2), 0);

  CHECK_EQ_I(arbint_xor(r1, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_xor_u32(r2, a, 0xFFu), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r1, r2), 0);

  /*  i32 negative b with multi-limb a.  */
  CHECK_EQ_I(arbint_set_i32(b, -2), ARBINT_OK);
  CHECK_EQ_I(arbint_and(r1, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_and_i32(r2, a, -2), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r1, r2), 0);

  CHECK_EQ_I(arbint_or(r1, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_or_i32(r2, a, -2), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r1, r2), 0);

  CHECK_EQ_I(arbint_xor(r1, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_xor_i32(r2, a, -2), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r1, r2), 0);

  arbint_clear(r2);
  arbint_clear(r1);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

int main(void) {
  test_testbit_zero();
  test_testbit_positive();
  test_testbit_negative();
  test_testbit_null();
  test_setbit_positive();
  test_setbit_negative();
  test_setbit_extend();
  test_clrbit_positive();
  test_clrbit_negative();
  test_ctz_basic();
  test_ctz_multilimb();
  test_clz_basic();
  test_popcount_basic();
  test_popcount_negative();
  test_popcount_large();
  test_popcount_complement();
  test_hammingdist_basic();
  test_hammingdist_negative();
  test_hammingdist_mixed_sign();
  test_hammingdist_symmetry();
  test_not_basic();
  test_not_aliasing();
  test_or_basic();
  test_or_negative();
  test_and_basic();
  test_and_negative();
  test_xor_basic();
  test_xor_negative();
  test_bitwise_aliasing();
  test_bitwise_identities();
  test_bitwise_multilimb();
  test_bitwise_multilimb_signed();
  test_bitwise_stochastic();
  test_stochastic_identity();
  test_and_u32();
  test_and_i32();
  test_or_u32();
  test_or_i32();
  test_xor_u32();
  test_xor_i32();
  test_imm_cross_validation();
  test_imm_multilimb();
  test_null_checks();

  ARBINT_TEST_FINISH("test_bitops");
}
