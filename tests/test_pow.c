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

static void test_pow_edge_cases(void) {
  arbint_ctx_t ctx;
  arbint_t base, rop;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(base, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(rop, &ctx), ARBINT_OK);

  /*  x^0 = 1 for any x.  */
  CHECK_EQ_I(arbint_set_i32(base, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rop, base, 0u), ARBINT_OK);
  check_i32_value(rop, 1);

  /*  0^0 = 1 (GMP convention).  */
  CHECK_EQ_I(arbint_set_i32(base, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rop, base, 0u), ARBINT_OK);
  check_i32_value(rop, 1);

  /*  0^1 = 0.  */
  CHECK_EQ_I(arbint_set_i32(base, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rop, base, 1u), ARBINT_OK);
  check_i32_value(rop, 0);

  /*  0^n = 0 for n > 0.  */
  CHECK_EQ_I(arbint_set_i32(base, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rop, base, 100u), ARBINT_OK);
  check_i32_value(rop, 0);

  /*  x^1 = x.  */
  CHECK_EQ_I(arbint_set_i32(base, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rop, base, 1u), ARBINT_OK);
  check_i32_value(rop, 42);

  /*  1^n = 1 for any n.  */
  CHECK_EQ_I(arbint_set_i32(base, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rop, base, 1000u), ARBINT_OK);
  check_i32_value(rop, 1);

  /*  (-1)^even = 1.  */
  CHECK_EQ_I(arbint_set_i32(base, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rop, base, 100u), ARBINT_OK);
  check_i32_value(rop, 1);

  /*  (-1)^odd = -1.  */
  CHECK_EQ_I(arbint_set_i32(base, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rop, base, 99u), ARBINT_OK);
  check_i32_value(rop, -1);

  /*  NULL inputs.  */
  CHECK_EQ_I(arbint_pow_u32(NULL, base, 2u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_pow_u32(rop, NULL, 2u), ARBINT_EINVAL);

  arbint_clear(rop);
  arbint_clear(base);
  arbint_ctx_clear(&ctx);
}

static void test_pow_small_values(void) {
  arbint_ctx_t ctx;
  arbint_t base, rop;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(base, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(rop, &ctx), ARBINT_OK);

  /*  2^10 = 1024.  */
  CHECK_EQ_I(arbint_set_i32(base, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rop, base, 10u), ARBINT_OK);
  check_i32_value(rop, 1024);

  /*  2^16 = 65536.  */
  CHECK_EQ_I(arbint_set_i32(base, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rop, base, 16u), ARBINT_OK);
  check_i32_value(rop, 65536);

  /*  3^7 = 2187.  */
  CHECK_EQ_I(arbint_set_i32(base, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rop, base, 7u), ARBINT_OK);
  check_i32_value(rop, 2187);

  /*  5^5 = 3125.  */
  CHECK_EQ_I(arbint_set_i32(base, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rop, base, 5u), ARBINT_OK);
  check_i32_value(rop, 3125);

  /*  10^9 = 1000000000.  */
  CHECK_EQ_I(arbint_set_i32(base, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rop, base, 9u), ARBINT_OK);
  check_i32_value(rop, 1000000000);

  /*  (-2)^10 = 1024.  */
  CHECK_EQ_I(arbint_set_i32(base, -2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rop, base, 10u), ARBINT_OK);
  check_i32_value(rop, 1024);

  /*  (-2)^11 = -2048.  */
  CHECK_EQ_I(arbint_set_i32(base, -2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rop, base, 11u), ARBINT_OK);
  check_i32_value(rop, -2048);

  /*  (-3)^5 = -243.  */
  CHECK_EQ_I(arbint_set_i32(base, -3), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rop, base, 5u), ARBINT_OK);
  check_i32_value(rop, -243);

  arbint_clear(rop);
  arbint_clear(base);
  arbint_ctx_clear(&ctx);
}

static void test_pow_aliasing(void) {
  arbint_ctx_t ctx;
  arbint_t a;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);

  /*  rop == base aliasing: a = a^10, where a starts as 2.  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 10u), ARBINT_OK);
  check_i32_value(a, 1024);

  /*  rop == base aliasing with exp == 0.  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 0u), ARBINT_OK);
  check_i32_value(a, 1);

  /*  rop == base aliasing with exp == 1.  */
  CHECK_EQ_I(arbint_set_i32(a, 99), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 1u), ARBINT_OK);
  check_i32_value(a, 99);

  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_pow_large(void) {
  arbint_ctx_t ctx;
  arbint_t base, rop, expected;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(base, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(rop, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  /*  2^64: should produce a 2-limb result on 64-bit, 3-limb on 32-bit.
      Verify by computing 2^32 * 2^32 and comparing.  */
  CHECK_EQ_I(arbint_set_i32(base, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rop, base, 64u), ARBINT_OK);

  /*  Build expected = 2^64 via arbint_set_i32(1) then shl 64.  */
  CHECK_EQ_I(arbint_set_i32(expected, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(expected, expected, 64u), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(rop, expected), 0);

  /*  2^128.  */
  CHECK_EQ_I(arbint_pow_u32(rop, base, 128u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(expected, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(expected, expected, 128u), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(rop, expected), 0);

  /*  2^256.  */
  CHECK_EQ_I(arbint_pow_u32(rop, base, 256u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(expected, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(expected, expected, 256u), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(rop, expected), 0);

  /*  Verify pow via repeated multiplication:
      3^20 computed two ways.  */
  CHECK_EQ_I(arbint_set_i32(base, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rop, base, 20u), ARBINT_OK);

  CHECK_EQ_I(arbint_set_i32(expected, 1), ARBINT_OK);
  {
    int i;
    for (i = 0; i < 20; ++i) {
      CHECK_EQ_I(arbint_mul_i32(expected, expected, 3), ARBINT_OK);
    }
  }
  CHECK_EQ_I(arbint_cmp(rop, expected), 0);

  arbint_clear(expected);
  arbint_clear(rop);
  arbint_clear(base);
  arbint_ctx_clear(&ctx);
}

static void test_pow_large_base(void) {
  arbint_ctx_t ctx;
  arbint_t base, rop, expected;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(base, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(rop, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  /*  (2^64 + 1)^2 = 2^128 + 2^65 + 1.
      Build base = 2^64 + 1.  */
  CHECK_EQ_I(arbint_set_i32(base, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(base, base, 64u), ARBINT_OK);
  CHECK_EQ_I(arbint_add_i32(base, base, 1), ARBINT_OK);

  CHECK_EQ_I(arbint_pow_u32(rop, base, 2u), ARBINT_OK);

  /*  expected = base * base.  */
  CHECK_EQ_I(arbint_mul(expected, base, base), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(rop, expected), 0);

  /*  (2^64 + 1)^3 = base^2 * base.  */
  CHECK_EQ_I(arbint_pow_u32(rop, base, 3u), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(expected, expected, base), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(rop, expected), 0);

  arbint_clear(expected);
  arbint_clear(rop);
  arbint_clear(base);
  arbint_ctx_clear(&ctx);
}

int main(void) {
  test_pow_edge_cases();
  test_pow_small_values();
  test_pow_aliasing();
  test_pow_large();
  test_pow_large_base();
  ARBINT_TEST_FINISH("test_pow");
}
