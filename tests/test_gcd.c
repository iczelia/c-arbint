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
#include <string.h>

ARBINT_TEST_DECLARE_FAILURES();

/*  Basic GCD tests.  */

static void test_gcd_basic(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);

  /* gcd(12, 8) = 4 */
  CHECK_EQ_I(arbint_set_i32(a, 12), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 8), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 4u);

  /* gcd(17, 13) = 1 (coprime) */
  CHECK_EQ_I(arbint_set_i32(a, 17), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 13), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 1u);

  /* gcd(100, 25) = 25 */
  CHECK_EQ_I(arbint_set_i32(a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 25), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 25u);

  /* gcd(48, 18) = 6 */
  CHECK_EQ_I(arbint_set_i32(a, 48), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 18), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 6u);

  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Edge cases: zeros.  */

static void test_gcd_zeros(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);

  /* gcd(0, 0) = 0 */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  CHECK(arbint_is_zero(g));

  /* gcd(0, 42) = 42 */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 42u);

  /* gcd(42, 0) = 42 */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 42u);

  /* gcd(0, -42) = 42 (absolute value) */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -42), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 42u);

  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Negative input tests.  */

static void test_gcd_negative(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);

  /* gcd(-36, 24) = 12 */
  CHECK_EQ_I(arbint_set_i32(a, -36), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 24), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 12u);

  /* gcd(36, -24) = 12 */
  CHECK_EQ_I(arbint_set_i32(a, 36), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -24), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 12u);

  /* gcd(-36, -24) = 12 */
  CHECK_EQ_I(arbint_set_i32(a, -36), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -24), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 12u);

  /* Result should always be non-negative */
  CHECK(g[0]._sz >= 0);

  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Symmetry test: gcd(a, b) == gcd(b, a).  */

static void test_gcd_symmetry(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g1, g2;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g2, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_set_i32(a, 123456), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 789012), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g1, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g2, b, a), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g1, g2), 0);

  arbint_clear(g2);
  arbint_clear(g1);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Aliasing tests.  */

static void test_gcd_aliasing(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);

  /* g = a case */
  CHECK_EQ_I(arbint_set_i32(a, 48), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 18), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(a, a, b), ARBINT_OK);
  check_u32_value(a, 6u);

  /* g = b case */
  CHECK_EQ_I(arbint_set_i32(a, 48), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 18), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(b, a, b), ARBINT_OK);
  check_u32_value(b, 6u);

  /* gcd(a, a, a) - all three alias */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(a, a, a), ARBINT_OK);
  check_u32_value(a, 42u);

  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  GCD u32 variant tests.  */

static void test_gcd_u32(void) {
  arbint_ctx_t ctx;
  arbint_t a, g;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);

  /* gcd(48, 18) = 6 */
  CHECK_EQ_I(arbint_set_i32(a, 48), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd_u32(g, a, 18u), ARBINT_OK);
  check_u32_value(g, 6u);

  /* gcd(0, 42) = 42 */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd_u32(g, a, 42u), ARBINT_OK);
  check_u32_value(g, 42u);

  /* gcd(42, 0) = 42 */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd_u32(g, a, 0u), ARBINT_OK);
  check_u32_value(g, 42u);

  /* gcd(-48, 18) = 6 (negative input) */
  CHECK_EQ_I(arbint_set_i32(a, -48), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd_u32(g, a, 18u), ARBINT_OK);
  check_u32_value(g, 6u);

  /* gcd(a, 1) = 1 */
  CHECK_EQ_I(arbint_set_i32(a, 123456789), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd_u32(g, a, 1u), ARBINT_OK);
  check_u32_value(g, 1u);

  arbint_clear(g);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Basic LCM tests.  */

static void test_lcm_basic(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, l;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(l, &ctx), ARBINT_OK);

  /* lcm(4, 6) = 12 */
  CHECK_EQ_I(arbint_set_i32(a, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 6), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm(l, a, b), ARBINT_OK);
  check_u32_value(l, 12u);

  /* lcm(12, 8) = 24 */
  CHECK_EQ_I(arbint_set_i32(a, 12), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 8), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm(l, a, b), ARBINT_OK);
  check_u32_value(l, 24u);

  /* lcm(a, a) = |a| */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm(l, a, a), ARBINT_OK);
  check_u32_value(l, 42u);

  /* lcm(1, n) = |n| */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm(l, a, b), ARBINT_OK);
  check_u32_value(l, 100u);

  arbint_clear(l);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  LCM zero cases.  */

static void test_lcm_zeros(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, l;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(l, &ctx), ARBINT_OK);

  /* lcm(0, x) = 0 */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm(l, a, b), ARBINT_OK);
  CHECK(arbint_is_zero(l));

  /* lcm(x, 0) = 0 */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm(l, a, b), ARBINT_OK);
  CHECK(arbint_is_zero(l));

  /* lcm(0, 0) = 0 */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm(l, a, b), ARBINT_OK);
  CHECK(arbint_is_zero(l));

  arbint_clear(l);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  LCM negative inputs.  */

static void test_lcm_negative(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, l;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(l, &ctx), ARBINT_OK);

  /* lcm(-4, 6) = 12 */
  CHECK_EQ_I(arbint_set_i32(a, -4), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 6), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm(l, a, b), ARBINT_OK);
  check_u32_value(l, 12u);

  /* lcm(4, -6) = 12 */
  CHECK_EQ_I(arbint_set_i32(a, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -6), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm(l, a, b), ARBINT_OK);
  check_u32_value(l, 12u);

  /* lcm(-4, -6) = 12 */
  CHECK_EQ_I(arbint_set_i32(a, -4), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -6), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm(l, a, b), ARBINT_OK);
  check_u32_value(l, 12u);

  /* Result should always be non-negative */
  CHECK(l[0]._sz >= 0);

  arbint_clear(l);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  LCM aliasing tests.  */

static void test_lcm_aliasing(void) {
  arbint_ctx_t ctx;
  arbint_t a, b;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);

  /* l = a case */
  CHECK_EQ_I(arbint_set_i32(a, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 6), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm(a, a, b), ARBINT_OK);
  check_u32_value(a, 12u);

  /* l = b case */
  CHECK_EQ_I(arbint_set_i32(a, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 6), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm(b, a, b), ARBINT_OK);
  check_u32_value(b, 12u);

  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  LCM u32 variant tests.  */

static void test_lcm_u32(void) {
  arbint_ctx_t ctx;
  arbint_t a, l;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(l, &ctx), ARBINT_OK);

  /* lcm(4, 6) = 12 */
  CHECK_EQ_I(arbint_set_i32(a, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm_u32(l, a, 6u), ARBINT_OK);
  check_u32_value(l, 12u);

  /* lcm(0, 6) = 0 */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm_u32(l, a, 6u), ARBINT_OK);
  CHECK(arbint_is_zero(l));

  /* lcm(4, 0) = 0 */
  CHECK_EQ_I(arbint_set_i32(a, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm_u32(l, a, 0u), ARBINT_OK);
  CHECK(arbint_is_zero(l));

  /* lcm(-4, 6) = 12 */
  CHECK_EQ_I(arbint_set_i32(a, -4), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm_u32(l, a, 6u), ARBINT_OK);
  check_u32_value(l, 12u);

  arbint_clear(l);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Multi-limb GCD test: powers of 2.  */

static void test_gcd_powers_of_two(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g, expected;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  /* gcd(2^128, 2^64) = 2^64 */
  /* Build 2^64 by repeated squaring from 2^16 */
  CHECK_EQ_I(arbint_set_u32(a, 65536u), ARBINT_OK); /* 2^16 */
  CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);          /* 2^32 */
  CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);          /* 2^64 */
  CHECK_EQ_I(arbint_set(expected, a), ARBINT_OK);   /* save 2^64 */
  CHECK_EQ_I(arbint_sqr(b, a), ARBINT_OK);          /* 2^128 */

  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, expected), 0);

  /* Symmetry: gcd(2^128, 2^64) = gcd(2^64, 2^128) */
  CHECK_EQ_I(arbint_gcd(g, b, a), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, expected), 0);

  arbint_clear(expected);
  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Multi-limb GCD test: shared factors.  */

static void test_gcd_shared_factors(void) {
  arbint_ctx_t ctx;
  arbint_t base, a, b, g;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(base, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);

  /* Build base = 3^10 = 59049 */
  CHECK_EQ_I(arbint_set_u32(base, 59049u), ARBINT_OK);

  /* a = base^2 = 3^20 */
  CHECK_EQ_I(arbint_sqr(a, base), ARBINT_OK);

  /* b = base^4 = 3^40 */
  CHECK_EQ_I(arbint_sqr(b, a), ARBINT_OK);

  /* gcd(3^20, 3^40) = 3^20 */
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, a), 0);

  /* gcd(3^40, 3^20) = 3^20 (symmetry) */
  CHECK_EQ_I(arbint_gcd(g, b, a), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, a), 0);

  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_clear(base);
  arbint_ctx_clear(&ctx);
}

/*  Divisibility property: gcd(a, b) divides both a and b.  */

static void test_gcd_divisibility(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /* Test with relatively prime numbers */
  CHECK_EQ_I(arbint_set_i32(a, 123456), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 789012), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);

  /* g should divide a */
  CHECK_EQ_I(arbint_tdiv_r(r, a, g), ARBINT_OK);
  CHECK(arbint_is_zero(r));

  /* g should divide b */
  CHECK_EQ_I(arbint_tdiv_r(r, b, g), ARBINT_OK);
  CHECK(arbint_is_zero(r));

  arbint_clear(r);
  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Identity: gcd(a, b) * lcm(a, b) = |a| * |b| (for small values).  */

static void test_gcd_lcm_identity(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g, l, prod1, prod2;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(l, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(prod1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(prod2, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_set_i32(a, 48), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 18), ARBINT_OK);

  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm(l, a, b), ARBINT_OK);

  /* prod1 = gcd * lcm */
  CHECK_EQ_I(arbint_mul(prod1, g, l), ARBINT_OK);

  /* prod2 = |a| * |b| */
  CHECK_EQ_I(arbint_abs(prod2, a), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(prod2, prod2, b), ARBINT_OK);
  CHECK_EQ_I(arbint_abs(prod2, prod2), ARBINT_OK);

  CHECK_EQ_I(arbint_cmp(prod1, prod2), 0);

  arbint_clear(prod2);
  arbint_clear(prod1);
  arbint_clear(l);
  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  LCM exactness: lcm(a, b) is divisible by both |a| and |b|.  */

static void test_lcm_exactness(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, l, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(l, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_set_i32(a, 123), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 456), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm(l, a, b), ARBINT_OK);

  /* l should be divisible by |a| */
  CHECK_EQ_I(arbint_tdiv_r(r, l, a), ARBINT_OK);
  CHECK(arbint_is_zero(r));

  /* l should be divisible by |b| */
  CHECK_EQ_I(arbint_tdiv_r(r, l, b), ARBINT_OK);
  CHECK(arbint_is_zero(r));

  arbint_clear(r);
  arbint_clear(l);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Size-disparate operands (tests the Euclidean path).  */

static void test_gcd_size_disparate(void) {
  arbint_ctx_t ctx;
  arbint_t big, small, g;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(big, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(small, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);

  /* Build big = 2^256 (large multi-limb number) */
  CHECK_EQ_I(arbint_set_u32(big, 2u), ARBINT_OK);
  for (int i = 0; i < 8; ++i)
    CHECK_EQ_I(arbint_sqr(big, big), ARBINT_OK);

  /* gcd(2^256, 6) = 2 */
  CHECK_EQ_I(arbint_set_u32(small, 6u), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, big, small), ARBINT_OK);
  check_u32_value(g, 2u);

  /* gcd(2^256, 5) = 1 (coprime) */
  CHECK_EQ_I(arbint_set_u32(small, 5u), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, big, small), ARBINT_OK);
  check_u32_value(g, 1u);

  arbint_clear(g);
  arbint_clear(small);
  arbint_clear(big);
  arbint_ctx_clear(&ctx);
}

/*  Large operand test using binary GCD path (both > 4 limbs).  */

static void test_gcd_large_binary(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g, factor;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(factor, &ctx), ARBINT_OK);

  /* Build large numbers with known GCD.
     a = 7 * 2^512
     b = 21 * 2^512
     gcd(a, b) = 7 * 2^512 */

  /* factor = 2^512 */
  CHECK_EQ_I(arbint_set_u32(factor, 2u), ARBINT_OK);
  for (int i = 0; i < 9; ++i)
    CHECK_EQ_I(arbint_sqr(factor, factor), ARBINT_OK);

  /* a = 7 * factor */
  CHECK_EQ_I(arbint_mul_u32(a, factor, 7u), ARBINT_OK);

  /* b = 21 * factor */
  CHECK_EQ_I(arbint_mul_u32(b, factor, 21u), ARBINT_OK);

  /* gcd(7*factor, 21*factor) should be 7*factor */
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, a), 0);

  arbint_clear(factor);
  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Size-disparate binary GCD (tests the mod optimization).
    Both operands > 4 limbs but one much larger than the other.  */

static void test_gcd_binary_mod_opt(void) {
  arbint_ctx_t ctx;
  arbint_t big, medium, g, expected;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(big, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(medium, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  /*  big = 2^2048 (32 limbs on 64-bit), medium = 7 * 2^320 (5 limbs).
      gcd should be 2^320 since 2^2048 = 2^320 * 2^1728.  */
  CHECK_EQ_I(arbint_set_u32(big, 2u), ARBINT_OK);
  for (int i = 0; i < 11; ++i)
    CHECK_EQ_I(arbint_sqr(big, big), ARBINT_OK); /* 2^2048 */

  CHECK_EQ_I(arbint_set_u32(expected, 2u), ARBINT_OK);
  for (int i = 0; i < 5; ++i)
    CHECK_EQ_I(arbint_sqr(expected, expected), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(expected, expected), ARBINT_OK); /* 2^320 */

  CHECK_EQ_I(arbint_mul_u32(medium, expected, 7u), ARBINT_OK); /* 7 * 2^320 */

  CHECK_EQ_I(arbint_gcd(g, big, medium), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, expected), 0);

  /*  Also test reverse order.  */
  CHECK_EQ_I(arbint_gcd(g, medium, big), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, expected), 0);

  /*  Test coprime case: gcd(2^2048, 3 * 2^320 + 1) = 1.  */
  CHECK_EQ_I(arbint_mul_u32(medium, expected, 3u), ARBINT_OK);
  CHECK_EQ_I(arbint_add_u32(medium, medium, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, big, medium), ARBINT_OK);
  check_u32_value(g, 1u);

  arbint_clear(expected);
  arbint_clear(g);
  arbint_clear(medium);
  arbint_clear(big);
  arbint_ctx_clear(&ctx);
}

/*  Normalization test: result has no leading zeros and _sz >= 0.  */

static void test_gcd_normalization(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g;
  size_t used;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_set_i32(a, 12345), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 54321), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);

  /* Result must be non-negative */
  CHECK(g[0]._sz >= 0);

  /* If nonzero, top limb must be nonzero (no leading zeros) */
  used = (size_t) g[0]._sz;
  if (used > 0u) {
    const arbint_limb_t * gp = ARBINT_CLIMBS(g);
    CHECK(gp[used - 1u] != 0u);
  }

  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Single-limb values > UINT32_MAX (64-bit regression test).
    Ensures the single-limb fast path doesn't truncate to 32 bits.  */
static void test_gcd_large_single_limb(void) {
#if ARBINT_LIMB_BITS == 64
  arbint_ctx_t ctx;
  arbint_t a, b, g, l;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(l, &ctx), ARBINT_OK);

  /*  gcd(2^32, 2^32) = 2^32 (not 0 from truncation).  */
  CHECK_EQ_I(arbint_set_u32(a, 65536u), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK); /* a = 2^32 */
  CHECK_EQ_I(arbint_set(b, a), ARBINT_OK); /* b = 2^32 */

  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, a), 0); /* gcd should equal a */

  /*  lcm(2^32, 2^32) = 2^32 (must not fail with EZERO).  */
  CHECK_EQ_I(arbint_lcm(l, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(l, a), 0); /* lcm should equal a */

  /*  gcd(2^33, 2^32) = 2^32.  */
  CHECK_EQ_I(arbint_mul_u32(a, a, 2u), ARBINT_OK); /* a = 2^33 */
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, b), 0); /* gcd should equal b = 2^32 */

  /*  gcd with large coprime single-limb values.  */
  CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 40u), ARBINT_OK); /* a = 2^40 */
  CHECK_EQ_I(arbint_set_u32(b, 3u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(b, b, 35u), ARBINT_OK); /* b = 3 * 2^35 */
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  /* gcd(2^40, 3*2^35) = 2^35 */
  CHECK_EQ_I(arbint_set_u32(l, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(l, l, 35u), ARBINT_OK); /* expected = 2^35 */
  CHECK_EQ_I(arbint_cmp(g, l), 0);

  arbint_clear(l);
  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
#endif /* ARBINT_LIMB_BITS == 64 */
}

/*  NULL pointer tests.  */
static void test_gcd_null_pointers(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_gcd(NULL, a, b), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_gcd(g, NULL, b), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_gcd(g, a, NULL), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_gcd_u32(NULL, a, 5u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_gcd_u32(g, NULL, 5u), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_lcm(NULL, a, b), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_lcm(g, NULL, b), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_lcm(g, a, NULL), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_lcm_u32(NULL, a, 5u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_lcm_u32(g, NULL, 5u), ARBINT_EINVAL);

  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Fibonacci pairs: gcd(F_n, F_{n-1}) = 1.
    Worst-case for Euclidean algorithm (maximum iterations).  */
static void test_gcd_fibonacci(void) {
  arbint_ctx_t ctx;
  arbint_t f_prev, f_curr, f_next, g;
  int i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(f_prev, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(f_curr, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(f_next, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);

  /*  F_0 = 0, F_1 = 1.  */
  CHECK_EQ_I(arbint_set_u32(f_prev, 0u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(f_curr, 1u), ARBINT_OK);

  /*  Compute up to F_100 (multi-limb).  */
  for (i = 2; i <= 100; ++i) {
    CHECK_EQ_I(arbint_add(f_next, f_curr, f_prev), ARBINT_OK);
    CHECK_EQ_I(arbint_set(f_prev, f_curr), ARBINT_OK);
    CHECK_EQ_I(arbint_set(f_curr, f_next), ARBINT_OK);
  }

  /*  gcd(F_100, F_99) = 1.  */
  CHECK_EQ_I(arbint_gcd(g, f_curr, f_prev), ARBINT_OK);
  check_u32_value(g, 1u);

  /*  gcd(F_99, F_100) = 1 (symmetry).  */
  CHECK_EQ_I(arbint_gcd(g, f_prev, f_curr), ARBINT_OK);
  check_u32_value(g, 1u);

  arbint_clear(g);
  arbint_clear(f_next);
  arbint_clear(f_curr);
  arbint_clear(f_prev);
  arbint_ctx_clear(&ctx);
}

/*  Consecutive integers are always coprime: gcd(n, n+1) = 1.  */
static void test_gcd_consecutive(void) {
  arbint_ctx_t ctx;
  arbint_t n, n_plus_1, g;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n_plus_1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);

  /*  Small consecutive.  */
  CHECK_EQ_I(arbint_set_u32(n, 999u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(n_plus_1, 1000u), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, n, n_plus_1), ARBINT_OK);
  check_u32_value(g, 1u);

  /*  Large consecutive (multi-limb): n = 2^256, n+1 = 2^256 + 1.  */
  CHECK_EQ_I(arbint_set_u32(n, 2u), ARBINT_OK);
  for (int i = 0; i < 8; ++i)
    CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 2^256 */
  CHECK_EQ_I(arbint_add_u32(n_plus_1, n, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, n, n_plus_1), ARBINT_OK);
  check_u32_value(g, 1u);

  arbint_clear(g);
  arbint_clear(n_plus_1);
  arbint_clear(n);
  arbint_ctx_clear(&ctx);
}

/*  Scaling property: gcd(k*a, k*b) = k * gcd(a, b).  */
static void test_gcd_scaling(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, ka, kb, g_ab, g_kakb, expected;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(ka, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(kb, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g_ab, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g_kakb, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_set_u32(a, 48u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 18u), ARBINT_OK);

  /*  gcd(48, 18) = 6.  */
  CHECK_EQ_I(arbint_gcd(g_ab, a, b), ARBINT_OK);
  check_u32_value(g_ab, 6u);

  /*  gcd(48*1000, 18*1000) = 6*1000 = 6000.  */
  CHECK_EQ_I(arbint_mul_u32(ka, a, 1000u), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(kb, b, 1000u), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g_kakb, ka, kb), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(expected, g_ab, 1000u), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g_kakb, expected), 0);

  /*  Large scale factor: k = 2^128.  */
  CHECK_EQ_I(arbint_set_u32(expected, 2u), ARBINT_OK);
  for (int i = 0; i < 7; ++i)
    CHECK_EQ_I(arbint_sqr(expected, expected), ARBINT_OK); /* k = 2^128 */

  CHECK_EQ_I(arbint_mul(ka, a, expected), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(kb, b, expected), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g_kakb, ka, kb), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(expected, g_ab, expected), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g_kakb, expected), 0);

  arbint_clear(expected);
  arbint_clear(g_kakb);
  arbint_clear(g_ab);
  arbint_clear(kb);
  arbint_clear(ka);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Threshold boundary tests: test around ARBINT_GCD_EUCLID_THRESHOLD (4).
    Ensures correct dispatch between Euclidean and binary GCD paths.  */
static void test_gcd_threshold_boundary(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g, r;
  int i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Build numbers of exactly 3, 4, 5, 6 limbs (on 64-bit: 192, 256, 320, 384
     bits). Use 2^(64*n) - 1 which has exactly n limbs.  */
  for (i = 3; i <= 6; ++i) {
    /*  a = 2^(64*i) - 1 (all 1-bits, exactly i limbs on 64-bit).  */
    CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(a, a, (unsigned) (ARBINT_LIMB_BITS * i)), ARBINT_OK);
    CHECK_EQ_I(arbint_sub_u32(a, a, 1u), ARBINT_OK);

    /*  b = a - 2 (same size, different value).  */
    CHECK_EQ_I(arbint_sub_u32(b, a, 2u), ARBINT_OK);

    CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);

    /*  Verify gcd divides both operands.  */
    CHECK_EQ_I(arbint_tdiv_r(r, a, g), ARBINT_OK);
    CHECK(arbint_is_zero(r));
    CHECK_EQ_I(arbint_tdiv_r(r, b, g), ARBINT_OK);
    CHECK(arbint_is_zero(r));

    /*  Result must be positive and normalized.  */
    CHECK(g[0]._sz > 0);
  }

  arbint_clear(r);
  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  GCD with prime numbers.  */
static void test_gcd_primes(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);

  /*  Two distinct primes are always coprime.  */
  CHECK_EQ_I(arbint_set_u32(a, 104729u), ARBINT_OK);  /* 10000th prime */
  CHECK_EQ_I(arbint_set_u32(b, 1299709u), ARBINT_OK); /* 100000th prime */
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 1u);

  /*  gcd(p, p^2) = p.  */
  CHECK_EQ_I(arbint_set_u32(a, 104729u), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(b, a), ARBINT_OK); /* b = p^2 */
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, a), 0);

  /*  gcd(p*q, p) = p for distinct primes p, q.  */
  CHECK_EQ_I(arbint_set_u32(a, 104729u), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(b, a, 1299709u), ARBINT_OK); /* b = p*q */
  CHECK_EQ_I(arbint_gcd(g, b, a), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, a), 0);

  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Stochastic GCD tests: verify properties over random-ish inputs.  */
static void test_gcd_stochastic(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g1, g2, r;
  int iter;
  uint32_t seed = 12345u;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g2, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  for (iter = 0; iter < 50; ++iter) {
    /*  Simple LCG for pseudo-random values.  */
    seed = seed * 1103515245u + 12345u;
    uint32_t av = (seed >> 16) | 1u; /* Ensure nonzero, odd */
    seed = seed * 1103515245u + 12345u;
    uint32_t bv = (seed >> 16) | 1u;

    CHECK_EQ_I(arbint_set_u32(a, av), ARBINT_OK);
    CHECK_EQ_I(arbint_set_u32(b, bv), ARBINT_OK);

    /*  Property 1: symmetry.  */
    CHECK_EQ_I(arbint_gcd(g1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_gcd(g2, b, a), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(g1, g2), 0);

    /*  Property 2: gcd divides both.  */
    CHECK_EQ_I(arbint_tdiv_r(r, a, g1), ARBINT_OK);
    CHECK(arbint_is_zero(r));
    CHECK_EQ_I(arbint_tdiv_r(r, b, g1), ARBINT_OK);
    CHECK(arbint_is_zero(r));

    /*  Property 3: result is positive.  */
    CHECK(g1[0]._sz >= 0);
  }

  arbint_clear(r);
  arbint_clear(g2);
  arbint_clear(g1);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Stochastic LCM tests.  */
static void test_lcm_stochastic(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, l, g, prod1, prod2, r;
  int iter;
  uint32_t seed = 54321u;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(l, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(prod1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(prod2, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  for (iter = 0; iter < 50; ++iter) {
    seed = seed * 1103515245u + 12345u;
    uint32_t av = ((seed >> 16) % 10000u) + 1u; /* 1 to 10000 */
    seed = seed * 1103515245u + 12345u;
    uint32_t bv = ((seed >> 16) % 10000u) + 1u;

    CHECK_EQ_I(arbint_set_u32(a, av), ARBINT_OK);
    CHECK_EQ_I(arbint_set_u32(b, bv), ARBINT_OK);

    CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_lcm(l, a, b), ARBINT_OK);

    /*  Property 1: lcm divides by both |a| and |b|.  */
    CHECK_EQ_I(arbint_tdiv_r(r, l, a), ARBINT_OK);
    CHECK(arbint_is_zero(r));
    CHECK_EQ_I(arbint_tdiv_r(r, l, b), ARBINT_OK);
    CHECK(arbint_is_zero(r));

    /*  Property 2: gcd * lcm = |a| * |b|.  */
    CHECK_EQ_I(arbint_mul(prod1, g, l), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(prod2, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(prod1, prod2), 0);

    /*  Property 3: result is positive.  */
    CHECK(l[0]._sz >= 0);
  }

  arbint_clear(r);
  arbint_clear(prod2);
  arbint_clear(prod1);
  arbint_clear(g);
  arbint_clear(l);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Large multi-limb stochastic test.  */
static void test_gcd_large_stochastic(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, factor, g, r;
  int iter;
  uint32_t seed = 99999u;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(factor, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Build base factor = 2^256.  */
  CHECK_EQ_I(arbint_set_u32(factor, 2u), ARBINT_OK);
  for (int i = 0; i < 8; ++i)
    CHECK_EQ_I(arbint_sqr(factor, factor), ARBINT_OK);

  for (iter = 0; iter < 20; ++iter) {
    seed = seed * 1103515245u + 12345u;
    uint32_t av = ((seed >> 16) % 1000u) + 1u;
    seed = seed * 1103515245u + 12345u;
    uint32_t bv = ((seed >> 16) % 1000u) + 1u;

    /*  a = av * factor, b = bv * factor.  */
    CHECK_EQ_I(arbint_mul_u32(a, factor, av), ARBINT_OK);
    CHECK_EQ_I(arbint_mul_u32(b, factor, bv), ARBINT_OK);

    CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);

    /*  gcd should divide both.  */
    CHECK_EQ_I(arbint_tdiv_r(r, a, g), ARBINT_OK);
    CHECK(arbint_is_zero(r));
    CHECK_EQ_I(arbint_tdiv_r(r, b, g), ARBINT_OK);
    CHECK(arbint_is_zero(r));

    /*  gcd(av * factor, bv * factor) >= factor (since factor divides both). */
    CHECK(arbint_cmpabs(g, factor) >= 0);
  }

  arbint_clear(r);
  arbint_clear(g);
  arbint_clear(factor);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  LCM u32 aliasing: l = a case.  */
static void test_lcm_u32_aliasing(void) {
  arbint_ctx_t ctx;
  arbint_t a;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_set_u32(a, 4u), ARBINT_OK);
  CHECK_EQ_I(arbint_lcm_u32(a, a, 6u), ARBINT_OK);
  check_u32_value(a, 12u);

  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  GCD u32 aliasing: g = a case.  */
static void test_gcd_u32_aliasing(void) {
  arbint_ctx_t ctx;
  arbint_t a;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_set_u32(a, 48u), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd_u32(a, a, 18u), ARBINT_OK);
  check_u32_value(a, 6u);

  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Edge case: gcd(1, n) = 1 and gcd(n, 1) = 1.  */
static void test_gcd_with_one(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);

  /*  Small case.  */
  CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 123456789u), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 1u);

  CHECK_EQ_I(arbint_gcd(g, b, a), ARBINT_OK);
  check_u32_value(g, 1u);

  /*  Large multi-limb case.  */
  CHECK_EQ_I(arbint_set_u32(b, 2u), ARBINT_OK);
  for (int i = 0; i < 10; ++i)
    CHECK_EQ_I(arbint_sqr(b, b), ARBINT_OK); /* b = 2^1024 */
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 1u);

  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Perfect powers: gcd(n^k, n^m) = n^min(k,m).  */
static void test_gcd_perfect_powers(void) {
  arbint_ctx_t ctx;
  arbint_t n2, n3, n5, g;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n2, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n3, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n5, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);

  /*  n = 7: compute 7^2, 7^3, 7^5.  */
  CHECK_EQ_I(arbint_set_u32(n2, 49u), ARBINT_OK);    /* 7^2 = 49 */
  CHECK_EQ_I(arbint_set_u32(n3, 343u), ARBINT_OK);   /* 7^3 = 343 */
  CHECK_EQ_I(arbint_set_u32(n5, 16807u), ARBINT_OK); /* 7^5 = 16807 */

  /*  gcd(7^2, 7^3) = 7^2.  */
  CHECK_EQ_I(arbint_gcd(g, n2, n3), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, n2), 0);

  /*  gcd(7^3, 7^5) = 7^3.  */
  CHECK_EQ_I(arbint_gcd(g, n3, n5), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, n3), 0);

  /*  gcd(7^2, 7^5) = 7^2.  */
  CHECK_EQ_I(arbint_gcd(g, n2, n5), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, n2), 0);

  arbint_clear(g);
  arbint_clear(n5);
  arbint_clear(n3);
  arbint_clear(n2);
  arbint_ctx_clear(&ctx);
}

/*  Lehmer GCD threshold boundary tests.
    ARBINT_LEHMER_THRESHOLD is 32 limbs; test at 31, 32, 33.  */
static void test_gcd_lehmer_threshold(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g, r;
  int nlimbs;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Test around Lehmer threshold (32 limbs).  */
  for (nlimbs = 31; nlimbs <= 35; ++nlimbs) {
    /*  a = 2^(LIMB_BITS * nlimbs) - 1 (exactly nlimbs limbs, all 1-bits).  */
    CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(a, a, (unsigned) (ARBINT_LIMB_BITS * nlimbs)),
               ARBINT_OK);
    CHECK_EQ_I(arbint_sub_u32(a, a, 1u), ARBINT_OK);

    /*  b = a - 2.  */
    CHECK_EQ_I(arbint_sub_u32(b, a, 2u), ARBINT_OK);

    CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);

    /*  Verify gcd divides both.  */
    CHECK_EQ_I(arbint_tdiv_r(r, a, g), ARBINT_OK);
    CHECK(arbint_is_zero(r));
    CHECK_EQ_I(arbint_tdiv_r(r, b, g), ARBINT_OK);
    CHECK(arbint_is_zero(r));

    /*  Result must be positive.  */
    CHECK(g[0]._sz > 0);
  }

  arbint_clear(r);
  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Large Fibonacci pairs test for Lehmer GCD.
    F_n and F_{n-1} are always coprime, and produce the maximum number
    of quotients (all 1s) -- worst case for Lehmer simulation.  */
static void test_gcd_lehmer_fibonacci(void) {
  arbint_ctx_t ctx;
  arbint_t f_prev, f_curr, f_next, g;
  int i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(f_prev, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(f_curr, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(f_next, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);

  /*  F_0 = 0, F_1 = 1.  */
  CHECK_EQ_I(arbint_set_u32(f_prev, 0u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(f_curr, 1u), ARBINT_OK);

  /*  Compute up to F_500 (multi-limb, well above Lehmer threshold).
      F_500 has about 104 decimal digits, roughly 346 bits (6+ limbs on 64-bit,
      but we want >32 limbs so go higher).  */
  for (i = 2; i <= 3000; ++i) {
    CHECK_EQ_I(arbint_add(f_next, f_curr, f_prev), ARBINT_OK);
    CHECK_EQ_I(arbint_set(f_prev, f_curr), ARBINT_OK);
    CHECK_EQ_I(arbint_set(f_curr, f_next), ARBINT_OK);
  }

  /*  gcd(F_3000, F_2999) = 1 (consecutive Fibonacci are coprime).  */
  CHECK_EQ_I(arbint_gcd(g, f_curr, f_prev), ARBINT_OK);
  check_u32_value(g, 1u);

  /*  gcd(F_2999, F_3000) = 1 (symmetry).  */
  CHECK_EQ_I(arbint_gcd(g, f_prev, f_curr), ARBINT_OK);
  check_u32_value(g, 1u);

  arbint_clear(g);
  arbint_clear(f_next);
  arbint_clear(f_curr);
  arbint_clear(f_prev);
  arbint_ctx_clear(&ctx);
}

/*  Test Lehmer GCD with shared large factors.
    gcd(k*a, k*b) = k * gcd(a, b) for large multi-limb k.  */
static void test_gcd_lehmer_shared_factor(void) {
  arbint_ctx_t ctx;
  arbint_t factor, a, b, g, expected;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(factor, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  /*  factor = 2^2048 (32 limbs on 64-bit, exactly at Lehmer threshold).  */
  CHECK_EQ_I(arbint_set_u32(factor, 2u), ARBINT_OK);
  for (int i = 0; i < 11; ++i)
    CHECK_EQ_I(arbint_sqr(factor, factor), ARBINT_OK);

  /*  a = 7 * factor, b = 21 * factor.
      gcd(a, b) = 7 * factor since gcd(7, 21) = 7.  */
  CHECK_EQ_I(arbint_mul_u32(a, factor, 7u), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(b, factor, 21u), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, a), 0);

  /*  a = 12 * factor, b = 8 * factor.
      gcd(a, b) = 4 * factor.  */
  CHECK_EQ_I(arbint_mul_u32(a, factor, 12u), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(b, factor, 8u), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(expected, factor, 4u), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, expected), 0);

  /*  a = 17 * factor, b = 13 * factor (coprime multipliers).
      gcd(a, b) = factor.  */
  CHECK_EQ_I(arbint_mul_u32(a, factor, 17u), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(b, factor, 13u), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, factor), 0);

  arbint_clear(expected);
  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_clear(factor);
  arbint_ctx_clear(&ctx);
}

/*  Test Lehmer GCD with very large operands (well above threshold).  */
static void test_gcd_lehmer_large(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  a = 2^4096 - 1 (64 limbs on 64-bit), b = 2^4096 - 3.
      These are coprime (consecutive odd numbers).  */
  CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 4096u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(a, a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(b, a, 2u), ARBINT_OK);

  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);

  /*  Verify gcd divides both.  */
  CHECK_EQ_I(arbint_tdiv_r(r, a, g), ARBINT_OK);
  CHECK(arbint_is_zero(r));
  CHECK_EQ_I(arbint_tdiv_r(r, b, g), ARBINT_OK);
  CHECK(arbint_is_zero(r));

  /*  Result is positive.  */
  CHECK(g[0]._sz > 0);

  /*  Test with known GCD: a = 3^100, b = 3^150.
      gcd = 3^100.  */
  CHECK_EQ_I(arbint_set_u32(a, 3u), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 100u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 3u), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(b, b, 150u), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, a), 0);

  arbint_clear(r);
  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test Lehmer GCD with quotient sequence edge cases.
    Powers of 2 +/- small values can trigger matrix overflow detection.  */
static void test_gcd_lehmer_quotient_edge(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g, expected;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  /*  Test: a = 2^2048, b = 7 (one operand >> other).
      gcd = 1 since 2^2048 is only divisible by 2.  */
  CHECK_EQ_I(arbint_set_u32(a, 2u), ARBINT_OK);
  for (int i = 0; i < 11; ++i)
    CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 7u), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 1u);

  /*  Test: a = 2^2048 + 1, b = 2^2048 - 1.
      gcd(2^n + 1, 2^n - 1) = gcd(2, 2^n - 1) = 1 (both are odd when n > 0). */
  CHECK_EQ_I(arbint_set_u32(a, 2u), ARBINT_OK);
  for (int i = 0; i < 11; ++i)
    CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(b, a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_add_u32(a, a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 1u);

  /*  Test: a = 2^2048 - 1, b = 2^1024 - 1.
      2^2048 - 1 = (2^1024 - 1)(2^1024 + 1), so gcd = 2^1024 - 1.  */
  CHECK_EQ_I(arbint_set_u32(a, 2u), ARBINT_OK);
  for (int i = 0; i < 11; ++i)
    CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(a, a, 1u), ARBINT_OK); /* 2^2048 - 1 */
  CHECK_EQ_I(arbint_set_u32(b, 2u), ARBINT_OK);
  for (int i = 0; i < 10; ++i)
    CHECK_EQ_I(arbint_sqr(b, b), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(b, b, 1u), ARBINT_OK); /* 2^1024 - 1 */
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, b), 0);

  arbint_clear(expected);
  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test Lehmer GCD with high-quotient case (regression test).
    When Lehmer step returns count=0, the fallback must use division
    rather than O(q) repeated subtraction. This test would take ~2 seconds
    with subtraction-based fallback but completes instantly with division.  */
static void test_gcd_lehmer_high_quotient(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g, expected;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  /*  Create a case with very high quotient: a = b * q + small_remainder.
      If the fallback uses subtraction, this takes O(q) time.
      With proper division, it should be O(1).

      b = 2^2048 (32 limbs, at Lehmer threshold)
      a = b * 100000000 + 1  (quotient = 10^8)
      gcd(a, b) = gcd(b, 1) = 1  */
  CHECK_EQ_I(arbint_set_u32(b, 2u), ARBINT_OK);
  for (int i = 0; i < 11; ++i)
    CHECK_EQ_I(arbint_sqr(b, b), ARBINT_OK); /* b = 2^2048 */

  CHECK_EQ_I(arbint_mul_u32(a, b, 100000000u), ARBINT_OK); /* a = b * 10^8 */
  CHECK_EQ_I(arbint_add_u32(a, a, 1u), ARBINT_OK); /* a = b * 10^8 + 1 */

  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 1u);

  /*  Another high-quotient case with non-trivial GCD.
      a = 7 * b, so gcd(a, b) = b.
      Then test gcd(a + small, b) where small < b to force reduction.  */
  CHECK_EQ_I(arbint_mul_u32(a, b, 7u), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, b), 0);

  /*  a = 1000000 * b + 3 (high quotient, gcd = gcd(b, 3)).  */
  CHECK_EQ_I(arbint_mul_u32(a, b, 1000000u), ARBINT_OK);
  CHECK_EQ_I(arbint_add_u32(a, a, 3u), ARBINT_OK);
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  /*  gcd(2^2048, 3) = 1 since 2^2048 mod 3 = 1 (2 = -1 mod 3, 2^even = 1).  */
  check_u32_value(g, 1u);

  arbint_clear(expected);
  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test Lehmer GCD stochastic with large operands.  */
static void test_gcd_lehmer_stochastic(void) {
  arbint_ctx_t ctx;
  arbint_t base, a, b, g, r;
  int iter;
  uint32_t seed = 77777u;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(base, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  base = 2^2048 (32 limbs on 64-bit).  */
  CHECK_EQ_I(arbint_set_u32(base, 2u), ARBINT_OK);
  for (int i = 0; i < 11; ++i)
    CHECK_EQ_I(arbint_sqr(base, base), ARBINT_OK);

  for (iter = 0; iter < 20; ++iter) {
    seed = seed * 1103515245u + 12345u;
    uint32_t av = ((seed >> 16) % 10000u) + 1u;
    seed = seed * 1103515245u + 12345u;
    uint32_t bv = ((seed >> 16) % 10000u) + 1u;

    /*  a = av * base + offset, b = bv * base + offset2.  */
    CHECK_EQ_I(arbint_mul_u32(a, base, av), ARBINT_OK);
    CHECK_EQ_I(arbint_add_u32(a, a, (seed >> 20) % 1000u), ARBINT_OK);
    CHECK_EQ_I(arbint_mul_u32(b, base, bv), ARBINT_OK);
    CHECK_EQ_I(arbint_add_u32(b, b, (seed >> 24) % 1000u), ARBINT_OK);

    CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);

    /*  gcd must divide both.  */
    CHECK_EQ_I(arbint_tdiv_r(r, a, g), ARBINT_OK);
    CHECK(arbint_is_zero(r));
    CHECK_EQ_I(arbint_tdiv_r(r, b, g), ARBINT_OK);
    CHECK(arbint_is_zero(r));

    /*  Result must be positive.  */
    CHECK(g[0]._sz > 0);
  }

  arbint_clear(r);
  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_clear(base);
  arbint_ctx_clear(&ctx);
}

/*  Main entry point.  */

int main(void) {
  ARBINT_TEST_START();
  /*  Basic functionality.  */
  test_gcd_basic();
  test_gcd_zeros();
  test_gcd_negative();
  test_gcd_symmetry();
  test_gcd_aliasing();
  test_gcd_u32();
  test_gcd_u32_aliasing();
  test_gcd_with_one();

  /*  LCM functionality.  */
  test_lcm_basic();
  test_lcm_zeros();
  test_lcm_negative();
  test_lcm_aliasing();
  test_lcm_u32();
  test_lcm_u32_aliasing();

  /*  Multi-limb tests.  */
  test_gcd_powers_of_two();
  test_gcd_shared_factors();
  test_gcd_perfect_powers();

  /*  Mathematical properties.  */
  test_gcd_divisibility();
  test_gcd_lcm_identity();
  test_lcm_exactness();
  test_gcd_scaling();
  test_gcd_primes();

  /*  Algorithm-specific tests.  */
  test_gcd_size_disparate();
  test_gcd_large_binary();
  test_gcd_binary_mod_opt();
  test_gcd_threshold_boundary();
  test_gcd_fibonacci();
  test_gcd_consecutive();

  /*  Edge cases and validation.  */
  test_gcd_normalization();
  test_gcd_large_single_limb();
  test_gcd_null_pointers();

  /*  Stochastic tests.  */
  test_gcd_stochastic();
  test_lcm_stochastic();
  test_gcd_large_stochastic();

  /*  Lehmer GCD specific tests.  */
  test_gcd_lehmer_threshold();
  test_gcd_lehmer_fibonacci();
  test_gcd_lehmer_shared_factor();
  test_gcd_lehmer_large();
  test_gcd_lehmer_quotient_edge();
  test_gcd_lehmer_high_quotient();
  test_gcd_lehmer_stochastic();

  ARBINT_TEST_FINISH("test_gcd");
}
