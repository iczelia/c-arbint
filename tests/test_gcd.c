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
    CHECK_EQ_I(arbint_sqr(big, big), ARBINT_OK);  /* 2^2048 */

  CHECK_EQ_I(arbint_set_u32(expected, 2u), ARBINT_OK);
  for (int i = 0; i < 5; ++i)
    CHECK_EQ_I(arbint_sqr(expected, expected), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(expected, expected), ARBINT_OK);  /* 2^320 */

  CHECK_EQ_I(arbint_mul_u32(medium, expected, 7u), ARBINT_OK);  /* 7 * 2^320 */

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
  CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);  /* a = 2^32 */
  CHECK_EQ_I(arbint_set(b, a), ARBINT_OK);  /* b = 2^32 */

  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, a), 0);  /* gcd should equal a */

  /*  lcm(2^32, 2^32) = 2^32 (must not fail with EZERO).  */
  CHECK_EQ_I(arbint_lcm(l, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(l, a), 0);  /* lcm should equal a */

  /*  gcd(2^33, 2^32) = 2^32.  */
  CHECK_EQ_I(arbint_mul_u32(a, a, 2u), ARBINT_OK);  /* a = 2^33 */
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, b), 0);  /* gcd should equal b = 2^32 */

  /*  gcd with large coprime single-limb values.  */
  CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 40u), ARBINT_OK);  /* a = 2^40 */
  CHECK_EQ_I(arbint_set_u32(b, 3u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(b, b, 35u), ARBINT_OK);  /* b = 3 * 2^35 */
  CHECK_EQ_I(arbint_gcd(g, a, b), ARBINT_OK);
  /* gcd(2^40, 3*2^35) = 2^35 */
  CHECK_EQ_I(arbint_set_u32(l, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(l, l, 35u), ARBINT_OK);  /* expected = 2^35 */
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

/*  Main entry point.  */

int main(void) {
  test_gcd_basic();
  test_gcd_zeros();
  test_gcd_negative();
  test_gcd_symmetry();
  test_gcd_aliasing();
  test_gcd_u32();
  test_lcm_basic();
  test_lcm_zeros();
  test_lcm_negative();
  test_lcm_aliasing();
  test_lcm_u32();
  test_gcd_powers_of_two();
  test_gcd_shared_factors();
  test_gcd_divisibility();
  test_gcd_lcm_identity();
  test_lcm_exactness();
  test_gcd_size_disparate();
  test_gcd_large_binary();
  test_gcd_binary_mod_opt();
  test_gcd_normalization();
  test_gcd_large_single_limb();
  test_gcd_null_pointers();

  ARBINT_TEST_FINISH("test_gcd");
}
