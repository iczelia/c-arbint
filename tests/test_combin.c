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

/*  Factorial tests.  */

static void test_fac_null(void) {
  CHECK_EQ_I(arbint_fac_u32(NULL, 5), ARBINT_EINVAL);
}

static void test_fac_basic(void) {
  arbint_ctx_t ctx;
  arbint_t r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  0! = 1.  */
  CHECK_EQ_I(arbint_fac_u32(r, 0), ARBINT_OK);
  check_i32_value(r, 1);

  /*  1! = 1.  */
  CHECK_EQ_I(arbint_fac_u32(r, 1), ARBINT_OK);
  check_i32_value(r, 1);

  /*  5! = 120.  */
  CHECK_EQ_I(arbint_fac_u32(r, 5), ARBINT_OK);
  check_i32_value(r, 120);

  /*  10! = 3628800.  */
  CHECK_EQ_I(arbint_fac_u32(r, 10), ARBINT_OK);
  check_i32_value(r, 3628800);

  arbint_clear(r);
  arbint_ctx_clear(&ctx);
}

static void test_fac_large(void) {
  arbint_ctx_t ctx;
  arbint_t r, expected, tmp;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(tmp, &ctx), ARBINT_OK);

  /*  20! = 2432902008176640000.
      Build via identity: 20! = 20 * 19!  */
  CHECK_EQ_I(arbint_fac_u32(r, 20), ARBINT_OK);
  CHECK_EQ_I(arbint_fac_u32(expected, 19), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(expected, expected, 20), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, expected), 0);

  /*  Verify 20! via a few multiplications from 10!.  */
  CHECK_EQ_I(arbint_fac_u32(expected, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(expected, expected, 11), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(expected, expected, 12), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(expected, expected, 13), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(expected, expected, 14), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(expected, expected, 15), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(expected, expected, 16), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(expected, expected, 17), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(expected, expected, 18), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(expected, expected, 19), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(expected, expected, 20), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, expected), 0);

  arbint_clear(tmp);
  arbint_clear(expected);
  arbint_clear(r);
  arbint_ctx_clear(&ctx);
}

/*  Binomial coefficient tests.  */

static void test_bin_null(void) {
  CHECK_EQ_I(arbint_bin_u32u32(NULL, 5, 2), ARBINT_EINVAL);
}

static void test_bin_edge(void) {
  arbint_ctx_t ctx;
  arbint_t r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  C(n, 0) = 1.  */
  CHECK_EQ_I(arbint_bin_u32u32(r, 5, 0), ARBINT_OK);
  check_i32_value(r, 1);
  CHECK_EQ_I(arbint_bin_u32u32(r, 100, 0), ARBINT_OK);
  check_i32_value(r, 1);

  /*  C(n, n) = 1.  */
  CHECK_EQ_I(arbint_bin_u32u32(r, 5, 5), ARBINT_OK);
  check_i32_value(r, 1);
  CHECK_EQ_I(arbint_bin_u32u32(r, 100, 100), ARBINT_OK);
  check_i32_value(r, 1);

  /*  C(n, k) = 0 when k > n.  */
  CHECK_EQ_I(arbint_bin_u32u32(r, 5, 6), ARBINT_OK);
  check_i32_value(r, 0);
  CHECK_EQ_I(arbint_bin_u32u32(r, 0, 1), ARBINT_OK);
  check_i32_value(r, 0);

  /*  C(0, 0) = 1.  */
  CHECK_EQ_I(arbint_bin_u32u32(r, 0, 0), ARBINT_OK);
  check_i32_value(r, 1);

  arbint_clear(r);
  arbint_ctx_clear(&ctx);
}

static void test_bin_basic(void) {
  arbint_ctx_t ctx;
  arbint_t r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  C(5, 2) = 10.  */
  CHECK_EQ_I(arbint_bin_u32u32(r, 5, 2), ARBINT_OK);
  check_i32_value(r, 10);

  /*  C(10, 5) = 252.  */
  CHECK_EQ_I(arbint_bin_u32u32(r, 10, 5), ARBINT_OK);
  check_i32_value(r, 252);

  /*  C(20, 10) = 184756.  */
  CHECK_EQ_I(arbint_bin_u32u32(r, 20, 10), ARBINT_OK);
  check_i32_value(r, 184756);

  /*  C(10, 3) = 120.  */
  CHECK_EQ_I(arbint_bin_u32u32(r, 10, 3), ARBINT_OK);
  check_i32_value(r, 120);

  arbint_clear(r);
  arbint_ctx_clear(&ctx);
}

static void test_bin_symmetry(void) {
  arbint_ctx_t ctx;
  arbint_t r1, r2;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r2, &ctx), ARBINT_OK);

  /*  C(100, 3) = C(100, 97).  */
  CHECK_EQ_I(arbint_bin_u32u32(r1, 100, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_bin_u32u32(r2, 100, 97), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r1, r2), 0);

  /*  C(50, 20) = C(50, 30).  */
  CHECK_EQ_I(arbint_bin_u32u32(r1, 50, 20), ARBINT_OK);
  CHECK_EQ_I(arbint_bin_u32u32(r2, 50, 30), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r1, r2), 0);

  arbint_clear(r2);
  arbint_clear(r1);
  arbint_ctx_clear(&ctx);
}

/*  Fibonacci tests.  */

static void test_fib_null(void) {
  CHECK_EQ_I(arbint_fib_u32(NULL, 10), ARBINT_EINVAL);
}

static void test_fib_basic(void) {
  arbint_ctx_t ctx;
  arbint_t r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  F(0) = 0.  */
  CHECK_EQ_I(arbint_fib_u32(r, 0), ARBINT_OK);
  check_i32_value(r, 0);

  /*  F(1) = 1.  */
  CHECK_EQ_I(arbint_fib_u32(r, 1), ARBINT_OK);
  check_i32_value(r, 1);

  /*  F(2) = 1.  */
  CHECK_EQ_I(arbint_fib_u32(r, 2), ARBINT_OK);
  check_i32_value(r, 1);

  /*  F(3) = 2.  */
  CHECK_EQ_I(arbint_fib_u32(r, 3), ARBINT_OK);
  check_i32_value(r, 2);

  /*  F(10) = 55.  */
  CHECK_EQ_I(arbint_fib_u32(r, 10), ARBINT_OK);
  check_i32_value(r, 55);

  /*  F(20) = 6765.  */
  CHECK_EQ_I(arbint_fib_u32(r, 20), ARBINT_OK);
  check_i32_value(r, 6765);

  arbint_clear(r);
  arbint_ctx_clear(&ctx);
}

static void test_fib_identity(void) {
  arbint_ctx_t ctx;
  arbint_t fn, fn1, fn_sq, fn1_sq, f2n1, lhs;
  uint32_t n;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(fn, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(fn1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(fn_sq, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(fn1_sq, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(f2n1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(lhs, &ctx), ARBINT_OK);

  /*  Verify identity: F(n)^2 + F(n+1)^2 = F(2n+1) for several n.  */
  for (n = 0; n <= 20 && g_failures == 0; ++n) {
    CHECK_EQ_I(arbint_fib_u32(fn, n), ARBINT_OK);
    CHECK_EQ_I(arbint_fib_u32(fn1, n + 1), ARBINT_OK);
    CHECK_EQ_I(arbint_fib_u32(f2n1, 2 * n + 1), ARBINT_OK);

    CHECK_EQ_I(arbint_sqr(fn_sq, fn), ARBINT_OK);
    CHECK_EQ_I(arbint_sqr(fn1_sq, fn1), ARBINT_OK);
    CHECK_EQ_I(arbint_add(lhs, fn_sq, fn1_sq), ARBINT_OK);

    CHECK_EQ_I(arbint_cmp(lhs, f2n1), 0);
  }

  arbint_clear(lhs);
  arbint_clear(f2n1);
  arbint_clear(fn1_sq);
  arbint_clear(fn_sq);
  arbint_clear(fn1);
  arbint_clear(fn);
  arbint_ctx_clear(&ctx);
}

static void test_fib_large(void) {
  arbint_ctx_t ctx;
  arbint_t f50, f49, f48, expected;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(f50, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(f49, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(f48, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  /*  Verify F(50) = F(49) + F(48).  */
  CHECK_EQ_I(arbint_fib_u32(f50, 50), ARBINT_OK);
  CHECK_EQ_I(arbint_fib_u32(f49, 49), ARBINT_OK);
  CHECK_EQ_I(arbint_fib_u32(f48, 48), ARBINT_OK);
  CHECK_EQ_I(arbint_add(expected, f49, f48), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(f50, expected), 0);

  arbint_clear(expected);
  arbint_clear(f48);
  arbint_clear(f49);
  arbint_clear(f50);
  arbint_ctx_clear(&ctx);
}

static void test_fib_high_n(void) {
  arbint_ctx_t ctx;
  arbint_t fn, fn1, fn2, expected;
  uint32_t n;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(fn, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(fn1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(fn2, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  /*  Test a few high-bit values to check for overflow in MSB search.
      Verify recurrence: F(n) = F(n-1) + F(n-2).  */
  n = 1000;
  CHECK_EQ_I(arbint_fib_u32(fn, n), ARBINT_OK);
  CHECK_EQ_I(arbint_fib_u32(fn1, n - 1), ARBINT_OK);
  CHECK_EQ_I(arbint_fib_u32(fn2, n - 2), ARBINT_OK);
  CHECK_EQ_I(arbint_add(expected, fn1, fn2), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(fn, expected), 0);

  /*  Test near 2^16 boundary.  */
  n = 65536u;
  CHECK_EQ_I(arbint_fib_u32(fn, n), ARBINT_OK);
  CHECK_EQ_I(arbint_fib_u32(fn1, n - 1), ARBINT_OK);
  CHECK_EQ_I(arbint_fib_u32(fn2, n - 2), ARBINT_OK);
  CHECK_EQ_I(arbint_add(expected, fn1, fn2), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(fn, expected), 0);

  arbint_clear(expected);
  arbint_clear(fn2);
  arbint_clear(fn1);
  arbint_clear(fn);
  arbint_ctx_clear(&ctx);
}

/*  Lucas tests.  */

static void test_lucas_null(void) {
  CHECK_EQ_I(arbint_lucas_u32(NULL, 10), ARBINT_EINVAL);
}

static void test_lucas_basic(void) {
  arbint_ctx_t ctx;
  arbint_t r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  L(0) = 2.  */
  CHECK_EQ_I(arbint_lucas_u32(r, 0), ARBINT_OK);
  check_i32_value(r, 2);

  /*  L(1) = 1.  */
  CHECK_EQ_I(arbint_lucas_u32(r, 1), ARBINT_OK);
  check_i32_value(r, 1);

  /*  L(2) = 3.  */
  CHECK_EQ_I(arbint_lucas_u32(r, 2), ARBINT_OK);
  check_i32_value(r, 3);

  /*  L(3) = 4.  */
  CHECK_EQ_I(arbint_lucas_u32(r, 3), ARBINT_OK);
  check_i32_value(r, 4);

  /*  L(10) = 123.  */
  CHECK_EQ_I(arbint_lucas_u32(r, 10), ARBINT_OK);
  check_i32_value(r, 123);

  /*  L(20) = 15127.  */
  CHECK_EQ_I(arbint_lucas_u32(r, 20), ARBINT_OK);
  check_i32_value(r, 15127);

  arbint_clear(r);
  arbint_ctx_clear(&ctx);
}

static void test_lucas_identity(void) {
  arbint_ctx_t ctx;
  arbint_t ln, fn_m1, fn_p1, expected;
  uint32_t n;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(ln, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(fn_m1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(fn_p1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  /*  Verify identity: L(n) = F(n-1) + F(n+1) for n >= 1.  */
  for (n = 1; n <= 20 && g_failures == 0; ++n) {
    CHECK_EQ_I(arbint_lucas_u32(ln, n), ARBINT_OK);
    CHECK_EQ_I(arbint_fib_u32(fn_m1, n - 1), ARBINT_OK);
    CHECK_EQ_I(arbint_fib_u32(fn_p1, n + 1), ARBINT_OK);
    CHECK_EQ_I(arbint_add(expected, fn_m1, fn_p1), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(ln, expected), 0);
  }

  arbint_clear(expected);
  arbint_clear(fn_p1);
  arbint_clear(fn_m1);
  arbint_clear(ln);
  arbint_ctx_clear(&ctx);
}

static void test_lucas_consistency(void) {
  arbint_ctx_t ctx;
  arbint_t ln, fn, fn1, expected;
  uint32_t n;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(ln, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(fn, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(fn1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  /*  Verify: L(n) = 2*F(n+1) - F(n).  */
  for (n = 0; n <= 20 && g_failures == 0; ++n) {
    CHECK_EQ_I(arbint_lucas_u32(ln, n), ARBINT_OK);
    CHECK_EQ_I(arbint_fib_u32(fn, n), ARBINT_OK);
    CHECK_EQ_I(arbint_fib_u32(fn1, n + 1), ARBINT_OK);
    CHECK_EQ_I(arbint_mul_i32(expected, fn1, 2), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(expected, expected, fn), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(ln, expected), 0);
  }

  arbint_clear(expected);
  arbint_clear(fn1);
  arbint_clear(fn);
  arbint_clear(ln);
  arbint_ctx_clear(&ctx);
}

/*  is_square tests.  */

static void test_is_square_null(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  int out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(a, 4), ARBINT_OK);

  CHECK_EQ_I(arbint_is_square(NULL, &out), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_is_square(a, NULL), ARBINT_EINVAL);

  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_is_square_negative(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  int out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_set_i32(a, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(a, &out), ARBINT_EDOM);

  CHECK_EQ_I(arbint_set_i32(a, -100), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(a, &out), ARBINT_EDOM);

  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_is_square_basic(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  int out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);

  /*  Perfect squares.  */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(a, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(a, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(a, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(a, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(a, 9), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(a, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(a, 16), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(a, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(a, 25), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(a, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(a, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_u32(a, 1000000), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(a, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  Non-squares.  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(a, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(a, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(a, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  CHECK_EQ_I(arbint_set_i32(a, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(a, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  CHECK_EQ_I(arbint_set_i32(a, 8), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(a, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  CHECK_EQ_I(arbint_set_i32(a, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(a, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_is_square_large(void) {
  arbint_ctx_t ctx;
  arbint_t x, sq;
  int out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(sq, &ctx), ARBINT_OK);

  /*  x = 2^100, sq = x^2 is a perfect square.  */
  CHECK_EQ_I(arbint_set_i32(x, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(x, x, 100u), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(sq, x), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(sq, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  sq + 1 is not a perfect square.  */
  CHECK_EQ_I(arbint_add_i32(sq, sq, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(sq, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  /*  x = 2^512 + 7, sq = x^2 is a perfect square.  */
  CHECK_EQ_I(arbint_set_i32(x, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(x, x, 512u), ARBINT_OK);
  CHECK_EQ_I(arbint_add_i32(x, x, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(sq, x), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(sq, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  sq - 1 is not a perfect square.  */
  CHECK_EQ_I(arbint_sub_i32(sq, sq, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_is_square(sq, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  arbint_clear(sq);
  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

/*  isprime tests.  */

static void test_isprime_basic(void) {
  arbint_ctx_t ctx;
  arbint_t n, p, c;
  int out = -1;
  size_t i;
  static const struct {
    int32_t n;
    int is_prime;
  } cases[] = {{-17, 0}, {-1, 0}, {0, 0}, {1, 0}, {2, 1}, {3, 1}, {4, 0},
               {5, 1},  {9, 0},  {17, 1}, {19, 1}, {21, 0}, {97, 1}};
  static const uint32_t carmichael[] = {561u, 1105u, 1729u, 2465u, 6601u};

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(p, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_isprime(NULL, 0, &out), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_isprime(n, 0, NULL), ARBINT_EINVAL);

  for (i = 0u; i < sizeof(cases) / sizeof(cases[0]); ++i) {
    CHECK_EQ_I(arbint_set_i32(n, cases[i].n), ARBINT_OK);
    CHECK_EQ_I(arbint_isprime(n, 0, &out), ARBINT_OK);
    CHECK_EQ_I(out, cases[i].is_prime);
  }

  for (i = 0u; i < sizeof(carmichael) / sizeof(carmichael[0]); ++i) {
    CHECK_EQ_I(arbint_set_u32(n, carmichael[i]), ARBINT_OK);
    CHECK_EQ_I(arbint_isprime(n, 8, &out), ARBINT_OK);
    CHECK_EQ_I(out, 0);
  }

  CHECK_EQ_I(arbint_set_u32(n, 2147483647u), ARBINT_OK);
  CHECK_EQ_I(arbint_isprime(n, 12, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_u32(n, 2147483645u), ARBINT_OK);
  CHECK_EQ_I(arbint_isprime(n, 12, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  /*  Large prime: 2^127 - 1 (Mersenne prime).  */
  CHECK_EQ_I(arbint_set_i32(p, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(p, p, 127u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_i32(p, p, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_isprime(p, 24, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_isprime(p, -1, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  Large composite derived from that prime.  */
  CHECK_EQ_I(arbint_mul_u32(c, p, 17u), ARBINT_OK);
  CHECK_EQ_I(arbint_isprime(c, 24, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  arbint_clear(c);
  arbint_clear(p);
  arbint_clear(n);
  arbint_ctx_clear(&ctx);
}

/*  Main.  */

int main(void) {
  /*  Factorial tests.  */
  test_fac_null();
  test_fac_basic();
  test_fac_large();

  /*  Binomial tests.  */
  test_bin_null();
  test_bin_edge();
  test_bin_basic();
  test_bin_symmetry();

  /*  Fibonacci tests.  */
  test_fib_null();
  test_fib_basic();
  test_fib_identity();
  test_fib_large();
  test_fib_high_n();

  /*  Lucas tests.  */
  test_lucas_null();
  test_lucas_basic();
  test_lucas_identity();
  test_lucas_consistency();

  /*  is_square tests.  */
  test_is_square_null();
  test_is_square_negative();
  test_is_square_basic();
  test_is_square_large();
  test_isprime_basic();

  ARBINT_TEST_FINISH("test_combin");
}
