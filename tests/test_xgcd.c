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

/*  Helper: verify Bezout identity a*x + b*y == g.  */
static void verify_bezout(const arbint_t a, const arbint_t b, const arbint_t g,
                          const arbint_t x, const arbint_t y, arbint_t tmp1,
                          arbint_t tmp2, const char * desc) {
  arbint_err_t rc;

  /*  tmp1 = a * x.  */
  rc = arbint_mul(tmp1, a, x);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "FAIL %s: a*x failed with %d\n", desc, (int) rc);
    ++g_failures;
    return;
  }

  /*  tmp2 = b * y.  */
  rc = arbint_mul(tmp2, b, y);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "FAIL %s: b*y failed with %d\n", desc, (int) rc);
    ++g_failures;
    return;
  }

  /*  tmp1 = a*x + b*y.  */
  rc = arbint_add(tmp1, tmp1, tmp2);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "FAIL %s: a*x + b*y failed with %d\n", desc, (int) rc);
    ++g_failures;
    return;
  }

  /*  Verify tmp1 == g.  */
  if (!arbint_eq(tmp1, g)) {
    fprintf(stderr, "FAIL %s: Bezout identity a*x + b*y != g\n", desc);
    ++g_failures;
  }
}

/*  Helper: verify g == gcd(a, b) via arbint_gcd.  */
static void verify_gcd(const arbint_t a, const arbint_t b, const arbint_t g,
                       arbint_t expected_g, const char * desc) {
  arbint_err_t rc;

  rc = arbint_gcd(expected_g, a, b);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "FAIL %s: arbint_gcd failed with %d\n", desc, (int) rc);
    ++g_failures;
    return;
  }

  if (!arbint_eq(g, expected_g)) {
    fprintf(stderr, "FAIL %s: g != gcd(a, b)\n", desc);
    ++g_failures;
  }
}

/*  Helper: test xgcd with i32 values and verify results.  */
static void test_xgcd_i32(arbint_t a, arbint_t b, arbint_t g, arbint_t x,
                          arbint_t y, arbint_t tmp1, arbint_t tmp2, int32_t av,
                          int32_t bv, const char * desc) {
  arbint_err_t rc;

  rc = arbint_set_i32(a, av);
  CHECK_EQ_I(rc, ARBINT_OK);
  rc = arbint_set_i32(b, bv);
  CHECK_EQ_I(rc, ARBINT_OK);

  rc = arbint_xgcd(g, x, y, a, b);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "FAIL %s: arbint_xgcd returned %d\n", desc, (int) rc);
    ++g_failures;
    return;
  }

  /*  g must be non-negative.  */
  if (g[0]._sz < 0) {
    fprintf(stderr, "FAIL %s: g is negative\n", desc);
    ++g_failures;
    return;
  }

  verify_bezout(a, b, g, x, y, tmp1, tmp2, desc);
  verify_gcd(a, b, g, tmp1, desc);
}

int main(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g, x, y, tmp1, tmp2;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(g, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(y, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(tmp1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(tmp2, &ctx), ARBINT_OK);

  /*  Test NULL inputs.  */
  {
    arbint_err_t rc;
    CHECK_EQ_I(arbint_set_i32(a, 12), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, 8), ARBINT_OK);

    rc = arbint_xgcd(NULL, x, y, a, b);
    CHECK_EQ_I(rc, ARBINT_EINVAL);
    rc = arbint_xgcd(g, NULL, y, a, b);
    CHECK_EQ_I(rc, ARBINT_EINVAL);
    rc = arbint_xgcd(g, x, NULL, a, b);
    CHECK_EQ_I(rc, ARBINT_EINVAL);
    rc = arbint_xgcd(g, x, y, NULL, b);
    CHECK_EQ_I(rc, ARBINT_EINVAL);
    rc = arbint_xgcd(g, x, y, a, NULL);
    CHECK_EQ_I(rc, ARBINT_EINVAL);
  }

  /*  Test gcd(0, 0) = 0.  */
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 0, 0, "gcd(0, 0)");
  check_i32_value(g, 0);
  check_i32_value(x, 0);
  check_i32_value(y, 0);

  /*  Test gcd(0, b) cases.  */
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 0, 5, "gcd(0, 5)");
  check_i32_value(g, 5);
  check_i32_value(x, 0);
  check_i32_value(y, 1);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 0, -5, "gcd(0, -5)");
  check_i32_value(g, 5);
  check_i32_value(x, 0);
  check_i32_value(y, -1);

  /*  Test gcd(a, 0) cases.  */
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 7, 0, "gcd(7, 0)");
  check_i32_value(g, 7);
  check_i32_value(x, 1);
  check_i32_value(y, 0);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, -7, 0, "gcd(-7, 0)");
  check_i32_value(g, 7);
  check_i32_value(x, -1);
  check_i32_value(y, 0);

  /*  Test coprime pairs.  */
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 17, 13, "gcd(17, 13)");
  check_i32_value(g, 1);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 100, 37, "gcd(100, 37)");
  check_i32_value(g, 1);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 35, 64, "gcd(35, 64)");
  check_i32_value(g, 1);

  /*  Test non-coprime pairs.  */
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 12, 8, "gcd(12, 8)");
  check_i32_value(g, 4);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 48, 18, "gcd(48, 18)");
  check_i32_value(g, 6);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 120, 45, "gcd(120, 45)");
  check_i32_value(g, 15);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 84, 126, "gcd(84, 126)");
  check_i32_value(g, 42);

  /*  Test all sign combinations for gcd(12, 8).  */
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 12, 8, "gcd(+12, +8)");
  check_i32_value(g, 4);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, -12, 8, "gcd(-12, +8)");
  check_i32_value(g, 4);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 12, -8, "gcd(+12, -8)");
  check_i32_value(g, 4);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, -12, -8, "gcd(-12, -8)");
  check_i32_value(g, 4);

  /*  Test output aliasing: g == x (should get x's value).  */
  {
    arbint_err_t rc;
    CHECK_EQ_I(arbint_set_i32(a, 12), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, 8), ARBINT_OK);

    /*  Use g as both g and x output.  */
    rc = arbint_xgcd(g, g, y, a, b);
    CHECK_EQ_I(rc, ARBINT_OK);

    /*  Per aliasing contract, g should have x's value (x written after g).
        Verify Bezout: 12*x + 8*y should equal gcd(12,8)=4.  */
    CHECK_EQ_I(arbint_gcd(tmp1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(tmp2, a, g), ARBINT_OK); /*  a * x (x is in g).  */
    CHECK_EQ_I(arbint_mul(x, b, y), ARBINT_OK);    /*  b * y.  */
    CHECK_EQ_I(arbint_add(tmp2, tmp2, x), ARBINT_OK);
    CHECK(arbint_eq(tmp2, tmp1)); /*  a*x + b*y == gcd.  */
  }

  /*  Test output aliasing: g == y (should get y's value).  */
  {
    arbint_err_t rc;
    CHECK_EQ_I(arbint_set_i32(a, 12), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, 8), ARBINT_OK);

    rc = arbint_xgcd(g, x, g, a, b);
    CHECK_EQ_I(rc, ARBINT_OK);

    /*  g should have y's value.  */
    CHECK_EQ_I(arbint_gcd(tmp1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(tmp2, a, x), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(y, b, g), ARBINT_OK); /*  b * y (y is in g).  */
    CHECK_EQ_I(arbint_add(tmp2, tmp2, y), ARBINT_OK);
    CHECK(arbint_eq(tmp2, tmp1));
  }

  /*  Test output aliasing: x == y (should get y's value).  */
  {
    arbint_err_t rc;
    CHECK_EQ_I(arbint_set_i32(a, 12), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, 8), ARBINT_OK);

    rc = arbint_xgcd(g, x, x, a, b);
    CHECK_EQ_I(rc, ARBINT_OK);

    /*  x should have y's value. Verify g == 4.  */
    check_i32_value(g, 4);
    /*  Bezout: 12*x + 8*y = 4 where y is stored in x.  */
  }

  /*  Test output aliasing: g == x == y (should get y's value).  */
  {
    arbint_err_t rc;
    CHECK_EQ_I(arbint_set_i32(a, 12), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, 8), ARBINT_OK);

    rc = arbint_xgcd(g, g, g, a, b);
    CHECK_EQ_I(rc, ARBINT_OK);

    /*  g should have y's value (last write wins).  */
  }

  /*  Test input-output aliasing: g == a.  */
  {
    arbint_err_t rc;
    CHECK_EQ_I(arbint_set_i32(g, 48), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, 18), ARBINT_OK);

    rc = arbint_xgcd(g, x, y, g, b);
    CHECK_EQ_I(rc, ARBINT_OK);

    check_i32_value(g, 6);
    /*  Verify Bezout with original a=48.  */
    CHECK_EQ_I(arbint_set_i32(a, 48), ARBINT_OK);
    verify_bezout(a, b, g, x, y, tmp1, tmp2, "g==a aliasing");
  }

  /*  Test input-output aliasing: g == b.  */
  {
    arbint_err_t rc;
    CHECK_EQ_I(arbint_set_i32(a, 48), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(g, 18), ARBINT_OK);

    rc = arbint_xgcd(g, x, y, a, g);
    CHECK_EQ_I(rc, ARBINT_OK);

    check_i32_value(g, 6);
    CHECK_EQ_I(arbint_set_i32(b, 18), ARBINT_OK);
    verify_bezout(a, b, g, x, y, tmp1, tmp2, "g==b aliasing");
  }

  /*  Test input-output aliasing: x == a.  */
  {
    arbint_err_t rc;
    CHECK_EQ_I(arbint_set_i32(x, 48), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, 18), ARBINT_OK);

    rc = arbint_xgcd(g, x, y, x, b);
    CHECK_EQ_I(rc, ARBINT_OK);

    check_i32_value(g, 6);
    CHECK_EQ_I(arbint_set_i32(a, 48), ARBINT_OK);
    verify_bezout(a, b, g, x, y, tmp1, tmp2, "x==a aliasing");
  }

  /*  Test input-output aliasing: y == b.  */
  {
    arbint_err_t rc;
    CHECK_EQ_I(arbint_set_i32(a, 48), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(y, 18), ARBINT_OK);

    rc = arbint_xgcd(g, x, y, a, y);
    CHECK_EQ_I(rc, ARBINT_OK);

    check_i32_value(g, 6);
    CHECK_EQ_I(arbint_set_i32(b, 18), ARBINT_OK);
    verify_bezout(a, b, g, x, y, tmp1, tmp2, "y==b aliasing");
  }

  /*  Test zero-path aliasing: gcd(0, b) with g==x.  */
  {
    arbint_err_t rc;
    CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, 17), ARBINT_OK);

    rc = arbint_xgcd(g, g, y, a, b);
    CHECK_EQ_I(rc, ARBINT_OK);

    /*  g == x, so g has x's value (0).  */
    check_i32_value(g, 0);
    check_i32_value(y, 1);
  }

  /*  Test large multi-limb inputs.  */
  {
    arbint_err_t rc;

    /*  a = 2^128, b = 2^64.  */
    CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
    CHECK_EQ_I(arbint_pow_u32(a, a, 128), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, 2), ARBINT_OK);
    CHECK_EQ_I(arbint_pow_u32(b, b, 64), ARBINT_OK);

    rc = arbint_xgcd(g, x, y, a, b);
    CHECK_EQ_I(rc, ARBINT_OK);

    /*  gcd(2^128, 2^64) = 2^64.  */
    CHECK_EQ_I(arbint_set_i32(tmp1, 2), ARBINT_OK);
    CHECK_EQ_I(arbint_pow_u32(tmp1, tmp1, 64), ARBINT_OK);
    CHECK(arbint_eq(g, tmp1));

    verify_bezout(a, b, g, x, y, tmp1, tmp2, "large: gcd(2^128, 2^64)");
  }

  /*  Test large coprime multi-limb inputs.  */
  {
    arbint_err_t rc;

    /*  a = 2^127 - 1 (Mersenne prime), b = 2^63 + 1.  */
    CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
    CHECK_EQ_I(arbint_pow_u32(a, a, 127), ARBINT_OK);
    CHECK_EQ_I(arbint_sub_i32(a, a, 1), ARBINT_OK);

    CHECK_EQ_I(arbint_set_i32(b, 2), ARBINT_OK);
    CHECK_EQ_I(arbint_pow_u32(b, b, 63), ARBINT_OK);
    CHECK_EQ_I(arbint_add_i32(b, b, 1), ARBINT_OK);

    rc = arbint_xgcd(g, x, y, a, b);
    CHECK_EQ_I(rc, ARBINT_OK);

    /*  These should be coprime, gcd = 1.  */
    check_i32_value(g, 1);
    verify_bezout(a, b, g, x, y, tmp1, tmp2, "large coprime");
  }

  /*  Test consecutive integers (always coprime).  */
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 99, 100, "gcd(99, 100)");
  check_i32_value(g, 1);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 1000, 1001, "gcd(1000, 1001)");
  check_i32_value(g, 1);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 9999, 10000, "gcd(9999, 10000)");
  check_i32_value(g, 1);

  /*  Test Fibonacci-like pairs (worst case for Euclidean).  */
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 89, 55, "gcd(89, 55) Fib");
  check_i32_value(g, 1);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 144, 89, "gcd(144, 89) Fib");
  check_i32_value(g, 1);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 233, 144, "gcd(233, 144) Fib");
  check_i32_value(g, 1);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 610, 377, "gcd(610, 377) Fib");
  check_i32_value(g, 1);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 987, 610, "gcd(987, 610) Fib");
  check_i32_value(g, 1);

  /*  Large Fibonacci numbers.  */
  {
    arbint_err_t rc;
    CHECK_EQ_I(arbint_fib_u32(a, 100), ARBINT_OK);
    CHECK_EQ_I(arbint_fib_u32(b, 99), ARBINT_OK);

    rc = arbint_xgcd(g, x, y, a, b);
    CHECK_EQ_I(rc, ARBINT_OK);
    check_i32_value(g, 1);
    verify_bezout(a, b, g, x, y, tmp1, tmp2, "gcd(Fib(100), Fib(99))");
  }

  /*  Test one divides the other.  */
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 100, 25, "gcd(100, 25)");
  check_i32_value(g, 25);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 1000, 8, "gcd(1000, 8)");
  check_i32_value(g, 8);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 144, 12, "gcd(144, 12)");
  check_i32_value(g, 12);

  /*  Test powers of 2.  */
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 1024, 256, "gcd(1024, 256)");
  check_i32_value(g, 256);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 65536, 1024, "gcd(65536, 1024)");
  check_i32_value(g, 1024);

  /*  Test large power-of-2 inputs.  */
  {
    arbint_err_t rc;

    CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
    CHECK_EQ_I(arbint_pow_u32(a, a, 256), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, 2), ARBINT_OK);
    CHECK_EQ_I(arbint_pow_u32(b, b, 128), ARBINT_OK);

    rc = arbint_xgcd(g, x, y, a, b);
    CHECK_EQ_I(rc, ARBINT_OK);

    CHECK_EQ_I(arbint_set_i32(tmp1, 2), ARBINT_OK);
    CHECK_EQ_I(arbint_pow_u32(tmp1, tmp1, 128), ARBINT_OK);
    CHECK(arbint_eq(g, tmp1));
    verify_bezout(a, b, g, x, y, tmp1, tmp2, "gcd(2^256, 2^128)");
  }

  /*  Test equal inputs.  */
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 42, 42, "gcd(42, 42)");
  check_i32_value(g, 42);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 1000, 1000, "gcd(1000, 1000)");
  check_i32_value(g, 1000);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, -42, -42, "gcd(-42, -42)");
  check_i32_value(g, 42);

  /*  Test with 1.  */
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 1, 1000, "gcd(1, 1000)");
  check_i32_value(g, 1);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 1000, 1, "gcd(1000, 1)");
  check_i32_value(g, 1);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 1, 1, "gcd(1, 1)");
  check_i32_value(g, 1);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, -1, 1000, "gcd(-1, 1000)");
  check_i32_value(g, 1);

  /*  Test with large primes.  */
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 104729, 104743,
                "gcd(104729, 104743) primes");
  check_i32_value(g, 1);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 999961, 999979,
                "gcd(999961, 999979) primes");
  check_i32_value(g, 1);

  /*  Test highly composite numbers.  */
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 720, 840, "gcd(720, 840)");
  check_i32_value(g, 120);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 5040, 2520, "gcd(5040, 2520)");
  check_i32_value(g, 2520);

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 27720, 55440, "gcd(27720, 55440)");
  check_i32_value(g, 27720);

  /*  Test symmetric pairs.  */
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 123, 321, "gcd(123, 321)");
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 321, 123, "gcd(321, 123) swapped");

  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 1234, 4321, "gcd(1234, 4321)");
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 4321, 1234,
                "gcd(4321, 1234) swapped");

  /*  Test all output aliasing in zero path: gcd(a, 0) with g==y.  */
  {
    arbint_err_t rc;
    CHECK_EQ_I(arbint_set_i32(a, 17), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, 0), ARBINT_OK);

    rc = arbint_xgcd(g, x, g, a, b);
    CHECK_EQ_I(rc, ARBINT_OK);

    /*  g == y, so g has y's value (0).  */
    check_i32_value(g, 0);
    check_i32_value(x, 1);
  }

  /*  Test zero path with x==y.  */
  {
    arbint_err_t rc;
    CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, 23), ARBINT_OK);

    rc = arbint_xgcd(g, x, x, a, b);
    CHECK_EQ_I(rc, ARBINT_OK);

    check_i32_value(g, 23);
    /*  x == y, x has y's value (1).  */
    check_i32_value(x, 1);
  }

  /*  Test complex aliasing: g==a, x==b.  */
  {
    arbint_err_t rc;
    CHECK_EQ_I(arbint_set_i32(g, 48), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(x, 18), ARBINT_OK);

    rc = arbint_xgcd(g, x, y, g, x);
    CHECK_EQ_I(rc, ARBINT_OK);

    check_i32_value(g, 6);
    CHECK_EQ_I(arbint_set_i32(a, 48), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, 18), ARBINT_OK);
    verify_bezout(a, b, g, x, y, tmp1, tmp2, "g==a, x==b aliasing");
  }

  /*  Test negative multi-limb inputs.  */
  {
    arbint_err_t rc;

    CHECK_EQ_I(arbint_set_i32(a, -2), ARBINT_OK);
    CHECK_EQ_I(arbint_pow_u32(a, a, 100), ARBINT_OK); /*  -2^100.  */
    CHECK_EQ_I(arbint_set_i32(b, 2), ARBINT_OK);
    CHECK_EQ_I(arbint_pow_u32(b, b, 50), ARBINT_OK);

    rc = arbint_xgcd(g, x, y, a, b);
    CHECK_EQ_I(rc, ARBINT_OK);

    CHECK_EQ_I(arbint_set_i32(tmp1, 2), ARBINT_OK);
    CHECK_EQ_I(arbint_pow_u32(tmp1, tmp1, 50), ARBINT_OK);
    CHECK(arbint_eq(g, tmp1));
    verify_bezout(a, b, g, x, y, tmp1, tmp2, "gcd(-2^100, 2^50)");
  }

  /*  Both negative multi-limb.  */
  {
    arbint_err_t rc;

    CHECK_EQ_I(arbint_set_i32(a, -3), ARBINT_OK);
    CHECK_EQ_I(arbint_pow_u32(a, a, 80), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, -3), ARBINT_OK);
    CHECK_EQ_I(arbint_pow_u32(b, b, 60), ARBINT_OK);

    rc = arbint_xgcd(g, x, y, a, b);
    CHECK_EQ_I(rc, ARBINT_OK);

    CHECK_EQ_I(arbint_set_i32(tmp1, 3), ARBINT_OK);
    CHECK_EQ_I(arbint_pow_u32(tmp1, tmp1, 60), ARBINT_OK);
    CHECK(arbint_eq(g, tmp1));
    verify_bezout(a, b, g, x, y, tmp1, tmp2, "gcd(-3^80, -3^60)");
  }

  /*  Test very unequal sizes.  */
  {
    arbint_err_t rc;

    CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
    CHECK_EQ_I(arbint_pow_u32(a, a, 500), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, 7), ARBINT_OK);

    rc = arbint_xgcd(g, x, y, a, b);
    CHECK_EQ_I(rc, ARBINT_OK);
    check_i32_value(g, 1);
    verify_bezout(a, b, g, x, y, tmp1, tmp2, "gcd(2^500, 7)");
  }

  /*  Test products of distinct primes.  */
  /*  gcd(2*3*5, 3*5*7) = 15.  */
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 30, 105, "gcd(30, 105)");
  check_i32_value(g, 15);

  /*  gcd(2*3*5*7, 3*5*7*11) = 105.  */
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 210, 1155, "gcd(210, 1155)");
  check_i32_value(g, 105);

  /*  gcd(2*3*5*7*11, 5*7*11*13) = 385.  */
  test_xgcd_i32(a, b, g, x, y, tmp1, tmp2, 2310, 5005, "gcd(2310, 5005)");
  check_i32_value(g, 385);

  arbint_clear(tmp2);
  arbint_clear(tmp1);
  arbint_clear(y);
  arbint_clear(x);
  arbint_clear(g);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);

  ARBINT_TEST_FINISH("test_xgcd");
}
