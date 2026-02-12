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

/*  Helper: compute a += b*c using arbint_mul + arbint_add and compare
    against result from arbint_addmul.  */
static void verify_addmul(arbint_t a, const arbint_t b, const arbint_t c,
                          arbint_ctx_t * ctx) {
  arbint_t expected;
  arbint_t prod;
  arbint_t a_copy;

  CHECK_EQ_I(arbint_init(expected, ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(prod, ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a_copy, ctx), ARBINT_OK);

  /*  Save original a.  */
  CHECK_EQ_I(arbint_set(a_copy, a), ARBINT_OK);

  /*  expected = a + b*c (naive).  */
  CHECK_EQ_I(arbint_mul(prod, b, c), ARBINT_OK);
  CHECK_EQ_I(arbint_add(expected, a, prod), ARBINT_OK);

  /*  Use addmul on a_copy.  */
  CHECK_EQ_I(arbint_addmul(a_copy, b, c), ARBINT_OK);

  /*  Results must match.  */
  CHECK(arbint_eq(a_copy, expected));

  arbint_clear(a_copy);
  arbint_clear(prod);
  arbint_clear(expected);
}

/*  Helper: compute a -= b*c using arbint_mul + arbint_sub and compare
    against result from arbint_submul.  */
static void verify_submul(arbint_t a, const arbint_t b, const arbint_t c,
                          arbint_ctx_t * ctx) {
  arbint_t expected;
  arbint_t prod;
  arbint_t a_copy;

  CHECK_EQ_I(arbint_init(expected, ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(prod, ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a_copy, ctx), ARBINT_OK);

  /*  Save original a.  */
  CHECK_EQ_I(arbint_set(a_copy, a), ARBINT_OK);

  /*  expected = a - b*c (naive).  */
  CHECK_EQ_I(arbint_mul(prod, b, c), ARBINT_OK);
  CHECK_EQ_I(arbint_sub(expected, a, prod), ARBINT_OK);

  /*  Use submul on a_copy.  */
  CHECK_EQ_I(arbint_submul(a_copy, b, c), ARBINT_OK);

  /*  Results must match.  */
  CHECK(arbint_eq(a_copy, expected));

  arbint_clear(a_copy);
  arbint_clear(prod);
  arbint_clear(expected);
}

static void test_addmul_null_checks(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_addmul(NULL, a, b), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_addmul(a, NULL, b), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_addmul(a, a, NULL), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_submul(NULL, a, b), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_submul(a, NULL, b), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_submul(a, a, NULL), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_addmul_u32(NULL, a, 5u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_addmul_u32(a, NULL, 5u), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_submul_u32(NULL, a, 5u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_submul_u32(a, NULL, 5u), ARBINT_EINVAL);

  arbint_clear(b);
  arbint_clear(a);
}

static void test_addmul_zero_cases(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t c;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);

  /*  a += 0*c = a (unchanged).  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_addmul(a, b, c), ARBINT_OK);
  check_i32_value(a, 42);

  /*  a += b*0 = a (unchanged).  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_addmul(a, b, c), ARBINT_OK);
  check_i32_value(a, 42);

  /*  0 += b*c = b*c.  */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 6), ARBINT_OK);
  CHECK_EQ_I(arbint_addmul(a, b, c), ARBINT_OK);
  check_i32_value(a, 42);

  /*  a -= 0*c = a (unchanged).  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_submul(a, b, c), ARBINT_OK);
  check_i32_value(a, 42);

  /*  0 -= b*c = -(b*c).  */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 6), ARBINT_OK);
  CHECK_EQ_I(arbint_submul(a, b, c), ARBINT_OK);
  check_i32_value(a, -42);

  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
}

static void test_addmul_basic(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t c;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);

  /*  10 += 3*4 = 22.  */
  CHECK_EQ_I(arbint_set_i32(a, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_addmul(a, b, c), ARBINT_OK);
  check_i32_value(a, 22);

  /*  10 -= 3*4 = -2.  */
  CHECK_EQ_I(arbint_set_i32(a, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_submul(a, b, c), ARBINT_OK);
  check_i32_value(a, -2);

  /*  -10 += (-3)*(-4) = 2 (positive + positive).  */
  CHECK_EQ_I(arbint_set_i32(a, -10), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -3), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, -4), ARBINT_OK);
  CHECK_EQ_I(arbint_addmul(a, b, c), ARBINT_OK);
  check_i32_value(a, 2);

  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
}

static void test_addmul_sign_combinations(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t c;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);

  /*  All 8 sign combinations for addmul.  */

  /*  (+a) += (+b)*(+c) = positive + positive.  */
  CHECK_EQ_I(arbint_set_i32(a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 8), ARBINT_OK);
  verify_addmul(a, b, c, &ctx);

  /*  (+a) += (+b)*(-c) = positive + negative.  */
  CHECK_EQ_I(arbint_set_i32(a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, -8), ARBINT_OK);
  verify_addmul(a, b, c, &ctx);

  /*  (+a) += (-b)*(+c) = positive + negative.  */
  CHECK_EQ_I(arbint_set_i32(a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 8), ARBINT_OK);
  verify_addmul(a, b, c, &ctx);

  /*  (+a) += (-b)*(-c) = positive + positive.  */
  CHECK_EQ_I(arbint_set_i32(a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, -8), ARBINT_OK);
  verify_addmul(a, b, c, &ctx);

  /*  (-a) += (+b)*(+c) = negative + positive.  */
  CHECK_EQ_I(arbint_set_i32(a, -100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 8), ARBINT_OK);
  verify_addmul(a, b, c, &ctx);

  /*  (-a) += (+b)*(-c) = negative + negative.  */
  CHECK_EQ_I(arbint_set_i32(a, -100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, -8), ARBINT_OK);
  verify_addmul(a, b, c, &ctx);

  /*  (-a) += (-b)*(+c) = negative + negative.  */
  CHECK_EQ_I(arbint_set_i32(a, -100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 8), ARBINT_OK);
  verify_addmul(a, b, c, &ctx);

  /*  (-a) += (-b)*(-c) = negative + positive.  */
  CHECK_EQ_I(arbint_set_i32(a, -100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, -8), ARBINT_OK);
  verify_addmul(a, b, c, &ctx);

  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
}

static void test_submul_sign_combinations(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t c;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);

  /*  All 8 sign combinations for submul.  */

  /*  (+a) -= (+b)*(+c) = positive - positive.  */
  CHECK_EQ_I(arbint_set_i32(a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 8), ARBINT_OK);
  verify_submul(a, b, c, &ctx);

  /*  (+a) -= (+b)*(-c) = positive - negative.  */
  CHECK_EQ_I(arbint_set_i32(a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, -8), ARBINT_OK);
  verify_submul(a, b, c, &ctx);

  /*  (+a) -= (-b)*(+c) = positive - negative.  */
  CHECK_EQ_I(arbint_set_i32(a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 8), ARBINT_OK);
  verify_submul(a, b, c, &ctx);

  /*  (+a) -= (-b)*(-c) = positive - positive.  */
  CHECK_EQ_I(arbint_set_i32(a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, -8), ARBINT_OK);
  verify_submul(a, b, c, &ctx);

  /*  (-a) -= (+b)*(+c) = negative - positive.  */
  CHECK_EQ_I(arbint_set_i32(a, -100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 8), ARBINT_OK);
  verify_submul(a, b, c, &ctx);

  /*  (-a) -= (+b)*(-c) = negative - negative.  */
  CHECK_EQ_I(arbint_set_i32(a, -100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, -8), ARBINT_OK);
  verify_submul(a, b, c, &ctx);

  /*  (-a) -= (-b)*(+c) = negative - negative.  */
  CHECK_EQ_I(arbint_set_i32(a, -100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 8), ARBINT_OK);
  verify_submul(a, b, c, &ctx);

  /*  (-a) -= (-b)*(-c) = negative - positive.  */
  CHECK_EQ_I(arbint_set_i32(a, -100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, -8), ARBINT_OK);
  verify_submul(a, b, c, &ctx);

  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
}

static void test_addmul_aliasing(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t c;
  arbint_t expected;
  arbint_t prod;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(prod, &ctx), ARBINT_OK);

  /*  a += a*c (a == b).  */
  CHECK_EQ_I(arbint_set_i32(a, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 3), ARBINT_OK);

  CHECK_EQ_I(arbint_mul(prod, a, c), ARBINT_OK);
  CHECK_EQ_I(arbint_add(expected, a, prod), ARBINT_OK);
  CHECK_EQ_I(arbint_addmul(a, a, c), ARBINT_OK);
  CHECK(arbint_eq(a, expected));

  /*  a += b*a (a == c).  */
  CHECK_EQ_I(arbint_set_i32(a, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 3), ARBINT_OK);

  CHECK_EQ_I(arbint_mul(prod, b, a), ARBINT_OK);
  CHECK_EQ_I(arbint_add(expected, a, prod), ARBINT_OK);
  CHECK_EQ_I(arbint_addmul(a, b, a), ARBINT_OK);
  CHECK(arbint_eq(a, expected));

  /*  a += a*a (a == b == c).  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);

  CHECK_EQ_I(arbint_mul(prod, a, a), ARBINT_OK);
  CHECK_EQ_I(arbint_add(expected, a, prod), ARBINT_OK);
  CHECK_EQ_I(arbint_addmul(a, a, a), ARBINT_OK);
  CHECK(arbint_eq(a, expected));

  /*  submul aliasing: a -= a*c.  */
  CHECK_EQ_I(arbint_set_i32(a, 50), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 2), ARBINT_OK);

  CHECK_EQ_I(arbint_mul(prod, a, c), ARBINT_OK);
  CHECK_EQ_I(arbint_sub(expected, a, prod), ARBINT_OK);
  CHECK_EQ_I(arbint_submul(a, a, c), ARBINT_OK);
  CHECK(arbint_eq(a, expected));

  /*  submul aliasing: a -= a*a.  */
  CHECK_EQ_I(arbint_set_i32(a, 10), ARBINT_OK);

  CHECK_EQ_I(arbint_mul(prod, a, a), ARBINT_OK);
  CHECK_EQ_I(arbint_sub(expected, a, prod), ARBINT_OK);
  CHECK_EQ_I(arbint_submul(a, a, a), ARBINT_OK);
  CHECK(arbint_eq(a, expected));

  arbint_clear(prod);
  arbint_clear(expected);
  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
}

static void test_addmul_u32(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t expected;
  arbint_t prod;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(prod, &ctx), ARBINT_OK);

  /*  a += b*0 = a (unchanged).  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_addmul_u32(a, b, 0u), ARBINT_OK);
  check_i32_value(a, 42);

  /*  a += b*1 = a + b.  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_addmul_u32(a, b, 1u), ARBINT_OK);
  check_i32_value(a, 52);

  /*  10 += 3*4 = 22.  */
  CHECK_EQ_I(arbint_set_i32(a, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_addmul_u32(a, b, 4u), ARBINT_OK);
  check_i32_value(a, 22);

  /*  10 += (-3)*4 = -2.  */
  CHECK_EQ_I(arbint_set_i32(a, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -3), ARBINT_OK);
  CHECK_EQ_I(arbint_addmul_u32(a, b, 4u), ARBINT_OK);
  check_i32_value(a, -2);

  /*  submul: 10 -= 3*4 = -2.  */
  CHECK_EQ_I(arbint_set_i32(a, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_submul_u32(a, b, 4u), ARBINT_OK);
  check_i32_value(a, -2);

  /*  submul: 0 -= 7*6 = -42.  */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_submul_u32(a, b, 6u), ARBINT_OK);
  check_i32_value(a, -42);

  /*  aliasing: a += a*5.  */
  CHECK_EQ_I(arbint_set_i32(a, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(prod, a, 5u), ARBINT_OK);
  CHECK_EQ_I(arbint_add(expected, a, prod), ARBINT_OK);
  CHECK_EQ_I(arbint_addmul_u32(a, a, 5u), ARBINT_OK);
  CHECK(arbint_eq(a, expected));

  arbint_clear(prod);
  arbint_clear(expected);
  arbint_clear(b);
  arbint_clear(a);
}

static void test_addmul_i32(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);

  /*  a += b*0 = a (unchanged).  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_addmul_i32(a, b, 0), ARBINT_OK);
  check_i32_value(a, 42);

  /*  10 += 3*4 = 22.  */
  CHECK_EQ_I(arbint_set_i32(a, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_addmul_i32(a, b, 4), ARBINT_OK);
  check_i32_value(a, 22);

  /*  10 += 3*(-4) = -2 (negative multiplier).  */
  CHECK_EQ_I(arbint_set_i32(a, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_addmul_i32(a, b, -4), ARBINT_OK);
  check_i32_value(a, -2);

  /*  submul: 10 -= 3*4 = -2.  */
  CHECK_EQ_I(arbint_set_i32(a, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_submul_i32(a, b, 4), ARBINT_OK);
  check_i32_value(a, -2);

  /*  submul with negative c: 10 -= 3*(-4) = 22.  */
  CHECK_EQ_I(arbint_set_i32(a, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_submul_i32(a, b, -4), ARBINT_OK);
  check_i32_value(a, 22);

  /*  INT32_MIN edge case.  */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_addmul_i32(a, b, INT32_MIN), ARBINT_OK);
  check_i32_value(a, INT32_MIN);

  arbint_clear(b);
  arbint_clear(a);
}

/*  Build a large non-trivial value with at least 'min_limbs' limbs.
    Uses 3^N via repeated squaring and multiply-by-3, producing
    non-power-of-two values with interesting bit patterns.

    On 64-bit, nbits(3^N) ~ N * 1.585, so limbs ~ N * 1.585 / 64.
    For 32 limbs (Karatsuba): N ~ 1300 -> use 3^1300 via pow_u32.
    For 96 limbs (Toom-3): N ~ 3900 -> build iteratively.

    Since pow_u32 is O(n^2 * log(exp)) and we need large exponents,
    we build incrementally: start with 3^16, then square repeatedly.  */
static void build_large_value(arbint_t x, size_t min_limbs,
                              arbint_ctx_t * ctx) {
  size_t i;

  (void) ctx;

  /*  Start with x = 3.  */
  CHECK_EQ_I(arbint_set_i32(x, 3), ARBINT_OK);

  /*  Square until we reach at least min_limbs.
      Each squaring roughly doubles the number of limbs.  */
  for (i = 0u; i < 20u; ++i) {
    size_t n = arbint_abs_sz(x[0]._sz);
    if (n >= min_limbs)
      break;
    CHECK_EQ_I(arbint_sqr(x, x), ARBINT_OK);
  }

  /*  Add 1 to avoid a perfect square (more interesting test value).  */
  CHECK_EQ_I(arbint_add_i32(x, x, 1), ARBINT_OK);
}

static void test_addmul_large(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t c;
  size_t i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);

  /*  Build large numbers via repeated squaring.  */
  CHECK_EQ_I(arbint_set_i32(b, 2), ARBINT_OK);
  for (i = 0u; i < 10u; ++i) {
    CHECK_EQ_I(arbint_sqr(b, b), ARBINT_OK);
  }
  /*  b = 2^1024.  */

  CHECK_EQ_I(arbint_set(c, b), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);

  /*  Verify addmul/submul with large operands.  */
  verify_addmul(a, b, c, &ctx);

  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(b, b), ARBINT_OK);
  verify_addmul(a, b, c, &ctx);

  CHECK_EQ_I(arbint_neg(a, a), ARBINT_OK);
  verify_submul(a, b, c, &ctx);

  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
}

/*  Test addmul/submul with Karatsuba-sized operands (32+ limbs).
    Exercises the threshold where fused schoolbook mulacc gives way
    to recursive multiplication + add.  */
static void test_addmul_karatsuba(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t c;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);

  /*  Build values with ~40 limbs (above Karatsuba threshold of 32).  */
  build_large_value(b, 40u, &ctx);
  build_large_value(c, 40u, &ctx);

  /*  Make c different from b.  */
  CHECK_EQ_I(arbint_add_i32(c, c, 17), ARBINT_OK);

  /*  Same sign: positive a, positive product.  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  verify_addmul(a, b, c, &ctx);

  /*  Same sign: large a, positive product (both large).  */
  CHECK_EQ_I(arbint_set(a, b), ARBINT_OK);
  verify_addmul(a, b, c, &ctx);

  /*  Different sign: positive a, negative product.  */
  CHECK_EQ_I(arbint_set(a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(c, c), ARBINT_OK);
  verify_addmul(a, b, c, &ctx);

  /*  submul with Karatsuba-sized operands.  */
  CHECK_EQ_I(arbint_neg(c, c), ARBINT_OK);
  CHECK_EQ_I(arbint_set(a, b), ARBINT_OK);
  verify_submul(a, b, c, &ctx);

  /*  submul: negative a, positive product.  */
  CHECK_EQ_I(arbint_neg(a, a), ARBINT_OK);
  verify_submul(a, b, c, &ctx);

  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
}

/*  Test addmul/submul with Toom-3-sized operands (96+ limbs).
    Exercises the full Toom-Cook-3 multiplication path.  */
static void test_addmul_toom3(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t c;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);

  /*  Build values with ~100 limbs (above Toom-3 threshold of 96).  */
  build_large_value(b, 100u, &ctx);
  build_large_value(c, 100u, &ctx);

  /*  Make c different from b.  */
  CHECK_EQ_I(arbint_add_i32(c, c, 31), ARBINT_OK);

  /*  Same sign addmul: small a + large product.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  verify_addmul(a, b, c, &ctx);

  /*  Same sign addmul: large a + large product.  */
  CHECK_EQ_I(arbint_set(a, b), ARBINT_OK);
  verify_addmul(a, b, c, &ctx);

  /*  Different sign addmul: large positive a + large negative product.  */
  CHECK_EQ_I(arbint_set(a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(c, c), ARBINT_OK);
  verify_addmul(a, b, c, &ctx);

  /*  Different sign addmul: negative a + positive product.  */
  CHECK_EQ_I(arbint_neg(a, a), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(c, c), ARBINT_OK);
  verify_addmul(a, b, c, &ctx);

  /*  Same sign submul: large a - large product.  */
  CHECK_EQ_I(arbint_set(a, b), ARBINT_OK);
  verify_submul(a, b, c, &ctx);

  /*  Different sign submul.  */
  CHECK_EQ_I(arbint_neg(a, a), ARBINT_OK);
  verify_submul(a, b, c, &ctx);

  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
}

/*  Test addmul at exact threshold boundaries (31, 32, 33 limbs).
    Catches off-by-one errors in the Karatsuba threshold check.  */
static void test_addmul_threshold_boundary(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t c;
  static const size_t edges[] = {31u, 32u, 33u, 95u, 96u, 97u};
  size_t ei;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);

  for (ei = 0u; ei < sizeof(edges) / sizeof(edges[0]); ++ei) {
    size_t target = edges[ei];

    build_large_value(b, target, &ctx);
    build_large_value(c, target, &ctx);
    CHECK_EQ_I(arbint_add_i32(c, c, 7), ARBINT_OK);

    /*  addmul same sign.  */
    CHECK_EQ_I(arbint_set(a, b), ARBINT_OK);
    verify_addmul(a, b, c, &ctx);

    /*  addmul different sign.  */
    CHECK_EQ_I(arbint_set(a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_neg(c, c), ARBINT_OK);
    verify_addmul(a, b, c, &ctx);

    /*  submul.  */
    CHECK_EQ_I(arbint_neg(c, c), ARBINT_OK);
    CHECK_EQ_I(arbint_set(a, b), ARBINT_OK);
    verify_submul(a, b, c, &ctx);
  }

  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
}

static void test_addmul_exact_cancellation(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t c;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);

  /*  12 -= 3*4 = 0 (exact cancellation).  */
  CHECK_EQ_I(arbint_set_i32(a, 12), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_submul(a, b, c), ARBINT_OK);
  CHECK(arbint_is_zero(a));

  /*  (-12) += 3*4 = 0 (exact cancellation).  */
  CHECK_EQ_I(arbint_set_i32(a, -12), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(c, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_addmul(a, b, c), ARBINT_OK);
  CHECK(arbint_is_zero(a));

  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
}

int main(void) {
  ARBINT_TEST_START();
  test_addmul_null_checks();
  test_addmul_zero_cases();
  test_addmul_basic();
  test_addmul_sign_combinations();
  test_submul_sign_combinations();
  test_addmul_aliasing();
  test_addmul_u32();
  test_addmul_i32();
  test_addmul_large();
  test_addmul_karatsuba();
  test_addmul_toom3();
  test_addmul_threshold_boundary();
  test_addmul_exact_cancellation();
  ARBINT_TEST_FINISH("test_addmul");
}
