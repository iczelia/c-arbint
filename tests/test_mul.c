/* arbint - portable arbitrary-precision computation library
 *
 * Copyright (C) 2026 Kamila Szewczyk (k@iczelia.net)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include "test_framework.h"

#include <stdint.h>
#include <string.h>

ARBINT_TEST_DECLARE_FAILURES();

static void test_mul_basic_and_alias(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_mul(NULL, a, b), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_mul(r, NULL, b), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_mul(r, a, NULL), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_sqr(NULL, a), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_sqr(r, NULL), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_set_i32(a, 1234), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -567), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  check_i32_value(r, -699678);

  CHECK_EQ_I(arbint_mul(a, a, b), ARBINT_OK);
  check_i32_value(a, -699678);

  CHECK_EQ_I(arbint_set_i32(a, 1234), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -567), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(b, a, b), ARBINT_OK);
  check_i32_value(b, -699678);

  CHECK_EQ_I(arbint_set_i32(a, 13), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(a, a, a), ARBINT_OK);
  check_i32_value(a, 169);

  CHECK_EQ_I(arbint_set_i32(a, -19), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
  check_i32_value(r, 361);

  CHECK_EQ_I(arbint_set_i32(a, -19), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
  check_i32_value(a, 361);

  CHECK_EQ_I(arbint_set_u32(a, 0u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -12345), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK(arbint_is_zero(r));

  arbint_clear(a);
  arbint_clear(b);
  arbint_clear(r);
}

static void test_mul_large_patterns(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t r;
  arbint_t t;
  arbint_limb_t am[64];
  arbint_limb_t bm[64];
  arbint_limb_t em[128];
  size_t n = 40u;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t, &ctx), ARBINT_OK);

  memset(am, 0, sizeof(am));
  memset(bm, 0, sizeof(bm));
  memset(em, 0, sizeof(em));

  am[0] = 1u;
  am[n - 1u] = 1u;
  set_mag_limbs(a, 1, am, n);

  CHECK_EQ_I(arbint_mul(r, a, a), ARBINT_OK);
  em[0] = 1u;
  em[n - 1u] = 2u;
  em[2u * n - 2u] = 1u;
  CHECK(limbs_equal(r, 1, em, 2u * n - 1u));

  CHECK_EQ_I(arbint_sqr(t, a), ARBINT_OK);
  CHECK(arbint_eq(r, t));

  memset(em, 0, sizeof(em));
  bm[0] = 3u;
  bm[n - 1u] = 1u;
  set_mag_limbs(b, 1, bm, n);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  em[0] = 3u;
  em[n - 1u] = 4u;
  em[2u * n - 2u] = 1u;
  CHECK(limbs_equal(r, 1, em, 2u * n - 1u));

  CHECK_EQ_I(arbint_set(a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(t, a), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(a, a, a), ARBINT_OK);
  CHECK(arbint_eq(a, t));

  set_mag_limbs(a, -1, am, n);
  CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
  CHECK(limbs_equal(r, 1, em, 2u * n - 1u) == 0);

  memset(em, 0, sizeof(em));
  em[0] = 1u;
  em[n - 1u] = 2u;
  em[2u * n - 2u] = 1u;
  CHECK(limbs_equal(r, 1, em, 2u * n - 1u));

  arbint_clear(a);
  arbint_clear(b);
  arbint_clear(r);
  arbint_clear(t);
}

/* Regression for carry accounting in muladd_limb:
 * (B + 2) * (B - 1) = B^2 + B - 2 */
static void test_mul_carry_regression(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t r;
  arbint_limb_t am[2];
  arbint_limb_t bm[1];
  arbint_limb_t em[3];
  arbint_limb_t maxv = (arbint_limb_t) ~((arbint_limb_t) 0u);

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  am[0] = (arbint_limb_t) 2u;
  am[1] = (arbint_limb_t) 1u;
  bm[0] = maxv;

  em[0] = maxv - (arbint_limb_t) 1u;
  em[1] = (arbint_limb_t) 0u;
  em[2] = (arbint_limb_t) 1u;

  set_mag_limbs(a, 1, am, 2u);
  set_mag_limbs(b, 1, bm, 1u);

  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK(limbs_equal(r, 1, em, 3u));

  CHECK_EQ_I(arbint_mul(a, a, b), ARBINT_OK);
  CHECK(limbs_equal(a, 1, em, 3u));

  CHECK_EQ_I(arbint_set(a, r), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(a, a, a), ARBINT_OK);
  CHECK(arbint_eq(a, r));

  arbint_clear(a);
  arbint_clear(b);
  arbint_clear(r);
}

static void test_mul_u32(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /* NULL checks */
  CHECK_EQ_I(arbint_mul_u32(NULL, a, 5u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_mul_u32(r, NULL, 5u), ARBINT_EINVAL);

  /* 0 * scalar = 0 */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(r, a, 42u), ARBINT_OK);
  CHECK(arbint_is_zero(r));

  /* a * 0 = 0 */
  CHECK_EQ_I(arbint_set_i32(a, 12345), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(r, a, 0u), ARBINT_OK);
  CHECK(arbint_is_zero(r));

  /* positive * positive */
  CHECK_EQ_I(arbint_set_i32(a, 1234), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(r, a, 567u), ARBINT_OK);
  check_i32_value(r, 699678);

  /* negative * positive */
  CHECK_EQ_I(arbint_set_i32(a, -1234), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(r, a, 567u), ARBINT_OK);
  check_i32_value(r, -699678);

  /* multiply by 1 */
  CHECK_EQ_I(arbint_set_i32(a, -42), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(r, a, 1u), ARBINT_OK);
  check_i32_value(r, -42);

  /* aliasing: rop == a */
  CHECK_EQ_I(arbint_set_i32(a, 1234), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(a, a, 567u), ARBINT_OK);
  check_i32_value(a, 699678);

  /* large scalar (max u32) */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(r, a, UINT32_MAX), ARBINT_OK);
  {
    /* 2 * 0xFFFFFFFF = 0x1FFFFFFFE */
    arbint_t expected;
    CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);
    CHECK_EQ_I(arbint_set_u32(expected, UINT32_MAX), ARBINT_OK);
    CHECK_EQ_I(arbint_add(expected, expected, expected), ARBINT_OK);
    CHECK(arbint_eq(r, expected));
    arbint_clear(expected);
  }

  /* multi-limb * scalar: cross-check against arbint_mul */
  {
    arbint_t b;
    arbint_t expected;
    arbint_limb_t am[4];

    CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
    CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

    am[0] = (arbint_limb_t) ~((arbint_limb_t) 0u);
    am[1] = (arbint_limb_t) ~((arbint_limb_t) 0u);
    am[2] = (arbint_limb_t) 1u;
    am[3] = (arbint_limb_t) 42u;
    set_mag_limbs(a, 1, am, 4u);

    CHECK_EQ_I(arbint_set_u32(b, 12345u), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(expected, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_mul_u32(r, a, 12345u), ARBINT_OK);
    CHECK(arbint_eq(r, expected));

    /* negative multi-limb */
    set_mag_limbs(a, -1, am, 4u);
    CHECK_EQ_I(arbint_mul(expected, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_mul_u32(r, a, 12345u), ARBINT_OK);
    CHECK(arbint_eq(r, expected));

    arbint_clear(expected);
    arbint_clear(b);
  }

  arbint_clear(r);
  arbint_clear(a);
}

static void test_mul_i32(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /* NULL checks */
  CHECK_EQ_I(arbint_mul_i32(NULL, a, 5), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_mul_i32(r, NULL, 5), ARBINT_EINVAL);

  /* 0 * scalar = 0 */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_i32(r, a, 42), ARBINT_OK);
  CHECK(arbint_is_zero(r));

  /* a * 0 = 0 */
  CHECK_EQ_I(arbint_set_i32(a, 12345), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_i32(r, a, 0), ARBINT_OK);
  CHECK(arbint_is_zero(r));

  /* positive * positive */
  CHECK_EQ_I(arbint_set_i32(a, 1234), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_i32(r, a, 567), ARBINT_OK);
  check_i32_value(r, 699678);

  /* positive * negative */
  CHECK_EQ_I(arbint_set_i32(a, 1234), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_i32(r, a, -567), ARBINT_OK);
  check_i32_value(r, -699678);

  /* negative * negative */
  CHECK_EQ_I(arbint_set_i32(a, -1234), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_i32(r, a, -567), ARBINT_OK);
  check_i32_value(r, 699678);

  /* negative * positive */
  CHECK_EQ_I(arbint_set_i32(a, -1234), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_i32(r, a, 567), ARBINT_OK);
  check_i32_value(r, -699678);

  /* INT32_MIN edge case */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_i32(r, a, INT32_MIN), ARBINT_OK);
  check_i32_value(r, INT32_MIN);

  /* multiply by -1 */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_i32(r, a, -1), ARBINT_OK);
  check_i32_value(r, -42);

  /* aliasing: rop == a */
  CHECK_EQ_I(arbint_set_i32(a, 1234), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_i32(a, a, -567), ARBINT_OK);
  check_i32_value(a, -699678);

  arbint_clear(r);
  arbint_clear(a);
}

int main(void) {
  test_mul_basic_and_alias();
  test_mul_large_patterns();
  test_mul_carry_regression();
  test_mul_u32();
  test_mul_i32();
  ARBINT_TEST_FINISH("test_mul");
}
