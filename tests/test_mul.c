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

int main(void) {
  test_mul_basic_and_alias();
  test_mul_large_patterns();
  test_mul_carry_regression();
  ARBINT_TEST_FINISH("test_mul");
}
