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

static void check_tdiv_i32_case(arbint_t n, arbint_t d, arbint_t q, arbint_t r,
                                int32_t nv, int32_t dv, int32_t qv,
                                int32_t rv) {
  CHECK_EQ_I(arbint_set_i32(n, nv), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(d, dv), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  check_i32_value(q, qv);
  check_i32_value(r, rv);
}

static void test_tdiv_small_and_signs(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t d;
  arbint_t q;
  arbint_t r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_tdiv_qr(NULL, r, n, d), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_tdiv_qr(q, NULL, n, d), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, NULL, d), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, NULL), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_set_i32(n, 123), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(d, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_EZERO);

  check_tdiv_i32_case(n, d, q, r, 17, 5, 3, 2);
  check_tdiv_i32_case(n, d, q, r, -17, 5, -3, -2);
  check_tdiv_i32_case(n, d, q, r, 17, -5, -3, 2);
  check_tdiv_i32_case(n, d, q, r, -17, -5, 3, -2);
  check_tdiv_i32_case(n, d, q, r, 2, 5, 0, 2);
  check_tdiv_i32_case(n, d, q, r, -2, 5, 0, -2);

  CHECK_EQ_I(arbint_set_i32(n, -17), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 5u), ARBINT_OK);
  check_i32_value(q, -3);
  check_i32_value(r, -2);

  CHECK_EQ_I(arbint_tdiv_qr_i32(q, r, n, -5), ARBINT_OK);
  check_i32_value(q, 3);
  check_i32_value(r, -2);

  CHECK_EQ_I(arbint_tdiv_q_u32(q, n, 5u), ARBINT_OK);
  check_i32_value(q, -3);
  CHECK_EQ_I(arbint_tdiv_r_u32(r, n, 5u), ARBINT_OK);
  check_i32_value(r, -2);

  CHECK_EQ_I(arbint_tdiv_q_i32(q, n, -5), ARBINT_OK);
  check_i32_value(q, 3);
  CHECK_EQ_I(arbint_tdiv_r_i32(r, n, -5), ARBINT_OK);
  check_i32_value(r, -2);

  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 0u), ARBINT_EZERO);
  CHECK_EQ_I(arbint_tdiv_qr_i32(q, r, n, 0), ARBINT_EZERO);

  arbint_clear(n);
  arbint_clear(d);
  arbint_clear(q);
  arbint_clear(r);
}

static void test_tdiv_multilimb(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t d;
  arbint_t q;
  arbint_t r;
  arbint_limb_t nmag[3];
  arbint_limb_t dmag[2];
  arbint_limb_t qmag[2];
  arbint_limb_t rmag[1];

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  nmag[0] = (arbint_limb_t) 7u;
  nmag[1] = (arbint_limb_t) 5u;
  nmag[2] = (arbint_limb_t) 1u;
  dmag[0] = (arbint_limb_t) 1u;
  dmag[1] = (arbint_limb_t) 1u;

  qmag[0] = (arbint_limb_t) 4u;
  qmag[1] = (arbint_limb_t) 1u;
  rmag[0] = (arbint_limb_t) 3u;

  set_mag_limbs(n, 1, nmag, 3u);
  set_mag_limbs(d, 1, dmag, 2u);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  CHECK(limbs_equal(q, 1, qmag, 2u));
  CHECK(limbs_equal(r, 1, rmag, 1u));

  set_mag_limbs(n, -1, nmag, 3u);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  CHECK(limbs_equal(q, -1, qmag, 2u));
  CHECK(limbs_equal(r, -1, rmag, 1u));

  set_mag_limbs(n, 1, nmag, 3u);
  set_mag_limbs(d, -1, dmag, 2u);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  CHECK(limbs_equal(q, -1, qmag, 2u));
  CHECK(limbs_equal(r, 1, rmag, 1u));

  arbint_clear(n);
  arbint_clear(d);
  arbint_clear(q);
  arbint_clear(r);
}

int main(void) {
  test_tdiv_small_and_signs();
  test_tdiv_multilimb();
  ARBINT_TEST_FINISH("test_tdiv");
}
