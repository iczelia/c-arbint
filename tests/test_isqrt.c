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

/*  Helper: verify that r = isqrt(a), i.e. r*r <= a < (r+1)*(r+1).  */
static void check_isqrt_invariant(const arbint_t r, const arbint_t a) {
  arbint_ctx_t * ctx = arbint_get_ctx(r);
  arbint_t sq, rp1, sq_rp1;

  CHECK_EQ_I(arbint_init(sq, ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(rp1, ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(sq_rp1, ctx), ARBINT_OK);

  /*  r*r <= a  */
  CHECK_EQ_I(arbint_sqr(sq, r), ARBINT_OK);
  CHECK(arbint_cmp(sq, a) <= 0);

  /*  (r+1)*(r+1) > a  */
  CHECK_EQ_I(arbint_add_i32(rp1, r, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(sq_rp1, rp1), ARBINT_OK);
  CHECK(arbint_cmp(sq_rp1, a) > 0);

  arbint_clear(sq_rp1);
  arbint_clear(rp1);
  arbint_clear(sq);
}

static void test_isqrt_edge_cases(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  isqrt(0) = 0.  */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  check_i32_value(r, 0);

  /*  isqrt(1) = 1.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  check_i32_value(r, 1);

  /*  isqrt(2) = 1.  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  check_i32_value(r, 1);

  /*  isqrt(3) = 1.  */
  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  check_i32_value(r, 1);

  /*  Negative input returns EDOM.  */
  CHECK_EQ_I(arbint_set_i32(a, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_EDOM);

  CHECK_EQ_I(arbint_set_i32(a, -100), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_EDOM);

  /*  NULL inputs.  */
  CHECK_EQ_I(arbint_isqrt(NULL, a), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_isqrt(r, NULL), ARBINT_EINVAL);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_isqrt_perfect_squares(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  isqrt(4) = 2.  */
  CHECK_EQ_I(arbint_set_i32(a, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  check_i32_value(r, 2);

  /*  isqrt(9) = 3.  */
  CHECK_EQ_I(arbint_set_i32(a, 9), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  check_i32_value(r, 3);

  /*  isqrt(16) = 4.  */
  CHECK_EQ_I(arbint_set_i32(a, 16), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  check_i32_value(r, 4);

  /*  isqrt(25) = 5.  */
  CHECK_EQ_I(arbint_set_i32(a, 25), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  check_i32_value(r, 5);

  /*  isqrt(100) = 10.  */
  CHECK_EQ_I(arbint_set_i32(a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  check_i32_value(r, 10);

  /*  isqrt(10000) = 100.  */
  CHECK_EQ_I(arbint_set_i32(a, 10000), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  check_i32_value(r, 100);

  /*  isqrt(1000000) = 1000.  */
  CHECK_EQ_I(arbint_set_u32(a, 1000000u), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  check_i32_value(r, 1000);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_isqrt_non_perfect(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  isqrt(5) = 2.  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  check_i32_value(r, 2);
  check_isqrt_invariant(r, a);

  /*  isqrt(8) = 2.  */
  CHECK_EQ_I(arbint_set_i32(a, 8), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  check_i32_value(r, 2);
  check_isqrt_invariant(r, a);

  /*  isqrt(99) = 9.  */
  CHECK_EQ_I(arbint_set_i32(a, 99), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  check_i32_value(r, 9);
  check_isqrt_invariant(r, a);

  /*  isqrt(101) = 10.  */
  CHECK_EQ_I(arbint_set_i32(a, 101), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  check_i32_value(r, 10);
  check_isqrt_invariant(r, a);

  /*  isqrt(2^31 - 1) = 46340.  */
  CHECK_EQ_I(arbint_set_i32(a, 2147483647), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  check_i32_value(r, 46340);
  check_isqrt_invariant(r, a);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_isqrt_aliasing(void) {
  arbint_ctx_t ctx;
  arbint_t a;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);

  /*  rop == a aliasing: isqrt(100) in place.  */
  CHECK_EQ_I(arbint_set_i32(a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(a, a), ARBINT_OK);
  check_i32_value(a, 10);

  /*  rop == a aliasing with non-perfect square.  */
  CHECK_EQ_I(arbint_set_i32(a, 99), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(a, a), ARBINT_OK);
  check_i32_value(a, 9);

  /*  rop == a aliasing with 0.  */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(a, a), ARBINT_OK);
  check_i32_value(a, 0);

  /*  rop == a aliasing with 1.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(a, a), ARBINT_OK);
  check_i32_value(a, 1);

  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_isqrt_large(void) {
  arbint_ctx_t ctx;
  arbint_t a, r, root;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(root, &ctx), ARBINT_OK);

  /*  isqrt(2^128) = 2^64.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 128u), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);

  CHECK_EQ_I(arbint_set_i32(root, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(root, root, 64u), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, root), 0);

  /*  isqrt(2^256) = 2^128.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 256u), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);

  CHECK_EQ_I(arbint_set_i32(root, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(root, root, 128u), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, root), 0);

  /*  Cross-check: build a large number via squaring, verify isqrt
      recovers the root. root = 2^512 + 7, a = root^2, isqrt(a) = root.  */
  CHECK_EQ_I(arbint_set_i32(root, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(root, root, 512u), ARBINT_OK);
  CHECK_EQ_I(arbint_add_i32(root, root, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(a, root), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, root), 0);

  /*  Non-perfect large square: isqrt(root^2 + 1) = root.  */
  CHECK_EQ_I(arbint_add_i32(a, a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, root), 0);
  check_isqrt_invariant(r, a);

  /*  Non-perfect large square: isqrt(root^2 - 1) = root - 1.  */
  CHECK_EQ_I(arbint_sub_i32(a, a, 2), ARBINT_OK); /*  a = root^2 - 1  */
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_i32(root, root, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, root), 0);
  check_isqrt_invariant(r, a);

  /*  isqrt(2^1024 - 1): verify invariant.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 1024u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_i32(a, a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);
  check_isqrt_invariant(r, a);

  arbint_clear(root);
  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_isqrt_consecutive(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;
  int32_t i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Exhaustive check for 0..1000: verify invariant for every value.  */
  for (i = 0; i <= 1000 && g_failures == 0; ++i) {
    CHECK_EQ_I(arbint_set_i32(a, i), ARBINT_OK);
    CHECK_EQ_I(arbint_isqrt(r, a), ARBINT_OK);

    if (i > 0)
      check_isqrt_invariant(r, a);
    else
      check_i32_value(r, 0);
  }

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

int main(void) {
  test_isqrt_edge_cases();
  test_isqrt_perfect_squares();
  test_isqrt_non_perfect();
  test_isqrt_aliasing();
  test_isqrt_large();
  test_isqrt_consecutive();
  ARBINT_TEST_FINISH("test_isqrt");
}
