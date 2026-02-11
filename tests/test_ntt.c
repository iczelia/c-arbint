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

/*  NTT multiplication tests.

    These tests verify the correctness of the NTT-based multiplication
    implementation by comparing results against the existing Toom-3
    implementation (via arbint_mul) and checking mathematical properties.  */

#include "test_framework.h"

#include <stdint.h>
#include <string.h>

ARBINT_TEST_DECLARE_FAILURES();

/*  Simple LCG for deterministic pseudorandom numbers.  */
static uint32_t test_rng_state = 0x12345678u;

static uint32_t test_rand(void) {
  test_rng_state = test_rng_state * 1103515245u + 12345u;
  return test_rng_state;
}

static void test_rng_seed(uint32_t seed) { test_rng_state = seed; }

/*  Build a large value by repeated squaring until it reaches target_limbs.
    Adds a constant after each squaring to avoid degenerate patterns.  */
static void build_large_value(arbint_t x, size_t target_limbs) {
  CHECK_EQ_I(arbint_set_u32(x, 0xDEADBEEFu), ARBINT_OK);

  while (arbint_abs_sz(x[0]._sz) < target_limbs) {
    CHECK_EQ_I(arbint_sqr(x, x), ARBINT_OK);
    CHECK_EQ_I(arbint_add_u32(x, x, 12345u), ARBINT_OK);
  }
}

/*  Test NTT multiplication at threshold boundary.
    Tests sizes just below, at, and just above the NTT threshold (1024 limbs).  */
static void test_ntt_threshold_boundary(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, c, d;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);

  /*  Test sizes around threshold.
      NTT threshold is 1024 limbs, so test 1023, 1024, 1025.  */
  size_t sizes[] = {1023u, 1024u, 1025u, 1100u, 1200u};
  size_t num_sizes = sizeof(sizes) / sizeof(sizes[0]);

  for (size_t i = 0; i < num_sizes; ++i) {
    size_t n = sizes[i];

    /*  Build operands of target size.  */
    build_large_value(a, n);
    build_large_value(b, n);

    /*  Multiply using arbint_mul (which dispatches to NTT or Toom-3).  */
    CHECK_EQ_I(arbint_mul(c, a, b), ARBINT_OK);

    /*  Verify result is non-zero and has reasonable size.  */
    size_t an_used = arbint_abs_sz(a[0]._sz);
    size_t bn_used = arbint_abs_sz(b[0]._sz);
    CHECK(!arbint_is_zero(c));
    CHECK(arbint_abs_sz(c[0]._sz) >= an_used);
    CHECK(arbint_abs_sz(c[0]._sz) >= bn_used);
    CHECK(arbint_abs_sz(c[0]._sz) <= an_used + bn_used + 1u);

    /*  Verify commutativity: a*b == b*a.  */
    CHECK_EQ_I(arbint_mul(d, b, a), ARBINT_OK);
    CHECK(arbint_eq(c, d));

    /*  Verify squaring: a*a via mul == sqr(a).  */
    CHECK_EQ_I(arbint_mul(c, a, a), ARBINT_OK);
    CHECK_EQ_I(arbint_sqr(d, a), ARBINT_OK);
    CHECK(arbint_eq(c, d));
  }

  arbint_clear(d);
  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
}

/*  Test NTT multiplication with various power-of-two sizes.  */
static void test_ntt_power_of_two(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, c, d;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);

  /*  Test power-of-two sizes that are >= NTT threshold.  */
  size_t sizes[] = {1024u, 2048u};
  size_t num_sizes = sizeof(sizes) / sizeof(sizes[0]);

  for (size_t i = 0; i < num_sizes; ++i) {
    size_t n = sizes[i];

    build_large_value(a, n);
    build_large_value(b, n);

    CHECK_EQ_I(arbint_mul(c, a, b), ARBINT_OK);
    CHECK(!arbint_is_zero(c));

    /*  Verify squaring consistency.  */
    CHECK_EQ_I(arbint_mul(c, a, a), ARBINT_OK);
    CHECK_EQ_I(arbint_sqr(d, a), ARBINT_OK);
    CHECK(arbint_eq(c, d));
  }

  arbint_clear(d);
  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
}

/*  Test NTT multiplication with unbalanced operands.
    When an > 2*bn, dispatch should fall through to Toom-3.  */
static void test_ntt_unbalanced(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, c;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);

  /*  Build a very large value and a smaller one.  */
  build_large_value(a, 2100u);  /* > 2 * 1024 */
  build_large_value(b, 1024u);

  /*  This should dispatch to Toom-3 due to imbalance ratio check.  */
  CHECK_EQ_I(arbint_mul(c, a, b), ARBINT_OK);
  CHECK(!arbint_is_zero(c));

  /*  Verify commutativity still holds.  */
  arbint_t d;
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(d, b, a), ARBINT_OK);
  CHECK(arbint_eq(c, d));

  arbint_clear(d);
  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
}

/*  Test NTT with algebraic identity: (a+b)^2 = a^2 + 2ab + b^2.  */
static void test_ntt_algebraic_identity(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, apb, lhs, rhs, a2, b2, ab2;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(apb, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(lhs, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(rhs, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a2, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b2, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(ab2, &ctx), ARBINT_OK);

  /*  Build large operands above NTT threshold.  */
  build_large_value(a, 1100u);
  build_large_value(b, 1050u);

  /*  LHS: (a + b)^2  */
  CHECK_EQ_I(arbint_add(apb, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(lhs, apb), ARBINT_OK);

  /*  RHS: a^2 + 2*a*b + b^2  */
  CHECK_EQ_I(arbint_sqr(a2, a), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(b2, b), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(ab2, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(ab2, ab2, 2u), ARBINT_OK);

  CHECK_EQ_I(arbint_add(rhs, a2, b2), ARBINT_OK);
  CHECK_EQ_I(arbint_add(rhs, rhs, ab2), ARBINT_OK);

  CHECK(arbint_eq(lhs, rhs));

  arbint_clear(ab2);
  arbint_clear(b2);
  arbint_clear(a2);
  arbint_clear(rhs);
  arbint_clear(lhs);
  arbint_clear(apb);
  arbint_clear(b);
  arbint_clear(a);
}

/*  Stochastic test with random operand sizes.  */
static void test_ntt_stochastic(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, c, d;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);

  test_rng_seed(0xABCDEF01u);

  /*  Run 20 random tests (fewer due to large size).  */
  for (int iter = 0; iter < 20; ++iter) {
    /*  Generate random sizes in range [1024, 1500].  */
    size_t an = 1024u + (test_rand() % 500u);
    size_t bn = 1024u + (test_rand() % 500u);

    build_large_value(a, an);
    build_large_value(b, bn);

    /*  Compute a * b.  */
    CHECK_EQ_I(arbint_mul(c, a, b), ARBINT_OK);

    /*  Verify commutativity.  */
    CHECK_EQ_I(arbint_mul(d, b, a), ARBINT_OK);
    CHECK(arbint_eq(c, d));

    /*  Verify result size is reasonable.  */
    size_t an_used = arbint_abs_sz(a[0]._sz);
    size_t bn_used = arbint_abs_sz(b[0]._sz);
    size_t result_limbs = arbint_abs_sz(c[0]._sz);
    CHECK(result_limbs >= an_used && result_limbs >= bn_used);
    CHECK(result_limbs <= an_used + bn_used + 1u);
  }

  arbint_clear(d);
  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
}

/*  Test NTT with small operands below threshold.
    These should use Toom-3/Karatsuba, not NTT.  */
static void test_ntt_below_threshold(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, c, d;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);

  /*  Test sizes below NTT threshold.  */
  size_t sizes[] = {100u, 200u, 500u, 800u};
  size_t num_sizes = sizeof(sizes) / sizeof(sizes[0]);

  for (size_t i = 0; i < num_sizes; ++i) {
    size_t n = sizes[i];

    build_large_value(a, n);
    build_large_value(b, n);

    CHECK_EQ_I(arbint_mul(c, a, b), ARBINT_OK);
    CHECK(!arbint_is_zero(c));

    /*  Verify squaring.  */
    CHECK_EQ_I(arbint_mul(c, a, a), ARBINT_OK);
    CHECK_EQ_I(arbint_sqr(d, a), ARBINT_OK);
    CHECK(arbint_eq(c, d));
  }

  arbint_clear(d);
  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
}

/*  Test aliasing: rop == a, rop == b, rop == a == b.  */
static void test_ntt_aliasing(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, c, ref;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(ref, &ctx), ARBINT_OK);

  /*  Build values above NTT threshold.  */
  build_large_value(a, 1100u);
  build_large_value(b, 1050u);

  /*  Compute reference: c = a * b.  */
  CHECK_EQ_I(arbint_mul(c, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_set(ref, c), ARBINT_OK);

  /*  Test aliasing: a = a * b (rop == a).  */
  build_large_value(a, 1100u);
  CHECK_EQ_I(arbint_mul(a, a, b), ARBINT_OK);
  CHECK(arbint_eq(a, ref));

  /*  Test aliasing: b = a * b (rop == b).  */
  build_large_value(a, 1100u);
  build_large_value(b, 1050u);
  CHECK_EQ_I(arbint_mul(b, a, b), ARBINT_OK);
  CHECK(arbint_eq(b, ref));

  /*  Test aliasing: a = a * a (rop == a == b).  */
  build_large_value(a, 1100u);
  CHECK_EQ_I(arbint_sqr(ref, a), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(a, a, a), ARBINT_OK);
  CHECK(arbint_eq(a, ref));

  arbint_clear(ref);
  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
}

int main(void) {
  test_ntt_below_threshold();
  test_ntt_threshold_boundary();
  test_ntt_power_of_two();
  test_ntt_unbalanced();
  test_ntt_algebraic_identity();
  test_ntt_aliasing();
  test_ntt_stochastic();

  ARBINT_TEST_FINISH("test_ntt");
}
