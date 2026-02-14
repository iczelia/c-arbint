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

/*  Test suite for Burnikel-Ziegler division.

    Tests division correctness across the BZ threshold boundaries and
    various operand sizes. Verifies that n = q*d + r and 0 <= r < d.  */

#include "test_framework.h"

#include <stdint.h>

ARBINT_TEST_DECLARE_FAILURES();

/*  Helper: verify n = q*d + r and 0 <= r < |d|.  */
static int verify_division(arbint_ctx_t * ctx, const arbint_t n,
                            const arbint_t d, const arbint_t q,
                            const arbint_t r) {
  arbint_t qd, check;
  int ok = 1;

  arbint_init(qd, ctx);
  arbint_init(check, ctx);

  /*  check = q * d + r  */
  arbint_mul(qd, q, d);
  arbint_add(check, qd, r);

  /*  Verify n == check.  */
  if (arbint_cmp(check, n) != 0)
    ok = 0;

  /*  Verify |r| < |d| (for truncated division).  */
  if (arbint_cmpabs(r, d) >= 0)
    ok = 0;

  arbint_clear(check);
  arbint_clear(qd);
  return ok;
}

/*  Build a large number by repeated squaring until target size reached.  */
static void build_large(arbint_t x, uint32_t seed, size_t target_limbs,
                         arbint_ctx_t * ctx) {
  (void) ctx;
  arbint_set_u32(x, seed);
  while (arbint_abs_sz(x[0]._sz) < target_limbs) {
    arbint_mul(x, x, x);
    arbint_add_u32(x, x, seed);
  }
}

/*  Test 1: Small divisions (below BZ threshold).
    These should still work via fallback to Knuth.  */
static void test_small_divisions(void) {
  arbint_ctx_t ctx;
  arbint_t n, d, q, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Test with ~10 limbs (below 40-limb threshold).  */
  build_large(n, 12345u, 20u, &ctx);
  build_large(d, 67890u, 10u, &ctx);

  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  CHECK(verify_division(&ctx, n, d, q, r));

  arbint_clear(r);
  arbint_clear(q);
  arbint_clear(d);
  arbint_clear(n);
}

/*  Test 2: Threshold boundary tests.
    Test at dn = 39, 40, 41 to verify dispatch works correctly.  */
static void test_threshold_boundaries(void) {
  arbint_ctx_t ctx;
  arbint_t n, d, q, r;
  size_t dn_vals[] = {39u, 40u, 41u, 63u, 64u, 65u};
  size_t i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  for (i = 0u; i < sizeof(dn_vals) / sizeof(dn_vals[0]); ++i) {
    size_t dn = dn_vals[i];
    size_t nn = 2u * dn;

    build_large(d, 12345u + (uint32_t) i, dn, &ctx);
    build_large(n, 67890u + (uint32_t) i, nn, &ctx);

    CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
    CHECK(verify_division(&ctx, n, d, q, r));
  }

  arbint_clear(r);
  arbint_clear(q);
  arbint_clear(d);
  arbint_clear(n);
}

/*  Test 3: Balanced divisions (nn = 2*dn).  */
static void test_balanced_divisions(void) {
  arbint_ctx_t ctx;
  arbint_t n, d, q, r;
  size_t dn_vals[] = {40u, 50u, 64u, 80u, 100u, 128u};
  size_t i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  for (i = 0u; i < sizeof(dn_vals) / sizeof(dn_vals[0]); ++i) {
    size_t dn = dn_vals[i];

    build_large(d, 11111u + (uint32_t) i, dn, &ctx);
    build_large(n, 22222u + (uint32_t) i, 2u * dn, &ctx);

    CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
    CHECK(verify_division(&ctx, n, d, q, r));
  }

  arbint_clear(r);
  arbint_clear(q);
  arbint_clear(d);
  arbint_clear(n);
}

/*  Test 4: Unbalanced divisions (nn > 2*dn).  */
static void test_unbalanced_divisions(void) {
  arbint_ctx_t ctx;
  arbint_t n, d, q, r;
  size_t ratios[] = {3u, 4u, 5u};
  size_t dn = 50u;
  size_t i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  build_large(d, 33333u, dn, &ctx);

  for (i = 0u; i < sizeof(ratios) / sizeof(ratios[0]); ++i) {
    size_t nn = ratios[i] * dn;

    build_large(n, 44444u + (uint32_t) i, nn, &ctx);

    CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
    CHECK(verify_division(&ctx, n, d, q, r));
  }

  arbint_clear(r);
  arbint_clear(q);
  arbint_clear(d);
  arbint_clear(n);
}

/*  Test 5: Exact multiples (remainder = 0).  */
static void test_exact_multiples(void) {
  arbint_ctx_t ctx;
  arbint_t n, d, q, r, expected_q;
  size_t dn = 50u;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected_q, &ctx), ARBINT_OK);

  build_large(d, 55555u, dn, &ctx);
  build_large(expected_q, 66666u, dn, &ctx);

  /*  n = expected_q * d (exact multiple).  */
  arbint_mul(n, expected_q, d);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(q, expected_q), 0);
  CHECK_EQ_I(arbint_signum(r), 0);

  /*  n = expected_q * d + 1 (remainder = 1).  */
  arbint_add_u32(n, n, 1u);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(q, expected_q), 0);
  check_u32_value(r, 1u);

  arbint_clear(expected_q);
  arbint_clear(r);
  arbint_clear(q);
  arbint_clear(d);
  arbint_clear(n);
}

/*  Test 6: Stochastic testing with random operands.  */
static void test_stochastic(void) {
  arbint_ctx_t ctx;
  arbint_t n, d, q, r;
  uint32_t seed = 77777u;
  size_t i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  for (i = 0u; i < 50u; ++i) {
    size_t dn = 40u + (seed % 100u);
    size_t nn = dn + (seed % (2u * dn));
    seed = seed * 1103515245u + 12345u;

    build_large(d, seed, dn, &ctx);
    seed = seed * 1103515245u + 12345u;
    build_large(n, seed, nn, &ctx);
    seed = seed * 1103515245u + 12345u;

    CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
    CHECK(verify_division(&ctx, n, d, q, r));
  }

  arbint_clear(r);
  arbint_clear(q);
  arbint_clear(d);
  arbint_clear(n);
}

/*  Test 7: Odd divisor sizes.  */
static void test_odd_divisor_sizes(void) {
  arbint_ctx_t ctx;
  arbint_t n, d, q, r;
  size_t dn_vals[] = {41u, 43u, 45u, 47u, 49u, 51u, 63u, 65u};
  size_t i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  for (i = 0u; i < sizeof(dn_vals) / sizeof(dn_vals[0]); ++i) {
    size_t dn = dn_vals[i];

    build_large(d, 88888u + (uint32_t) i, dn, &ctx);
    build_large(n, 99999u + (uint32_t) i, 2u * dn, &ctx);

    CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
    CHECK(verify_division(&ctx, n, d, q, r));
  }

  arbint_clear(r);
  arbint_clear(q);
  arbint_clear(d);
  arbint_clear(n);
}

int main(void) {
  ARBINT_TEST_START();
  test_small_divisions();
  test_threshold_boundaries();
  test_balanced_divisions();
  test_unbalanced_divisions();
  test_exact_multiples();
  test_stochastic();
  test_odd_divisor_sizes();
  ARBINT_TEST_FINISH("test_div_bz");
}
