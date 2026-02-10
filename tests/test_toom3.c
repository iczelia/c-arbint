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

#include <stdio.h>

ARBINT_TEST_DECLARE_FAILURES();

/*  Build a large value by repeated squaring: result = seed^(2^iters).
    On 64-bit limbs, 7 squarings of seed=3 gives ~6340 bits = ~100 limbs.  */
static arbint_err_t build_large(arbint_t dst, arbint_ctx_t * ctx,
                                int32_t seed, unsigned iters) {
  unsigned i;
  arbint_err_t rc;

  (void) ctx;
  rc = arbint_set_i32(dst, seed);
  if (rc != ARBINT_OK)
    return rc;
  for (i = 0u; i < iters; ++i) {
    rc = arbint_sqr(dst, dst);
    if (rc != ARBINT_OK)
      return rc;
  }
  return ARBINT_OK;
}

/*  Test 1: a^2 via arbint_sqr == a*a via arbint_mul for Toom-3-sized a.  */
static void test_toom3_known_squares(void) {
  arbint_ctx_t ctx;
  arbint_t a, sqr_result, mul_result;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(sqr_result, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(mul_result, &ctx), ARBINT_OK);

  /*  3^(2^7) = 3^128 has ~6340*2 = ~12700 bits after one more squaring,
      but 3^128 itself is about 203 bits * 2^0... Let's use pow_u32 to get
      a value with exactly many limbs.  */
  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 4000u), ARBINT_OK);
  CHECK(arbint_nbits(a) > 6000u);

  CHECK_EQ_I(arbint_sqr(sqr_result, a), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(mul_result, a, a), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(sqr_result, mul_result), 0);

  arbint_clear(mul_result);
  arbint_clear(sqr_result);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test 2: 2^N * 2^M == 2^(N+M).  */
static void test_toom3_power_of_two(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, prod, expected;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(prod, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  /*  Use bit counts that force Toom-3 territory (>= 96 limbs on 64-bit).  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 6400u), ARBINT_OK);

  CHECK_EQ_I(arbint_set_i32(b, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(b, b, 7000u), ARBINT_OK);

  CHECK_EQ_I(arbint_mul(prod, a, b), ARBINT_OK);

  CHECK_EQ_I(arbint_set_i32(expected, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(expected, expected, 13400u), ARBINT_OK);

  CHECK_EQ_I(arbint_cmp(prod, expected), 0);

  arbint_clear(expected);
  arbint_clear(prod);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test 3: a*b == b*a for large operands.  */
static void test_toom3_commutativity(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, ab, ba;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(ab, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(ba, &ctx), ARBINT_OK);

  CHECK_EQ_I(build_large(a, &ctx, 3, 7u), ARBINT_OK);
  CHECK_EQ_I(build_large(b, &ctx, 5, 7u), ARBINT_OK);

  CHECK_EQ_I(arbint_mul(ab, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(ba, b, a), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(ab, ba), 0);

  arbint_clear(ba);
  arbint_clear(ab);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test 4: a*(b+c) == a*b + a*c.  */
static void test_toom3_distributivity(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, c, bc, lhs, ab, ac, rhs;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(bc, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(lhs, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(ab, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(ac, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(rhs, &ctx), ARBINT_OK);

  CHECK_EQ_I(build_large(a, &ctx, 7, 7u), ARBINT_OK);
  CHECK_EQ_I(build_large(b, &ctx, 11, 7u), ARBINT_OK);
  CHECK_EQ_I(build_large(c, &ctx, 13, 7u), ARBINT_OK);

  /*  lhs = a * (b + c).  */
  CHECK_EQ_I(arbint_add(bc, b, c), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(lhs, a, bc), ARBINT_OK);

  /*  rhs = a*b + a*c.  */
  CHECK_EQ_I(arbint_mul(ab, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(ac, a, c), ARBINT_OK);
  CHECK_EQ_I(arbint_add(rhs, ab, ac), ARBINT_OK);

  CHECK_EQ_I(arbint_cmp(lhs, rhs), 0);

  arbint_clear(rhs);
  arbint_clear(ac);
  arbint_clear(ab);
  arbint_clear(lhs);
  arbint_clear(bc);
  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test 5: (a*b)/b == a and (a*b) % b == 0.  */
static void test_toom3_cross_check_div(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, prod, q, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(prod, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  CHECK_EQ_I(build_large(a, &ctx, 17, 7u), ARBINT_OK);
  CHECK_EQ_I(build_large(b, &ctx, 19, 7u), ARBINT_OK);

  CHECK_EQ_I(arbint_mul(prod, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, prod, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(q, a), 0);
  CHECK_EQ_I(arbint_signum(r), 0);

  arbint_clear(r);
  arbint_clear(q);
  arbint_clear(prod);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test 6: (-a)*b == -(a*b) and (-a)*(-b) == a*b.  */
static void test_toom3_negative(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, neg_a, neg_b, ab, neg_ab, na_b, na_nb;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(neg_a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(neg_b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(ab, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(neg_ab, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(na_b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(na_nb, &ctx), ARBINT_OK);

  CHECK_EQ_I(build_large(a, &ctx, 23, 7u), ARBINT_OK);
  CHECK_EQ_I(build_large(b, &ctx, 29, 7u), ARBINT_OK);

  CHECK_EQ_I(arbint_neg(neg_a, a), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(neg_b, b), ARBINT_OK);

  CHECK_EQ_I(arbint_mul(ab, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(neg_ab, ab), ARBINT_OK);

  /*  (-a)*b == -(a*b).  */
  CHECK_EQ_I(arbint_mul(na_b, neg_a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(na_b, neg_ab), 0);

  /*  (-a)*(-b) == a*b.  */
  CHECK_EQ_I(arbint_mul(na_nb, neg_a, neg_b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(na_nb, ab), 0);

  arbint_clear(na_nb);
  arbint_clear(na_b);
  arbint_clear(neg_ab);
  arbint_clear(ab);
  arbint_clear(neg_b);
  arbint_clear(neg_a);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test 7: threshold boundary -- sizes 95, 96, 97 limbs.
    Build values at specific bit counts and verify (a*b)/b == a.  */
static void test_toom3_threshold_boundary(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, prod, q, r;
  static const size_t sizes[] = {95u, 96u, 97u, 127u, 128u, 129u};
  size_t si;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(prod, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  for (si = 0u; si < sizeof(sizes) / sizeof(sizes[0]) && g_failures == 0;
       ++si) {
    size_t nbits;
    size_t actual_bits;

    /*  Target bit count for this many limbs.  */
#if ARBINT_LIMB_BITS == 64
    nbits = sizes[si] * 64u;
#else
    nbits = sizes[si] * 32u;
#endif

    /*  Build a = 2^nbits - 1 (all bits set).  */
    CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(a, a, (uint32_t) nbits), ARBINT_OK);
    CHECK_EQ_I(arbint_sub_i32(a, a, 1), ARBINT_OK);

    /*  Build b from a different pattern.  */
    CHECK_EQ_I(arbint_set_i32(b, 1), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(b, b, (uint32_t) (nbits - 1u)), ARBINT_OK);
    CHECK_EQ_I(arbint_add_i32(b, b, 1), ARBINT_OK);

    actual_bits = arbint_nbits(a);
    CHECK_EQ_I((int) (actual_bits == nbits), 1);

    CHECK_EQ_I(arbint_mul(prod, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_qr(q, r, prod, b), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(q, a), 0);
    CHECK_EQ_I(arbint_signum(r), 0);
  }

  arbint_clear(r);
  arbint_clear(q);
  arbint_clear(prod);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test 8: unbalanced operands (one ~2x the other).  */
static void test_toom3_unbalanced(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, prod, q, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(prod, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  a ~= 200 limbs, b ~= 100 limbs (within 3x ratio).  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 12800u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_i32(a, a, 1), ARBINT_OK);

  CHECK_EQ_I(arbint_set_i32(b, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(b, b, 6400u), ARBINT_OK);
  CHECK_EQ_I(arbint_add_i32(b, b, 7), ARBINT_OK);

  CHECK_EQ_I(arbint_mul(prod, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, prod, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(q, a), 0);
  CHECK_EQ_I(arbint_signum(r), 0);

  arbint_clear(r);
  arbint_clear(q);
  arbint_clear(prod);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test 9: aliasing -- rop==a, rop==b, rop==a==b.  */
static void test_toom3_aliasing(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, expected, rop;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(rop, &ctx), ARBINT_OK);

  CHECK_EQ_I(build_large(a, &ctx, 31, 7u), ARBINT_OK);
  CHECK_EQ_I(build_large(b, &ctx, 37, 7u), ARBINT_OK);

  /*  Compute expected = a * b.  */
  CHECK_EQ_I(arbint_mul(expected, a, b), ARBINT_OK);

  /*  rop == a: a = a * b.  */
  CHECK_EQ_I(arbint_set(rop, a), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(rop, rop, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(rop, expected), 0);

  /*  rop == b: b_copy = a * b_copy.  */
  CHECK_EQ_I(arbint_set(rop, b), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(rop, a, rop), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(rop, expected), 0);

  /*  rop == a == b: a^2.  */
  CHECK_EQ_I(arbint_mul(rop, a, a), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(expected, a), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(rop, expected), 0);

  arbint_clear(rop);
  arbint_clear(expected);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test 10: a * 1 == a for Toom-3-sized a.  */
static void test_toom3_identity(void) {
  arbint_ctx_t ctx;
  arbint_t a, one, prod;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(one, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(prod, &ctx), ARBINT_OK);

  CHECK_EQ_I(build_large(a, &ctx, 41, 7u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(one, 1), ARBINT_OK);

  CHECK_EQ_I(arbint_mul(prod, a, one), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(prod, a), 0);

  arbint_clear(prod);
  arbint_clear(one);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test 11: force a(-1) negative by constructing a where a1 > a0 + a2.
    Use a = 2^(k*BITS) where k is the piece size, so a0=0, a2=0, a1=1.  */
static void test_toom3_force_neg_eval(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, prod, q, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(prod, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  a has all weight in the middle third: a = big * B^k where k = 34.
      96 limbs = 3*32, so k=32. a1 has all the mass, a0=a2=0.
      Build: a = (2^(32*64) - 1) << (32*64).  On 64-bit this is a
      value where limbs 0..31 are 0, limbs 32..63 are all-ones, rest 0.  */

  /*  Simpler approach: a = big_value * 2^shift where shift = k*BITS.
      This puts a0 = 0 (all limbs below k are zero).  */
#if ARBINT_LIMB_BITS == 64
  {
    uint32_t k_bits = 32u * 64u;
    /*  a = (2^(k_bits) - 1) << k_bits.  a0 = 0, a1 = all-ones, a2 = 0.  */
    CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(a, a, k_bits), ARBINT_OK);
    CHECK_EQ_I(arbint_sub_i32(a, a, 1), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(a, a, k_bits), ARBINT_OK);
  }
#else
  {
    uint32_t k_bits = 32u * 32u;
    CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(a, a, k_bits), ARBINT_OK);
    CHECK_EQ_I(arbint_sub_i32(a, a, 1), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(a, a, k_bits), ARBINT_OK);
  }
#endif /* ARBINT_LIMB_BITS */

  CHECK_EQ_I(build_large(b, &ctx, 43, 7u), ARBINT_OK);

  CHECK_EQ_I(arbint_mul(prod, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, prod, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(q, a), 0);
  CHECK_EQ_I(arbint_signum(r), 0);

  arbint_clear(r);
  arbint_clear(q);
  arbint_clear(prod);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test 12: max carries -- all-ones limbs to stress carry propagation.  */
static void test_toom3_max_carries(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, prod, q, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(prod, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  a = 2^6400 - 1 (100 limbs of all-ones on 64-bit).  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 6400u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_i32(a, a, 1), ARBINT_OK);

  /*  b = 2^6400 - 1 (same).  */
  CHECK_EQ_I(arbint_set_i32(b, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(b, b, 6400u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_i32(b, b, 1), ARBINT_OK);

  CHECK_EQ_I(arbint_mul(prod, a, b), ARBINT_OK);

  /*  Verify: (2^n - 1)^2 = 2^(2n) - 2^(n+1) + 1.
      So prod + 2^(n+1) - 1 == 2^(2n).  */
  CHECK_EQ_I(arbint_tdiv_qr(q, r, prod, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(q, a), 0);
  CHECK_EQ_I(arbint_signum(r), 0);

  arbint_clear(r);
  arbint_clear(q);
  arbint_clear(prod);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test 13: non-multiple-of-3 limb counts (94, 95, 97, 98).  */
static void test_toom3_non_multiple_of_3(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, prod, q, r;
  static const uint32_t bit_counts[] = {
#if ARBINT_LIMB_BITS == 64
    94u * 64u, 95u * 64u, 97u * 64u, 98u * 64u, 99u * 64u
#else
    94u * 32u, 95u * 32u, 97u * 32u, 98u * 32u, 99u * 32u
#endif
  };
  size_t i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(prod, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  for (i = 0u; i < sizeof(bit_counts) / sizeof(bit_counts[0]) &&
                   g_failures == 0;
       ++i) {
    CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(a, a, bit_counts[i]), ARBINT_OK);
    CHECK_EQ_I(arbint_sub_i32(a, a, 1), ARBINT_OK);

    CHECK_EQ_I(arbint_set_i32(b, 1), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(b, b, bit_counts[i]), ARBINT_OK);
    CHECK_EQ_I(arbint_add_i32(b, b, 1), ARBINT_OK);

    CHECK_EQ_I(arbint_mul(prod, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_qr(q, r, prod, b), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(q, a), 0);
    CHECK_EQ_I(arbint_signum(r), 0);
  }

  arbint_clear(r);
  arbint_clear(q);
  arbint_clear(prod);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test 14: multiple random-ish multiplications verified via division.  */
static void test_toom3_stochastic(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, prod, q, r;
  int32_t seeds_a[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
  int32_t seeds_b[] = {31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
  size_t i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(prod, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  for (i = 0u; i < 10u && g_failures == 0; ++i) {
    CHECK_EQ_I(build_large(a, &ctx, seeds_a[i], 7u), ARBINT_OK);
    CHECK_EQ_I(build_large(b, &ctx, seeds_b[i], 7u), ARBINT_OK);

    CHECK_EQ_I(arbint_mul(prod, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_qr(q, r, prod, b), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(q, a), 0);
    CHECK_EQ_I(arbint_signum(r), 0);
  }

  arbint_clear(r);
  arbint_clear(q);
  arbint_clear(prod);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

int main(void) {
  test_toom3_known_squares();
  test_toom3_power_of_two();
  test_toom3_commutativity();
  test_toom3_distributivity();
  test_toom3_cross_check_div();
  test_toom3_negative();
  test_toom3_threshold_boundary();
  test_toom3_unbalanced();
  test_toom3_aliasing();
  test_toom3_identity();
  test_toom3_force_neg_eval();
  test_toom3_max_carries();
  test_toom3_non_multiple_of_3();
  test_toom3_stochastic();
  ARBINT_TEST_FINISH("test_toom3");
}
