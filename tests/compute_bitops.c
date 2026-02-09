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

/*  Comprehensive bitwise stress test.

    Algorithms implemented entirely in terms of bitwise operations,
    plus exhaustive identity verification across all sign quadrants.

    1. Bitwise ripple-carry adder: sum = a + b computed purely with
       xor/and/shl.  Validated against arbint_add for growing operands.

    2. Bitwise subtractor (two's complement): diff = a - b computed via
       not/add/inc.  Validated against arbint_sub.

    3. Binary GCD (Stein's algorithm): gcd(a, b) computed using only
       ctz/shr/sub/cmp/is_zero/swap.  Validated against Fibonacci identities
       and known GCD values.

    4. Identity gauntlet: 15 boolean algebra identities verified over all
       N*N operand pairs (both orderings) spanning all 4 sign quadrants
       with small, multi-limb, and asymmetric-size operands.

    5. Popcount / hamming cross-checks: hamming(a,b) = popcount(a^b),
       popcount(x) + popcount(not(x)) = W(x), triangle inequality,
       hamming invariance under XOR.

    6. Bit manipulation stress: setbit/clrbit/testbit for positive and
       negative, multi-limb, limb-boundary positions.

    7. CTZ / nbits / sizeinbase consistency.

    8. Shift-bitwise cross-checks: a << 1 = a + a, mask-via-shift.

    9. Aliasing stress for multi-limb signed operands.

   10. Immediate variant (_u32/_i32) cross-validation for multi-limb.  */

#include "test_framework.h"

#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

ARBINT_TEST_DECLARE_FAILURES();

/* ========== Local helper: pow_u32 via repeated squaring ========== */

static arbint_err_t local_pow_u32(arbint_t rop, const arbint_t base,
                                  uint32_t exp) {
  arbint_err_t rc;
  arbint_ctx_t * ctx = arbint_get_ctx(base);
  arbint_t acc, b;

  rc = arbint_init(acc, ctx);
  if (rc != ARBINT_OK)
    return rc;
  rc = arbint_init(b, ctx);
  if (rc != ARBINT_OK) {
    arbint_clear(acc);
    return rc;
  }

  rc = arbint_set_u32(acc, 1u);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_set(b, base);
  if (rc != ARBINT_OK)
    goto cleanup;

  while (exp != 0u) {
    if (exp & 1u) {
      rc = arbint_mul(acc, acc, b);
      if (rc != ARBINT_OK)
        goto cleanup;
    }
    exp >>= 1u;
    if (exp != 0u) {
      rc = arbint_sqr(b, b);
      if (rc != ARBINT_OK)
        goto cleanup;
    }
  }

  rc = arbint_set(rop, acc);

cleanup:
  arbint_clear(b);
  arbint_clear(acc);
  return rc;
}

/* ========== Section 1: Bitwise ripple-carry adder ========== */

/*  Add two non-negative magnitudes using only xor, and, shl.  */
static arbint_err_t bitwise_add_mag(arbint_t sum, const arbint_t a,
                                    const arbint_t b) {
  arbint_err_t rc;
  arbint_ctx_t * ctx = arbint_get_ctx(a);
  arbint_t carry, tmp;
  unsigned guard = 0u;

  rc = arbint_init(carry, ctx);
  if (rc != ARBINT_OK)
    return rc;
  rc = arbint_init(tmp, ctx);
  if (rc != ARBINT_OK) {
    arbint_clear(carry);
    return rc;
  }

  rc = arbint_abs(sum, a);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_abs(carry, b);
  if (rc != ARBINT_OK)
    goto cleanup;

  while (!arbint_is_zero(carry)) {
    rc = arbint_xor(tmp, sum, carry);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_and(carry, sum, carry);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_shl(carry, carry, 1u);
    if (rc != ARBINT_OK)
      goto cleanup;
    arbint_swap(sum, tmp);

    if (++guard > 100000u) {
      rc = ARBINT_EDOM;
      goto cleanup;
    }
  }

cleanup:
  arbint_clear(tmp);
  arbint_clear(carry);
  return rc;
}

static void test_bitwise_adder(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, t, sum_bit, sum_ref, base;
  size_t i;

  fprintf(stderr, "  bitwise adder: ");

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init_all(&ctx, a, b, t, sum_bit, sum_ref, base,
                              (arbint_t *) NULL),
             ARBINT_OK);

  /*  a = 2^256 - 1, b = 1.  Full carry chain.  */
  CHECK_EQ_I(arbint_set_u32(base, 2u), ARBINT_OK);
  CHECK_EQ_I(local_pow_u32(a, base, 256u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(a, a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_add(sum_ref, a, b), ARBINT_OK);
  CHECK_EQ_I(bitwise_add_mag(sum_bit, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(sum_ref, sum_bit), 0);

  /*  Both multi-limb: 3^100 + 5^80.  */
  CHECK_EQ_I(arbint_set_u32(base, 3u), ARBINT_OK);
  CHECK_EQ_I(local_pow_u32(a, base, 100u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(base, 5u), ARBINT_OK);
  CHECK_EQ_I(local_pow_u32(b, base, 80u), ARBINT_OK);
  CHECK_EQ_I(arbint_add(sum_ref, a, b), ARBINT_OK);
  CHECK_EQ_I(bitwise_add_mag(sum_bit, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(sum_ref, sum_bit), 0);

  /*  a + a = 2*a.  */
  CHECK_EQ_I(arbint_add(sum_ref, a, a), ARBINT_OK);
  CHECK_EQ_I(bitwise_add_mag(sum_bit, a, a), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(sum_ref, sum_bit), 0);

  /*  a + 0 = a.  */
  CHECK_EQ_I(arbint_set_u32(t, 0u), ARBINT_OK);
  CHECK_EQ_I(bitwise_add_mag(sum_bit, a, t), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(a, sum_bit), 0);

  /*  0 + b = b.  */
  CHECK_EQ_I(bitwise_add_mag(sum_bit, t, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(b, sum_bit), 0);

  /*  Small + large (very asymmetric).  */
  CHECK_EQ_I(arbint_set_u32(t, 42u), ARBINT_OK);
  CHECK_EQ_I(arbint_add(sum_ref, a, t), ARBINT_OK);
  CHECK_EQ_I(bitwise_add_mag(sum_bit, a, t), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(sum_ref, sum_bit), 0);

  /*  Iterated doubling: 40 iterations growing t.  */
  CHECK_EQ_I(arbint_set_u32(t, 1u), ARBINT_OK);
  for (i = 0u; i < 40u && g_failures == 0; ++i) {
    CHECK_EQ_I(arbint_add(sum_ref, t, b), ARBINT_OK);
    CHECK_EQ_I(bitwise_add_mag(sum_bit, t, b), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(sum_ref, sum_bit), 0);
    CHECK_EQ_I(arbint_add(t, t, t), ARBINT_OK);
  }

  /*  All-ones + all-ones: (2^512-1) + (2^512-1).  */
  CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 512u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(a, a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_add(sum_ref, a, a), ARBINT_OK);
  CHECK_EQ_I(bitwise_add_mag(sum_bit, a, a), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(sum_ref, sum_bit), 0);

  fprintf(stderr, "passed\n");

  arbint_clear_all(base, t, sum_bit, sum_ref, b, a, (arbint_t *) NULL);
  arbint_ctx_clear(&ctx);
}

/* ========== Section 2: Bitwise subtractor ========== */

/*  Compute |a| - |b| using bitwise NOT + adder.
    sub(a, b) = a + not(b) + 1 = a + (~b + 1).
    Precondition: |a| >= |b|.  */
static arbint_err_t bitwise_sub_mag(arbint_t diff, const arbint_t a,
                                    const arbint_t b) {
  arbint_err_t rc;
  arbint_ctx_t * ctx = arbint_get_ctx(a);
  arbint_t neg_b, one;

  rc = arbint_init(neg_b, ctx);
  if (rc != ARBINT_OK)
    return rc;
  rc = arbint_init(one, ctx);
  if (rc != ARBINT_OK) {
    arbint_clear(neg_b);
    return rc;
  }

  /*  neg_b = not(|b|) -- this gives -(|b|+1) via TC.  */
  rc = arbint_abs(neg_b, b);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_not(neg_b, neg_b);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  diff = |a| + not(|b|) + 1 = |a| + (-|b|-1) + 1 = |a| - |b|.  */
  rc = arbint_abs(diff, a);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_add(diff, diff, neg_b);
  if (rc != ARBINT_OK)
    goto cleanup;
  CHECK_EQ_I(arbint_set_u32(one, 1u), ARBINT_OK);
  rc = arbint_add(diff, diff, one);

cleanup:
  arbint_clear(one);
  arbint_clear(neg_b);
  return rc;
}

static void test_bitwise_subtractor(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, diff_bit, diff_ref, base;

  fprintf(stderr, "  bitwise subtractor: ");

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init_all(&ctx, a, b, diff_bit, diff_ref, base,
                              (arbint_t *) NULL),
             ARBINT_OK);

  /*  2^256 - 1.  */
  CHECK_EQ_I(arbint_set_u32(base, 2u), ARBINT_OK);
  CHECK_EQ_I(local_pow_u32(a, base, 256u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub(diff_ref, a, b), ARBINT_OK);
  CHECK_EQ_I(bitwise_sub_mag(diff_bit, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(diff_ref, diff_bit), 0);

  /*  3^100 - 5^80 (3^100 > 5^80).  */
  CHECK_EQ_I(arbint_set_u32(base, 3u), ARBINT_OK);
  CHECK_EQ_I(local_pow_u32(a, base, 100u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(base, 5u), ARBINT_OK);
  CHECK_EQ_I(local_pow_u32(b, base, 80u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub(diff_ref, a, b), ARBINT_OK);
  CHECK_EQ_I(bitwise_sub_mag(diff_bit, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(diff_ref, diff_bit), 0);

  /*  a - 0 = a.  */
  CHECK_EQ_I(arbint_set_u32(b, 0u), ARBINT_OK);
  CHECK_EQ_I(bitwise_sub_mag(diff_bit, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(a, diff_bit), 0);

  /*  a - a = 0.  */
  CHECK_EQ_I(bitwise_sub_mag(diff_bit, a, a), ARBINT_OK);
  CHECK(arbint_is_zero(diff_bit));

  /*  (2^1000 - 1) - (2^500 - 1) = 2^1000 - 2^500.  */
  CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 1000u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(a, a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(b, b, 500u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(b, b, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub(diff_ref, a, b), ARBINT_OK);
  CHECK_EQ_I(bitwise_sub_mag(diff_bit, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(diff_ref, diff_bit), 0);

  fprintf(stderr, "passed\n");

  arbint_clear_all(base, diff_bit, diff_ref, b, a, (arbint_t *) NULL);
  arbint_ctx_clear(&ctx);
}

/* ========== Section 3: Binary GCD (Stein's algorithm) ========== */

static arbint_err_t binary_gcd(arbint_t g, const arbint_t a,
                                const arbint_t b) {
  arbint_t u, v;
  arbint_err_t rc;
  arbint_ctx_t * ctx = arbint_get_ctx(a);
  size_t ctz_u, ctz_v, shift;
  unsigned guard = 0u;

  if (arbint_is_zero(a))
    return arbint_abs(g, b);
  if (arbint_is_zero(b))
    return arbint_abs(g, a);

  rc = arbint_init(u, ctx);
  if (rc != ARBINT_OK)
    return rc;
  rc = arbint_init(v, ctx);
  if (rc != ARBINT_OK) {
    arbint_clear(u);
    return rc;
  }

  rc = arbint_abs(u, a);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_abs(v, b);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_ctz(u, &ctz_u);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_ctz(v, &ctz_v);
  if (rc != ARBINT_OK)
    goto cleanup;
  shift = (ctz_u < ctz_v) ? ctz_u : ctz_v;

  rc = arbint_shr(u, u, (uint32_t) ctz_u);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_shr(v, v, (uint32_t) ctz_v);
  if (rc != ARBINT_OK)
    goto cleanup;

  while (!arbint_is_zero(v)) {
    int c = arbint_cmp(u, v);
    if (c == 0)
      break;
    if (c < 0)
      arbint_swap(u, v);

    rc = arbint_sub(u, u, v);
    if (rc != ARBINT_OK)
      goto cleanup;

    if (!arbint_is_zero(u)) {
      rc = arbint_ctz(u, &ctz_u);
      if (rc != ARBINT_OK)
        goto cleanup;
      rc = arbint_shr(u, u, (uint32_t) ctz_u);
      if (rc != ARBINT_OK)
        goto cleanup;
    }

    if (++guard > 1000000u) {
      rc = ARBINT_EDOM;
      goto cleanup;
    }
  }

  rc = arbint_shl(g, u, (uint32_t) shift);

cleanup:
  arbint_clear(v);
  arbint_clear(u);
  return rc;
}

static arbint_err_t fibonacci(arbint_t out, uint32_t n, arbint_ctx_t * ctx) {
  arbint_t a, b, tmp;
  arbint_err_t rc;
  uint32_t i;

  rc = arbint_init(a, ctx);
  if (rc != ARBINT_OK)
    return rc;
  rc = arbint_init(b, ctx);
  if (rc != ARBINT_OK) {
    arbint_clear(a);
    return rc;
  }
  rc = arbint_init(tmp, ctx);
  if (rc != ARBINT_OK) {
    arbint_clear(b);
    arbint_clear(a);
    return rc;
  }

  CHECK_EQ_I(arbint_set_u32(a, 0u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 1u), ARBINT_OK);

  if (n == 0u) {
    rc = arbint_set(out, a);
    goto cleanup;
  }

  for (i = 1u; i < n; ++i) {
    rc = arbint_add(tmp, a, b);
    if (rc != ARBINT_OK)
      goto cleanup;
    arbint_swap(a, b);
    arbint_swap(b, tmp);
  }

  rc = arbint_set(out, b);

cleanup:
  arbint_clear(tmp);
  arbint_clear(b);
  arbint_clear(a);
  return rc;
}

static void test_binary_gcd(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, g, expected, base;

  fprintf(stderr, "  binary GCD: ");

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init_all(&ctx, a, b, g, expected, base,
                              (arbint_t *) NULL),
             ARBINT_OK);

  /*  gcd(0, 0) = 0.  */
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  CHECK(arbint_is_zero(g));

  /*  gcd(12, 8) = 4.  */
  CHECK_EQ_I(arbint_set_u32(a, 12u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 8u), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 4u);

  /*  gcd(n, 0) = n, gcd(0, n) = n.  */
  CHECK_EQ_I(arbint_set_u32(a, 42u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 0u), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 42u);
  CHECK_EQ_I(binary_gcd(g, b, a), ARBINT_OK);
  check_u32_value(g, 42u);

  /*  gcd(n, n) = n.  */
  CHECK_EQ_I(arbint_set_u32(a, 100u), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, a), ARBINT_OK);
  check_u32_value(g, 100u);

  /*  gcd(2^128, 2^64) = 2^64.  */
  CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 128u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(b, b, 64u), ARBINT_OK);
  CHECK_EQ_I(arbint_set(expected, b), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, expected), 0);

  /*  Consecutive Fibonacci coprimality: gcd(F(100), F(99)) = 1.  */
  CHECK_EQ_I(fibonacci(a, 100u, &ctx), ARBINT_OK);
  CHECK_EQ_I(fibonacci(b, 99u, &ctx), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 1u);

  /*  Fibonacci GCD identity: gcd(F(200), F(100)) = F(100).  */
  CHECK_EQ_I(fibonacci(a, 200u, &ctx), ARBINT_OK);
  CHECK_EQ_I(fibonacci(b, 100u, &ctx), ARBINT_OK);
  CHECK_EQ_I(fibonacci(expected, 100u, &ctx), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, expected), 0);

  /*  Fibonacci GCD identity: gcd(F(300), F(150)) = F(150).  */
  CHECK_EQ_I(fibonacci(a, 300u, &ctx), ARBINT_OK);
  CHECK_EQ_I(fibonacci(b, 150u, &ctx), ARBINT_OK);
  CHECK_EQ_I(fibonacci(expected, 150u, &ctx), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, expected), 0);

  /*  Negative inputs: gcd(-36, 24) = 12.  */
  CHECK_EQ_I(arbint_set_i32(a, -36), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 24), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 12u);

  /*  Large coprime pair: 2^127-1 and 2^64.  */
  CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 127u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(a, a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(b, b, 64u), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 1u);

  /*  gcd(6 * 2^500, 10 * 2^500) = 2 * 2^500.  */
  CHECK_EQ_I(arbint_set_u32(a, 6u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 500u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 10u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(b, b, 500u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(expected, 2u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(expected, expected, 500u), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, expected), 0);

  /*  gcd(3^200, 3^100) = 3^100.  */
  CHECK_EQ_I(arbint_set_u32(base, 3u), ARBINT_OK);
  CHECK_EQ_I(local_pow_u32(a, base, 200u), ARBINT_OK);
  CHECK_EQ_I(local_pow_u32(b, base, 100u), ARBINT_OK);
  CHECK_EQ_I(local_pow_u32(expected, base, 100u), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, expected), 0);

  fprintf(stderr, "passed\n");

  arbint_clear_all(base, expected, g, b, a, (arbint_t *) NULL);
  arbint_ctx_clear(&ctx);
}

/* ========== Section 4: Identity gauntlet ========== */

#define N_OPERANDS 14

static void build_operands(arbint_t * ops, arbint_ctx_t * ctx) {
  arbint_t base;

  CHECK_EQ_I(arbint_init(base, ctx), ARBINT_OK);

  /*  ops[0]  = 0  */
  /*  ops[1]  = 1  */
  CHECK_EQ_I(arbint_set_u32(ops[1], 1u), ARBINT_OK);
  /*  ops[2]  = -1  */
  CHECK_EQ_I(arbint_set_i32(ops[2], -1), ARBINT_OK);
  /*  ops[3]  = 42  */
  CHECK_EQ_I(arbint_set_i32(ops[3], 42), ARBINT_OK);
  /*  ops[4]  = -42  */
  CHECK_EQ_I(arbint_set_i32(ops[4], -42), ARBINT_OK);
  /*  ops[5]  = 0x7FFFFFFF  */
  CHECK_EQ_I(arbint_set_u32(ops[5], 0x7FFFFFFFu), ARBINT_OK);
  /*  ops[6]  = 2^500  (multi-limb, single set bit)  */
  CHECK_EQ_I(arbint_set_u32(ops[6], 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(ops[6], ops[6], 500u), ARBINT_OK);
  /*  ops[7]  = -(2^500)  */
  CHECK_EQ_I(arbint_set(ops[7], ops[6]), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(ops[7], ops[7]), ARBINT_OK);
  /*  ops[8]  = 2^1000 - 1  (all-ones in magnitude, 16 limbs)  */
  CHECK_EQ_I(arbint_set_u32(ops[8], 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(ops[8], ops[8], 1000u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(ops[8], ops[8], 1u), ARBINT_OK);
  /*  ops[9]  = -(2^1000 - 1)  */
  CHECK_EQ_I(arbint_set(ops[9], ops[8]), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(ops[9], ops[9]), ARBINT_OK);
  /*  ops[10] = 2^500 - 1  (all-ones, 8 limbs -- asymmetric vs ops[8])  */
  CHECK_EQ_I(arbint_set_u32(ops[10], 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(ops[10], ops[10], 500u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(ops[10], ops[10], 1u), ARBINT_OK);
  /*  ops[11] = -(2^500 - 1)  */
  CHECK_EQ_I(arbint_set(ops[11], ops[10]), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(ops[11], ops[11]), ARBINT_OK);
  /*  ops[12] = 3^200  (dense multi-limb positive)  */
  CHECK_EQ_I(arbint_set_u32(base, 3u), ARBINT_OK);
  CHECK_EQ_I(local_pow_u32(ops[12], base, 200u), ARBINT_OK);
  /*  ops[13] = -(3^200)  (dense multi-limb negative)  */
  CHECK_EQ_I(arbint_set(ops[13], ops[12]), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(ops[13], ops[13]), ARBINT_OK);

  arbint_clear(base);
}

static void test_identity_gauntlet(void) {
  arbint_ctx_t ctx;
  arbint_t ops[N_OPERANDS];
  arbint_t r1, r2, t1, t2, zero;
  int i, j;
  unsigned pair_count = 0u;

  fprintf(stderr, "  identity gauntlet: ");

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);

  for (i = 0; i < N_OPERANDS; ++i)
    CHECK_EQ_I(arbint_init(ops[i], &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init_all(&ctx, r1, r2, t1, t2, zero,
                              (arbint_t *) NULL),
             ARBINT_OK);

  build_operands(ops, &ctx);

  /*  Test ALL (i, j) pairs -- both orderings for non-commutative checks.  */
  for (i = 0; i < N_OPERANDS && g_failures == 0; ++i) {
    for (j = 0; j < N_OPERANDS && g_failures == 0; ++j) {
      ++pair_count;

      /*  1. a | a = a  */
      CHECK_EQ_I(arbint_or(r1, ops[i], ops[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, ops[i]), 0);

      /*  2. a & a = a  */
      CHECK_EQ_I(arbint_and(r1, ops[i], ops[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, ops[i]), 0);

      /*  3. a ^ a = 0  */
      CHECK_EQ_I(arbint_xor(r1, ops[i], ops[i]), ARBINT_OK);
      CHECK(arbint_is_zero(r1));

      /*  4. a | 0 = a  */
      CHECK_EQ_I(arbint_or(r1, ops[i], zero), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, ops[i]), 0);

      /*  5. a & 0 = 0  */
      CHECK_EQ_I(arbint_and(r1, ops[i], zero), ARBINT_OK);
      CHECK(arbint_is_zero(r1));

      /*  6. not(not(a)) = a  */
      CHECK_EQ_I(arbint_not(r1, ops[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_not(r2, r1), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r2, ops[i]), 0);

      /*  7. a ^ b ^ b = a  */
      CHECK_EQ_I(arbint_xor(r1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_xor(r1, r1, ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, ops[i]), 0);

      /*  8. (a & b) | (a ^ b) = a | b  */
      CHECK_EQ_I(arbint_and(t1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_xor(t2, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(r1, t1, t2), ARBINT_OK);
      CHECK_EQ_I(arbint_or(r2, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      /*  9. not(a & b) = not(a) | not(b)  */
      CHECK_EQ_I(arbint_and(t1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_not(r1, t1), ARBINT_OK);
      CHECK_EQ_I(arbint_not(t1, ops[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_not(t2, ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(r2, t1, t2), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      /* 10. not(a | b) = not(a) & not(b)  */
      CHECK_EQ_I(arbint_or(t1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_not(r1, t1), ARBINT_OK);
      CHECK_EQ_I(arbint_not(t1, ops[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_not(t2, ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_and(r2, t1, t2), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      /* 11. a & (a | b) = a  */
      CHECK_EQ_I(arbint_or(t1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_and(r1, ops[i], t1), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, ops[i]), 0);

      /* 12. a | (a & b) = a  */
      CHECK_EQ_I(arbint_and(t1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(r1, ops[i], t1), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, ops[i]), 0);

      /* 13. (a ^ b) = (a | b) & not(a & b)  */
      CHECK_EQ_I(arbint_xor(r1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(t1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_and(t2, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_not(t2, t2), ARBINT_OK);
      CHECK_EQ_I(arbint_and(r2, t1, t2), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      /* 14. Commutativity: a op b = b op a.  */
      CHECK_EQ_I(arbint_or(r1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(r2, ops[j], ops[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);
      CHECK_EQ_I(arbint_and(r1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_and(r2, ops[j], ops[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);
      CHECK_EQ_I(arbint_xor(r1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_xor(r2, ops[j], ops[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      /* 15. not(a) = -(a + 1)  (algebraic definition)  */
      CHECK_EQ_I(arbint_not(r1, ops[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_add_i32(r2, ops[i], 1), ARBINT_OK);
      CHECK_EQ_I(arbint_neg(r2, r2), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);
    }
  }

  /*  Distributive laws with triples.  */
  {
    static const int triples[][3] = {
      {1, 3, 6},  {2, 4, 7},  {0, 8, 9},  {3, 6, 10},
      {1, 7, 9},  {4, 8, 10}, {2, 6, 8},  {0, 1, 2},
      {5, 11, 13}, {3, 12, 7}, {4, 10, 13}, {1, 9, 12}
    };
    int t;
    for (t = 0;
         t < (int) (sizeof(triples) / sizeof(triples[0])) && g_failures == 0;
         ++t) {
      int ai = triples[t][0];
      int bi = triples[t][1];
      int ci = triples[t][2];

      /*  a & (b | c) = (a & b) | (a & c)  */
      CHECK_EQ_I(arbint_or(t1, ops[bi], ops[ci]), ARBINT_OK);
      CHECK_EQ_I(arbint_and(r1, ops[ai], t1), ARBINT_OK);
      CHECK_EQ_I(arbint_and(t1, ops[ai], ops[bi]), ARBINT_OK);
      CHECK_EQ_I(arbint_and(t2, ops[ai], ops[ci]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(r2, t1, t2), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      /*  a | (b & c) = (a | b) & (a | c)  */
      CHECK_EQ_I(arbint_and(t1, ops[bi], ops[ci]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(r1, ops[ai], t1), ARBINT_OK);
      CHECK_EQ_I(arbint_or(t1, ops[ai], ops[bi]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(t2, ops[ai], ops[ci]), ARBINT_OK);
      CHECK_EQ_I(arbint_and(r2, t1, t2), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      /*  Associativity: (a & b) & c = a & (b & c)  */
      CHECK_EQ_I(arbint_and(t1, ops[ai], ops[bi]), ARBINT_OK);
      CHECK_EQ_I(arbint_and(r1, t1, ops[ci]), ARBINT_OK);
      CHECK_EQ_I(arbint_and(t1, ops[bi], ops[ci]), ARBINT_OK);
      CHECK_EQ_I(arbint_and(r2, ops[ai], t1), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      /*  Associativity: (a | b) | c = a | (b | c)  */
      CHECK_EQ_I(arbint_or(t1, ops[ai], ops[bi]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(r1, t1, ops[ci]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(t1, ops[bi], ops[ci]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(r2, ops[ai], t1), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      /*  Associativity: (a ^ b) ^ c = a ^ (b ^ c)  */
      CHECK_EQ_I(arbint_xor(t1, ops[ai], ops[bi]), ARBINT_OK);
      CHECK_EQ_I(arbint_xor(r1, t1, ops[ci]), ARBINT_OK);
      CHECK_EQ_I(arbint_xor(t1, ops[bi], ops[ci]), ARBINT_OK);
      CHECK_EQ_I(arbint_xor(r2, ops[ai], t1), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);
    }
  }

  fprintf(stderr, "%u pairs, passed\n", pair_count);

  arbint_clear_all(r1, r2, t1, t2, zero, (arbint_t *) NULL);
  for (i = N_OPERANDS - 1; i >= 0; --i)
    arbint_clear(ops[i]);
  arbint_ctx_clear(&ctx);
}

/* ========== Section 5: Popcount / hamming cross-checks ========== */

static void test_popcount_hamming_cross(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, axb, na;
  size_t ham, pc_xor, pc_a, pc_na, nb_a, nb_na, w;
  int i, j;

  fprintf(stderr, "  popcount/hamming cross: ");

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init_all(&ctx, a, b, axb, na, (arbint_t *) NULL),
             ARBINT_OK);

  /*  hamming(a, b) = popcount(a ^ b) for non-negative a, b.  */
  {
    static const uint32_t vals[] = {0u, 1u, 42u, 0xFFu, 0xFFFFu,
                                    0x7FFFFFFFu, 0xFFFFFFFFu};
    int nvals = (int) (sizeof(vals) / sizeof(vals[0]));

    for (i = 0; i < nvals && g_failures == 0; ++i) {
      for (j = i; j < nvals && g_failures == 0; ++j) {
        CHECK_EQ_I(arbint_set_u32(a, vals[i]), ARBINT_OK);
        CHECK_EQ_I(arbint_set_u32(b, vals[j]), ARBINT_OK);
        CHECK_EQ_I(arbint_xor(axb, a, b), ARBINT_OK);
        CHECK_EQ_I(arbint_popcount(axb, &pc_xor), ARBINT_OK);
        CHECK_EQ_I(arbint_hammingdist(a, b, &ham), ARBINT_OK);
        CHECK_EQ_I(ham, pc_xor);
      }
    }
  }

  /*  popcount(x) + popcount(not(x)) = W(x) for various values.  */
  {
    static const int32_t pcvals[] = {1, -1, 42, -42, 127, -128, 0x7FFFFFFF};
    int npcvals = (int) (sizeof(pcvals) / sizeof(pcvals[0]));

    for (i = 0; i < npcvals && g_failures == 0; ++i) {
      CHECK_EQ_I(arbint_set_i32(a, pcvals[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_not(na, a), ARBINT_OK);
      CHECK_EQ_I(arbint_popcount(a, &pc_a), ARBINT_OK);
      CHECK_EQ_I(arbint_popcount(na, &pc_na), ARBINT_OK);
      nb_a = arbint_nbits(a);
      nb_na = arbint_nbits(na);
      w = ((((nb_a > nb_na) ? nb_a : nb_na) + ARBINT_LIMB_BITS - 1u)
            / ARBINT_LIMB_BITS)
          * ARBINT_LIMB_BITS;
      CHECK_EQ_I(pc_a + pc_na, w);
    }
  }

  /*  Multi-limb popcount complement: 3^128.  */
  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  {
    int k;
    for (k = 0; k < 7; ++k)
      CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
  }
  CHECK_EQ_I(arbint_not(na, a), ARBINT_OK);
  CHECK_EQ_I(arbint_popcount(a, &pc_a), ARBINT_OK);
  CHECK_EQ_I(arbint_popcount(na, &pc_na), ARBINT_OK);
  nb_a = arbint_nbits(a);
  nb_na = arbint_nbits(na);
  w = ((((nb_a > nb_na) ? nb_a : nb_na) + ARBINT_LIMB_BITS - 1u)
        / ARBINT_LIMB_BITS)
      * ARBINT_LIMB_BITS;
  CHECK_EQ_I(pc_a + pc_na, w);

  /*  hamming(2^k, 2^k + 1) = 1.  */
  CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 256u), ARBINT_OK);
  CHECK_EQ_I(arbint_add_u32(b, a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_hammingdist(a, b, &ham), ARBINT_OK);
  CHECK_EQ_I(ham, 1u);

  /*  hamming(a, a) = 0.  */
  CHECK_EQ_I(arbint_hammingdist(a, a, &ham), ARBINT_OK);
  CHECK_EQ_I(ham, 0u);

  /*  hamming symmetry for mixed-sign multi-limb.  */
  {
    size_t h1, h2;
    CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(a, a, 500u), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, -42), ARBINT_OK);
    CHECK_EQ_I(arbint_hammingdist(a, b, &h1), ARBINT_OK);
    CHECK_EQ_I(arbint_hammingdist(b, a, &h2), ARBINT_OK);
    CHECK_EQ_I(h1, h2);
    CHECK(h1 > 0u);
  }

  /*  hamming symmetry for both-negative multi-limb.  */
  {
    size_t h1, h2;
    CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(a, a, 500u), ARBINT_OK);
    CHECK_EQ_I(arbint_neg(a, a), ARBINT_OK);
    CHECK_EQ_I(arbint_set_u32(b, 1u), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(b, b, 1000u), ARBINT_OK);
    CHECK_EQ_I(arbint_neg(b, b), ARBINT_OK);
    CHECK_EQ_I(arbint_hammingdist(a, b, &h1), ARBINT_OK);
    CHECK_EQ_I(arbint_hammingdist(b, a, &h2), ARBINT_OK);
    CHECK_EQ_I(h1, h2);
    CHECK(h1 > 0u);
  }

  /*  Hamming invariance under XOR: hamming(a, b) = hamming(a^c, b^c).  */
  {
    size_t h1, h2;
    arbint_t ac, bc, c;
    CHECK_EQ_I(arbint_init_all(&ctx, ac, bc, c, (arbint_t *) NULL),
               ARBINT_OK);
    CHECK_EQ_I(arbint_set_u32(a, 123u), ARBINT_OK);
    CHECK_EQ_I(arbint_set_u32(b, 456u), ARBINT_OK);
    CHECK_EQ_I(arbint_set_u32(c, 789u), ARBINT_OK);
    CHECK_EQ_I(arbint_hammingdist(a, b, &h1), ARBINT_OK);
    CHECK_EQ_I(arbint_xor(ac, a, c), ARBINT_OK);
    CHECK_EQ_I(arbint_xor(bc, b, c), ARBINT_OK);
    CHECK_EQ_I(arbint_hammingdist(ac, bc, &h2), ARBINT_OK);
    CHECK_EQ_I(h1, h2);
    arbint_clear_all(ac, bc, c, (arbint_t *) NULL);
  }

  fprintf(stderr, "passed\n");

  arbint_clear_all(na, axb, b, a, (arbint_t *) NULL);
  arbint_ctx_clear(&ctx);
}

/* ========== Section 6: Bit manipulation stress ========== */

static void test_bit_manipulation_stress(void) {
  arbint_ctx_t ctx;
  arbint_t x, y;
  int bit;
  size_t i;

  fprintf(stderr, "  bit manipulation stress: ");

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init_all(&ctx, x, y, (arbint_t *) NULL), ARBINT_OK);

  /*  Build x by setting individual bits at limb boundaries and beyond.  */
  {
    static const unsigned bits[] = {0u, 1u, 7u, 63u, 64u, 65u, 127u,
                                    128u, 255u, 256u, 500u, 511u, 512u};
    int nbits = (int) (sizeof(bits) / sizeof(bits[0]));

    for (i = 0; i < (size_t) nbits; ++i)
      CHECK_EQ_I(arbint_setbit(x, bits[i]), ARBINT_OK);

    /*  Build the same value via shl + or.  */
    for (i = 0; i < (size_t) nbits; ++i) {
      arbint_t tmp;
      CHECK_EQ_I(arbint_init(tmp, &ctx), ARBINT_OK);
      CHECK_EQ_I(arbint_set_u32(tmp, 1u), ARBINT_OK);
      CHECK_EQ_I(arbint_shl(tmp, tmp, bits[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(y, y, tmp), ARBINT_OK);
      arbint_clear(tmp);
    }

    CHECK_EQ_I(arbint_cmp(x, y), 0);

    /*  Verify each bit is set.  */
    for (i = 0; i < (size_t) nbits; ++i) {
      CHECK_EQ_I(arbint_testbit(x, bits[i], &bit), ARBINT_OK);
      CHECK_EQ_I(bit, 1);
    }

    /*  Verify some non-set bits are 0.  */
    {
      static const unsigned non_bits[] = {2u, 3u, 62u, 66u, 126u,
                                          129u, 254u, 257u, 499u};
      int nnon = (int) (sizeof(non_bits) / sizeof(non_bits[0]));
      for (i = 0; i < (size_t) nnon; ++i) {
        CHECK_EQ_I(arbint_testbit(x, non_bits[i], &bit), ARBINT_OK);
        CHECK_EQ_I(bit, 0);
      }
    }

    /*  Clear bits one by one in reverse and verify.  */
    for (i = (size_t) nbits; i > 0u; --i) {
      CHECK_EQ_I(arbint_clrbit(x, bits[i - 1u]), ARBINT_OK);
      CHECK_EQ_I(arbint_testbit(x, bits[i - 1u], &bit), ARBINT_OK);
      CHECK_EQ_I(bit, 0);
    }
    CHECK(arbint_is_zero(x));
  }

  /*  Negative: set bit 0 of -16.  TC(-16) = ...11110000.  Result = -15.  */
  CHECK_EQ_I(arbint_set_i32(x, -16), ARBINT_OK);
  CHECK_EQ_I(arbint_setbit(x, 0u), ARBINT_OK);
  check_i32_value(x, -15);

  /*  Clear bit 1 of -1.  TC(-1) = ...1111.  Result = -3.  */
  CHECK_EQ_I(arbint_set_i32(x, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_clrbit(x, 1u), ARBINT_OK);
  check_i32_value(x, -3);

  /*  Clear bit 0 of -1.  Result = -2.  */
  CHECK_EQ_I(arbint_set_i32(x, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_clrbit(x, 0u), ARBINT_OK);
  check_i32_value(x, -2);

  /*  Set bit on negative multi-limb: -(2^500).
      TC(-(2^500)) = ...110...0 (bit 500 set, all below = 0).
      Setting bit 0: TC becomes ...110...001, which is -(2^500 - 1).  */
  CHECK_EQ_I(arbint_set_u32(x, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(x, x, 500u), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(x, x), ARBINT_OK);
  CHECK_EQ_I(arbint_setbit(x, 0u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(y, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(y, y, 500u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(y, y, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(y, y), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(x, y), 0);

  /*  Clear bit 500 on multi-limb positive: 2^500 + 2^100.
      Clear bit 500 -> 2^100.  */
  CHECK_EQ_I(arbint_set_u32(x, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(x, x, 500u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(y, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(y, y, 100u), ARBINT_OK);
  CHECK_EQ_I(arbint_or(x, x, y), ARBINT_OK);
  CHECK_EQ_I(arbint_clrbit(x, 500u), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(x, y), 0);

  fprintf(stderr, "passed\n");

  arbint_clear_all(y, x, (arbint_t *) NULL);
  arbint_ctx_clear(&ctx);
}

/* ========== Section 7: CTZ / nbits / sizeinbase consistency ========== */

static void test_ctz_nbits_consistency(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  size_t ctz_val, nb;
  uint32_t i;

  fprintf(stderr, "  ctz/nbits consistency: ");

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  Powers of 2: ctz(2^k) = k, nbits(2^k) = k + 1.  */
  for (i = 0u; i < 20u; ++i) {
    CHECK_EQ_I(arbint_set_u32(x, 1u), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(x, x, i), ARBINT_OK);
    CHECK_EQ_I(arbint_ctz(x, &ctz_val), ARBINT_OK);
    CHECK_EQ_I(ctz_val, (size_t) i);
    nb = arbint_nbits(x);
    CHECK_EQ_I(nb, (size_t) i + 1u);
  }

  /*  Large: ctz(2^500) = 500, nbits = 501.  */
  CHECK_EQ_I(arbint_set_u32(x, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(x, x, 500u), ARBINT_OK);
  CHECK_EQ_I(arbint_ctz(x, &ctz_val), ARBINT_OK);
  CHECK_EQ_I(ctz_val, 500u);
  nb = arbint_nbits(x);
  CHECK_EQ_I(nb, 501u);

  /*  sizeinbase(2^500, 2) = 501.  */
  CHECK_EQ_I(arbint_sizeinbase(x, 2), 501u);

  /*  sizeinbase(2^500, 16) = ceil(501/4) = 126.  */
  CHECK_EQ_I(arbint_sizeinbase(x, 16), 126u);

  /*  Negative: ctz(-2^500) = 500 (operates on magnitude).  */
  CHECK_EQ_I(arbint_neg(x, x), ARBINT_OK);
  CHECK_EQ_I(arbint_ctz(x, &ctz_val), ARBINT_OK);
  CHECK_EQ_I(ctz_val, 500u);

  /*  nbits operates on magnitude too.  */
  nb = arbint_nbits(x);
  CHECK_EQ_I(nb, 501u);

  /*  ctz(0) = EDOM.  */
  arbint_zero(x);
  CHECK_EQ_I(arbint_ctz(x, &ctz_val), ARBINT_EDOM);

  /*  nbits(0) = 0.  */
  CHECK_EQ_I(arbint_nbits(x), 0u);

  /*  sizeinbase(0, 10) = 1 (convention).  */
  CHECK_EQ_I(arbint_sizeinbase(x, 10), 1u);

  fprintf(stderr, "passed\n");

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

/* ========== Section 8: Shift-bitwise cross-checks ========== */

static void test_shift_bitwise_cross(void) {
  arbint_ctx_t ctx;
  arbint_t x, doubled, shifted, mask, masked_ref, masked;

  fprintf(stderr, "  shift-bitwise cross: ");

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init_all(&ctx, x, doubled, shifted, mask, masked_ref,
                              masked, (arbint_t *) NULL),
             ARBINT_OK);

  /*  a << 1 = a + a for various values.  */
  {
    static const int32_t vals[] = {1, 42, -1, -42, 0x7FFFFFFF};
    int nvals = (int) (sizeof(vals) / sizeof(vals[0]));
    int i;
    for (i = 0; i < nvals && g_failures == 0; ++i) {
      CHECK_EQ_I(arbint_set_i32(x, vals[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_shl(shifted, x, 1u), ARBINT_OK);
      CHECK_EQ_I(arbint_add(doubled, x, x), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(shifted, doubled), 0);
    }
  }

  /*  Multi-limb: 3^100 << 1 = 3^100 + 3^100.  */
  {
    arbint_t base;
    CHECK_EQ_I(arbint_init(base, &ctx), ARBINT_OK);
    CHECK_EQ_I(arbint_set_u32(base, 3u), ARBINT_OK);
    CHECK_EQ_I(local_pow_u32(x, base, 100u), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(shifted, x, 1u), ARBINT_OK);
    CHECK_EQ_I(arbint_add(doubled, x, x), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(shifted, doubled), 0);
    arbint_clear(base);
  }

  /*  Mask-via-shift: (a >> k) << k = a & ~(2^k - 1).
      Test with a = 2^1000 - 1, k = 64.  */
  CHECK_EQ_I(arbint_set_u32(x, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(x, x, 1000u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(x, x, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shr(shifted, x, 64u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(shifted, shifted, 64u), ARBINT_OK);
  /*  Build mask = ~(2^64 - 1) = -(2^64).  */
  CHECK_EQ_I(arbint_set_u32(mask, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(mask, mask, 64u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(mask, mask, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_not(mask, mask), ARBINT_OK);
  CHECK_EQ_I(arbint_and(masked, x, mask), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(shifted, masked), 0);

  /*  Same with k = 500.  */
  CHECK_EQ_I(arbint_shr(shifted, x, 500u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(shifted, shifted, 500u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(mask, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(mask, mask, 500u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(mask, mask, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_not(mask, mask), ARBINT_OK);
  CHECK_EQ_I(arbint_and(masked, x, mask), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(shifted, masked), 0);

  fprintf(stderr, "passed\n");

  arbint_clear_all(masked, masked_ref, mask, shifted, doubled, x,
                    (arbint_t *) NULL);
  arbint_ctx_clear(&ctx);
}

/* ========== Section 9: Aliasing stress ========== */

static void test_aliasing_stress(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r, expected;

  fprintf(stderr, "  aliasing stress: ");

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init_all(&ctx, a, b, r, expected, (arbint_t *) NULL),
             ARBINT_OK);

  /*  Build multi-limb signed values.  */
  CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 1000u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(a, a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(a, a), ARBINT_OK);   /*  -(2^1000 - 1)  */

  CHECK_EQ_I(arbint_set_u32(b, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(b, b, 500u), ARBINT_OK);   /*  2^500  */

  /*  OR: rop == a.  */
  CHECK_EQ_I(arbint_or(expected, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_set(r, a), ARBINT_OK);
  CHECK_EQ_I(arbint_or(r, r, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, expected), 0);

  /*  OR: rop == b.  */
  CHECK_EQ_I(arbint_set(r, b), ARBINT_OK);
  CHECK_EQ_I(arbint_or(r, a, r), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, expected), 0);

  /*  AND: rop == a.  */
  CHECK_EQ_I(arbint_and(expected, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_set(r, a), ARBINT_OK);
  CHECK_EQ_I(arbint_and(r, r, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, expected), 0);

  /*  AND: rop == b.  */
  CHECK_EQ_I(arbint_set(r, b), ARBINT_OK);
  CHECK_EQ_I(arbint_and(r, a, r), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, expected), 0);

  /*  XOR: rop == a.  */
  CHECK_EQ_I(arbint_xor(expected, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_set(r, a), ARBINT_OK);
  CHECK_EQ_I(arbint_xor(r, r, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, expected), 0);

  /*  XOR: rop == b.  */
  CHECK_EQ_I(arbint_set(r, b), ARBINT_OK);
  CHECK_EQ_I(arbint_xor(r, a, r), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, expected), 0);

  /*  NOT: rop == a (in-place).  */
  CHECK_EQ_I(arbint_not(expected, a), ARBINT_OK);
  CHECK_EQ_I(arbint_set(r, a), ARBINT_OK);
  CHECK_EQ_I(arbint_not(r, r), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, expected), 0);

  /*  OR with both-negative multi-limb, rop == a.  */
  CHECK_EQ_I(arbint_neg(b, b), ARBINT_OK);   /*  b = -(2^500)  */
  CHECK_EQ_I(arbint_or(expected, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_set(r, a), ARBINT_OK);
  CHECK_EQ_I(arbint_or(r, r, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, expected), 0);

  /*  AND with both-negative multi-limb, rop == b.  */
  CHECK_EQ_I(arbint_and(expected, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_set(r, b), ARBINT_OK);
  CHECK_EQ_I(arbint_and(r, a, r), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, expected), 0);

  /*  XOR with both-negative multi-limb, rop == a.  */
  CHECK_EQ_I(arbint_xor(expected, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_set(r, a), ARBINT_OK);
  CHECK_EQ_I(arbint_xor(r, r, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, expected), 0);

  fprintf(stderr, "passed\n");

  arbint_clear_all(expected, r, b, a, (arbint_t *) NULL);
  arbint_ctx_clear(&ctx);
}

/* ========== Section 10: Immediate variant cross-validation ========== */

static void test_immediate_cross_validation(void) {
  arbint_ctx_t ctx;
  arbint_t a, bval, r1, r2, base;
  int i;

  fprintf(stderr, "  immediate cross-validation: ");

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init_all(&ctx, a, bval, r1, r2, base,
                              (arbint_t *) NULL),
             ARBINT_OK);

  /*  Multi-limb positive: 3^128.  */
  CHECK_EQ_I(arbint_set_u32(base, 3u), ARBINT_OK);
  CHECK_EQ_I(local_pow_u32(a, base, 128u), ARBINT_OK);

  /*  Test u32 variants.  */
  {
    static const uint32_t u32vals[] = {0u, 1u, 0xFFu, 0xFFFFFFFFu};
    int nvals = (int) (sizeof(u32vals) / sizeof(u32vals[0]));
    for (i = 0; i < nvals && g_failures == 0; ++i) {
      CHECK_EQ_I(arbint_set_u32(bval, u32vals[i]), ARBINT_OK);

      CHECK_EQ_I(arbint_and(r1, a, bval), ARBINT_OK);
      CHECK_EQ_I(arbint_and_u32(r2, a, u32vals[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      CHECK_EQ_I(arbint_or(r1, a, bval), ARBINT_OK);
      CHECK_EQ_I(arbint_or_u32(r2, a, u32vals[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      CHECK_EQ_I(arbint_xor(r1, a, bval), ARBINT_OK);
      CHECK_EQ_I(arbint_xor_u32(r2, a, u32vals[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);
    }
  }

  /*  Multi-limb negative: -(3^128).  */
  CHECK_EQ_I(arbint_neg(a, a), ARBINT_OK);

  /*  Test u32 variants with negative a.  */
  {
    static const uint32_t u32vals[] = {0u, 1u, 0xFFu, 0xFFFFFFFFu};
    int nvals = (int) (sizeof(u32vals) / sizeof(u32vals[0]));
    for (i = 0; i < nvals && g_failures == 0; ++i) {
      CHECK_EQ_I(arbint_set_u32(bval, u32vals[i]), ARBINT_OK);

      CHECK_EQ_I(arbint_and(r1, a, bval), ARBINT_OK);
      CHECK_EQ_I(arbint_and_u32(r2, a, u32vals[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      CHECK_EQ_I(arbint_or(r1, a, bval), ARBINT_OK);
      CHECK_EQ_I(arbint_or_u32(r2, a, u32vals[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      CHECK_EQ_I(arbint_xor(r1, a, bval), ARBINT_OK);
      CHECK_EQ_I(arbint_xor_u32(r2, a, u32vals[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);
    }
  }

  /*  Test i32 variants with positive and negative a.  */
  {
    static const int32_t i32vals[] = {0, 1, -1, -2, 42, -42, 0x7FFFFFFF,
                                      -0x7FFFFFFF};
    int nvals = (int) (sizeof(i32vals) / sizeof(i32vals[0]));

    /*  Positive a.  */
    CHECK_EQ_I(arbint_neg(a, a), ARBINT_OK);   /*  back to positive  */
    for (i = 0; i < nvals && g_failures == 0; ++i) {
      CHECK_EQ_I(arbint_set_i32(bval, i32vals[i]), ARBINT_OK);

      CHECK_EQ_I(arbint_and(r1, a, bval), ARBINT_OK);
      CHECK_EQ_I(arbint_and_i32(r2, a, i32vals[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      CHECK_EQ_I(arbint_or(r1, a, bval), ARBINT_OK);
      CHECK_EQ_I(arbint_or_i32(r2, a, i32vals[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      CHECK_EQ_I(arbint_xor(r1, a, bval), ARBINT_OK);
      CHECK_EQ_I(arbint_xor_i32(r2, a, i32vals[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);
    }

    /*  Negative a.  */
    CHECK_EQ_I(arbint_neg(a, a), ARBINT_OK);
    for (i = 0; i < nvals && g_failures == 0; ++i) {
      CHECK_EQ_I(arbint_set_i32(bval, i32vals[i]), ARBINT_OK);

      CHECK_EQ_I(arbint_and(r1, a, bval), ARBINT_OK);
      CHECK_EQ_I(arbint_and_i32(r2, a, i32vals[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      CHECK_EQ_I(arbint_or(r1, a, bval), ARBINT_OK);
      CHECK_EQ_I(arbint_or_i32(r2, a, i32vals[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      CHECK_EQ_I(arbint_xor(r1, a, bval), ARBINT_OK);
      CHECK_EQ_I(arbint_xor_i32(r2, a, i32vals[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);
    }
  }

  fprintf(stderr, "passed\n");

  arbint_clear_all(base, r2, r1, bval, a, (arbint_t *) NULL);
  arbint_ctx_clear(&ctx);
}

/* ========== main ========== */

int main(void) {
  fprintf(stderr, "compute_bitops: comprehensive bitwise stress test\n");

  test_bitwise_adder();
  test_bitwise_subtractor();
  test_binary_gcd();
  test_identity_gauntlet();
  test_popcount_hamming_cross();
  test_bit_manipulation_stress();
  test_ctz_nbits_consistency();
  test_shift_bitwise_cross();
  test_aliasing_stress();
  test_immediate_cross_validation();

  ARBINT_TEST_FINISH("compute_bitops");
}
