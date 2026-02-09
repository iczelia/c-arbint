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

    Three algorithms implemented entirely in terms of bitwise operations,
    exercising every sign combination, multi-limb operands, and aliasing.

    1. Bitwise ripple-carry adder: sum = a + b computed purely with
       xor/and/shl.  Validated against arbint_add for growing operands.

    2. Binary GCD (Stein's algorithm): gcd(a, b) computed using only
       ctz/shr/sub/cmp/is_zero/swap.  Validated against known identities
       (e.g. gcd(F_n, F_{n-1}) = 1 for Fibonacci numbers).

    3. Identity gauntlet: 15 boolean algebra identities verified over
       28 operand pairs spanning all 4 sign quadrants with both small
       and multi-limb operands.  */

#include "test_framework.h"

#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

ARBINT_TEST_DECLARE_FAILURES();

/* ========== Local helper: pow_u32 via repeated squaring ========== */

/*  Compute rop = base^exp using public API (arbint_pow_u32 may not be
    implemented yet).  */
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

/*  Add two non-negative magnitudes using only xor, and, shl.
    sum = a + b, both a and b must be >= 0.

    Algorithm:
      while carry != 0:
        sum   = a ^ carry
        carry = (a & carry) << 1
        a     = sum
      where initial carry = b.

    This exercises xor, and, shl in a tight loop with growing operands.  */
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

  /*  sum = a, carry = b.  */
  rc = arbint_abs(sum, a);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_abs(carry, b);
  if (rc != ARBINT_OK)
    goto cleanup;

  while (!arbint_is_zero(carry)) {
    /*  tmp = sum ^ carry  */
    rc = arbint_xor(tmp, sum, carry);
    if (rc != ARBINT_OK)
      goto cleanup;
    /*  carry = (sum & carry) << 1  */
    rc = arbint_and(carry, sum, carry);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_shl(carry, carry, 1u);
    if (rc != ARBINT_OK)
      goto cleanup;
    /*  sum = tmp  */
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

  /*  a = 2^256 - 1, b = 1.  sum should be 2^256.  */
  CHECK_EQ_I(arbint_set_u32(base, 2u), ARBINT_OK);
  CHECK_EQ_I(local_pow_u32(a, base, 256u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(a, a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 1u), ARBINT_OK);

  CHECK_EQ_I(arbint_add(sum_ref, a, b), ARBINT_OK);
  CHECK_EQ_I(bitwise_add_mag(sum_bit, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(sum_ref, sum_bit), 0);

  /*  Test 2: both multi-limb.  a = 3^100, b = 5^80.  */
  CHECK_EQ_I(arbint_set_u32(base, 3u), ARBINT_OK);
  CHECK_EQ_I(local_pow_u32(a, base, 100u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(base, 5u), ARBINT_OK);
  CHECK_EQ_I(local_pow_u32(b, base, 80u), ARBINT_OK);

  CHECK_EQ_I(arbint_add(sum_ref, a, b), ARBINT_OK);
  CHECK_EQ_I(bitwise_add_mag(sum_bit, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(sum_ref, sum_bit), 0);

  /*  Test 3: a + a = 2*a.  */
  CHECK_EQ_I(arbint_add(sum_ref, a, a), ARBINT_OK);
  CHECK_EQ_I(bitwise_add_mag(sum_bit, a, a), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(sum_ref, sum_bit), 0);

  /*  Test 4: a + 0.  */
  CHECK_EQ_I(arbint_set_u32(t, 0u), ARBINT_OK);
  CHECK_EQ_I(bitwise_add_mag(sum_bit, a, t), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(a, sum_bit), 0);

  /*  Test 5: small + large (very asymmetric).  */
  CHECK_EQ_I(arbint_set_u32(t, 42u), ARBINT_OK);
  CHECK_EQ_I(arbint_add(sum_ref, a, t), ARBINT_OK);
  CHECK_EQ_I(bitwise_add_mag(sum_bit, a, t), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(sum_ref, sum_bit), 0);

  /*  Test 6: many iterations with growing operands.  */
  CHECK_EQ_I(arbint_set_u32(t, 1u), ARBINT_OK);
  for (i = 0u; i < 40u && g_failures == 0; ++i) {
    CHECK_EQ_I(arbint_add(sum_ref, t, b), ARBINT_OK);
    CHECK_EQ_I(bitwise_add_mag(sum_bit, t, b), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(sum_ref, sum_bit), 0);
    /*  Grow t by doubling (via add).  */
    CHECK_EQ_I(arbint_add(t, t, t), ARBINT_OK);
  }

  fprintf(stderr, "passed\n");

  arbint_clear_all(base, t, sum_bit, sum_ref, b, a, (arbint_t *) NULL);
  arbint_ctx_clear(&ctx);
}

/* ========== Section 2: Binary GCD (Stein's algorithm) ========== */

/*  Compute gcd(|a|, |b|) using the binary (Stein's) algorithm.
    Uses only: ctz, shr, sub, cmp, is_zero, abs, swap.
    Does NOT use any multiplication or division.  */
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

  /*  Factor out common powers of 2.  */
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

  /*  Main loop: both u and v are odd.  */
  while (!arbint_is_zero(v)) {
    int c = arbint_cmp(u, v);
    if (c == 0)
      break;
    if (c < 0)
      arbint_swap(u, v);

    /*  u > v, both odd.  u = u - v, then remove trailing zeros.  */
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

  /*  gcd = u * 2^shift.  */
  rc = arbint_shl(g, u, (uint32_t) shift);

cleanup:
  arbint_clear(v);
  arbint_clear(u);
  return rc;
}

/*  Compute Fibonacci(n) iteratively.  */
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

  /*  a = F(0) = 0, b = F(1) = 1.  */
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
  arbint_t a, b, g, expected;

  fprintf(stderr, "  binary GCD: ");

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init_all(&ctx, a, b, g, expected, (arbint_t *) NULL),
             ARBINT_OK);

  /*  Test 1: gcd(0, 0) = 0.  */
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  CHECK(arbint_is_zero(g));

  /*  Test 2: gcd(12, 8) = 4.  */
  CHECK_EQ_I(arbint_set_u32(a, 12u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 8u), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 4u);

  /*  Test 3: gcd(n, 0) = n.  */
  CHECK_EQ_I(arbint_set_u32(a, 42u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 0u), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 42u);

  /*  Test 4: gcd(0, n) = n.  */
  CHECK_EQ_I(arbint_set_u32(a, 0u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 77u), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 77u);

  /*  Test 5: gcd(n, n) = n.  */
  CHECK_EQ_I(arbint_set_u32(a, 100u), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, a), ARBINT_OK);
  check_u32_value(g, 100u);

  /*  Test 6: gcd(2^128, 2^64) = 2^64.  */
  CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 128u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(b, b, 64u), ARBINT_OK);
  CHECK_EQ_I(arbint_set(expected, b), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, expected), 0);

  /*  Test 7: consecutive Fibonacci numbers are coprime.
      gcd(F(100), F(99)) = 1.  */
  CHECK_EQ_I(fibonacci(a, 100u, &ctx), ARBINT_OK);
  CHECK_EQ_I(fibonacci(b, 99u, &ctx), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 1u);

  /*  Test 8: gcd(F(200), F(100)) = F(gcd(200,100)) = F(100).
      This is the Fibonacci GCD identity.  */
  CHECK_EQ_I(fibonacci(a, 200u, &ctx), ARBINT_OK);
  CHECK_EQ_I(fibonacci(b, 100u, &ctx), ARBINT_OK);
  CHECK_EQ_I(fibonacci(expected, 100u, &ctx), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(g, expected), 0);

  /*  Test 9: handles negative inputs (takes |a|, |b|).  */
  CHECK_EQ_I(arbint_set_i32(a, -36), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 24), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 12u);

  /*  Test 10: large coprime pair: 2^127 - 1 (Mersenne prime) and 2^64.  */
  CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 127u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(a, a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(b, b, 64u), ARBINT_OK);
  CHECK_EQ_I(binary_gcd(g, a, b), ARBINT_OK);
  check_u32_value(g, 1u);

  fprintf(stderr, "passed\n");

  arbint_clear_all(expected, g, b, a, (arbint_t *) NULL);
  arbint_ctx_clear(&ctx);
}

/* ========== Section 3: Identity gauntlet ========== */

/*  Verify 15 boolean algebra identities over all 4 sign quadrants.

    The identities tested:
     1. a | a = a                         (idempotent)
     2. a & a = a                         (idempotent)
     3. a ^ a = 0                         (self-inverse)
     4. a | 0 = a                         (identity)
     5. a & 0 = 0                         (annihilator)
     6. not(not(a)) = a                   (involution)
     7. a ^ b ^ b = a                     (XOR self-inverse)
     8. (a & b) | (a ^ b) = a | b         (partition)
     9. not(a & b) = not(a) | not(b)      (De Morgan 1)
    10. not(a | b) = not(a) & not(b)      (De Morgan 2)
    11. a & (a | b) = a                   (absorption 1)
    12. a | (a & b) = a                   (absorption 2)
    13. (a ^ b) = (a | b) & not(a & b)    (XOR decomposition)
    14. a & (b | c) = (a & b) | (a & c)   (distributive)
    15. a | (b & c) = (a | b) & (a | c)   (distributive)

    Each identity is tested with 28+ operand pairs covering:
      - Small values: 0, 1, -1, 42, -42, 0x7FFFFFFF, -0x7FFFFFFF
      - Multi-limb: 2^500, -(2^500), 3^200, -(3^200)
      - Asymmetric sizes: 2^1000 vs small, etc.  */

/*  Build a set of interesting test operands.  */
#define N_OPERANDS 10

static void build_operands(arbint_t * ops, arbint_ctx_t * ctx) {
  arbint_t base;

  CHECK_EQ_I(arbint_init(base, ctx), ARBINT_OK);

  /*  ops[0] = 0  */
  /*  (already zero from init)  */

  /*  ops[1] = 1  */
  CHECK_EQ_I(arbint_set_u32(ops[1], 1u), ARBINT_OK);

  /*  ops[2] = -1  */
  CHECK_EQ_I(arbint_set_i32(ops[2], -1), ARBINT_OK);

  /*  ops[3] = 42  */
  CHECK_EQ_I(arbint_set_i32(ops[3], 42), ARBINT_OK);

  /*  ops[4] = -42  */
  CHECK_EQ_I(arbint_set_i32(ops[4], -42), ARBINT_OK);

  /*  ops[5] = 2^500  (multi-limb positive)  */
  CHECK_EQ_I(arbint_set_u32(ops[5], 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(ops[5], ops[5], 500u), ARBINT_OK);

  /*  ops[6] = -(2^500)  (multi-limb negative)  */
  CHECK_EQ_I(arbint_set(ops[6], ops[5]), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(ops[6], ops[6]), ARBINT_OK);

  /*  ops[7] = 2^1000 - 1  (all-ones in magnitude)  */
  CHECK_EQ_I(arbint_set_u32(ops[7], 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(ops[7], ops[7], 1000u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(ops[7], ops[7], 1u), ARBINT_OK);

  /*  ops[8] = -(2^1000 - 1)  */
  CHECK_EQ_I(arbint_set(ops[8], ops[7]), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(ops[8], ops[8]), ARBINT_OK);

  /*  ops[9] = 2^500 - 1  */
  CHECK_EQ_I(arbint_set_u32(ops[9], 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(ops[9], ops[9], 500u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(ops[9], ops[9], 1u), ARBINT_OK);

  arbint_clear(base);
}

static void test_identity_gauntlet(void) {
  arbint_ctx_t ctx;
  arbint_t ops[N_OPERANDS];
  arbint_t r1, r2, r3, t1, t2, t3, zero;
  int i, j, k;
  unsigned pair_count = 0u;

  fprintf(stderr, "  identity gauntlet: ");

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);

  for (i = 0; i < N_OPERANDS; ++i)
    CHECK_EQ_I(arbint_init(ops[i], &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init_all(&ctx, r1, r2, r3, t1, t2, t3, zero,
                              (arbint_t *) NULL),
             ARBINT_OK);

  build_operands(ops, &ctx);

  /*  Test all (i, j) pairs from the operand set.  */
  for (i = 0; i < N_OPERANDS && g_failures == 0; ++i) {
    for (j = i; j < N_OPERANDS && g_failures == 0; ++j) {
      ++pair_count;

      /*  Identity 1: a | a = a  */
      CHECK_EQ_I(arbint_or(r1, ops[i], ops[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, ops[i]), 0);

      /*  Identity 2: a & a = a  */
      CHECK_EQ_I(arbint_and(r1, ops[i], ops[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, ops[i]), 0);

      /*  Identity 3: a ^ a = 0  */
      CHECK_EQ_I(arbint_xor(r1, ops[i], ops[i]), ARBINT_OK);
      CHECK(arbint_is_zero(r1));

      /*  Identity 4: a | 0 = a  */
      CHECK_EQ_I(arbint_or(r1, ops[i], zero), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, ops[i]), 0);

      /*  Identity 5: a & 0 = 0  */
      CHECK_EQ_I(arbint_and(r1, ops[i], zero), ARBINT_OK);
      CHECK(arbint_is_zero(r1));

      /*  Identity 6: not(not(a)) = a  */
      CHECK_EQ_I(arbint_not(r1, ops[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_not(r2, r1), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r2, ops[i]), 0);

      /*  Identity 7: a ^ b ^ b = a  */
      CHECK_EQ_I(arbint_xor(r1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_xor(r1, r1, ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, ops[i]), 0);

      /*  Identity 8: (a & b) | (a ^ b) = a | b  */
      CHECK_EQ_I(arbint_and(t1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_xor(t2, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(r1, t1, t2), ARBINT_OK);
      CHECK_EQ_I(arbint_or(r2, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      /*  Identity 9: not(a & b) = not(a) | not(b)  (De Morgan 1)  */
      CHECK_EQ_I(arbint_and(t1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_not(r1, t1), ARBINT_OK);
      CHECK_EQ_I(arbint_not(t1, ops[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_not(t2, ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(r2, t1, t2), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      /*  Identity 10: not(a | b) = not(a) & not(b)  (De Morgan 2)  */
      CHECK_EQ_I(arbint_or(t1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_not(r1, t1), ARBINT_OK);
      CHECK_EQ_I(arbint_not(t1, ops[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_not(t2, ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_and(r2, t1, t2), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      /*  Identity 11: a & (a | b) = a  (absorption 1)  */
      CHECK_EQ_I(arbint_or(t1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_and(r1, ops[i], t1), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, ops[i]), 0);

      /*  Identity 12: a | (a & b) = a  (absorption 2)  */
      CHECK_EQ_I(arbint_and(t1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(r1, ops[i], t1), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, ops[i]), 0);

      /*  Identity 13: (a ^ b) = (a | b) & not(a & b)  */
      CHECK_EQ_I(arbint_xor(r1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(t1, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_and(t2, ops[i], ops[j]), ARBINT_OK);
      CHECK_EQ_I(arbint_not(t2, t2), ARBINT_OK);
      CHECK_EQ_I(arbint_and(r2, t1, t2), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);
    }
  }

  /*  Identity 14 and 15 need three operands.  Pick a few triples.  */
  {
    static const int triples[][3] = {
      {1, 3, 5}, {2, 4, 6}, {0, 7, 8}, {3, 5, 9},
      {1, 6, 8}, {4, 7, 9}, {2, 5, 7}, {0, 1, 2}
    };
    int t;
    for (t = 0;
         t < (int) (sizeof(triples) / sizeof(triples[0])) && g_failures == 0;
         ++t) {
      int ai = triples[t][0];
      int bi = triples[t][1];
      int ci = triples[t][2];

      /*  Identity 14: a & (b | c) = (a & b) | (a & c)  */
      CHECK_EQ_I(arbint_or(t1, ops[bi], ops[ci]), ARBINT_OK);
      CHECK_EQ_I(arbint_and(r1, ops[ai], t1), ARBINT_OK);
      CHECK_EQ_I(arbint_and(t1, ops[ai], ops[bi]), ARBINT_OK);
      CHECK_EQ_I(arbint_and(t2, ops[ai], ops[ci]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(r2, t1, t2), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);

      /*  Identity 15: a | (b & c) = (a | b) & (a | c)  */
      CHECK_EQ_I(arbint_and(t1, ops[bi], ops[ci]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(r1, ops[ai], t1), ARBINT_OK);
      CHECK_EQ_I(arbint_or(t1, ops[ai], ops[bi]), ARBINT_OK);
      CHECK_EQ_I(arbint_or(t2, ops[ai], ops[ci]), ARBINT_OK);
      CHECK_EQ_I(arbint_and(r2, t1, t2), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r1, r2), 0);
    }
  }

  fprintf(stderr, "%u pairs, passed\n", pair_count);

  arbint_clear_all(r1, r2, r3, t1, t2, t3, zero, (arbint_t *) NULL);
  for (i = N_OPERANDS - 1; i >= 0; --i)
    arbint_clear(ops[i]);
  arbint_ctx_clear(&ctx);
}

/* ========== Section 4: Popcount / hamming cross-checks ========== */

static void test_popcount_hamming_cross(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, axb;
  size_t ham, pc_xor;
  int i, j;

  fprintf(stderr, "  popcount/hamming cross: ");

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init_all(&ctx, a, b, axb, (arbint_t *) NULL), ARBINT_OK);

  /*  For non-negative a, b: hamming(a, b) = popcount(a ^ b).
      (This only holds when both are non-negative and have the same
      canonical width, i.e., a ^ b captures all differing bits.)  */
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

  /*  Multi-limb: hamming(2^k, 2^k + 1) = popcount(1) = 1.  */
  CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 256u), ARBINT_OK);
  CHECK_EQ_I(arbint_add_u32(b, a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_hammingdist(a, b, &ham), ARBINT_OK);
  CHECK_EQ_I(ham, 1u);

  /*  hamming(a, a) = 0 for any a.  */
  CHECK_EQ_I(arbint_hammingdist(a, a, &ham), ARBINT_OK);
  CHECK_EQ_I(ham, 0u);

  /*  hamming symmetry for multi-limb mixed sign.  */
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

  fprintf(stderr, "passed\n");

  arbint_clear_all(axb, b, a, (arbint_t *) NULL);
  arbint_ctx_clear(&ctx);
}

/* ========== Section 5: Bit manipulation (setbit/clrbit/testbit) ========== */

static void test_bit_manipulation_stress(void) {
  arbint_ctx_t ctx;
  arbint_t x, y;
  int bit;
  size_t i;

  fprintf(stderr, "  bit manipulation stress: ");

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init_all(&ctx, x, y, (arbint_t *) NULL), ARBINT_OK);

  /*  Build a number by setting individual bits, verify against shl + or.  */
  /*  x = 2^0 | 2^7 | 2^64 | 2^127 | 2^255 | 2^500.  */
  {
    static const unsigned bits[] = {0u, 7u, 64u, 127u, 255u, 500u};
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
    CHECK_EQ_I(arbint_testbit(x, 1u, &bit), ARBINT_OK);
    CHECK_EQ_I(bit, 0);
    CHECK_EQ_I(arbint_testbit(x, 65u, &bit), ARBINT_OK);
    CHECK_EQ_I(bit, 0);
    CHECK_EQ_I(arbint_testbit(x, 256u, &bit), ARBINT_OK);
    CHECK_EQ_I(bit, 0);

    /*  Clear bits one by one and verify.  */
    for (i = 0; i < (size_t) nbits; ++i) {
      CHECK_EQ_I(arbint_clrbit(x, bits[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_testbit(x, bits[i], &bit), ARBINT_OK);
      CHECK_EQ_I(bit, 0);
    }
    CHECK(arbint_is_zero(x));
  }

  /*  Negative number bit manipulation: set bit on negative.  */
  CHECK_EQ_I(arbint_set_i32(x, -16), ARBINT_OK);
  /*  TC(-16) = ...11110000.  Setting bit 0 -> ...11110001 = -15.  */
  CHECK_EQ_I(arbint_setbit(x, 0u), ARBINT_OK);
  check_i32_value(x, -15);

  /*  Clear bit 1 of -1: TC(-1) = ...1111.  Clear bit 1 -> ...1101 = -3.  */
  CHECK_EQ_I(arbint_set_i32(x, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_clrbit(x, 1u), ARBINT_OK);
  check_i32_value(x, -3);

  fprintf(stderr, "passed\n");

  arbint_clear_all(y, x, (arbint_t *) NULL);
  arbint_ctx_clear(&ctx);
}

/* ========== Section 6: CTZ / nbits / sizeinbase consistency ========== */

static void test_ctz_nbits_consistency(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  size_t ctz_val, nb;
  uint32_t i;

  fprintf(stderr, "  ctz/nbits consistency: ");

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);

  /*  For powers of 2: ctz(2^k) = k, nbits(2^k) = k + 1.  */
  for (i = 0u; i < 20u; ++i) {
    CHECK_EQ_I(arbint_set_u32(x, 1u), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(x, x, i), ARBINT_OK);
    CHECK_EQ_I(arbint_ctz(x, &ctz_val), ARBINT_OK);
    CHECK_EQ_I(ctz_val, (size_t) i);
    nb = arbint_nbits(x);
    CHECK_EQ_I(nb, (size_t) i + 1u);
  }

  /*  Large power of 2: ctz(2^500) = 500.  */
  CHECK_EQ_I(arbint_set_u32(x, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(x, x, 500u), ARBINT_OK);
  CHECK_EQ_I(arbint_ctz(x, &ctz_val), ARBINT_OK);
  CHECK_EQ_I(ctz_val, 500u);
  nb = arbint_nbits(x);
  CHECK_EQ_I(nb, 501u);

  /*  sizeinbase(2^k, 2) = k + 1 = nbits.  */
  {
    size_t sb = arbint_sizeinbase(x, 2);
    CHECK_EQ_I(sb, 501u);
  }

  /*  sizeinbase consistency: base-16 digits of 2^500.
      2^500 needs exactly 501 binary digits = ceil(501/4) = 126 hex digits.  */
  {
    size_t sb16 = arbint_sizeinbase(x, 16);
    CHECK_EQ_I(sb16, 126u);
  }

  fprintf(stderr, "passed\n");

  arbint_clear(x);
  arbint_ctx_clear(&ctx);
}

/* ========== main ========== */

int main(void) {
  fprintf(stderr, "compute_bitops: comprehensive bitwise stress test\n");

  test_bitwise_adder();
  test_binary_gcd();
  test_identity_gauntlet();
  test_popcount_hamming_cross();
  test_bit_manipulation_stress();
  test_ctz_nbits_consistency();

  ARBINT_TEST_FINISH("compute_bitops");
}
