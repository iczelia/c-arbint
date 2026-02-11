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

#include "arbint_base.h"
#include "arbint_mul.h"
#include "arbint_div.h"
#include "config.h"

/*  Compute (F(n), F(n+1)) simultaneously via fast doubling.

    Uses the identities:
      F(2k)   = F(k) * [2*F(k+1) - F(k)]
      F(2k+1) = F(k)^2 + F(k+1)^2

    Complexity: O(log n) multiplications.

    Algorithm: Start with (a, b) = (F(1), F(2)) = (1, 1). Process bits of n
    from bit position (floor(log2(n)) - 1) down to 0. The MSB is implicitly 1
    and sets the initial state.  */
static arbint_err_t fib_pair(arbint_t fn, arbint_t fn1, uint32_t n,
                             arbint_ctx_t * ctx) {
  arbint_t a, b, c, d, tmp;
  arbint_err_t rc;
  uint32_t mask;

  /*  Base case: F(0) = 0, F(1) = 1.  */
  if (n == 0) {
    rc = arbint_set_i32(fn, 0);
    if (rc != ARBINT_OK)
      return rc;
    return arbint_set_i32(fn1, 1);
  }

  /*  Base case: F(1) = 1, F(2) = 1.  */
  if (n == 1) {
    rc = arbint_set_i32(fn, 1);
    if (rc != ARBINT_OK)
      return rc;
    return arbint_set_i32(fn1, 1);
  }

  /*  Initialize temporaries.  */
  rc = arbint_init(a, ctx);
  if (rc != ARBINT_OK)
    return rc;
  rc = arbint_init(b, ctx);
  if (rc != ARBINT_OK)
    goto cleanup_a;
  rc = arbint_init(c, ctx);
  if (rc != ARBINT_OK)
    goto cleanup_b;
  rc = arbint_init(d, ctx);
  if (rc != ARBINT_OK)
    goto cleanup_c;
  rc = arbint_init(tmp, ctx);
  if (rc != ARBINT_OK)
    goto cleanup_d;

  /*  Start with (a, b) = (F(1), F(2)) = (1, 1).
      The MSB of n is implicitly 1, which sets this initial state.  */
  rc = arbint_set_i32(a, 1);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_set_i32(b, 1);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Find the mask for the second-highest bit.
      For n >= 2, we process bits from (MSB - 1) down to 0.
      The MSB is implicitly 1 and gives us the initial state (F(1), F(2)).
      Use shift-based check to avoid overflow when mask approaches 2^31.  */
  mask = 1u;
  while (mask <= n >> 1u)
    mask <<= 1u;
  /*  Now mask is at the MSB. Move to second-highest bit.  */
  mask >>= 1u;

  /*  Process bits from (MSB - 1) down to LSB.  */
  while (mask != 0) {
    /*  Double step: compute (F(2k), F(2k+1)) from (F(k), F(k+1)).
        c = F(2k)   = a * (2*b - a)
        d = F(2k+1) = a^2 + b^2  */
    rc = arbint_mul_i32(c, b, 2);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_sub(c, c, a);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_mul(c, a, c);
    if (rc != ARBINT_OK)
      goto cleanup;

    rc = arbint_sqr(d, a);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_sqr(tmp, b);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_add(d, d, tmp);
    if (rc != ARBINT_OK)
      goto cleanup;

    if (n & mask) {
      /*  Bit is 1: (a, b) = (F(2k+1), F(2k+2)) = (d, c+d).  */
      rc = arbint_set(a, d);
      if (rc != ARBINT_OK)
        goto cleanup;
      rc = arbint_add(b, c, d);
      if (rc != ARBINT_OK)
        goto cleanup;
    } else {
      /*  Bit is 0: (a, b) = (F(2k), F(2k+1)) = (c, d).  */
      rc = arbint_set(a, c);
      if (rc != ARBINT_OK)
        goto cleanup;
      rc = arbint_set(b, d);
      if (rc != ARBINT_OK)
        goto cleanup;
    }

    mask >>= 1u;
  }

  /*  Copy results to output.  */
  rc = arbint_set(fn, a);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_set(fn1, b);

cleanup:
  arbint_clear(tmp);
cleanup_d:
  arbint_clear(d);
cleanup_c:
  arbint_clear(c);
cleanup_b:
  arbint_clear(b);
cleanup_a:
  arbint_clear(a);
  return rc;
}

/*  Compute rop = n! (factorial).

    Uses simple iterative multiplication: rop = 1 * 2 * 3 * ... * n.
    0! = 1! = 1.  */
arbint_err_t arbint_fac_u32(arbint_t rop, uint32_t n) {
  arbint_err_t rc;
  uint32_t i;

  if (rop == NULL)
    return ARBINT_EINVAL;

  rc = arbint_set_i32(rop, 1);
  if (rc != ARBINT_OK)
    return rc;

  for (i = 2; i <= n; ++i) {
    rc = arbint_mul_u32(rop, rop, i);
    if (rc != ARBINT_OK)
      return rc;
  }

  return ARBINT_OK;
}

/*  Compute rop = C(n, k) (binomial coefficient).

    Uses iterative multiplication with exact division:
      C(n, k) = (n * (n-1) * ... * (n-k+1)) / (k * (k-1) * ... * 1)

    Each division is exact because the partial product is always
    divisible by the corresponding denominator.

    Exploits symmetry: C(n, k) = C(n, n-k), so we use k = min(k, n-k).  */
arbint_err_t arbint_bin_u32u32(arbint_t rop, uint32_t n, uint32_t k) {
  arbint_err_t rc;
  uint32_t i;

  if (rop == NULL)
    return ARBINT_EINVAL;

  /*  Edge case: k > n means C(n, k) = 0.  */
  if (k > n) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  /*  Edge case: C(n, 0) = C(n, n) = 1.  */
  if (k == 0 || k == n)
    return arbint_set_i32(rop, 1);

  /*  Exploit symmetry: C(n, k) = C(n, n-k).  */
  if (k > n - k)
    k = n - k;

  /*  Compute iteratively: rop = n/1 * (n-1)/2 * (n-2)/3 * ... * (n-k+1)/k.  */
  rc = arbint_set_u32(rop, n);
  if (rc != ARBINT_OK)
    return rc;

  for (i = 1; i < k; ++i) {
    rc = arbint_mul_u32(rop, rop, n - i);
    if (rc != ARBINT_OK)
      return rc;
    rc = arbint_tdiv_q_u32(rop, rop, i + 1);
    if (rc != ARBINT_OK)
      return rc;
  }

  return ARBINT_OK;
}

/*  Compute rop = F(n) (n-th Fibonacci number).

    F(0) = 0, F(1) = 1, F(n) = F(n-1) + F(n-2).
    Uses fast doubling for O(log n) multiplications.  */
arbint_err_t arbint_fib_u32(arbint_t rop, uint32_t n) {
  arbint_ctx_t * ctx;
  arbint_t fn1;
  arbint_err_t rc;

  if (rop == NULL)
    return ARBINT_EINVAL;

  ctx = rop[0]._ctx;

  rc = arbint_init(fn1, ctx);
  if (rc != ARBINT_OK)
    return rc;

  rc = fib_pair(rop, fn1, n, ctx);

  arbint_clear(fn1);
  return rc;
}

/*  Compute rop = L(n) (n-th Lucas number).

    L(0) = 2, L(1) = 1, L(n) = L(n-1) + L(n-2).
    Uses the identity: L(n) = 2*F(n+1) - F(n).  */
arbint_err_t arbint_lucas_u32(arbint_t rop, uint32_t n) {
  arbint_ctx_t * ctx;
  arbint_t fn, fn1;
  arbint_err_t rc;

  if (rop == NULL)
    return ARBINT_EINVAL;

  ctx = rop[0]._ctx;

  rc = arbint_init(fn, ctx);
  if (rc != ARBINT_OK)
    return rc;

  rc = arbint_init(fn1, ctx);
  if (rc != ARBINT_OK) {
    arbint_clear(fn);
    return rc;
  }

  /*  Compute F(n) and F(n+1).  */
  rc = fib_pair(fn, fn1, n, ctx);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  L(n) = 2*F(n+1) - F(n).  */
  rc = arbint_mul_i32(rop, fn1, 2);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_sub(rop, rop, fn);

cleanup:
  arbint_clear(fn1);
  arbint_clear(fn);
  return rc;
}

/*  Test if a is a perfect square.

    Sets *out = 1 if a = k^2 for some integer k, 0 otherwise.
    Returns ARBINT_EDOM if a < 0.  */
arbint_err_t arbint_is_square(const arbint_t a, int * out) {
  arbint_ctx_t * ctx;
  arbint_t root, sq;
  arbint_err_t rc;

  if (a == NULL || out == NULL)
    return ARBINT_EINVAL;

  /*  Negative numbers are not perfect squares.  */
  if (a[0]._sz < 0)
    return ARBINT_EDOM;

  /*  Zero is a perfect square (0 = 0^2).  */
  if (a[0]._sz == 0) {
    *out = 1;
    return ARBINT_OK;
  }

  ctx = a[0]._ctx;

  rc = arbint_init(root, ctx);
  if (rc != ARBINT_OK)
    return rc;

  rc = arbint_init(sq, ctx);
  if (rc != ARBINT_OK) {
    arbint_clear(root);
    return rc;
  }

  /*  Compute root = floor(sqrt(a)).  */
  rc = arbint_isqrt(root, a);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Check if root^2 == a.  */
  rc = arbint_sqr(sq, root);
  if (rc != ARBINT_OK)
    goto cleanup;

  *out = arbint_eq(sq, a) ? 1 : 0;
  rc = ARBINT_OK;

cleanup:
  arbint_clear(sq);
  arbint_clear(root);
  return rc;
}
