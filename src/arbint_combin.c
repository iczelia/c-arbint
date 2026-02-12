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
#include "arbint_cache.h"
#include "arbint_div.h"
#include "arbint_mul.h"
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

/*  Small-value lookup tables to avoid setup overhead for tiny n.  */
static const uint32_t arbint_fib_small_u32[] = {0u, 1u,  1u,  2u,  3u, 5u,
                                                8u, 13u, 21u, 34u, 55u};
static const uint32_t arbint_lucas_small_u32[] = {2u,  1u,  3u,  4u,  7u,  11u,
                                                  18u, 29u, 47u, 76u, 123u};
#define ARBINT_FIB_SMALL_U32_COUNT                                            \
  (sizeof(arbint_fib_small_u32) / sizeof(arbint_fib_small_u32[0]))
#define ARBINT_LUCAS_SMALL_U32_COUNT                                          \
  (sizeof(arbint_lucas_small_u32) / sizeof(arbint_lucas_small_u32[0]))

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

  if ((size_t) n < ARBINT_FIB_SMALL_U32_COUNT)
    return arbint_set_u32(rop, arbint_fib_small_u32[n]);

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

  if ((size_t) n < ARBINT_LUCAS_SMALL_U32_COUNT)
    return arbint_set_u32(rop, arbint_lucas_small_u32[n]);

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

/*  Compute Euler's totient phi(n).

    Domain: n > 0. Returns ARBINT_EDOM for n <= 0.
    Uses a fast uint32 path and a generic big-int factorization fallback:
      phi(n) = n * product_{p|n} (1 - 1/p).  */
arbint_err_t arbint_totient(arbint_t rop, const arbint_t n) {
  arbint_ctx_t * ctx;
  arbint_err_t rc;
  uint32_t n_u32;
  arbint_t m, phi, p, p2, q, r, tmp;

  if (rop == NULL || n == NULL)
    return ARBINT_EINVAL;
  if (n[0]._sz <= 0)
    return ARBINT_EDOM;

  /*  phi(1) = 1.  */
  if (arbint_cmp_u32(n, 1u) == 0)
    return arbint_set_u32(rop, 1u);

  /*  Fast path for small inputs using pure u32 arithmetic.  */
  if (arbint_get_u32(n, &n_u32) == ARBINT_OK) {
    uint32_t m = n_u32;
    uint32_t phi = n_u32;
    uint32_t p;

    if ((m & 1u) == 0u) {
      phi -= phi / 2u;
      do {
        m >>= 1u;
      } while ((m & 1u) == 0u);
    }

    for (p = 3u; p <= m / p; p += 2u) {
      if (m % p != 0u)
        continue;
      phi -= phi / p;
      do {
        m /= p;
      } while (m % p == 0u);
    }

    if (m > 1u)
      phi -= phi / m;

    return arbint_set_u32(rop, phi);
  }

  /*  Generic big-int fallback.  */
  ctx = rop[0]._ctx;
  if (ctx == NULL)
    ctx = n[0]._ctx;

  rc = arbint_init_all(ctx, m, phi, p, p2, q, r, tmp, (arbint_t *) NULL);
  if (rc != ARBINT_OK)
    return rc;

  rc = arbint_set(m, n);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_set(phi, n);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Factor 2 quickly.  */
  if ((ARBINT_CLIMBS(m)[0] & 1u) == 0u) {
    rc = arbint_tdiv_q_u32(tmp, phi, 2u);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_sub(phi, phi, tmp);
    if (rc != ARBINT_OK)
      goto cleanup;

    do {
      rc = arbint_tdiv_q_u32(m, m, 2u);
      if (rc != ARBINT_OK)
        goto cleanup;
    } while (!arbint_is_zero(m) && (ARBINT_CLIMBS(m)[0] & 1u) == 0u);
  }

  rc = arbint_set_u32(p, 3u);
  if (rc != ARBINT_OK)
    goto cleanup;

  while (!arbint_is_zero(m)) {
    rc = arbint_sqr(p2, p);
    if (rc != ARBINT_OK)
      goto cleanup;
    if (arbint_cmp(p2, m) > 0)
      break;

    rc = arbint_tdiv_qr(q, r, m, p);
    if (rc != ARBINT_OK)
      goto cleanup;

    if (!arbint_is_zero(r)) {
      rc = arbint_nextprime(p, p);
      if (rc != ARBINT_OK)
        goto cleanup;
      continue;
    }

    /*  Distinct prime factor contribution: phi -= phi / p.  */
    rc = arbint_tdiv_q(tmp, phi, p);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_sub(phi, phi, tmp);
    if (rc != ARBINT_OK)
      goto cleanup;

    rc = arbint_set(m, q);
    if (rc != ARBINT_OK)
      goto cleanup;

    /*  Remove full p-adic power.  */
    for (;;) {
      rc = arbint_tdiv_qr(q, r, m, p);
      if (rc != ARBINT_OK)
        goto cleanup;
      if (!arbint_is_zero(r))
        break;
      rc = arbint_set(m, q);
      if (rc != ARBINT_OK)
        goto cleanup;
    }
  }

  if (arbint_cmp_u32(m, 1u) > 0) {
    rc = arbint_tdiv_q(tmp, phi, m);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_sub(phi, phi, tmp);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

  rc = arbint_set(rop, phi);

cleanup:
  arbint_clear_all(m, phi, p, p2, q, r, tmp, (arbint_t *) NULL);

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

/*  Static small-prime table shared by is_power and isprime.
    Primes up to 521 (97 entries).  */
static const uint32_t arbint_is_power_primes[] = {
    2,   3,   5,   7,   11,  13,  17,  19,  23,  29,  31,  37,  41,  43,
    47,  53,  59,  61,  67,  71,  73,  79,  83,  89,  97,  101, 103, 107,
    109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
    191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263,
    269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
    353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433,
    439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521};
#define ARBINT_IS_POWER_NPRIMES                                               \
  (sizeof(arbint_is_power_primes) / sizeof(arbint_is_power_primes[0]))

#define ARBINT_ISPRIME_DEFAULT_REPS 16u

/*  x = base^exp mod mod, where base is a u32 and exp/mod are arbint.  */
static arbint_err_t arbint_isprime_pow_u32_tmod(arbint_t x, uint32_t base,
                                                const arbint_t exp,
                                                const arbint_t mod,
                                                arbint_ctx_t * ctx) {
  arbint_t acc, pow_base, e, tmp;
  arbint_err_t rc;

  rc = arbint_init(acc, ctx);
  if (rc != ARBINT_OK)
    return rc;

  rc = arbint_init(pow_base, ctx);
  if (rc != ARBINT_OK) {
    arbint_clear(acc);
    return rc;
  }

  rc = arbint_init(e, ctx);
  if (rc != ARBINT_OK) {
    arbint_clear(pow_base);
    arbint_clear(acc);
    return rc;
  }

  rc = arbint_init(tmp, ctx);
  if (rc != ARBINT_OK) {
    arbint_clear(e);
    arbint_clear(pow_base);
    arbint_clear(acc);
    return rc;
  }

  rc = arbint_set_i32(acc, 1);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_set_u32(pow_base, base);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_tdiv_r(pow_base, pow_base, mod);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_set(e, exp);
  if (rc != ARBINT_OK)
    goto cleanup;

  while (!arbint_is_zero(e)) {
    if ((ARBINT_CLIMBS(e)[0] & 1u) != 0u) {
      rc = arbint_mul(tmp, acc, pow_base);
      if (rc != ARBINT_OK)
        goto cleanup;
      rc = arbint_tdiv_r(acc, tmp, mod);
      if (rc != ARBINT_OK)
        goto cleanup;
    }

    rc = arbint_shr(e, e, 1u);
    if (rc != ARBINT_OK)
      goto cleanup;
    if (arbint_is_zero(e))
      break;

    rc = arbint_sqr(tmp, pow_base);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_tdiv_r(pow_base, tmp, mod);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

  rc = arbint_set(x, acc);

cleanup:
  arbint_clear(tmp);
  arbint_clear(e);
  arbint_clear(pow_base);
  arbint_clear(acc);
  return rc;
}

/*  One strong Miller-Rabin round. out_composite=1 means definitely
    composite.  */
static arbint_err_t arbint_isprime_mr_round(const arbint_t n, const arbint_t d,
                                            size_t s, const arbint_t n_minus_1,
                                            uint32_t base, int * out_composite,
                                            arbint_ctx_t * ctx) {
  arbint_t x, tmp;
  arbint_err_t rc;
  size_t i;

  rc = arbint_init(x, ctx);
  if (rc != ARBINT_OK)
    return rc;

  rc = arbint_init(tmp, ctx);
  if (rc != ARBINT_OK) {
    arbint_clear(x);
    return rc;
  }

  rc = arbint_isprime_pow_u32_tmod(x, base, d, n, ctx);
  if (rc != ARBINT_OK)
    goto cleanup;

  if (arbint_cmp_u32(x, 1u) == 0 || arbint_eq(x, n_minus_1)) {
    *out_composite = 0;
    rc = ARBINT_OK;
    goto cleanup;
  }

  for (i = 1u; i < s; ++i) {
    rc = arbint_sqr(tmp, x);
    if (rc != ARBINT_OK)
      goto cleanup;

    rc = arbint_tdiv_r(x, tmp, n);
    if (rc != ARBINT_OK)
      goto cleanup;

    if (arbint_eq(x, n_minus_1)) {
      *out_composite = 0;
      rc = ARBINT_OK;
      goto cleanup;
    }

    if (arbint_cmp_u32(x, 1u) == 0) {
      *out_composite = 1;
      rc = ARBINT_OK;
      goto cleanup;
    }
  }

  *out_composite = 1;
  rc = ARBINT_OK;

cleanup:
  arbint_clear(tmp);
  arbint_clear(x);
  return rc;
}

/*  Probable-prime test using trial division + strong Miller-Rabin.  */
arbint_err_t arbint_isprime(const arbint_t n, int reps, int * out) {
  arbint_ctx_t * ctx;
  arbint_t n_minus_1, d;
  const arbint_limb_t * np;
  size_t nn;
  size_t i;
  size_t s = 0u;
  size_t rounds;
  int init_n_minus_1 = 0;
  int init_d = 0;
  arbint_err_t rc;

  if (n == NULL || out == NULL)
    return ARBINT_EINVAL;

  /*  Primes are positive integers >= 2.  */
  if (n[0]._sz <= 0) {
    *out = 0;
    return ARBINT_OK;
  }

  np = ARBINT_CLIMBS(n);
  nn = arbint_abs_sz(n[0]._sz);

  if (nn == 1u && np[0] < 2u) {
    *out = 0;
    return ARBINT_OK;
  }

  /*  Handle even numbers quickly.  */
  if ((np[0] & 1u) == 0u) {
    *out = (nn == 1u && np[0] == 2u) ? 1 : 0;
    return ARBINT_OK;
  }

  /*  Trial divide by a fixed small-prime table.  */
  for (i = 0u; i < ARBINT_IS_POWER_NPRIMES; ++i) {
    uint32_t p = arbint_is_power_primes[i];
    arbint_limb_t rem = 0u;

    if (nn == 1u && np[0] == (arbint_limb_t) p) {
      *out = 1;
      return ARBINT_OK;
    }

    rc = arbint_mod_mag_single_limb(np, nn, (arbint_limb_t) p, &rem);
    if (rc != ARBINT_OK)
      return rc;
    if (rem == 0u) {
      *out = 0;
      return ARBINT_OK;
    }
  }

  /*  If n <= 521^2 and survived all primes <= 521, it is prime.  */
  if (nn == 1u && np[0] <= (arbint_limb_t) (521u * 521u)) {
    *out = 1;
    return ARBINT_OK;
  }

  ctx = n[0]._ctx;
  rc = arbint_init(n_minus_1, ctx);
  if (rc != ARBINT_OK)
    return rc;
  init_n_minus_1 = 1;

  rc = arbint_init(d, ctx);
  if (rc != ARBINT_OK)
    goto cleanup;
  init_d = 1;

  rc = arbint_sub_i32(n_minus_1, n, 1);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_set(d, n_minus_1);
  if (rc != ARBINT_OK)
    goto cleanup;

  while ((ARBINT_CLIMBS(d)[0] & 1u) == 0u) {
    rc = arbint_shr(d, d, 1u);
    if (rc != ARBINT_OK)
      goto cleanup;
    ++s;
  }

  rounds = (reps > 0) ? (size_t) reps : ARBINT_ISPRIME_DEFAULT_REPS;
  if (rounds > ARBINT_IS_POWER_NPRIMES)
    rounds = ARBINT_IS_POWER_NPRIMES;
  if (rounds == 0u)
    rounds = 1u;

  for (i = 0u; i < rounds; ++i) {
    int composite = 0;
    rc = arbint_isprime_mr_round(n, d, s, n_minus_1, arbint_is_power_primes[i],
                                 &composite, ctx);
    if (rc != ARBINT_OK)
      goto cleanup;
    if (composite) {
      *out = 0;
      rc = ARBINT_OK;
      goto cleanup;
    }
  }

  *out = 1;
  rc = ARBINT_OK;

cleanup:
  if (init_d)
    arbint_clear(d);
  if (init_n_minus_1)
    arbint_clear(n_minus_1);
  return rc;
}

/*  Find the smallest prime strictly greater than n.  */
arbint_err_t arbint_nextprime(arbint_t rop, const arbint_t n) {
  arbint_ctx_t * ctx;
  arbint_t cand;
  arbint_err_t rc;
  int is_prime = 0;

  if (rop == NULL || n == NULL)
    return ARBINT_EINVAL;

  ctx = rop[0]._ctx;
  if (ctx == NULL)
    ctx = n[0]._ctx;

  rc = arbint_init(cand, ctx);
  if (rc != ARBINT_OK)
    return rc;

  rc = arbint_add_i32(cand, n, 1);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  For n < 2, next prime is 2.  */
  if (cand[0]._sz <= 0 || arbint_cmp_u32(cand, 2u) <= 0) {
    rc = arbint_set_u32(rop, 2u);
    goto cleanup;
  }

  /*  Keep candidate odd (all odd primes > 2).  */
  if ((ARBINT_CLIMBS(cand)[0] & 1u) == 0u) {
    rc = arbint_add_i32(cand, cand, 1);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

  for (;;) {
    rc = arbint_isprime(cand, 0, &is_prime);
    if (rc != ARBINT_OK)
      goto cleanup;

    if (is_prime) {
      rc = arbint_set(rop, cand);
      goto cleanup;
    }

    rc = arbint_add_i32(cand, cand, 2);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

cleanup:
  arbint_clear(cand);
  return rc;
}

/*  Find the largest prime strictly smaller than n.
    Returns ARBINT_EDOM when no such prime exists (n <= 2).  */
arbint_err_t arbint_prevprime(arbint_t rop, const arbint_t n) {
  arbint_ctx_t * ctx;
  arbint_t cand;
  arbint_err_t rc;
  int is_prime = 0;

  if (rop == NULL || n == NULL)
    return ARBINT_EINVAL;

  /*  No prime exists below 2.  */
  if (n[0]._sz <= 0 || arbint_cmp_u32(n, 2u) <= 0)
    return ARBINT_EDOM;

  ctx = rop[0]._ctx;
  if (ctx == NULL)
    ctx = n[0]._ctx;

  rc = arbint_init(cand, ctx);
  if (rc != ARBINT_OK)
    return rc;

  rc = arbint_sub_i32(cand, n, 1);
  if (rc != ARBINT_OK)
    goto cleanup;

  if (arbint_cmp_u32(cand, 2u) == 0) {
    rc = arbint_set_u32(rop, 2u);
    goto cleanup;
  }

  if (cand[0]._sz <= 0 || arbint_cmp_u32(cand, 2u) < 0) {
    rc = ARBINT_EDOM;
    goto cleanup;
  }

  /*  Keep candidate odd (all odd primes > 2).  */
  if ((ARBINT_CLIMBS(cand)[0] & 1u) == 0u) {
    rc = arbint_sub_i32(cand, cand, 1);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

  for (;;) {
    if (arbint_cmp_u32(cand, 2u) < 0) {
      rc = ARBINT_EDOM;
      goto cleanup;
    }

    rc = arbint_isprime(cand, 0, &is_prime);
    if (rc != ARBINT_OK)
      goto cleanup;

    if (is_prime) {
      rc = arbint_set(rop, cand);
      goto cleanup;
    }

    rc = arbint_sub_i32(cand, cand, 2);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

cleanup:
  arbint_clear(cand);
  return rc;
}

/*  primorial(n) = product of all primes p <= n.

    Domain: n >= 0. Returns ARBINT_EDOM for n < 0.
    For n <= 1, result is 1.  */
arbint_err_t arbint_primorial(arbint_t rop, const arbint_t n) {
  arbint_err_t rc;
  uint32_t n_u32;
  arbint_ctx_t * ctx;
  arbint_t bound, acc, p, tmp;

  if (rop == NULL || n == NULL)
    return ARBINT_EINVAL;
  if (n[0]._sz < 0)
    return ARBINT_EDOM;

  if (n[0]._sz == 0 || arbint_cmp_u32(n, 1u) <= 0)
    return arbint_set_u32(rop, 1u);

  /*  Fast path for u32 bounds via sieve cache.  */
  if (arbint_get_u32(n, &n_u32) == ARBINT_OK) {
    uint32_t k;

    rc = arbint_set_u32(rop, 1u);
    if (rc != ARBINT_OK)
      return rc;

    if (n_u32 >= 2u) {
      rc = arbint_cache_prime_sieve_ensure(n_u32);
      if (rc != ARBINT_OK)
        return rc;
    }

    for (k = 2u; k <= n_u32; ++k) {
      if (arbint_cache_prime_sieve_is_prime(k)) {
        rc = arbint_mul_u32(rop, rop, k);
        if (rc != ARBINT_OK)
          return rc;
      }
    }
    return ARBINT_OK;
  }

  /*  Generic big-int bound path.  */
  ctx = rop[0]._ctx;
  if (ctx == NULL)
    ctx = n[0]._ctx;

  rc = arbint_init_all(ctx, bound, acc, p, tmp, (arbint_t *) NULL);
  if (rc != ARBINT_OK)
    return rc;

  rc = arbint_set(bound, n);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_set_u32(acc, 1u);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_set_u32(p, 2u);
  if (rc != ARBINT_OK)
    goto cleanup;

  while (arbint_cmp(p, bound) <= 0) {
    rc = arbint_mul(tmp, acc, p);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_set(acc, tmp);
    if (rc != ARBINT_OK)
      goto cleanup;

    rc = arbint_nextprime(p, p);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

  rc = arbint_set(rop, acc);

cleanup:
  arbint_clear_all(bound, acc, p, tmp, (arbint_t *) NULL);

  return rc;
}

/*  Test if a is a perfect power (a = b^k for some integers b, k >= 2).

    Sets *out = 1 if such b, k exist, 0 otherwise.

    GMP-compatible semantics:
    - 0: out=1 (0 = 0^k for any k >= 2)
    - 1: out=1 (1 = 1^k for any k >= 2)
    - -1: out=1 (-1 = (-1)^3 = (-1)^5 = ...)
    - Negative odd powers: out=1 (e.g., -8 = (-2)^3)
    - Negative even powers: out=0 (no real even root of negative)

    Algorithm: Only check prime exponents k (composite k = p*q means
    a = (b^q)^p is also a p-th power). Check all primes up to log2(|a|).

    Note: arbint_root and arbint_pow_u32 use uint32_t exponents.
    For inputs with |a| > 2^UINT32_MAX (requiring exponent checks beyond
    UINT32_MAX), returns ARBINT_EOVERFLOW. Such inputs would have ~500M
    limbs on 64-bit systems -- far beyond practical use.  */
arbint_err_t arbint_is_power(const arbint_t a, int * out) {
  arbint_ctx_t * ctx;
  arbint_t abs_a, r, r_pow;
  arbint_err_t rc;
  size_t max_k;
  size_t nbits;
  size_t i;
  int a_neg;
  int is_sq;
  int init_abs = 0;
  int init_r = 0;
  int init_rpow = 0;

  if (a == NULL || out == NULL)
    return ARBINT_EINVAL;

  /*  Handle 0: 0 = 0^k for any k >= 2.  */
  if (a[0]._sz == 0) {
    *out = 1;
    return ARBINT_OK;
  }

  a_neg = (a[0]._sz < 0);

  /*  Handle 1 and -1.  */
  {
    size_t an = arbint_abs_sz(a[0]._sz);
    if (an == 1u && ARBINT_CLIMBS(a)[0] == 1u) {
      /*  1 = 1^k, -1 = (-1)^3.  */
      *out = 1;
      return ARBINT_OK;
    }
  }

  ctx = a[0]._ctx;

  /*  Allocate temporaries.  */
  rc = arbint_init(abs_a, ctx);
  if (rc != ARBINT_OK)
    return rc;
  init_abs = 1;

  rc = arbint_init(r, ctx);
  if (rc != ARBINT_OK)
    goto cleanup;
  init_r = 1;

  rc = arbint_init(r_pow, ctx);
  if (rc != ARBINT_OK)
    goto cleanup;
  init_rpow = 1;

  /*  Get |a|.  */
  rc = arbint_abs(abs_a, a);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Compute max exponent to check: floor(log2(|a|)).  */
  nbits = arbint_nbits(abs_a);
  max_k = nbits - 1u;

  /*  Check if max_k exceeds API limit.  */
  if (max_k > (size_t) UINT32_MAX) {
    rc = ARBINT_EOVERFLOW;
    goto cleanup;
  }

  /*  Fast path for positive: check if perfect square.  */
  if (!a_neg) {
    rc = arbint_is_square(a, &is_sq);
    if (rc != ARBINT_OK)
      goto cleanup;
    if (is_sq) {
      *out = 1;
      rc = ARBINT_OK;
      goto cleanup;
    }
  }

  /*  Check prime exponents from static table.
      For positive a: start at index 1 (k=3) since squares already checked.
      For negative a: start at index 1 (k=3) since even exponents skipped.  */
  for (i = 1u; i < ARBINT_IS_POWER_NPRIMES; ++i) {
    uint32_t k = arbint_is_power_primes[i];
    if ((size_t) k > max_k)
      break;

    /*  Skip even exponents for negative a.  */
    if (a_neg && (k % 2u == 0u))
      continue;

    /*  Compute r = floor(|a|^(1/k)).  */
    rc = arbint_root(r, abs_a, k);
    if (rc != ARBINT_OK)
      goto cleanup;

    /*  Skip if r < 2 (only 0, 1 have k-th roots < 2, already handled).  */
    if (arbint_cmp_u32(r, 2u) < 0)
      continue;

    /*  Check if r^k == |a|.  */
    rc = arbint_pow_u32(r_pow, r, k);
    if (rc != ARBINT_OK)
      goto cleanup;

    if (arbint_eq(r_pow, abs_a)) {
      *out = 1;
      rc = ARBINT_OK;
      goto cleanup;
    }
  }

  /*  Check primes beyond static table using cached sieve.  */
  if (max_k > (size_t) arbint_is_power_primes[ARBINT_IS_POWER_NPRIMES - 1]) {
    uint32_t k;
    uint32_t max_k_u32 = (uint32_t) max_k;

    rc = arbint_cache_prime_sieve_ensure(max_k_u32);
    if (rc != ARBINT_OK)
      goto cleanup;

    for (k = arbint_is_power_primes[ARBINT_IS_POWER_NPRIMES - 1] + 1u;
         k <= max_k_u32; ++k) {
      if (!arbint_cache_prime_sieve_is_prime(k))
        continue;

      /*  Skip even exponents for negative a.  */
      if (a_neg && (k % 2u == 0u))
        continue;

      /*  Compute r = floor(|a|^(1/k)).  */
      rc = arbint_root(r, abs_a, k);
      if (rc != ARBINT_OK)
        goto cleanup;

      /*  Skip if r < 2.  */
      if (arbint_cmp_u32(r, 2u) < 0)
        continue;

      /*  Check if r^k == |a|.  */
      rc = arbint_pow_u32(r_pow, r, k);
      if (rc != ARBINT_OK)
        goto cleanup;

      if (arbint_eq(r_pow, abs_a)) {
        *out = 1;
        rc = ARBINT_OK;
        goto cleanup;
      }
    }
  }

  /*  No perfect power found.  */
  *out = 0;
  rc = ARBINT_OK;

cleanup:
  if (init_rpow)
    arbint_clear(r_pow);
  if (init_r)
    arbint_clear(r);
  if (init_abs)
    arbint_clear(abs_a);
  return rc;
}

/*  Remove all factors of p from a: write a = p^k * rest.
    Returns rest and k.  For a == 0, returns rest = 0, k = 0.  */
ARBINT_API arbint_err_t arbint_removefactor_u32(arbint_t rest,
                                                const arbint_t a, uint32_t p,
                                                uint32_t * k) {
  arbint_t cur, q, r;
  arbint_ctx_t * ctx;
  arbint_err_t rc;
  uint32_t count = 0;

  if (rest == NULL || a == NULL || k == NULL)
    return ARBINT_EINVAL;

  if (p == 0u)
    return ARBINT_EZERO;

  /*  p == 1 divides everything infinitely; treat as removing 0 factors.  */
  if (p == 1u) {
    *k = 0;
    return arbint_set(rest, a);
  }

  /*  a == 0: no factors to remove.  */
  if (arbint_is_zero(a)) {
    *k = 0;
    arbint_zero(rest);
    return ARBINT_OK;
  }

  ctx = rest[0]._ctx;
  if (ctx == NULL)
    ctx = a[0]._ctx;

  rc = arbint_init_all(ctx, cur, q, r, (arbint_t *) NULL);
  if (rc != ARBINT_OK)
    return rc;

  /*  cur = a (handles aliasing with rest).  */
  rc = arbint_set(cur, a);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Repeatedly divide by p while remainder is zero.  */
  for (;;) {
    rc = arbint_tdiv_qr_u32(q, r, cur, p);
    if (rc != ARBINT_OK)
      goto cleanup;

    if (!arbint_is_zero(r))
      break;

    /*  cur <- q for next iteration.  */
    arbint_swap(cur, q);
    ++count;
  }

  rc = arbint_set(rest, cur);
  *k = count;

cleanup:
  arbint_clear_all(cur, q, r, (arbint_t *) NULL);
  return rc;
}

/*  Jacobi symbol for u32 operands.  Uses binary algorithm.  */
static int arbint_jacobi_u32u32(uint32_t a, uint32_t n) {
  int result = 1;
  uint32_t u = a % n;
  uint32_t v = n;
  uint32_t t;

  while (u != 0u) {
    /*  Factor out powers of 2 from u.  */
    while ((u & 1u) == 0u) {
      u >>= 1u;
      /*  (2/v) = (-1)^((v^2-1)/8) = -1 iff v == 3 or 5 (mod 8).  */
      t = v & 7u;
      if (t == 3u || t == 5u)
        result = -result;
    }

    /*  Quadratic reciprocity: (u/v)(v/u) = (-1)^((u-1)/2 * (v-1)/2).  */
    if ((u & 3u) == 3u && (v & 3u) == 3u)
      result = -result;

    /*  Swap u and v, then reduce u mod v.  */
    t = u;
    u = v % t;
    v = t;
  }

  return (v == 1u) ? result : 0;
}

/*  Jacobi symbol (a/n).  n must be odd and positive.  */
ARBINT_API arbint_err_t arbint_jacobi(const arbint_t a, const arbint_t n,
                                      int * out) {
  arbint_t u, v;
  arbint_ctx_t * ctx;
  arbint_err_t rc;
  int result = 1;
  size_t ctz_val;

  if (a == NULL || n == NULL || out == NULL)
    return ARBINT_EINVAL;

  /*  n must be positive and odd.  */
  if (n[0]._sz <= 0)
    return ARBINT_EDOM;
  if ((ARBINT_CLIMBS(n)[0] & 1u) == 0u)
    return ARBINT_EDOM;

  /*  (a/1) = 1 for all a.  */
  if (n[0]._sz == 1 && ARBINT_CLIMBS(n)[0] == 1u) {
    *out = 1;
    return ARBINT_OK;
  }

  /*  Fast u32 path: if n fits in u32, reduce a mod n and compute directly.  */
  {
    uint32_t n_u32;
    if (arbint_get_u32(n, &n_u32) == ARBINT_OK) {
      arbint_t tmp;
      uint32_t a_u32;

      ctx = a[0]._ctx;
      if (ctx == NULL)
        ctx = n[0]._ctx;

      rc = arbint_init(tmp, ctx);
      if (rc != ARBINT_OK)
        return rc;

      rc = arbint_fdiv_r(tmp, a, n);
      if (rc != ARBINT_OK) {
        arbint_clear(tmp);
        return rc;
      }

      /*  tmp is now in [0, n), which fits in u32.  */
      rc = arbint_get_u32(tmp, &a_u32);
      arbint_clear(tmp);
      if (rc != ARBINT_OK)
        return rc;

      *out = arbint_jacobi_u32u32(a_u32, n_u32);
      return ARBINT_OK;
    }
  }

  /*  General big-int path using binary Jacobi algorithm.  */
  ctx = a[0]._ctx;
  if (ctx == NULL)
    ctx = n[0]._ctx;

  rc = arbint_init_all(ctx, u, v, (arbint_t *) NULL);
  if (rc != ARBINT_OK)
    return rc;

  /*  u = a mod n (floor division for canonical [0, n) result).  */
  rc = arbint_fdiv_r(u, a, n);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_set(v, n);
  if (rc != ARBINT_OK)
    goto cleanup;

  while (!arbint_is_zero(u)) {
    /*  Factor out powers of 2 from u.  */
    rc = arbint_ctz(u, &ctz_val);
    if (rc != ARBINT_OK)
      goto cleanup;

    if (ctz_val > 0u) {
      rc = arbint_shr(u, u, (uint32_t) ctz_val);
      if (rc != ARBINT_OK)
        goto cleanup;

      /*  (2/v)^ctz_val: flip result for each odd power of 2.  */
      if (ctz_val & 1u) {
        arbint_limb_t v8 = ARBINT_CLIMBS(v)[0] & 7u;
        if (v8 == 3u || v8 == 5u)
          result = -result;
      }
    }

    /*  Quadratic reciprocity.  */
    {
      arbint_limb_t u4 = ARBINT_CLIMBS(u)[0] & 3u;
      arbint_limb_t v4 = ARBINT_CLIMBS(v)[0] & 3u;
      if (u4 == 3u && v4 == 3u)
        result = -result;
    }

    /*  Swap u and v, then reduce u mod v.  */
    arbint_swap(u, v);
    rc = arbint_fdiv_r(u, u, v);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

  /*  If v == 1, gcd(a,n) == 1; otherwise gcd > 1.  */
  *out = (v[0]._sz == 1 && ARBINT_CLIMBS(v)[0] == 1u) ? result : 0;
  rc = ARBINT_OK;

cleanup:
  arbint_clear_all(u, v, (arbint_t *) NULL);
  return rc;
}

/*  Legendre symbol (a/p) where p is an odd prime.
    This is just the Jacobi symbol; primality is caller's responsibility.  */
ARBINT_API arbint_err_t arbint_legendre(const arbint_t a, const arbint_t p,
                                        int * out) {
  return arbint_jacobi(a, p, out);
}

/*  Kronecker symbol (a/n).  Extends Jacobi to all integers n.  */
ARBINT_API arbint_err_t arbint_kronecker(const arbint_t a, const arbint_t n,
                                         int * out) {
  arbint_t abs_n, odd_part;
  arbint_ctx_t * ctx;
  arbint_err_t rc;
  int result = 1;
  size_t v2;

  if (a == NULL || n == NULL || out == NULL)
    return ARBINT_EINVAL;

  /*  (a/0) = 1 if |a| = 1, else 0.  */
  if (arbint_is_zero(n)) {
    size_t an = arbint_abs_sz(a[0]._sz);
    *out = (an == 1u && ARBINT_CLIMBS(a)[0] == 1u) ? 1 : 0;
    return ARBINT_OK;
  }

  ctx = a[0]._ctx;
  if (ctx == NULL)
    ctx = n[0]._ctx;

  rc = arbint_init_all(ctx, abs_n, odd_part, (arbint_t *) NULL);
  if (rc != ARBINT_OK)
    return rc;

  /*  Handle negative n: (a/-1) = -1 if a < 0.  */
  if (n[0]._sz < 0) {
    if (a[0]._sz < 0)
      result = -result;
    rc = arbint_abs(abs_n, n);
  } else {
    rc = arbint_set(abs_n, n);
  }
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Handle the case |n| == 1 early.  */
  if (abs_n[0]._sz == 1 && ARBINT_CLIMBS(abs_n)[0] == 1u) {
    *out = result;
    rc = ARBINT_OK;
    goto cleanup;
  }

  /*  Extract 2^v2 from |n|.  */
  rc = arbint_ctz(abs_n, &v2);
  if (rc != ARBINT_OK)
    goto cleanup;

  if (v2 > 0u) {
    rc = arbint_shr(odd_part, abs_n, (uint32_t) v2);
    if (rc != ARBINT_OK)
      goto cleanup;

    /*  (a/2) = 0 if a is even and a != 0.  */
    if (a[0]._sz != 0 && (ARBINT_CLIMBS(a)[0] & 1u) == 0u) {
      *out = 0;
      rc = ARBINT_OK;
      goto cleanup;
    }

    /*  (a/2^v2): for each 2, flip if a == 3 or 5 (mod 8).  */
    if (v2 & 1u) {
      arbint_limb_t a8 = (a[0]._sz != 0) ? (ARBINT_CLIMBS(a)[0] & 7u) : 0u;
      if (a8 == 3u || a8 == 5u)
        result = -result;
    }
  } else {
    rc = arbint_set(odd_part, abs_n);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

  /*  If odd_part == 1, we are done.  */
  if (odd_part[0]._sz == 1 && ARBINT_CLIMBS(odd_part)[0] == 1u) {
    *out = result;
    rc = ARBINT_OK;
    goto cleanup;
  }

  /*  Delegate odd part to Jacobi.  */
  {
    int jac;
    rc = arbint_jacobi(a, odd_part, &jac);
    if (rc != ARBINT_OK)
      goto cleanup;
    result *= jac;
  }

  *out = result;
  rc = ARBINT_OK;

cleanup:
  arbint_clear_all(abs_n, odd_part, (arbint_t *) NULL);
  return rc;
}

/*  Moebius function for u32.  */
static int arbint_moebius_u32(uint32_t n) {
  int omega = 0;
  uint32_t p;

  if (n == 1u)
    return 1;

  /*  Handle factor of 2.  */
  if ((n & 1u) == 0u) {
    if ((n & 3u) == 0u)
      return 0; /*  4 | n  */
    n >>= 1u;
    ++omega;
  }

  /*  Trial divide by odd primes.  */
  for (p = 3u; p <= n / p; p += 2u) {
    if (n % p == 0u) {
      n /= p;
      if (n % p == 0u)
        return 0; /*  p^2 | n  */
      ++omega;
    }
  }

  /*  If n > 1, it is a prime factor.  */
  if (n > 1u)
    ++omega;

  return (omega & 1) ? -1 : 1;
}

/*  Moebius function mu(n).
    mu(n) = 0 if n has a squared prime factor.
    mu(n) = (-1)^k if n is product of k distinct primes.  */
ARBINT_API arbint_err_t arbint_moebius(int * out, const arbint_t n) {
  arbint_t m;
  arbint_ctx_t * ctx;
  arbint_err_t rc;
  int omega = 0;
  size_t i;

  if (out == NULL || n == NULL)
    return ARBINT_EINVAL;

  /*  mu(n) is only defined for positive integers.  */
  if (n[0]._sz <= 0)
    return ARBINT_EDOM;

  /*  mu(1) = 1.  */
  if (n[0]._sz == 1 && ARBINT_CLIMBS(n)[0] == 1u) {
    *out = 1;
    return ARBINT_OK;
  }

  /*  Fast u32 path.  */
  {
    uint32_t n_u32;
    if (arbint_get_u32(n, &n_u32) == ARBINT_OK) {
      *out = arbint_moebius_u32(n_u32);
      return ARBINT_OK;
    }
  }

  /*  General path: trial division with early exit on squared factors.  */
  ctx = n[0]._ctx;
  rc = arbint_init(m, ctx);
  if (rc != ARBINT_OK)
    return rc;

  rc = arbint_set(m, n);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Trial divide by small primes from table.  */
  for (i = 0u; i < ARBINT_IS_POWER_NPRIMES; ++i) {
    uint32_t prime = arbint_is_power_primes[i];
    uint32_t k;

    rc = arbint_removefactor_u32(m, m, prime, &k);
    if (rc != ARBINT_OK)
      goto cleanup;

    if (k >= 2u) {
      *out = 0;
      rc = ARBINT_OK;
      goto cleanup;
    }
    if (k == 1u)
      ++omega;

    /*  Early exit if m == 1.  */
    if (m[0]._sz == 1 && ARBINT_CLIMBS(m)[0] == 1u)
      break;
  }

  /*  If m > 1 after trial division, it is either prime or has large factors.
   */
  if (arbint_cmp_u32(m, 1u) > 0) {
    arbint_t p, p2, q, r;
    uint32_t p_start;

    rc = arbint_init_all(ctx, p, p2, q, r, (arbint_t *) NULL);
    if (rc != ARBINT_OK)
      goto cleanup;

    p_start = arbint_is_power_primes[ARBINT_IS_POWER_NPRIMES - 1] + 2u;
    rc = arbint_set_u32(p, p_start);
    if (rc != ARBINT_OK) {
      arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
      goto cleanup;
    }

    while (arbint_cmp_u32(m, 1u) > 0) {
      rc = arbint_sqr(p2, p);
      if (rc != ARBINT_OK) {
        arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
        goto cleanup;
      }

      /*  If p^2 > m, then m is prime.  */
      if (arbint_cmp(p2, m) > 0) {
        ++omega;
        break;
      }

      rc = arbint_tdiv_qr(q, r, m, p);
      if (rc != ARBINT_OK) {
        arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
        goto cleanup;
      }

      if (arbint_is_zero(r)) {
        /*  p divides m; check for p^2.  */
        rc = arbint_set(m, q);
        if (rc != ARBINT_OK) {
          arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
          goto cleanup;
        }

        rc = arbint_tdiv_r(r, m, p);
        if (rc != ARBINT_OK) {
          arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
          goto cleanup;
        }

        if (arbint_is_zero(r)) {
          /*  p^2 | original n.  */
          *out = 0;
          arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
          rc = ARBINT_OK;
          goto cleanup;
        }
        ++omega;
      }

      rc = arbint_add_u32(p, p, 2u);
      if (rc != ARBINT_OK) {
        arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
        goto cleanup;
      }
    }

    arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
  }

  *out = (omega & 1) ? -1 : 1;
  rc = ARBINT_OK;

cleanup:
  arbint_clear(m);
  return rc;
}

/*  Carmichael function lambda(n).
    lambda(2) = 1, lambda(4) = 2, lambda(2^k) = 2^(k-2) for k >= 3.
    lambda(p^k) = p^(k-1) * (p-1) for odd prime p.
    lambda(n) = lcm of lambda(p_i^{k_i}) for prime factorization.  */
ARBINT_API arbint_err_t arbint_carmichael(arbint_t rop, const arbint_t n) {
  arbint_t m, lambda, contrib;
  arbint_ctx_t * ctx;
  arbint_err_t rc;
  size_t i;

  if (rop == NULL || n == NULL)
    return ARBINT_EINVAL;

  /*  Carmichael is only defined for positive integers.  */
  if (n[0]._sz <= 0)
    return ARBINT_EDOM;

  /*  lambda(1) = 1.  */
  if (n[0]._sz == 1 && ARBINT_CLIMBS(n)[0] == 1u)
    return arbint_set_u32(rop, 1u);

  ctx = rop[0]._ctx;
  if (ctx == NULL)
    ctx = n[0]._ctx;

  rc = arbint_init_all(ctx, m, lambda, contrib, (arbint_t *) NULL);
  if (rc != ARBINT_OK)
    return rc;

  rc = arbint_set(m, n);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_set_u32(lambda, 1u);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Handle factor of 2 specially.  */
  {
    uint32_t k2;
    rc = arbint_removefactor_u32(m, m, 2u, &k2);
    if (rc != ARBINT_OK)
      goto cleanup;

    if (k2 > 0u) {
      uint32_t lam2;
      /*  lambda(2) = 1, lambda(4) = 2, lambda(2^k) = 2^(k-2) for k >= 3.  */
      if (k2 == 1u)
        lam2 = 1u;
      else if (k2 == 2u)
        lam2 = 2u;
      else
        lam2 = 1u << (k2 - 2u);

      rc = arbint_set_u32(contrib, lam2);
      if (rc != ARBINT_OK)
        goto cleanup;

      rc = arbint_lcm(lambda, lambda, contrib);
      if (rc != ARBINT_OK)
        goto cleanup;
    }
  }

  /*  Handle odd primes from the small primes table.  */
  for (i = 1u; i < ARBINT_IS_POWER_NPRIMES; ++i) {
    uint32_t prime = arbint_is_power_primes[i];
    uint32_t k;
    uint32_t j;

    if (m[0]._sz == 1 && ARBINT_CLIMBS(m)[0] == 1u)
      break;

    rc = arbint_removefactor_u32(m, m, prime, &k);
    if (rc != ARBINT_OK)
      goto cleanup;

    if (k > 0u) {
      /*  lambda(p^k) = p^(k-1) * (p-1).  */
      rc = arbint_set_u32(contrib, prime - 1u);
      if (rc != ARBINT_OK)
        goto cleanup;

      for (j = 1u; j < k; ++j) {
        rc = arbint_mul_u32(contrib, contrib, prime);
        if (rc != ARBINT_OK)
          goto cleanup;
      }

      rc = arbint_lcm(lambda, lambda, contrib);
      if (rc != ARBINT_OK)
        goto cleanup;
    }
  }

  /*  Handle remaining large prime factors.  */
  if (arbint_cmp_u32(m, 1u) > 0) {
    arbint_t p, p2, q, r;
    uint32_t p_start;

    rc = arbint_init_all(ctx, p, p2, q, r, (arbint_t *) NULL);
    if (rc != ARBINT_OK)
      goto cleanup;

    p_start = arbint_is_power_primes[ARBINT_IS_POWER_NPRIMES - 1] + 2u;
    rc = arbint_set_u32(p, p_start);
    if (rc != ARBINT_OK) {
      arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
      goto cleanup;
    }

    while (arbint_cmp_u32(m, 1u) > 0) {
      rc = arbint_sqr(p2, p);
      if (rc != ARBINT_OK) {
        arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
        goto cleanup;
      }

      /*  If p^2 > m, then m is prime: lambda(m) = m - 1.  */
      if (arbint_cmp(p2, m) > 0) {
        rc = arbint_sub_i32(contrib, m, 1);
        if (rc != ARBINT_OK) {
          arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
          goto cleanup;
        }

        rc = arbint_lcm(lambda, lambda, contrib);
        arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
        if (rc != ARBINT_OK)
          goto cleanup;
        break;
      }

      rc = arbint_tdiv_qr(q, r, m, p);
      if (rc != ARBINT_OK) {
        arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
        goto cleanup;
      }

      if (arbint_is_zero(r)) {
        uint32_t k = 1u;
        rc = arbint_set(m, q);
        if (rc != ARBINT_OK) {
          arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
          goto cleanup;
        }

        /*  Count additional factors of p.  */
        for (;;) {
          rc = arbint_tdiv_qr(q, r, m, p);
          if (rc != ARBINT_OK) {
            arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
            goto cleanup;
          }
          if (!arbint_is_zero(r))
            break;
          rc = arbint_set(m, q);
          if (rc != ARBINT_OK) {
            arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
            goto cleanup;
          }
          ++k;
        }

        /*  lambda(p^k) = p^(k-1) * (p-1).  */
        rc = arbint_sub_i32(contrib, p, 1);
        if (rc != ARBINT_OK) {
          arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
          goto cleanup;
        }

        if (k > 1u) {
          rc = arbint_pow_u32(p2, p, k - 1u);
          if (rc != ARBINT_OK) {
            arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
            goto cleanup;
          }

          rc = arbint_mul(contrib, contrib, p2);
          if (rc != ARBINT_OK) {
            arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
            goto cleanup;
          }
        }

        rc = arbint_lcm(lambda, lambda, contrib);
        if (rc != ARBINT_OK) {
          arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
          goto cleanup;
        }
      }

      rc = arbint_add_u32(p, p, 2u);
      if (rc != ARBINT_OK) {
        arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
        goto cleanup;
      }
    }

    /*  Clean up if we broke out of the loop without cleaning.  */
    if (arbint_cmp_u32(m, 1u) <= 0)
      arbint_clear_all(p, p2, q, r, (arbint_t *) NULL);
  }

  rc = arbint_set(rop, lambda);

cleanup:
  arbint_clear_all(m, lambda, contrib, (arbint_t *) NULL);
  return rc;
}
