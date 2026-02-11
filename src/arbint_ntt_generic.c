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

#include "arbint_ntt.h"
#include "config.h"

#include <stdlib.h>
#include <string.h>

/* ========== Wide Multiplication ========== */

/*  Multiply two 64-bit integers producing full 128-bit result (hi:lo = a * b).

    Uses half-limb (32x32) multiplication to synthesize 128-bit product
    via Karatsuba-like decomposition (same pattern as arbint_div_generic.c).  */
static inline void ntt_umul(uint64_t * hi, uint64_t * lo, uint64_t a,
                            uint64_t b) {
  /*  Half-limb decomposition:
      a = a1*2^32 + a0, b = b1*2^32 + b0
      a*b = a1*b1*2^64 + (a1*b0 + a0*b1)*2^32 + a0*b0  */
  uint64_t a0 = a & 0xFFFFFFFFu;
  uint64_t a1 = a >> 32;
  uint64_t b0 = b & 0xFFFFFFFFu;
  uint64_t b1 = b >> 32;

  uint64_t w0 = a0 * b0;
  uint64_t t = a1 * b0 + (w0 >> 32);
  uint64_t w1 = t & 0xFFFFFFFFu;
  uint64_t w2 = t >> 32;

  w1 = a0 * b1 + w1;

  *hi = a1 * b1 + w2 + (w1 >> 32);
  *lo = (w1 << 32) | (w0 & 0xFFFFFFFFu);
}

/* ========== Montgomery Arithmetic ========== */

/*  Montgomery reduction: compute (t_hi:t_lo) * R^(-1) mod p.

    REDC algorithm:
    1. m = t_lo * p_neg_inv (mod 2^64)
    2. u = (t + m * p) / 2^64
    3. return (u >= p) ? u - p : u

    For primes close to 2^64 (like our NTT primes), the addition
    t_hi + mp_hi can overflow. When this happens, the true result is
    2^64 + sum. Since 2^64 > p, we must reduce by adding (2^64 - p) to
    the truncated sum, then do the final conditional subtraction.  */
static inline uint64_t mont_redc(uint64_t t_hi, uint64_t t_lo,
                                 const arbint_mont_params_t * m) {
  uint64_t mm = t_lo * m->p_neg_inv;
  uint64_t mp_hi, mp_lo;
  ntt_umul(&mp_hi, &mp_lo, mm, m->p);

  /*  Add t_lo + mp_lo. By construction of mm, this equals 0 mod 2^64.  */
  uint64_t carry1 = (t_lo + mp_lo < t_lo) ? 1u : 0u;

  /*  Add t_hi + mp_hi, tracking overflow.  */
  uint64_t sum = t_hi + mp_hi;
  uint64_t overflow = (sum < t_hi) ? 1u : 0u;

  /*  Add the carry from low word, which might cause another overflow.  */
  sum += carry1;
  overflow += (sum < carry1) ? 1u : 0u;

  /*  If overflow occurred, the true result is 2^64 + sum.
      Since 2^64 > p, we reduce: result = 2^64 + sum - p = sum + (2^64 - p).
      Note: 2^64 - p = (UINT64_MAX - p) + 1.  */
  if (overflow) {
    uint64_t adjustment = (UINT64_MAX - m->p) + 1u;
    sum += adjustment;
    /*  After adjustment, sum is in [0, 2p), may need final reduction.  */
  }

  /*  Final conditional subtraction (branchless).  */
  uint64_t mask = (sum >= m->p) ? UINT64_MAX : 0u;
  return sum - (mask & m->p);
}

/*  Montgomery multiplication: compute (a * b * R^(-1)) mod p.
    Both a and b must already be in Montgomery form.  */
static inline uint64_t mont_mul(uint64_t a, uint64_t b,
                                const arbint_mont_params_t * m) {
  uint64_t hi, lo;
  ntt_umul(&hi, &lo, a, b);
  return mont_redc(hi, lo, m);
}

/*  Convert to Montgomery form: compute (x * R) mod p
    = (x * R^2 * R^(-1)) mod p.  */
static inline uint64_t to_mont(uint64_t x, const arbint_mont_params_t * m) {
  return mont_mul(x, m->r_squared, m);
}

/*  Convert from Montgomery form: compute (x * R^(-1)) mod p.  */
static inline uint64_t from_mont(uint64_t x, const arbint_mont_params_t * m) {
  return mont_redc(0, x, m);
}

/*  Modular addition: (a + b) mod p, where a, b < p.
    Branchless implementation.  */
static inline uint64_t mod_add(uint64_t a, uint64_t b, uint64_t p) {
  /*  Avoid 64-bit overflow pitfalls for primes close to 2^64.
      Since a,b < p, compute:
        if a >= p-b: a+b-p
        else        : a+b
      This is exact in uint64_t arithmetic without wide intermediates.  */
  uint64_t p_minus_b = p - b;
  if (a >= p_minus_b)
    return a - p_minus_b;
  return a + b;
}

/*  Modular subtraction: (a - b) mod p, where a, b < p.
    Branchless implementation.  */
static inline uint64_t mod_sub(uint64_t a, uint64_t b, uint64_t p) {
  uint64_t diff = a - b;
  /*  If a < b (borrow), add p.  */
  uint64_t borrow = (a < b) ? UINT64_MAX : 0u;
  return diff + (borrow & p);
}

/* ========== Montgomery Parameter Computation ========== */

/*  Compute -p^(-1) mod 2^64 using Newton iteration.
    Uses the identity: p * p^(-1) = 1 (mod 2^k) =>
    p^(-1) = p^(-1) * (2 - p * p^(-1)) (mod 2^(2k)).

    Starting with x = 1 (valid mod 2^1 since p is odd), iterate to get
    inverse mod 2^64.  */
static uint64_t compute_p_neg_inv(uint64_t p) {
  uint64_t x = 1u;
  /*  Newton iteration: x = x * (2 - p * x) mod 2^64.
      Each iteration doubles the number of correct bits.
      6 iterations: 1 -> 2 -> 4 -> 8 -> 16 -> 32 -> 64 bits.  */
  x = x * (2u - p * x); /* mod 2^2 */
  x = x * (2u - p * x); /* mod 2^4 */
  x = x * (2u - p * x); /* mod 2^8 */
  x = x * (2u - p * x); /* mod 2^16 */
  x = x * (2u - p * x); /* mod 2^32 */
  x = x * (2u - p * x); /* mod 2^64 */
  return (uint64_t) 0u - x;  /* negate to get -p^(-1) */
}

/*  Compute R^2 mod p where R = 2^64.

    Since R mod p = 2^64 mod p, and R < 2p for our primes (all near 2^64),
    we compute R mod p = 2^64 - p, then square it modularly.

    For squaring mod p without overflow, we use repeated doubling:
    R mod p, then (R mod p)^2 mod p.  */
static uint64_t compute_r_squared(uint64_t p) {
  /*  R mod p = 2^64 mod p.
      Since p > 2^63, we have R mod p = 2^64 - p.  */
  uint64_t r_mod_p = (uint64_t) 0u - p;

  /*  Compute r_mod_p^2 mod p using 128-bit arithmetic.  */
  uint64_t hi, lo;
  ntt_umul(&hi, &lo, r_mod_p, r_mod_p);

  /*  Reduce hi:lo mod p.
      Since hi < p and lo < 2^64, we have hi:lo < p * 2^64.
      Compute hi:lo - floor(hi:lo / p) * p.

      For efficiency, we use: result = hi:lo - hi * floor(2^64 / p) * p
      and correct.

      Simpler: use the fact that p is close to 2^64.
      hi:lo mod p = (hi * (2^64 mod p) + lo) mod p
                  = (hi * r_mod_p + lo) mod p.  */
  uint64_t hi2, lo2;
  ntt_umul(&hi2, &lo2, hi, r_mod_p);

  /*  Add lo.  */
  uint64_t sum_lo = lo2 + lo;
  uint64_t carry = (sum_lo < lo2) ? 1u : 0u;
  uint64_t sum_hi = hi2 + carry;

  /*  Now reduce sum_hi:sum_lo mod p, recursively.  */
  if (sum_hi == 0u) {
    return (sum_lo >= p) ? (sum_lo - p) : sum_lo;
  }

  /*  sum_hi:sum_lo = sum_hi * 2^64 + sum_lo
      = sum_hi * (p + r_mod_p) + sum_lo (since 2^64 = p + r_mod_p)
      = sum_hi * p + sum_hi * r_mod_p + sum_lo
      = sum_hi * r_mod_p + sum_lo (mod p).  */
  uint64_t hi3, lo3;
  ntt_umul(&hi3, &lo3, sum_hi, r_mod_p);
  uint64_t final_lo = lo3 + sum_lo;
  uint64_t final_hi = hi3 + ((final_lo < lo3) ? 1u : 0u);

  /*  One more iteration if needed.  */
  if (final_hi > 0u) {
    uint64_t hi4, lo4;
    ntt_umul(&hi4, &lo4, final_hi, r_mod_p);
    final_lo = lo4 + final_lo;
    final_hi = hi4 + ((final_lo < lo4) ? 1u : 0u);
  }

  /*  At this point final_hi should be 0. Reduce final_lo.  */
  while (final_lo >= p)
    final_lo -= p;

  return final_lo;
}

/*  Initialize Montgomery parameters for a prime p.  */
static void init_mont_params(arbint_mont_params_t * m, uint64_t p) {
  m->p = p;
  m->p_neg_inv = compute_p_neg_inv(p);
  m->r_squared = compute_r_squared(p);
}

/* ========== Modular Exponentiation (for root computation) ========== */

/*  Modular exponentiation: compute base^exp mod p.
    Uses square-and-multiply algorithm.  */
static uint64_t powmod(uint64_t base, uint64_t exp, uint64_t p) {
  uint64_t result = 1u;
  base = base % p;

  while (exp > 0u) {
    if (exp & 1u) {
      uint64_t hi, lo;
      ntt_umul(&hi, &lo, result, base);
      /*  Reduce mod p using same technique as compute_r_squared.  */
      uint64_t r_mod_p = (uint64_t) 0u - p;
      uint64_t hi2, lo2;
      ntt_umul(&hi2, &lo2, hi, r_mod_p);
      uint64_t sum = lo2 + lo;
      uint64_t carry = (sum < lo2) ? 1u : 0u;
      uint64_t sum_hi = hi2 + carry;
      while (sum_hi > 0u) {
        uint64_t hi3, lo3;
        ntt_umul(&hi3, &lo3, sum_hi, r_mod_p);
        sum = lo3 + sum;
        sum_hi = hi3 + ((sum < lo3) ? 1u : 0u);
      }
      while (sum >= p)
        sum -= p;
      result = sum;
    }
    exp >>= 1;
    if (exp > 0u) {
      uint64_t hi, lo;
      ntt_umul(&hi, &lo, base, base);
      uint64_t r_mod_p = (uint64_t) 0u - p;
      uint64_t hi2, lo2;
      ntt_umul(&hi2, &lo2, hi, r_mod_p);
      uint64_t sum = lo2 + lo;
      uint64_t carry = (sum < lo2) ? 1u : 0u;
      uint64_t sum_hi = hi2 + carry;
      while (sum_hi > 0u) {
        uint64_t hi3, lo3;
        ntt_umul(&hi3, &lo3, sum_hi, r_mod_p);
        sum = lo3 + sum;
        sum_hi = hi3 + ((sum < lo3) ? 1u : 0u);
      }
      while (sum >= p)
        sum -= p;
      base = sum;
    }
  }
  return result;
}

/*  Compute modular inverse: a^(-1) mod p using Fermat's little theorem.
    a^(-1) = a^(p-2) mod p.  */
static uint64_t modinv(uint64_t a, uint64_t p) { return powmod(a, p - 2u, p); }

/* ========== 192-bit Wide Integer Arithmetic for CRT ========== */

/*  Add two 192-bit wide integers: r = a + b.  */
static inline void crt_wide_add(arbint_crt_wide_t * r,
                                const arbint_crt_wide_t * a,
                                const arbint_crt_wide_t * b) {
  uint64_t lo = a->lo + b->lo;
  uint64_t c1 = (lo < a->lo) ? 1u : 0u;
  uint64_t mid = a->mid + b->mid + c1;
  uint64_t c2 = (mid < a->mid || (c1 && mid == a->mid)) ? 1u : 0u;
  uint64_t hi = a->hi + b->hi + c2;
  r->lo = lo;
  r->mid = mid;
  r->hi = hi;
}

/*  Divide wide integer by limb base B = 2^64, return remainder.
    x = x / B (in place), returns x mod B.  */
static inline uint64_t crt_wide_divmod_base64(arbint_crt_wide_t * x) {
  uint64_t rem = x->lo;
  x->lo = x->mid;
  x->mid = x->hi;
  x->hi = 0u;
  return rem;
}

#if ARBINT_LIMB_BITS == 32
/*  Divide wide integer by limb base B = 2^32, return remainder.  */
static inline uint32_t crt_wide_divmod_base32(arbint_crt_wide_t * x) {
  uint32_t rem = (uint32_t) x->lo;
  x->lo = (x->lo >> 32) | (x->mid << 32);
  x->mid = (x->mid >> 32) | (x->hi << 32);
  x->hi = x->hi >> 32;
  return rem;
}
#endif

/*  Set wide integer from a single 64-bit value.  */
static inline void crt_wide_set_u64(arbint_crt_wide_t * r, uint64_t v) {
  r->lo = v;
  r->mid = 0u;
  r->hi = 0u;
}

/* ========== NTT Context Initialization ========== */

/*  Global NTT context - lazy initialized.  */
static arbint_ntt_ctx_t g_ntt_ctx;
static int g_ntt_initialized = 0;

/*  Initialize NTT context with precomputed parameters.
    Called once on first NTT use.  */
static arbint_err_t arbint_ntt_ctx_init_internal(void) {
  /*  Prime 1: Goldilocks prime.  */
  g_ntt_ctx.primes[0].p = ARBINT_NTT_P1;
  g_ntt_ctx.primes[0].g = ARBINT_NTT_G1;
  g_ntt_ctx.primes[0].max_log2 = ARBINT_NTT_P1_MAX_LOG2;
  init_mont_params(&g_ntt_ctx.primes[0].mont, ARBINT_NTT_P1);

  /*  Prime 2.  */
  g_ntt_ctx.primes[1].p = ARBINT_NTT_P2;
  g_ntt_ctx.primes[1].g = ARBINT_NTT_G2;
  g_ntt_ctx.primes[1].max_log2 = ARBINT_NTT_P2_MAX_LOG2;
  init_mont_params(&g_ntt_ctx.primes[1].mont, ARBINT_NTT_P2);

  /*  Prime 3.  */
  g_ntt_ctx.primes[2].p = ARBINT_NTT_P3;
  g_ntt_ctx.primes[2].g = ARBINT_NTT_G3;
  g_ntt_ctx.primes[2].max_log2 = ARBINT_NTT_P3_MAX_LOG2;
  init_mont_params(&g_ntt_ctx.primes[2].mont, ARBINT_NTT_P3);

  /*  Initialize root tables to NULL (allocated on demand).  */
  for (int i = 0; i < 3; ++i) {
    g_ntt_ctx.roots[i].omega = NULL;
    g_ntt_ctx.roots[i].omega_inv = NULL;
    g_ntt_ctx.roots[i].max_log2 = 0u;
  }

  /*  Compute CRT constants.  */
  g_ntt_ctx.crt.p1_inv_mod_p2 = modinv(ARBINT_NTT_P1, ARBINT_NTT_P2);

  /*  Compute p1 * p2 mod p3, then invert.  */
  uint64_t p1_mod_p3 = ARBINT_NTT_P1 % ARBINT_NTT_P3;
  uint64_t p2_mod_p3 = ARBINT_NTT_P2 % ARBINT_NTT_P3;
  uint64_t hi, lo;
  ntt_umul(&hi, &lo, p1_mod_p3, p2_mod_p3);
  /*  Reduce hi:lo mod p3.  */
  uint64_t r_mod_p3 = (uint64_t) 0u - ARBINT_NTT_P3;
  uint64_t hi2, lo2;
  ntt_umul(&hi2, &lo2, hi, r_mod_p3);
  uint64_t sum = lo2 + lo;
  uint64_t sum_hi = hi2 + ((sum < lo2) ? 1u : 0u);
  while (sum_hi > 0u) {
    uint64_t hi3, lo3;
    ntt_umul(&hi3, &lo3, sum_hi, r_mod_p3);
    uint64_t new_sum = lo3 + sum;
    sum_hi = hi3 + ((new_sum < lo3) ? 1u : 0u);
    sum = new_sum;
  }
  while (sum >= ARBINT_NTT_P3)
    sum -= ARBINT_NTT_P3;
  uint64_t p1p2_mod_p3 = sum;
  g_ntt_ctx.crt.p1p2_inv_mod_p3 = modinv(p1p2_mod_p3, ARBINT_NTT_P3);

  /*  Store p1 and p1*p2 for CRT reconstruction.  */
  g_ntt_ctx.crt.p1_lo = ARBINT_NTT_P1;
  ntt_umul(&g_ntt_ctx.crt.p1p2_mid, &g_ntt_ctx.crt.p1p2_lo, ARBINT_NTT_P1,
           ARBINT_NTT_P2);

  g_ntt_ctx.initialized = 1;

#ifndef NDEBUG
  /*  Debug validation: verify primitive roots.  */
  for (int i = 0; i < 3; ++i) {
    uint64_t p = g_ntt_ctx.primes[i].p;
    uint64_t g = g_ntt_ctx.primes[i].g;
    /*  g^((p-1)/2) should not be 1 (g is not a quadratic residue).  */
    uint64_t half_order = powmod(g, (p - 1u) / 2u, p);
    if (half_order == 1u) {
      /*  Invalid primitive root - this shouldn't happen with correct
       * constants.  */
      return ARBINT_EINVAL;
    }
    /*  g^(p-1) should be 1 (Fermat's little theorem).  */
    uint64_t full_order = powmod(g, p - 1u, p);
    if (full_order != 1u) {
      return ARBINT_EINVAL;
    }
    /*  Verify Montgomery roundtrip.  */
    uint64_t test_val = 0x123456789ABCDEFull % p;
    uint64_t mont_val = to_mont(test_val, &g_ntt_ctx.primes[i].mont);
    uint64_t back_val = from_mont(mont_val, &g_ntt_ctx.primes[i].mont);
    if (back_val != test_val) {
      return ARBINT_EINVAL;
    }
  }
#endif

  return ARBINT_OK;
}

/*  Ensure NTT context is initialized.  */
static arbint_err_t arbint_ntt_ensure_init(void) {
  if (g_ntt_initialized)
    return ARBINT_OK;
  arbint_err_t rc = arbint_ntt_ctx_init_internal();
  if (rc == ARBINT_OK)
    g_ntt_initialized = 1;
  return rc;
}

/* ========== Root Table Management ========== */

/*  Ensure root tables are allocated for at least the given log2 size.  */
static arbint_err_t ensure_roots(size_t log2_n, const arbint_alloc_t * alloc) {
  if (log2_n > ARBINT_NTT_MAX_LOG2)
    return ARBINT_EOVERFLOW;

  for (int i = 0; i < 3; ++i) {
    arbint_ntt_roots_t * r = &g_ntt_ctx.roots[i];
    /*  Need stages [0, log2_n-1]. max_log2 stores count of initialized
        stages (exclusive upper bound).  */
    if (r->max_log2 >= log2_n)
      continue;

    /*  Allocate or reallocate root tables.  */
    size_t new_size = log2_n;
    if (new_size == 0u)
      new_size = 1u;
    uint64_t * new_omega =
        (uint64_t *) alloc->realloc(alloc->ud, r->omega,
                                    new_size * sizeof(uint64_t));
    if (new_omega == NULL)
      return ARBINT_ENOMEM;

    uint64_t * new_omega_inv =
        (uint64_t *) alloc->realloc(alloc->ud, r->omega_inv,
                                    new_size * sizeof(uint64_t));
    if (new_omega_inv == NULL) {
      alloc->realloc(alloc->ud, new_omega, 0);
      return ARBINT_ENOMEM;
    }

    r->omega = new_omega;
    r->omega_inv = new_omega_inv;

    /*  Compute roots for new sizes.  */
    uint64_t p = g_ntt_ctx.primes[i].p;
    uint64_t g = g_ntt_ctx.primes[i].g;
    const arbint_mont_params_t * m = &g_ntt_ctx.primes[i].mont;

    for (size_t k = r->max_log2; k < log2_n; ++k) {
      /*  omega[k] = g^((p-1) / 2^(k+1)) mod p.
          This is a primitive 2^(k+1)-th root of unity.  */
      uint64_t exp = (p - 1u) >> (k + 1u);
      uint64_t root = powmod(g, exp, p);
      uint64_t root_inv = modinv(root, p);

      /*  Store in Montgomery form.  */
      r->omega[k] = to_mont(root, m);
      r->omega_inv[k] = to_mont(root_inv, m);
    }

    r->max_log2 = (unsigned) log2_n;
  }

  return ARBINT_OK;
}

/* ========== NTT Transforms ========== */

/*  Bit-reversal permutation index.  */
static inline size_t bit_reverse(size_t x, size_t log2_n) {
  size_t result = 0u;
  for (size_t i = 0u; i < log2_n; ++i) {
    result = (result << 1) | (x & 1u);
    x >>= 1;
  }
  return result;
}

/*  Forward NTT: Cooley-Tukey decimation-in-time.
    x must be in Montgomery form on entry; remains in Montgomery form on exit.
    Uses iterative algorithm with implicit bit-reversal.  */
static void ntt_forward(uint64_t * x, size_t log2_n,
                        const arbint_ntt_roots_t * roots,
                        const arbint_mont_params_t * mont) {
  size_t n = (size_t) 1u << log2_n;
  uint64_t p = mont->p;

  /*  Bit-reversal permutation.  */
  for (size_t i = 0u; i < n; ++i) {
    size_t j = bit_reverse(i, log2_n);
    if (i < j) {
      uint64_t tmp = x[i];
      x[i] = x[j];
      x[j] = tmp;
    }
  }

  /*  Cooley-Tukey butterfly iterations.  */
  for (size_t s = 0u; s < log2_n; ++s) {
    size_t m = (size_t) 1u << (s + 1u);
    size_t m_half = m >> 1;

    /*  Get primitive m-th root of unity (omega[s] is 2^(s+1)-th root).  */
    uint64_t omega_m = roots->omega[s];

    for (size_t k = 0u; k < n; k += m) {
      uint64_t omega = to_mont(1u, mont); /* omega^0 = 1 in Montgomery form */
      for (size_t j = 0u; j < m_half; ++j) {
        uint64_t t = mont_mul(omega, x[k + j + m_half], mont);
        uint64_t u = x[k + j];
        x[k + j] = mod_add(u, t, p);
        x[k + j + m_half] = mod_sub(u, t, p);
        omega = mont_mul(omega, omega_m, mont);
      }
    }
  }
}

/*  Inverse NTT: Gentleman-Sande decimation-in-frequency.
    x must be in Montgomery form on entry.
    On exit, x contains INTT(x) * n (not yet scaled by n^-1).  */
static void ntt_inverse(uint64_t * x, size_t log2_n,
                        const arbint_ntt_roots_t * roots,
                        const arbint_mont_params_t * mont) {
  size_t n = (size_t) 1u << log2_n;
  uint64_t p = mont->p;

  /*  Gentleman-Sande butterfly iterations (reverse order of forward).  */
  for (size_t s = log2_n; s > 0u; --s) {
    size_t m = (size_t) 1u << s;
    size_t m_half = m >> 1;

    /*  Get inverse of primitive m-th root of unity.  */
    uint64_t omega_m_inv = roots->omega_inv[s - 1u];

    for (size_t k = 0u; k < n; k += m) {
      uint64_t omega_inv = to_mont(1u, mont);
      for (size_t j = 0u; j < m_half; ++j) {
        uint64_t u = x[k + j];
        uint64_t v = x[k + j + m_half];
        x[k + j] = mod_add(u, v, p);
        x[k + j + m_half] = mont_mul(mod_sub(u, v, p), omega_inv, mont);
        omega_inv = mont_mul(omega_inv, omega_m_inv, mont);
      }
    }
  }

  /*  Bit-reversal permutation.  */
  for (size_t i = 0u; i < n; ++i) {
    size_t j = bit_reverse(i, log2_n);
    if (i < j) {
      uint64_t tmp = x[i];
      x[i] = x[j];
      x[j] = tmp;
    }
  }
}

/* ========== CRT Reconstruction ========== */

/*  Combine results from 3 NTTs via Garner's algorithm.
    r1, r2, r3 are residues mod p1, p2, p3 (in normal form, not Montgomery).
    Output is written to dst as arbint limbs.  */
static void ntt_crt_combine(arbint_limb_t * dst, size_t * dst_used,
                            const uint64_t * r1, const uint64_t * r2,
                            const uint64_t * r3, size_t n,
                            const arbint_ntt_crt_t * crt) {
  arbint_crt_wide_t carry;
  crt_wide_set_u64(&carry, 0u);

  size_t out_idx = 0u;

  for (size_t i = 0u; i < n; ++i) {
    /*  Garner's algorithm:
        v1 = r1[i]
        v2 = (r2[i] - v1) * p1_inv_mod_p2  mod p2
        v3 = ((r3[i] - v1) - v2 * p1) * p1p2_inv_mod_p3  mod p3

        result = v1 + v2 * p1 + v3 * p1 * p2  */

    uint64_t v1 = r1[i];

    /*  v2 = (r2[i] - v1) * p1_inv_mod_p2 mod p2.  */
    uint64_t diff2 = mod_sub(r2[i], v1 % ARBINT_NTT_P2, ARBINT_NTT_P2);
    uint64_t hi, lo;
    ntt_umul(&hi, &lo, diff2, crt->p1_inv_mod_p2);
    /*  Reduce mod p2.  */
    uint64_t r_mod_p2 = (uint64_t) 0u - ARBINT_NTT_P2;
    uint64_t hi2, lo2;
    ntt_umul(&hi2, &lo2, hi, r_mod_p2);
    uint64_t sum = lo2 + lo;
    uint64_t sum_hi = hi2 + ((sum < lo2) ? 1u : 0u);
    while (sum_hi > 0u) {
      uint64_t hi3, lo3;
      ntt_umul(&hi3, &lo3, sum_hi, r_mod_p2);
      uint64_t new_sum = lo3 + sum;
      sum_hi = hi3 + ((new_sum < lo3) ? 1u : 0u);
      sum = new_sum;
    }
    while (sum >= ARBINT_NTT_P2)
      sum -= ARBINT_NTT_P2;
    uint64_t v2 = sum;

    /*  v3 = ((r3[i] - v1) - v2 * p1) * p1p2_inv_mod_p3 mod p3.  */
    uint64_t v1_mod_p3 = v1 % ARBINT_NTT_P3;
    uint64_t diff3a = mod_sub(r3[i], v1_mod_p3, ARBINT_NTT_P3);

    /*  v2 * p1 mod p3.  */
    uint64_t p1_mod_p3 = ARBINT_NTT_P1 % ARBINT_NTT_P3;
    ntt_umul(&hi, &lo, v2, p1_mod_p3);
    uint64_t r_mod_p3 = (uint64_t) 0u - ARBINT_NTT_P3;
    ntt_umul(&hi2, &lo2, hi, r_mod_p3);
    sum = lo2 + lo;
    sum_hi = hi2 + ((sum < lo2) ? 1u : 0u);
    while (sum_hi > 0u) {
      uint64_t hi3, lo3;
      ntt_umul(&hi3, &lo3, sum_hi, r_mod_p3);
      uint64_t new_sum = lo3 + sum;
      sum_hi = hi3 + ((new_sum < lo3) ? 1u : 0u);
      sum = new_sum;
    }
    while (sum >= ARBINT_NTT_P3)
      sum -= ARBINT_NTT_P3;
    uint64_t v2p1_mod_p3 = sum;

    uint64_t diff3b = mod_sub(diff3a, v2p1_mod_p3, ARBINT_NTT_P3);

    /*  diff3b * p1p2_inv_mod_p3 mod p3.  */
    ntt_umul(&hi, &lo, diff3b, crt->p1p2_inv_mod_p3);
    ntt_umul(&hi2, &lo2, hi, r_mod_p3);
    sum = lo2 + lo;
    sum_hi = hi2 + ((sum < lo2) ? 1u : 0u);
    while (sum_hi > 0u) {
      uint64_t hi3, lo3;
      ntt_umul(&hi3, &lo3, sum_hi, r_mod_p3);
      uint64_t new_sum = lo3 + sum;
      sum_hi = hi3 + ((new_sum < lo3) ? 1u : 0u);
      sum = new_sum;
    }
    while (sum >= ARBINT_NTT_P3)
      sum -= ARBINT_NTT_P3;
    uint64_t v3 = sum;

    /*  Reconstruct: coeff = v1 + v2 * p1 + v3 * p1 * p2.
        Using 192-bit arithmetic to avoid overflow.  */
    arbint_crt_wide_t coeff;
    crt_wide_set_u64(&coeff, v1);

    /*  Add v2 * p1.  */
    arbint_crt_wide_t term2;
    ntt_umul(&term2.mid, &term2.lo, v2, crt->p1_lo);
    term2.hi = 0u;
    crt_wide_add(&coeff, &coeff, &term2);

    /*  Add v3 * p1 * p2.  */
    arbint_crt_wide_t term3;
    uint64_t v3_p1p2_lo_hi, v3_p1p2_lo_lo;
    ntt_umul(&v3_p1p2_lo_hi, &v3_p1p2_lo_lo, v3, crt->p1p2_lo);
    uint64_t v3_p1p2_mid_hi, v3_p1p2_mid_lo;
    ntt_umul(&v3_p1p2_mid_hi, &v3_p1p2_mid_lo, v3, crt->p1p2_mid);

    term3.lo = v3_p1p2_lo_lo;
    term3.mid = v3_p1p2_lo_hi + v3_p1p2_mid_lo;
    uint64_t c = (term3.mid < v3_p1p2_lo_hi) ? 1u : 0u;
    term3.hi = v3_p1p2_mid_hi + c;

    crt_wide_add(&coeff, &coeff, &term3);

    /*  Add carry from previous coefficient.  */
    crt_wide_add(&coeff, &coeff, &carry);

    /*  Extract limbs and compute new carry.  */
#if ARBINT_LIMB_BITS == 64
    dst[out_idx++] = (arbint_limb_t) crt_wide_divmod_base64(&coeff);
#elif ARBINT_LIMB_BITS == 32
    dst[out_idx++] = (arbint_limb_t) crt_wide_divmod_base32(&coeff);
    dst[out_idx++] = (arbint_limb_t) crt_wide_divmod_base32(&coeff);
#endif

    carry = coeff;
  }

  /*  Emit remaining carry digits.  */
#if ARBINT_LIMB_BITS == 64
  while (carry.lo != 0u || carry.mid != 0u || carry.hi != 0u) {
    dst[out_idx++] = (arbint_limb_t) crt_wide_divmod_base64(&carry);
  }
#elif ARBINT_LIMB_BITS == 32
  while (carry.lo != 0u || carry.mid != 0u || carry.hi != 0u) {
    dst[out_idx++] = (arbint_limb_t) crt_wide_divmod_base32(&carry);
  }
#endif

  *dst_used = out_idx;
}

/* ========== Main NTT Multiplication ========== */

/*  Compute next power of 2 >= n.  */
static inline size_t next_pow2(size_t n) {
  if (n == 0u)
    return 1u;
  n--;
  n |= n >> 1;
  n |= n >> 2;
  n |= n >> 4;
  n |= n >> 8;
  n |= n >> 16;
#if SIZE_MAX > 0xFFFFFFFFu
  n |= n >> 32;
#endif
  return n + 1u;
}

/*  Compute log2 of a power of 2.  */
static inline size_t log2_of_pow2(size_t n) {
  size_t log2_n = 0u;
  while ((n >> log2_n) > 1u)
    ++log2_n;
  return log2_n;
}

arbint_err_t arbint_mul_mag_ntt_generic(arbint_limb_t * dst, size_t * out_used,
                                        const arbint_limb_t * a, size_t an,
                                        const arbint_limb_t * b, size_t bn,
                                        const arbint_alloc_t * alloc) {
  arbint_err_t rc;
  uint64_t * work = NULL;

  if (dst == NULL || out_used == NULL || a == NULL || b == NULL)
    return ARBINT_EINVAL;

  if (an == 0u || bn == 0u) {
    dst[0] = 0u;
    *out_used = 0u;
    return ARBINT_OK;
  }

  /*  Ensure NTT context is initialized.  */
  rc = arbint_ntt_ensure_init();
  if (rc != ARBINT_OK)
    return rc;

  /*  Determine NTT size: next power of 2 >= an + bn - 1.  */
  size_t conv_len = an + bn - 1u;
  size_t ntt_size = next_pow2(conv_len);
  size_t log2_n = log2_of_pow2(ntt_size);

  if (log2_n > ARBINT_NTT_MAX_LOG2)
    return ARBINT_EOVERFLOW;

  /*  Ensure root tables are large enough.  */
  rc = ensure_roots(log2_n, alloc);
  if (rc != ARBINT_OK)
    return rc;

  /*  Allocate workspace.
      Layout: work_a[ntt_size], work_b[ntt_size],
              result_0[ntt_size], result_1[ntt_size], result_2[ntt_size],
              crt_tmp[conv_len + 3] (extra space for CRT carry overflow).  */
  size_t crt_tmp_size = conv_len + 3u;
  size_t work_size = 5u * ntt_size + crt_tmp_size;
  work = (uint64_t *) alloc->realloc(alloc->ud, NULL,
                                     work_size * sizeof(uint64_t));
  if (work == NULL)
    return ARBINT_ENOMEM;

  uint64_t * work_a = work;
  uint64_t * work_b = work_a + ntt_size;
  uint64_t * result[3];
  result[0] = work_b + ntt_size;
  result[1] = result[0] + ntt_size;
  result[2] = result[1] + ntt_size;
  arbint_limb_t * crt_tmp = (arbint_limb_t *) (result[2] + ntt_size);

  /*  Process each prime.  */
  for (int pi = 0; pi < 3; ++pi) {
    const arbint_ntt_prime_t * prime = &g_ntt_ctx.primes[pi];
    const arbint_mont_params_t * mont = &prime->mont;
    const arbint_ntt_roots_t * roots = &g_ntt_ctx.roots[pi];

    /*  Convert operand a to Montgomery form and zero-pad.  */
    for (size_t i = 0u; i < an; ++i) {
#if ARBINT_LIMB_BITS == 64
      work_a[i] = to_mont(a[i], mont);
#elif ARBINT_LIMB_BITS == 32
      /*  Combine two 32-bit limbs into one 64-bit value for NTT.  */
      uint64_t limb_val = a[i];
      if (i + 1u < an) {
        /*  This packing is only valid if we adjust the convolution
            interpretation. For now, treat 32-bit as 64-bit slots
            with high 32 bits zero.  */
      }
      work_a[i] = to_mont(limb_val, mont);
#endif
    }
    memset(&work_a[an], 0, (ntt_size - an) * sizeof(uint64_t));

    /*  Convert operand b to Montgomery form and zero-pad.  */
    for (size_t i = 0u; i < bn; ++i) {
#if ARBINT_LIMB_BITS == 64
      work_b[i] = to_mont(b[i], mont);
#elif ARBINT_LIMB_BITS == 32
      work_b[i] = to_mont((uint64_t) b[i], mont);
#endif
    }
    memset(&work_b[bn], 0, (ntt_size - bn) * sizeof(uint64_t));

    /*  Forward NTT.  */
    ntt_forward(work_a, log2_n, roots, mont);
    ntt_forward(work_b, log2_n, roots, mont);

    /*  Pointwise multiplication.  */
    for (size_t i = 0u; i < ntt_size; ++i) {
      result[pi][i] = mont_mul(work_a[i], work_b[i], mont);
    }

    /*  Inverse NTT.  */
    ntt_inverse(result[pi], log2_n, roots, mont);

    /*  Scale by n^(-1) and convert from Montgomery form.  */
    uint64_t n_inv = modinv(ntt_size, prime->p);
    uint64_t n_inv_mont = to_mont(n_inv, mont);
    for (size_t i = 0u; i < ntt_size; ++i) {
      result[pi][i] = from_mont(mont_mul(result[pi][i], n_inv_mont, mont), mont);
    }

  }

  /*  CRT reconstruction into temporary buffer.
      CRT can output up to conv_len + 2 limbs due to carry propagation.
      We use crt_tmp which has conv_len + 3 limbs of space.  */
  size_t crt_used = 0u;
  ntt_crt_combine(crt_tmp, &crt_used, result[0], result[1], result[2], conv_len,
                  &g_ntt_ctx.crt);

  /*  The actual product has at most an + bn limbs. Copy to dst.  */
  size_t max_out = an + bn;
  if (crt_used > max_out) {
    /*  This shouldn't happen for correct inputs, but truncate safely.  */
    crt_used = max_out;
  }

  memcpy(dst, crt_tmp, crt_used * sizeof(arbint_limb_t));

  /*  Normalize result (trim leading zeros).  */
  *out_used = arbint_norm_used(dst, crt_used);

  /*  Cleanup.  */
  alloc->realloc(alloc->ud, work, 0);

  return ARBINT_OK;
}
