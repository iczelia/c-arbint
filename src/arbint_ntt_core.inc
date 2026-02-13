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

/*  Shared NTT core included by arbint_ntt_generic.c, arbint_ntt_bmi2.c,
    and arbint_ntt_avx2.c.

    Before including this file, the translation unit must define:
      - static inline void ntt_umul(uint64_t * hi, uint64_t * lo,
                                    uint64_t a, uint64_t b);
      - ARBINT_NTT_MUL_MAG_FN  -- name for the NTT multiplication function  */

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

/*  Reduce a 128-bit value represented as hi:lo modulo p.
    Assumes p is an odd 64-bit prime near 2^64 (NTT primes).  */
static inline uint64_t ntt_reduce_u128_mod(uint64_t hi, uint64_t lo,
                                           uint64_t p) {
  uint64_t r_mod_p = (uint64_t) 0u - p;
  uint64_t hi2, lo2;
  uint64_t sum;
  uint64_t sum_hi;

  ntt_umul(&hi2, &lo2, hi, r_mod_p);
  sum = lo2 + lo;
  sum_hi = hi2 + ((sum < lo2) ? 1u : 0u);

  while (sum_hi > 0u) {
    uint64_t hi3, lo3;
    uint64_t new_sum;
    ntt_umul(&hi3, &lo3, sum_hi, r_mod_p);
    new_sum = lo3 + sum;
    sum_hi = hi3 + ((new_sum < lo3) ? 1u : 0u);
    sum = new_sum;
  }

  while (sum >= p)
    sum -= p;

  return sum;
}

/*  Compute (a * b) mod p without requiring compiler 128-bit integer type.  */
static inline uint64_t ntt_mul_mod_u64(uint64_t a, uint64_t b, uint64_t p) {
  uint64_t hi, lo;
  ntt_umul(&hi, &lo, a, b);
  return ntt_reduce_u128_mod(hi, lo, p);
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
  x = x * (2u - p * x);     /* mod 2^2 */
  x = x * (2u - p * x);     /* mod 2^4 */
  x = x * (2u - p * x);     /* mod 2^8 */
  x = x * (2u - p * x);     /* mod 2^16 */
  x = x * (2u - p * x);     /* mod 2^32 */
  x = x * (2u - p * x);     /* mod 2^64 */
  return (uint64_t) 0u - x; /* negate to get -p^(-1) */
}

/*  Compute R^2 mod p where R = 2^64.

    Since R mod p = 2^64 mod p, and R < 2p for our primes (all near 2^64),
    we compute R mod p = 2^64 - p, then square it modularly.

    For squaring mod p without overflow, we use repeated doubling:
    R mod p, then (R mod p)^2 mod p.  */
static uint64_t compute_r_squared(uint64_t p) {
  uint64_t r_mod_p = (uint64_t) 0u - p;
  return ntt_mul_mod_u64(r_mod_p, r_mod_p, p);
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
    if (exp & 1u)
      result = ntt_mul_mod_u64(result, base, p);
    exp >>= 1;
    if (exp > 0u)
      base = ntt_mul_mod_u64(base, base, p);
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
static arbint_alloc_t g_ntt_cache_alloc;
static int g_ntt_cache_alloc_set = 0;

static const arbint_alloc_t *
arbint_ntt_get_cache_alloc(const arbint_alloc_t * alloc) {
  if (!g_ntt_cache_alloc_set && alloc != NULL && alloc->realloc != NULL) {
    g_ntt_cache_alloc = *alloc;
    g_ntt_cache_alloc_set = 1;
  }

  if (!g_ntt_cache_alloc_set || g_ntt_cache_alloc.realloc == NULL)
    return NULL;

  return &g_ntt_cache_alloc;
}

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
    g_ntt_ctx.roots[i].omega_full = NULL;
    g_ntt_ctx.roots[i].omega_inv_full = NULL;
    g_ntt_ctx.roots[i].full_max_log2 = 0u;
  }

  /*  Compute CRT constants.  */
  g_ntt_ctx.crt.p1_inv_mod_p2 = modinv(ARBINT_NTT_P1, ARBINT_NTT_P2);

  /*  Compute p1 * p2 mod p3, then invert.  */
  uint64_t p1_mod_p3 = ARBINT_NTT_P1 % ARBINT_NTT_P3;
  uint64_t p2_mod_p3 = ARBINT_NTT_P2 % ARBINT_NTT_P3;
  uint64_t p1p2_mod_p3 = ntt_mul_mod_u64(p1_mod_p3, p2_mod_p3, ARBINT_NTT_P3);
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

#ifndef ARBINT_NTT_CACHE_CLEAR_FN
  #error "Including translation unit must define ARBINT_NTT_CACHE_CLEAR_FN"
#endif

void ARBINT_NTT_CACHE_CLEAR_FN(void) {
  const arbint_alloc_t * alloc = arbint_ntt_get_cache_alloc(NULL);
  size_t i;
  size_t s;

  if (alloc != NULL) {
    for (i = 0u; i < 3u; ++i) {
      arbint_ntt_roots_t * r = &g_ntt_ctx.roots[i];

      if (r->omega_full != NULL) {
        for (s = 0u; s < ARBINT_NTT_MAX_LOG2; ++s)
          alloc->realloc(alloc->ud, r->omega_full[s], 0u);
        alloc->realloc(alloc->ud, r->omega_full, 0u);
      }

      if (r->omega_inv_full != NULL) {
        for (s = 0u; s < ARBINT_NTT_MAX_LOG2; ++s)
          alloc->realloc(alloc->ud, r->omega_inv_full[s], 0u);
        alloc->realloc(alloc->ud, r->omega_inv_full, 0u);
      }

      alloc->realloc(alloc->ud, r->omega, 0u);
      alloc->realloc(alloc->ud, r->omega_inv, 0u);

      r->omega = NULL;
      r->omega_inv = NULL;
      r->omega_full = NULL;
      r->omega_inv_full = NULL;
      r->max_log2 = 0u;
      r->full_max_log2 = 0u;
    }
  }

  memset(&g_ntt_ctx, 0, sizeof(g_ntt_ctx));
  memset(&g_ntt_cache_alloc, 0, sizeof(g_ntt_cache_alloc));
  g_ntt_cache_alloc_set = 0;
  g_ntt_initialized = 0;
}

/* ========== Root Table Management ========== */

/*  Ensure root tables are allocated for at least the given log2 size.  */
static arbint_err_t ensure_roots(size_t log2_n, const arbint_alloc_t * alloc) {
  const arbint_alloc_t * cache_alloc;

  if (log2_n > ARBINT_NTT_MAX_LOG2)
    return ARBINT_EOVERFLOW;

  cache_alloc = arbint_ntt_get_cache_alloc(alloc);
  if (cache_alloc == NULL)
    return ARBINT_EINVAL;

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
    uint64_t * new_omega = (uint64_t *) cache_alloc->realloc(
        cache_alloc->ud, r->omega, new_size * sizeof(uint64_t));
    if (new_omega == NULL)
      return ARBINT_ENOMEM;

    uint64_t * new_omega_inv = (uint64_t *) cache_alloc->realloc(
        cache_alloc->ud, r->omega_inv, new_size * sizeof(uint64_t));
    if (new_omega_inv == NULL) {
      cache_alloc->realloc(cache_alloc->ud, new_omega, 0);
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

/*  Minimum stage size for which full omega tables are beneficial.
    Below this, the overhead of table lookup isn't worth it.
    Stage s has m_half = 2^s elements, so s >= 4 means >= 16 elements.  */
#ifndef NTT_FULL_OMEGA_MIN_STAGE
  #define NTT_FULL_OMEGA_MIN_STAGE 4u
#endif

/*  Ensure full omega tables are allocated for stages [min_stage, log2_n).
    This trades O(n) memory for eliminating omega accumulation in butterflies.

    For stage s, we store omega_m^j for j = 0, 1, ..., 2^s - 1
    where omega_m = omega[s] is the primitive 2^(s+1)-th root of unity.  */
#ifdef NTT_USE_FULL_OMEGA_TABLES
static arbint_err_t ensure_full_roots(size_t log2_n,
                                      const arbint_alloc_t * alloc) {
  const arbint_alloc_t * cache_alloc;

  if (log2_n > ARBINT_NTT_MAX_LOG2)
    return ARBINT_EOVERFLOW;

  /*  First ensure compact roots exist.  */
  arbint_err_t rc = ensure_roots(log2_n, alloc);
  if (rc != ARBINT_OK)
    return rc;
  cache_alloc = arbint_ntt_get_cache_alloc(alloc);
  if (cache_alloc == NULL)
    return ARBINT_EINVAL;

  for (int i = 0; i < 3; ++i) {
    arbint_ntt_roots_t * r = &g_ntt_ctx.roots[i];
    const arbint_mont_params_t * m = &g_ntt_ctx.primes[i].mont;

    if (r->full_max_log2 >= log2_n)
      continue;

    /*  Allocate array of pointers for full tables if needed.  */
    if (r->omega_full == NULL) {
      r->omega_full =
          (uint64_t **) cache_alloc->realloc(
              cache_alloc->ud, NULL,
              ARBINT_NTT_MAX_LOG2 * sizeof(uint64_t *));
      if (r->omega_full == NULL)
        return ARBINT_ENOMEM;
      for (size_t s = 0u; s < ARBINT_NTT_MAX_LOG2; ++s)
        r->omega_full[s] = NULL;
    }
    if (r->omega_inv_full == NULL) {
      r->omega_inv_full =
          (uint64_t **) cache_alloc->realloc(
              cache_alloc->ud, NULL,
              ARBINT_NTT_MAX_LOG2 * sizeof(uint64_t *));
      if (r->omega_inv_full == NULL)
        return ARBINT_ENOMEM;
      for (size_t s = 0u; s < ARBINT_NTT_MAX_LOG2; ++s)
        r->omega_inv_full[s] = NULL;
    }

    /*  Allocate and fill full tables for new stages.  */
    for (size_t s = r->full_max_log2; s < log2_n; ++s) {
      /*  Skip small stages where full tables aren't beneficial.  */
      if (s < NTT_FULL_OMEGA_MIN_STAGE)
        continue;

      size_t table_size = (size_t) 1u << s;  /* 2^s entries */

      /*  Allocate omega_full[s].  */
      if (r->omega_full[s] == NULL) {
        r->omega_full[s] =
            (uint64_t *) cache_alloc->realloc(
                cache_alloc->ud, NULL, table_size * sizeof(uint64_t));
        if (r->omega_full[s] == NULL)
          return ARBINT_ENOMEM;

        /*  Fill with omega_m^0, omega_m^1, ..., omega_m^(2^s - 1).  */
        uint64_t omega_m = r->omega[s];
        r->omega_full[s][0] = to_mont(1u, m);  /* omega_m^0 = 1 */
        for (size_t j = 1u; j < table_size; ++j) {
          r->omega_full[s][j] = mont_mul(r->omega_full[s][j - 1], omega_m, m);
        }
      }

      /*  Allocate omega_inv_full[s].  */
      if (r->omega_inv_full[s] == NULL) {
        r->omega_inv_full[s] =
            (uint64_t *) cache_alloc->realloc(
                cache_alloc->ud, NULL, table_size * sizeof(uint64_t));
        if (r->omega_inv_full[s] == NULL)
          return ARBINT_ENOMEM;

        /*  Fill with omega_m_inv^0, omega_m_inv^1, ..., omega_m_inv^(2^s - 1).  */
        uint64_t omega_m_inv = r->omega_inv[s];
        r->omega_inv_full[s][0] = to_mont(1u, m);  /* omega_m_inv^0 = 1 */
        for (size_t j = 1u; j < table_size; ++j) {
          r->omega_inv_full[s][j] =
              mont_mul(r->omega_inv_full[s][j - 1], omega_m_inv, m);
        }
      }
    }

    r->full_max_log2 = (unsigned) log2_n;
  }

  return ARBINT_OK;
}
#endif /* NTT_USE_FULL_OMEGA_TABLES */

/* ========== NTT Transforms ========== */

/*  If ARBINT_NTT_CUSTOM_TRANSFORMS is defined, the including file provides
    its own ntt_forward and ntt_inverse implementations (e.g., AVX2-optimized).
    Otherwise, use the generic scalar implementations below.  */
#ifndef ARBINT_NTT_CUSTOM_TRANSFORMS

/*  Forward NTT: Gentleman-Sande decimation-in-frequency.
    x must be in Montgomery form on entry; remains in Montgomery form on exit.
    Input is in natural order and output is in bit-reversed order.  */
static void ntt_forward(uint64_t * x, size_t log2_n,
                        const arbint_ntt_roots_t * roots,
                        const arbint_mont_params_t * mont) {
  size_t n = (size_t) 1u << log2_n;
  uint64_t p = mont->p;

  /*  Gentleman-Sande butterfly iterations.  */
  for (size_t s = log2_n; s > 0u; --s) {
    size_t m = (size_t) 1u << s;
    size_t m_half = m >> 1;
    size_t stage_idx = s - 1u; /* omega index for this stage */

    /*  Get primitive m-th root of unity.  */
    uint64_t omega_m = roots->omega[stage_idx];

    for (size_t k = 0u; k < n; k += m) {
      uint64_t omega = to_mont(1u, mont); /* omega^0 = 1 in Montgomery form */
      for (size_t j = 0u; j < m_half; ++j) {
        uint64_t u = x[k + j];
        uint64_t v = x[k + j + m_half];
        x[k + j] = mod_add(u, v, p);
        x[k + j + m_half] = mont_mul(mod_sub(u, v, p), omega, mont);
        omega = mont_mul(omega, omega_m, mont);
      }
    }
  }
}

/*  Inverse NTT: Cooley-Tukey decimation-in-time.
    x must be in Montgomery form on entry.
    Input is in bit-reversed order and output is in natural order.
    On exit, x contains INTT(x) * n (not yet scaled by n^-1).  */
static void ntt_inverse(uint64_t * x, size_t log2_n,
                        const arbint_ntt_roots_t * roots,
                        const arbint_mont_params_t * mont) {
  size_t n = (size_t) 1u << log2_n;
  uint64_t p = mont->p;

  /*  Cooley-Tukey butterfly iterations.  */
  for (size_t s = 0u; s < log2_n; ++s) {
    size_t m = (size_t) 1u << (s + 1u);
    size_t m_half = m >> 1;

    /*  Get inverse of primitive m-th root of unity (2^(s+1)-th root).  */
    uint64_t omega_m_inv = roots->omega_inv[s];

    for (size_t k = 0u; k < n; k += m) {
      uint64_t omega_inv = to_mont(1u, mont);
      for (size_t j = 0u; j < m_half; ++j) {
        uint64_t t = mont_mul(omega_inv, x[k + j + m_half], mont);
        uint64_t u = x[k + j];
        x[k + j] = mod_add(u, t, p);
        x[k + j + m_half] = mod_sub(u, t, p);
        omega_inv = mont_mul(omega_inv, omega_m_inv, mont);
      }
    }
  }
}

#endif /* ARBINT_NTT_CUSTOM_TRANSFORMS */

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
    uint64_t v2 = ntt_mul_mod_u64(diff2, crt->p1_inv_mod_p2, ARBINT_NTT_P2);

    /*  v3 = ((r3[i] - v1) - v2 * p1) * p1p2_inv_mod_p3 mod p3.  */
    uint64_t v1_mod_p3 = v1 % ARBINT_NTT_P3;
    uint64_t diff3a = mod_sub(r3[i], v1_mod_p3, ARBINT_NTT_P3);

    /*  v2 * p1 mod p3.  */
    uint64_t p1_mod_p3 = ARBINT_NTT_P1 % ARBINT_NTT_P3;
    uint64_t v2p1_mod_p3 = ntt_mul_mod_u64(v2, p1_mod_p3, ARBINT_NTT_P3);

    uint64_t diff3b = mod_sub(diff3a, v2p1_mod_p3, ARBINT_NTT_P3);

    /*  diff3b * p1p2_inv_mod_p3 mod p3.  */
    uint64_t v3 =
        ntt_mul_mod_u64(diff3b, crt->p1p2_inv_mod_p3, ARBINT_NTT_P3);

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
    /*  In the 32-bit configuration each convolution coefficient maps to
        one base-2^32 output limb.  */
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

#ifndef ARBINT_NTT_SQR_MAG_FN
  #error "Including translation unit must define ARBINT_NTT_SQR_MAG_FN"
#endif

arbint_err_t ARBINT_NTT_MUL_MAG_FN(arbint_limb_t * dst, size_t * out_used,
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

#ifdef NTT_USE_FULL_OMEGA_TABLES
  /*  Ensure full omega tables are allocated for AVX2 table-based butterflies.  */
  rc = ensure_full_roots(log2_n, alloc);
  if (rc != ARBINT_OK)
    return rc;
#endif

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

    /*  Convert operand a to Montgomery form and zero-pad.
        If NTT_TO_MONT_LOOP is defined, use custom implementation.  */
#ifdef NTT_TO_MONT_LOOP
    NTT_TO_MONT_LOOP(work_a, a, an, mont);
#else
    for (size_t i = 0u; i < an; ++i) {
#if ARBINT_LIMB_BITS == 64
      work_a[i] = to_mont(a[i], mont);
#elif ARBINT_LIMB_BITS == 32
      work_a[i] = to_mont((uint64_t) a[i], mont);
#endif
    }
#endif
    memset(&work_a[an], 0, (ntt_size - an) * sizeof(uint64_t));

    /*  Convert operand b to Montgomery form and zero-pad.  */
#ifdef NTT_TO_MONT_LOOP
    NTT_TO_MONT_LOOP(work_b, b, bn, mont);
#else
    for (size_t i = 0u; i < bn; ++i) {
#if ARBINT_LIMB_BITS == 64
      work_b[i] = to_mont(b[i], mont);
#elif ARBINT_LIMB_BITS == 32
      work_b[i] = to_mont((uint64_t) b[i], mont);
#endif
    }
#endif
    memset(&work_b[bn], 0, (ntt_size - bn) * sizeof(uint64_t));

    /*  Forward NTT.  */
    ntt_forward(work_a, log2_n, roots, mont);
    ntt_forward(work_b, log2_n, roots, mont);

    /*  Pointwise multiplication.
        If NTT_POINTWISE_MUL is defined, use custom implementation.  */
#ifdef NTT_POINTWISE_MUL
    NTT_POINTWISE_MUL(result[pi], work_a, work_b, ntt_size, mont);
#else
    for (size_t i = 0u; i < ntt_size; ++i) {
      result[pi][i] = mont_mul(work_a[i], work_b[i], mont);
    }
#endif

    /*  Inverse NTT.  */
    ntt_inverse(result[pi], log2_n, roots, mont);

    /*  Scale by n^(-1) and convert from Montgomery form.
        If NTT_SCALE_FROM_MONT is defined, use custom implementation.  */
    uint64_t n_inv = modinv(ntt_size, prime->p);
    uint64_t n_inv_mont = to_mont(n_inv, mont);
#ifdef NTT_SCALE_FROM_MONT
    NTT_SCALE_FROM_MONT(result[pi], ntt_size, n_inv_mont, mont);
#else
    for (size_t i = 0u; i < ntt_size; ++i) {
      result[pi][i] =
          from_mont(mont_mul(result[pi][i], n_inv_mont, mont), mont);
    }
#endif
  }

  /*  CRT reconstruction into temporary buffer.
      CRT can output up to conv_len + 2 limbs due to carry propagation.
      We use crt_tmp which has conv_len + 3 limbs of space.  */
  size_t crt_used = 0u;
  ntt_crt_combine(crt_tmp, &crt_used, result[0], result[1], result[2],
                  conv_len, &g_ntt_ctx.crt);

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

/*  NTT-based squaring: dst = a^2 for large operands.

    Squaring optimization over multiplication:
    1. Only one forward NTT (vs two for multiplication) - saves ~33% of work
    2. Pointwise squaring instead of multiplication

    The convolution length for squaring is 2*an - 1 (vs an + bn - 1 for mul).
    Transform size is next power of 2 >= 2*an - 1.  */
arbint_err_t ARBINT_NTT_SQR_MAG_FN(arbint_limb_t * dst, size_t * out_used,
                                   const arbint_limb_t * a, size_t an,
                                   const arbint_alloc_t * alloc) {
  arbint_err_t rc;
  uint64_t * work = NULL;

  if (dst == NULL || out_used == NULL || a == NULL)
    return ARBINT_EINVAL;

  if (an == 0u) {
    dst[0] = 0u;
    *out_used = 0u;
    return ARBINT_OK;
  }

  /*  Ensure NTT context is initialized.  */
  rc = arbint_ntt_ensure_init();
  if (rc != ARBINT_OK)
    return rc;

  /*  Determine NTT size: next power of 2 >= 2*an - 1.
      For squaring, conv_len = 2*an - 1.  */
  size_t conv_len = 2u * an - 1u;
  size_t ntt_size = next_pow2(conv_len);
  size_t log2_n = log2_of_pow2(ntt_size);

  if (log2_n > ARBINT_NTT_MAX_LOG2)
    return ARBINT_EOVERFLOW;

  /*  Ensure root tables are large enough.  */
  rc = ensure_roots(log2_n, alloc);
  if (rc != ARBINT_OK)
    return rc;

#ifdef NTT_USE_FULL_OMEGA_TABLES
  /*  Ensure full omega tables are allocated for AVX2 table-based butterflies.  */
  rc = ensure_full_roots(log2_n, alloc);
  if (rc != ARBINT_OK)
    return rc;
#endif

  /*  Allocate workspace.
      Layout: work_a[ntt_size],
              result_0[ntt_size], result_1[ntt_size], result_2[ntt_size],
              crt_tmp[conv_len + 3] (extra space for CRT carry overflow).

      Note: no work_b needed since we square work_a in place.  */
  size_t crt_tmp_size = conv_len + 3u;
  size_t work_size = 4u * ntt_size + crt_tmp_size;
  work = (uint64_t *) alloc->realloc(alloc->ud, NULL,
                                     work_size * sizeof(uint64_t));
  if (work == NULL)
    return ARBINT_ENOMEM;

  uint64_t * work_a = work;
  uint64_t * result[3];
  result[0] = work_a + ntt_size;
  result[1] = result[0] + ntt_size;
  result[2] = result[1] + ntt_size;
  arbint_limb_t * crt_tmp = (arbint_limb_t *) (result[2] + ntt_size);

  /*  Process each prime.  */
  for (int pi = 0; pi < 3; ++pi) {
    const arbint_ntt_prime_t * prime = &g_ntt_ctx.primes[pi];
    const arbint_mont_params_t * mont = &prime->mont;
    const arbint_ntt_roots_t * roots = &g_ntt_ctx.roots[pi];

    /*  Convert operand a to Montgomery form and zero-pad.  */
#ifdef NTT_TO_MONT_LOOP
    NTT_TO_MONT_LOOP(work_a, a, an, mont);
#else
    for (size_t i = 0u; i < an; ++i) {
#if ARBINT_LIMB_BITS == 64
      work_a[i] = to_mont(a[i], mont);
#elif ARBINT_LIMB_BITS == 32
      work_a[i] = to_mont((uint64_t) a[i], mont);
#endif
    }
#endif
    memset(&work_a[an], 0, (ntt_size - an) * sizeof(uint64_t));

    /*  Forward NTT - only one transform needed for squaring.  */
    ntt_forward(work_a, log2_n, roots, mont);

    /*  Pointwise squaring (vs multiplication for two distinct operands).
        If NTT_POINTWISE_SQR is defined, use custom implementation.  */
#ifdef NTT_POINTWISE_SQR
    NTT_POINTWISE_SQR(result[pi], work_a, ntt_size, mont);
#else
    for (size_t i = 0u; i < ntt_size; ++i) {
      result[pi][i] = mont_mul(work_a[i], work_a[i], mont);
    }
#endif

    /*  Inverse NTT.  */
    ntt_inverse(result[pi], log2_n, roots, mont);

    /*  Scale by n^(-1) and convert from Montgomery form.  */
    uint64_t n_inv = modinv(ntt_size, prime->p);
    uint64_t n_inv_mont = to_mont(n_inv, mont);
#ifdef NTT_SCALE_FROM_MONT
    NTT_SCALE_FROM_MONT(result[pi], ntt_size, n_inv_mont, mont);
#else
    for (size_t i = 0u; i < ntt_size; ++i) {
      result[pi][i] =
          from_mont(mont_mul(result[pi][i], n_inv_mont, mont), mont);
    }
#endif
  }

  /*  CRT reconstruction into temporary buffer.  */
  size_t crt_used = 0u;
  ntt_crt_combine(crt_tmp, &crt_used, result[0], result[1], result[2],
                  conv_len, &g_ntt_ctx.crt);

  /*  The actual product has at most 2*an limbs. Copy to dst.  */
  size_t max_out = 2u * an;
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
