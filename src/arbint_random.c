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

#include "arbint_random.h"
#include <string.h>

/* ========== MT19937 Core Algorithm ========== */

/*  Initialize MT19937 with a 32-bit seed.
    Parameters:
      rng - MT19937 state to initialize
      seed - 32-bit seed value
    Algorithm:
      Uses standard MT19937 initialization formula with the multiplier
      1812433253.  */
static void arbint_mt19937_init_u32(arbint_rng_t * rng, uint32_t seed) {
  unsigned int i;
  rng->mt[0] = seed;
  for (i = 1; i < ARBINT_MT_N; i++) {
    rng->mt[i] =
        (1812433253UL * (rng->mt[i - 1] ^ (rng->mt[i - 1] >> 30)) + i);
  }
  rng->mti = ARBINT_MT_N;
}

/*  Initialize MT19937 with an array of seeds.
    Allows full state initialization from arbitrary-length seed array.
    This provides better entropy distribution than single-word initialization.
    Parameters:
      rng - MT19937 state to initialize
      init_key - Array of uint32_t seeds
      key_length - Number of uint32_t elements in init_key
    Algorithm:
      Two-phase mixing:
      1. Forward mixing: Combine init_key into state with multiplier 1664525
      2. Backward mixing: Further mix state with multiplier 1566083941  */
static void arbint_mt19937_init_array(arbint_rng_t * rng,
                                      const uint32_t * init_key,
                                      size_t key_length) {
  unsigned int i, j;
  size_t k;

  arbint_mt19937_init_u32(rng, 19650218UL);

  i = 1;
  j = 0;
  k = (ARBINT_MT_N > key_length) ? ARBINT_MT_N : key_length;

  for (; k != 0; k--) {
    rng->mt[i] = (rng->mt[i] ^
                  ((rng->mt[i - 1] ^ (rng->mt[i - 1] >> 30)) * 1664525UL)) +
                 init_key[j] + j;
    i++;
    j++;
    if (i >= ARBINT_MT_N) {
      rng->mt[0] = rng->mt[ARBINT_MT_N - 1];
      i = 1;
    }
    if (j >= key_length)
      j = 0;
  }

  for (k = ARBINT_MT_N - 1; k != 0; k--) {
    rng->mt[i] = (rng->mt[i] ^
                  ((rng->mt[i - 1] ^ (rng->mt[i - 1] >> 30)) * 1566083941UL)) -
                 i;
    i++;
    if (i >= ARBINT_MT_N) {
      rng->mt[0] = rng->mt[ARBINT_MT_N - 1];
      i = 1;
    }
  }

  rng->mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/*  Generate next 32-bit random value from MT19937.
    Implements the standard tempered output transformation.
    Parameters:
      rng - MT19937 state
    Returns:
      32-bit uniformly distributed random value
    Complexity: Amortized O(1)  */
static uint32_t arbint_mt19937_u32(arbint_rng_t * rng) {
  uint32_t y;
  static const uint32_t mag01[2] = {0x0UL, ARBINT_MT_MATRIX_A};

  /* Generate N words at once */
  if (rng->mti >= ARBINT_MT_N) {
    int kk;

    for (kk = 0; kk < ARBINT_MT_N - ARBINT_MT_M; kk++) {
      y = (rng->mt[kk] & ARBINT_MT_UPPER_MASK) |
          (rng->mt[kk + 1] & ARBINT_MT_LOWER_MASK);
      rng->mt[kk] = rng->mt[kk + ARBINT_MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for (; kk < ARBINT_MT_N - 1; kk++) {
      y = (rng->mt[kk] & ARBINT_MT_UPPER_MASK) |
          (rng->mt[kk + 1] & ARBINT_MT_LOWER_MASK);
      rng->mt[kk] = rng->mt[kk + (ARBINT_MT_M - ARBINT_MT_N)] ^ (y >> 1) ^
                    mag01[y & 0x1UL];
    }
    y = (rng->mt[ARBINT_MT_N - 1] & ARBINT_MT_UPPER_MASK) |
        (rng->mt[0] & ARBINT_MT_LOWER_MASK);
    rng->mt[ARBINT_MT_N - 1] =
        rng->mt[ARBINT_MT_M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

    rng->mti = 0;
  }

  y = rng->mt[rng->mti++];

  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  return y;
}

/* ========== Platform Entropy Selection ========== */

/*  Select platform entropy source at runtime.
    Returns NULL if no platform source is available.
    Platform priority (in order):
      1. Windows: BCryptGenRandom or CryptGenRandom
      2. macOS/BSD: arc4random_buf
      3. Linux/POSIX: /dev/urandom  */
arbint_entropy_fn_t arbint_select_entropy_source(void) {
#if defined(ARBINT_HAS_WINAPI_ENTROPY)
  return arbint_entropy_winapi;
#elif defined(ARBINT_HAS_ARC4RANDOM)
  return arbint_entropy_arc4;
#elif defined(ARBINT_HAS_URANDOM)
  return arbint_entropy_urandom;
#else
  return NULL;
#endif
}

/* ========== Public API ========== */

/*  Initialize RNG state from platform entropy or provided seed.
    If seed is NULL, reads 2496 bytes from platform entropy source and
    uses it to fully initialize the MT19937 state array.
    If seed is non-NULL, interprets it as an array of uint32_t values
    and uses them to initialize MT19937.
    Parameters:
      rng - MT19937 state to initialize
      seed - Seed buffer (NULL for platform entropy) (may alias nothing)
      seed_len - Seed buffer length in bytes (must be multiple of 4 if
                 non-zero).
    Returns: ARBINT_OK on success ARBINT_EINVAL if rng is NULL, or
    seed_len invalid, or no platform entropy ARBINT_EDOM if platform
    entropy source fails.  */
arbint_err_t arbint_rng_init(arbint_rng_t * rng, const void * seed,
                             size_t seed_len) {
  if (rng == NULL)
    return ARBINT_EINVAL;

  if (seed == NULL) {
    /* Use platform entropy source */
    arbint_entropy_fn_t entropy_fn;
    uint32_t entropy_buf[ARBINT_MT_N];
    int rc;
    entropy_fn = arbint_select_entropy_source();
    if (entropy_fn == NULL)
      return ARBINT_EINVAL; /* No platform entropy source */
    rc = entropy_fn((uint8_t *) entropy_buf, sizeof(entropy_buf));
    if (rc != 0)
      return ARBINT_EDOM; /* Entropy source failed */
    arbint_mt19937_init_array(rng, entropy_buf, ARBINT_MT_N);
    /* Zero entropy buffer for security */
    arbint_secure_zero(entropy_buf, sizeof(entropy_buf));
  } else {
    /* Use provided seed, copying to aligned buffer to avoid UB from
       potentially misaligned caller-supplied pointer.  */
    size_t key_len;
    uint32_t seed_buf[ARBINT_MT_N];
    if (seed_len == 0 || (seed_len % sizeof(uint32_t)) != 0)
      return ARBINT_EINVAL;
    key_len = seed_len / sizeof(uint32_t);
    if (key_len > ARBINT_MT_N)
      key_len = ARBINT_MT_N;
    memcpy(seed_buf, seed, key_len * sizeof(uint32_t));
    arbint_mt19937_init_array(rng, seed_buf, key_len);
    arbint_secure_zero(seed_buf, key_len * sizeof(uint32_t));
  }

  return ARBINT_OK;
}

/*  Clear RNG state (zeros memory for security).
    This prevents state from leaking in memory dumps.
    Parameters:
      rng - MT19937 state to clear (may be NULL)  */
void arbint_rng_clear(arbint_rng_t * rng) {
  if (rng != NULL)
    arbint_secure_zero(rng, sizeof(arbint_rng_t));
}

/*  Generate uniform random integer in [0, 2^k).
    Fills rop with exactly k random bits from MT19937.
    Parameters:
      rop - Destination (may alias nothing)
      rng - MT19937 state
      k - Number of random bits to generate
    Returns:
      ARBINT_OK on success
      ARBINT_EINVAL if rop or rng is NULL
      ARBINT_ENOMEM if allocation fails
    Algorithm:
      1. Calculate limbs_needed = ceil(k / LIMB_BITS)
      2. Calculate bytes_needed = ceil(k / 8)
      3. Resize rop to hold limbs_needed
      4. Generate random bytes using MT19937
      5. Mask high bits if k not byte-aligned
      6. Normalize and set _sz
    Complexity: O(k)  */
arbint_err_t arbint_urandomb(arbint_t rop, arbint_rng_t * rng, size_t k) {
  size_t limbs_needed, i;
  size_t words_per_limb;
  arbint_limb_t * rp;
  arbint_err_t rc;

  if (rop == NULL || rng == NULL)
    return ARBINT_EINVAL;

  /* k == 0 -> result is 0 */
  if (k == 0u) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  /* Calculate capacity needed */
  limbs_needed = (k + ARBINT_LIMB_BITS - 1u) / ARBINT_LIMB_BITS;

  /* Resize rop */
  rc = arbint_resize(rop, limbs_needed);
  if (rc != ARBINT_OK)
    return rc;

  /* Fill limbs directly from MT19937 32-bit words (endian-safe). */
  rp = ARBINT_LIMBS(rop);
  words_per_limb = ARBINT_LIMB_BITS / 32u;

  for (i = 0; i < limbs_needed; i++) {
    arbint_limb_t limb = 0;
    size_t w;
    for (w = 0; w < words_per_limb; w++)
      limb |= (arbint_limb_t) arbint_mt19937_u32(rng) << (w * 32u);
    rp[i] = limb;
  }

  /* Mask unused high bits in the top limb */
  {
    size_t top_bits = k % ARBINT_LIMB_BITS;
    if (top_bits != 0u) {
      arbint_limb_t mask = ((arbint_limb_t) 1u << top_bits) - 1u;
      rp[limbs_needed - 1u] &= mask;
    }
  }

  /* Normalize and set _sz */
  limbs_needed = arbint_norm_used(rp, limbs_needed);
  rop[0]._sz = (ptrdiff_t) limbs_needed;

  return ARBINT_OK;
}

/*  Generate uniform random integer in [0, bound) via rejection sampling.
    Uses rejection sampling to ensure perfectly uniform distribution with
    no modulo bias.
    Parameters:
      rop - Destination (may alias bound)
      rng - MT19937 state
      bound - Upper bound (exclusive) (may alias rop)
    Returns:
      ARBINT_OK on success
      ARBINT_EINVAL if rop, rng, or bound is NULL
      ARBINT_EDOM if bound <= 0, or if rejection sampling fails after 256 tries
      ARBINT_ENOMEM if allocation fails
    Algorithm:
      1. Get bit count of bound
      2. Loop (max 256 times):
         a. Generate random value with bound_bits random bits
         b. If value < bound, accept and return
      3. If exceeded retry limit, return ARBINT_EDOM
    Expected iterations: ~1.5 (worst case ~2 for bound just above power of 2)
    Complexity: O(k) expected, where k = ceil(log2(bound))  */
arbint_err_t arbint_urandomm(arbint_t rop, arbint_rng_t * rng,
                             const arbint_t bound) {
  size_t bound_bits;
  int retry;
  arbint_err_t rc;
  arbint_t bound_copy;
  int need_copy;

  if (rop == NULL || rng == NULL || bound == NULL)
    return ARBINT_EINVAL;

  if (arbint_signum(bound) <= 0)
    return ARBINT_EDOM;

  /* bound == 1 → always return 0 */
  if (arbint_cmp_u32(bound, 1u) == 0) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  /* Handle aliasing: if rop == bound, we need to copy bound */
  need_copy = (rop == bound);
  if (need_copy) {
    rc = arbint_init(bound_copy, bound[0]._ctx);
    if (rc != ARBINT_OK)
      return rc;
    rc = arbint_set(bound_copy, bound);
    if (rc != ARBINT_OK) {
      arbint_clear(bound_copy);
      return rc;
    }
  }

  /* Get bit count of bound */
  bound_bits = arbint_nbits(need_copy ? bound_copy : bound);

  /* Rejection sampling loop (typically 1-2 iterations) */
  for (retry = 0; retry < 256; ++retry) {
    rc = arbint_urandomb(rop, rng, bound_bits);
    if (rc != ARBINT_OK) {
      if (need_copy)
        arbint_clear(bound_copy);
      return rc;
    }

    if (arbint_cmp(rop, need_copy ? bound_copy : bound) < 0) {
      if (need_copy)
        arbint_clear(bound_copy);
      return ARBINT_OK; /* Success: rop < bound */
    }
  }

  /* Exceeded retry limit (probability < 2^-256) */
  if (need_copy)
    arbint_clear(bound_copy);
  return ARBINT_EDOM;
}
