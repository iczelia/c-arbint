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

#ifndef ARBINT_NTT_H
#define ARBINT_NTT_H

#include "arbint_base.h"
#include "config.h"

#include <stddef.h>
#include <stdint.h>

/*  NTT multiplication threshold is defined in arbint_mul.h to allow
    centralized threshold management.  */
#ifndef ARBINT_NTT_THRESHOLD
  #define ARBINT_NTT_THRESHOLD 1024u
#endif

/*  Maximum supported NTT size as log2(size).
    2^24 = 16M limbs = 128MB per operand (64-bit limbs).
    Limited by p1's max_log2 = 32, but 24 is practical.  */
#define ARBINT_NTT_MAX_LOG2 24u
#define ARBINT_NTT_MAX_SIZE (1u << ARBINT_NTT_MAX_LOG2)

/*  NTT-friendly primes for 3-prime CRT reconstruction.

    Coefficient growth analysis (64-bit limbs, w=64):
    - Each coefficient: c[i] <= min(an,bn) * (2^64 - 1)^2
    - log2(c[i]) <= log2(n) + 128
    - At n = 2^20 (1M limbs): log2(c[i]) <= 148 bits

    Product of 3 primes must satisfy: P > 2 * n_max * (2^64 - 1)^2
    With these primes: P > 2^190 > 2^149 (safe margin).

    All primes are of form k * 2^n + 1, allowing radix-2 NTT.  */

/*  p1 = 2^64 - 2^32 + 1 (Goldilocks prime)
    Factorization of p1-1: 2^32 * 3 * 5 * 17 * 257 * 65537
    max_log2 = 32, primitive root g = 7  */
#define ARBINT_NTT_P1 UINT64_C(0xFFFFFFFF00000001)
#define ARBINT_NTT_G1 UINT64_C(7)
#define ARBINT_NTT_P1_MAX_LOG2 32u

/*  p2 = 2^64 - 2^34 + 1
    Factorization of p2-1: 2^34 * (2^30 - 1)
    2^30 - 1 = 3^2 * 7 * 11 * 31 * 151 * 331
    max_log2 = 34, primitive root g = 10  */
#define ARBINT_NTT_P2 UINT64_C(0xFFFFFFFC00000001)
#define ARBINT_NTT_G2 UINT64_C(10)
#define ARBINT_NTT_P2_MAX_LOG2 34u

/*  p3 = 2^64 - 2^24 + 1
    Factorization of p3-1: 2^24 * (2^40 - 1)
    2^40 - 1 = 3 * 5^2 * 11 * 17 * 31 * 41 * 61681
    max_log2 = 24, primitive root g = 43  */
#define ARBINT_NTT_P3 UINT64_C(0xFFFFFFFFFF000001)
#define ARBINT_NTT_G3 UINT64_C(43)
#define ARBINT_NTT_P3_MAX_LOG2 24u

/*  Montgomery parameters for a 64-bit prime p.
    R = 2^64 (implicit).  */
typedef struct arbint_mont_params {
  uint64_t p;         /* Prime modulus */
  uint64_t p_neg_inv; /* -p^(-1) mod 2^64 */
  uint64_t r_squared; /* R^2 mod p (for to-Montgomery conversion) */
} arbint_mont_params_t;

/*  NTT prime configuration including Montgomery parameters.  */
typedef struct arbint_ntt_prime {
  uint64_t p;                /* Prime modulus */
  uint64_t g;                /* Primitive root mod p */
  unsigned max_log2;         /* Max power of 2 dividing (p-1) */
  arbint_mont_params_t mont; /* Montgomery parameters */
} arbint_ntt_prime_t;

/*  Precomputed twiddle factors (roots of unity) for one prime.
    Stored in Montgomery form for fast modular multiply.  */
typedef struct arbint_ntt_roots {
  uint64_t * omega;     /* omega[k] = g^((p-1)/2^(k+1)) mod p (Montgomery) */
  uint64_t * omega_inv; /* Inverses for INTT (Montgomery) */
  unsigned max_log2;    /* Allocated table size */
} arbint_ntt_roots_t;

/*  CRT constants for 3-prime reconstruction via Garner's algorithm.  */
typedef struct arbint_ntt_crt {
  uint64_t p1_inv_mod_p2;   /* p1^(-1) mod p2 */
  uint64_t p1p2_inv_mod_p3; /* (p1*p2)^(-1) mod p3 */
  uint64_t p1_lo;           /* p1 (low 64 bits) */
  uint64_t p1p2_lo;         /* (p1*p2) mod 2^64 */
  uint64_t p1p2_mid;        /* (p1*p2) >> 64 mod 2^64 */
} arbint_ntt_crt_t;

/*  Global NTT context holding all precomputed data for the 3 primes.
    Lazy-initialized on first NTT call.  */
typedef struct arbint_ntt_ctx {
  arbint_ntt_prime_t primes[3];
  arbint_ntt_roots_t roots[3];
  arbint_ntt_crt_t crt;
  int initialized;
} arbint_ntt_ctx_t;

/*  Extended precision type for CRT reconstruction.
    Holds values up to p1*p2*p3 (~192 bits).  */
typedef struct arbint_crt_wide {
  uint64_t lo;  /* Bits 0-63 */
  uint64_t mid; /* Bits 64-127 */
  uint64_t hi;  /* Bits 128-191 */
} arbint_crt_wide_t;

/*  Core NTT magnitude multiplication.
    dst must have capacity for an + bn limbs.
    Returns result limb count via out_used.

    Preconditions:
    - an >= bn (caller swaps if needed)
    - bn >= ARBINT_NTT_THRESHOLD
    - an <= 2 * bn (balanced operands)  */
arbint_err_t arbint_mul_mag_ntt_generic(arbint_limb_t * dst, size_t * out_used,
                                        const arbint_limb_t * a, size_t an,
                                        const arbint_limb_t * b, size_t bn,
                                        const arbint_alloc_t * alloc);

#if HAS_BMI2
arbint_err_t arbint_mul_mag_ntt_bmi2(arbint_limb_t * dst, size_t * out_used,
                                     const arbint_limb_t * a, size_t an,
                                     const arbint_limb_t * b, size_t bn,
                                     const arbint_alloc_t * alloc);
#endif

#if HAS_AVX2
arbint_err_t arbint_mul_mag_ntt_avx2(arbint_limb_t * dst, size_t * out_used,
                                     const arbint_limb_t * a, size_t an,
                                     const arbint_limb_t * b, size_t bn,
                                     const arbint_alloc_t * alloc);
#endif

/*  Public dispatch function for NTT multiplication.
    Selects optimal implementation based on CPU features.  */
arbint_err_t arbint_mul_mag_ntt(arbint_limb_t * dst, size_t * out_used,
                                const arbint_limb_t * a, size_t an,
                                const arbint_limb_t * b, size_t bn,
                                const arbint_alloc_t * alloc);

#endif /* ARBINT_NTT_H */
