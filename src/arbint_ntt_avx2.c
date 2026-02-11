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

/*  AVX2-optimized NTT multiplication.

    This file provides AVX2-vectorized butterfly operations for the NTT.
    Key optimizations:
    1. Vectorized mod_add/mod_sub: process 4 elements in parallel using
       AVX2 256-bit operations with unsigned comparison via XOR trick.
    2. Vectorized butterfly memory access: _mm256_loadu_si256 for loading
       4 elements at once.
    3. Uses BMI2 _mulx_u64 for Montgomery multiplication (AVX2 implies BMI2
       on all modern CPUs: Intel Haswell, AMD Zen and later).
    4. Precomputed full omega tables: for stages >= 4, we precompute all
       omega^j values and use direct table lookups instead of computing
       omega values on the fly. This eliminates the omega accumulation
       chain which was a bottleneck in the inner loop.
    5. AVX2-optimized loops for to_mont, pointwise multiply, and from_mont
       conversions with prefetching.

    The Montgomery multiplication itself cannot be vectorized with AVX2
    because it requires 64x64->128 bit multiplication which AVX2 lacks.
    We use scalar BMI2 multiplies but vectorize the surrounding add/sub.

    Future optimization potential:
    - Radix-4 butterflies: combine two adjacent stages into a single pass,
      reducing loop overhead and improving cache locality. Would require
      restructuring the stage loop to process pairs of stages.
    - AVX-512: 8 elements per vector, native 64-bit unsigned compare.  */

#include "arbint_ntt.h"
#include "config.h"

#include <immintrin.h>
#include <stdlib.h>
#include <string.h>

/* ========== Wide Multiplication ========== */

/*  Multiply two 64-bit integers producing full 128-bit result (hi:lo = a * b).

    Uses BMI2 _mulx_u64 intrinsic. AVX2 CPUs always have BMI2 support, so
    this is safe to use unconditionally in AVX2-compiled code.  */
static inline void ntt_umul(uint64_t * hi, uint64_t * lo, uint64_t a,
                            uint64_t b) {
  unsigned long long hi64 = 0ull;
  unsigned long long lo64 =
      _mulx_u64((unsigned long long) a, (unsigned long long) b, &hi64);
  *lo = (uint64_t) lo64;
  *hi = (uint64_t) hi64;
}

/*  Tell the core include that we provide custom transforms.  */
#define ARBINT_NTT_CUSTOM_TRANSFORMS 1

/*  Enable precomputed full omega tables for faster butterfly operations.
    This trades O(n) memory for eliminating the omega accumulation chain.  */
#define NTT_USE_FULL_OMEGA_TABLES 1

/*  Forward declarations for custom NTT transforms (defined after the include).
    The core include calls these functions but we define them below.  */
static void ntt_forward(uint64_t * x, size_t log2_n,
                        const arbint_ntt_roots_t * roots,
                        const arbint_mont_params_t * mont);
static void ntt_inverse(uint64_t * x, size_t log2_n,
                        const arbint_ntt_roots_t * roots,
                        const arbint_mont_params_t * mont);

/* ========== AVX2 Optimized Loop Macros ========== */

/*  Threshold for using AVX2 loops. Below this, scalar loops have less overhead.
    16 elements = 4 AVX2 vectors = 128 bytes = 2 cache lines.  */
#define NTT_AVX2_LOOP_THRESHOLD 16u

/*  AVX2-optimized to_mont conversion loop.
    Uses vectorized loads/stores with scalar Montgomery multiply.
    The to_mont function is defined in the include, so we inline the loop body.  */
#define NTT_TO_MONT_LOOP(dst, src, n, mont)                                    \
  do {                                                                         \
    size_t _i = 0u;                                                            \
    if ((n) >= NTT_AVX2_LOOP_THRESHOLD) {                                      \
      for (; _i + 4u <= (n); _i += 4u) {                                       \
        if (_i + 16u < (n))                                                    \
          _mm_prefetch((const char *) &(src)[_i + 16u], _MM_HINT_T0);          \
        __m256i _v = _mm256_loadu_si256((const __m256i *) &(src)[_i]);         \
        uint64_t _s0 = (uint64_t) _mm256_extract_epi64(_v, 0);                 \
        uint64_t _s1 = (uint64_t) _mm256_extract_epi64(_v, 1);                 \
        uint64_t _s2 = (uint64_t) _mm256_extract_epi64(_v, 2);                 \
        uint64_t _s3 = (uint64_t) _mm256_extract_epi64(_v, 3);                 \
        uint64_t _r0 = to_mont(_s0, mont);                                     \
        uint64_t _r1 = to_mont(_s1, mont);                                     \
        uint64_t _r2 = to_mont(_s2, mont);                                     \
        uint64_t _r3 = to_mont(_s3, mont);                                     \
        __m256i _result = _mm256_set_epi64x(                                   \
            (long long) _r3, (long long) _r2, (long long) _r1, (long long) _r0); \
        _mm256_storeu_si256((__m256i *) &(dst)[_i], _result);                  \
      }                                                                        \
    }                                                                          \
    for (; _i < (n); ++_i) {                                                   \
      (dst)[_i] = to_mont((src)[_i], mont);                                    \
    }                                                                          \
  } while (0)

/*  AVX2-optimized pointwise multiplication loop.  */
#define NTT_POINTWISE_MUL(dst, a, b, n, mont)                                  \
  do {                                                                         \
    size_t _i = 0u;                                                            \
    if ((n) >= NTT_AVX2_LOOP_THRESHOLD) {                                      \
      for (; _i + 4u <= (n); _i += 4u) {                                       \
        if (_i + 16u < (n)) {                                                  \
          _mm_prefetch((const char *) &(a)[_i + 16u], _MM_HINT_T0);            \
          _mm_prefetch((const char *) &(b)[_i + 16u], _MM_HINT_T0);            \
        }                                                                      \
        __m256i _va = _mm256_loadu_si256((const __m256i *) &(a)[_i]);          \
        __m256i _vb = _mm256_loadu_si256((const __m256i *) &(b)[_i]);          \
        uint64_t _a0 = (uint64_t) _mm256_extract_epi64(_va, 0);                \
        uint64_t _a1 = (uint64_t) _mm256_extract_epi64(_va, 1);                \
        uint64_t _a2 = (uint64_t) _mm256_extract_epi64(_va, 2);                \
        uint64_t _a3 = (uint64_t) _mm256_extract_epi64(_va, 3);                \
        uint64_t _b0 = (uint64_t) _mm256_extract_epi64(_vb, 0);                \
        uint64_t _b1 = (uint64_t) _mm256_extract_epi64(_vb, 1);                \
        uint64_t _b2 = (uint64_t) _mm256_extract_epi64(_vb, 2);                \
        uint64_t _b3 = (uint64_t) _mm256_extract_epi64(_vb, 3);                \
        uint64_t _r0 = mont_mul(_a0, _b0, mont);                               \
        uint64_t _r1 = mont_mul(_a1, _b1, mont);                               \
        uint64_t _r2 = mont_mul(_a2, _b2, mont);                               \
        uint64_t _r3 = mont_mul(_a3, _b3, mont);                               \
        __m256i _result = _mm256_set_epi64x(                                   \
            (long long) _r3, (long long) _r2, (long long) _r1, (long long) _r0); \
        _mm256_storeu_si256((__m256i *) &(dst)[_i], _result);                  \
      }                                                                        \
    }                                                                          \
    for (; _i < (n); ++_i) {                                                   \
      (dst)[_i] = mont_mul((a)[_i], (b)[_i], mont);                            \
    }                                                                          \
  } while (0)

/*  AVX2-optimized scale and from_mont conversion loop.  */
#define NTT_SCALE_FROM_MONT(x, n, scale, mont)                                 \
  do {                                                                         \
    size_t _i = 0u;                                                            \
    if ((n) >= NTT_AVX2_LOOP_THRESHOLD) {                                      \
      for (; _i + 4u <= (n); _i += 4u) {                                       \
        if (_i + 16u < (n))                                                    \
          _mm_prefetch((const char *) &(x)[_i + 16u], _MM_HINT_T0);            \
        __m256i _v = _mm256_loadu_si256((const __m256i *) &(x)[_i]);           \
        uint64_t _v0 = (uint64_t) _mm256_extract_epi64(_v, 0);                 \
        uint64_t _v1 = (uint64_t) _mm256_extract_epi64(_v, 1);                 \
        uint64_t _v2 = (uint64_t) _mm256_extract_epi64(_v, 2);                 \
        uint64_t _v3 = (uint64_t) _mm256_extract_epi64(_v, 3);                 \
        uint64_t _r0 = from_mont(mont_mul(_v0, scale, mont), mont);            \
        uint64_t _r1 = from_mont(mont_mul(_v1, scale, mont), mont);            \
        uint64_t _r2 = from_mont(mont_mul(_v2, scale, mont), mont);            \
        uint64_t _r3 = from_mont(mont_mul(_v3, scale, mont), mont);            \
        __m256i _result = _mm256_set_epi64x(                                   \
            (long long) _r3, (long long) _r2, (long long) _r1, (long long) _r0); \
        _mm256_storeu_si256((__m256i *) &(x)[_i], _result);                    \
      }                                                                        \
    }                                                                          \
    for (; _i < (n); ++_i) {                                                   \
      (x)[_i] = from_mont(mont_mul((x)[_i], scale, mont), mont);               \
    }                                                                          \
  } while (0)

/*  Include the core for Montgomery arithmetic, CRT, etc.
    We need the helper functions before defining our custom transforms.  */
#define ARBINT_NTT_MUL_MAG_FN arbint_mul_mag_ntt_avx2
#include "arbint_ntt_core.inc"

/* ========== AVX2 Vectorized Modular Arithmetic ========== */

/*  Threshold for using AVX2 butterflies. When m_half < this value, use
    scalar code to avoid SIMD setup overhead.  */
#define NTT_AVX2_BUTTERFLY_THRESHOLD 4u

/*  Vectorized modular addition: (a + b) mod p for 4 elements.

    For primes p close to 2^64, we cannot simply compute a + b because it
    may overflow. Instead we use: result = (a >= p - b) ? (a - (p - b)) : (a + b)

    AVX2 has no unsigned 64-bit compare, so we use the XOR-with-sign-bit trick:
    unsigned a >= b iff signed (a ^ 0x8000...) >= (b ^ 0x8000...)  */
static inline __m256i ntt_mod_add_avx2(__m256i a, __m256i b, __m256i p) {
  /*  p_minus_b = p - b  */
  __m256i p_minus_b = _mm256_sub_epi64(p, b);

  /*  sum = a + b (wraps on overflow, but we select the correct result)  */
  __m256i sum = _mm256_add_epi64(a, b);

  /*  diff = a - (p - b) = a + b - p  */
  __m256i diff = _mm256_sub_epi64(a, p_minus_b);

  /*  Unsigned compare a >= p_minus_b using XOR trick.  */
  __m256i sign_bit = _mm256_set1_epi64x((long long) 0x8000000000000000ULL);
  __m256i a_signed = _mm256_xor_si256(a, sign_bit);
  __m256i pmb_signed = _mm256_xor_si256(p_minus_b, sign_bit);

  /*  cmp = (a_signed > pmb_signed) ? 0xFFFF... : 0  */
  __m256i cmp_gt = _mm256_cmpgt_epi64(a_signed, pmb_signed);

  /*  eq = (a_signed == pmb_signed) ? 0xFFFF... : 0  */
  __m256i cmp_eq = _mm256_cmpeq_epi64(a_signed, pmb_signed);

  /*  ge = (a >= p_minus_b) = (a > p_minus_b) || (a == p_minus_b)  */
  __m256i ge = _mm256_or_si256(cmp_gt, cmp_eq);

  /*  Select: ge ? diff : sum  */
  return _mm256_blendv_epi8(sum, diff, ge);
}

/*  Vectorized modular subtraction: (a - b) mod p for 4 elements.

    If a >= b, result = a - b.
    If a < b, result = a - b + p (add p to correct the underflow).  */
static inline __m256i ntt_mod_sub_avx2(__m256i a, __m256i b, __m256i p) {
  /*  diff = a - b (wraps on underflow)  */
  __m256i diff = _mm256_sub_epi64(a, b);

  /*  Unsigned compare a < b using XOR trick.  */
  __m256i sign_bit = _mm256_set1_epi64x((long long) 0x8000000000000000ULL);
  __m256i a_signed = _mm256_xor_si256(a, sign_bit);
  __m256i b_signed = _mm256_xor_si256(b, sign_bit);

  /*  borrow = (a < b) ? 0xFFFF... : 0
      a < b iff signed(a) < signed(b) iff signed(b) > signed(a)  */
  __m256i borrow = _mm256_cmpgt_epi64(b_signed, a_signed);

  /*  correction = borrow ? p : 0  */
  __m256i correction = _mm256_and_si256(borrow, p);

  /*  result = diff + correction  */
  return _mm256_add_epi64(diff, correction);
}

/* ========== AVX2 Butterfly Helpers ========== */

/*  Process 4 forward butterflies in parallel.

    For each j in [0, 4):
      t[j] = omega[j] * x[k + j + m_half]  (Montgomery multiply)
      u[j] = x[k + j]
      x[k + j] = (u[j] + t[j]) mod p
      x[k + j + m_half] = (u[j] - t[j]) mod p

    omega_vec contains [omega^0, omega^1, omega^2, omega^3] in Montgomery form.
    Returns omega^4 for continuing the sequence.  */
static inline uint64_t ntt_butterfly_forward_4(
    uint64_t * x, size_t k, size_t m_half, __m256i omega_vec,
    uint64_t omega_m, const arbint_mont_params_t * mont) {
  uint64_t p = mont->p;
  __m256i p_vec = _mm256_set1_epi64x((long long) p);

  /*  Load u = x[k:k+4] and v = x[k+m_half:k+m_half+4]  */
  __m256i u = _mm256_loadu_si256((const __m256i *) &x[k]);
  __m256i v = _mm256_loadu_si256((const __m256i *) &x[k + m_half]);

  /*  Extract omega values for scalar Montgomery multiplication.  */
  uint64_t omega0 = (uint64_t) _mm256_extract_epi64(omega_vec, 0);
  uint64_t omega1 = (uint64_t) _mm256_extract_epi64(omega_vec, 1);
  uint64_t omega2 = (uint64_t) _mm256_extract_epi64(omega_vec, 2);
  uint64_t omega3 = (uint64_t) _mm256_extract_epi64(omega_vec, 3);

  /*  Extract v values for scalar Montgomery multiplication.  */
  uint64_t v0 = (uint64_t) _mm256_extract_epi64(v, 0);
  uint64_t v1 = (uint64_t) _mm256_extract_epi64(v, 1);
  uint64_t v2 = (uint64_t) _mm256_extract_epi64(v, 2);
  uint64_t v3 = (uint64_t) _mm256_extract_epi64(v, 3);

  /*  Compute t[j] = omega[j] * v[j] using scalar Montgomery multiply.
      This is the bottleneck - AVX2 cannot do 64x64->128 multiply.  */
  uint64_t hi, lo;

  ntt_umul(&hi, &lo, omega0, v0);
  uint64_t t0 = mont_redc(hi, lo, mont);

  ntt_umul(&hi, &lo, omega1, v1);
  uint64_t t1 = mont_redc(hi, lo, mont);

  ntt_umul(&hi, &lo, omega2, v2);
  uint64_t t2 = mont_redc(hi, lo, mont);

  ntt_umul(&hi, &lo, omega3, v3);
  uint64_t t3 = mont_redc(hi, lo, mont);

  __m256i t = _mm256_set_epi64x((long long) t3, (long long) t2, (long long) t1,
                                (long long) t0);

  /*  Vectorized mod_add and mod_sub.  */
  __m256i result_add = ntt_mod_add_avx2(u, t, p_vec);
  __m256i result_sub = ntt_mod_sub_avx2(u, t, p_vec);

  /*  Store results.  */
  _mm256_storeu_si256((__m256i *) &x[k], result_add);
  _mm256_storeu_si256((__m256i *) &x[k + m_half], result_sub);

  /*  Compute omega^4 for the next iteration.
      omega^4 = omega^3 * omega_m  */
  ntt_umul(&hi, &lo, omega3, omega_m);
  return mont_redc(hi, lo, mont);
}

/*  Process 4 inverse butterflies in parallel.

    For each j in [0, 4):
      u[j] = x[k + j]
      v[j] = x[k + j + m_half]
      x[k + j] = (u[j] + v[j]) mod p
      x[k + j + m_half] = (u[j] - v[j]) * omega_inv[j] mod p

    omega_inv_vec contains [omega_inv^0, omega_inv^1, omega_inv^2, omega_inv^3].
    Returns omega_inv^4 for continuing the sequence.  */
static inline uint64_t ntt_butterfly_inverse_4(
    uint64_t * x, size_t k, size_t m_half, __m256i omega_inv_vec,
    uint64_t omega_m_inv, const arbint_mont_params_t * mont) {
  uint64_t p = mont->p;
  __m256i p_vec = _mm256_set1_epi64x((long long) p);

  /*  Load u = x[k:k+4] and v = x[k+m_half:k+m_half+4]  */
  __m256i u = _mm256_loadu_si256((const __m256i *) &x[k]);
  __m256i v = _mm256_loadu_si256((const __m256i *) &x[k + m_half]);

  /*  Vectorized mod_add for x[k + j] = (u + v) mod p.  */
  __m256i result_add = ntt_mod_add_avx2(u, v, p_vec);

  /*  Vectorized mod_sub for diff = (u - v) mod p.  */
  __m256i diff = ntt_mod_sub_avx2(u, v, p_vec);

  /*  Store the addition result immediately.  */
  _mm256_storeu_si256((__m256i *) &x[k], result_add);

  /*  Extract values for scalar Montgomery multiplication.  */
  uint64_t omega_inv0 = (uint64_t) _mm256_extract_epi64(omega_inv_vec, 0);
  uint64_t omega_inv1 = (uint64_t) _mm256_extract_epi64(omega_inv_vec, 1);
  uint64_t omega_inv2 = (uint64_t) _mm256_extract_epi64(omega_inv_vec, 2);
  uint64_t omega_inv3 = (uint64_t) _mm256_extract_epi64(omega_inv_vec, 3);

  uint64_t d0 = (uint64_t) _mm256_extract_epi64(diff, 0);
  uint64_t d1 = (uint64_t) _mm256_extract_epi64(diff, 1);
  uint64_t d2 = (uint64_t) _mm256_extract_epi64(diff, 2);
  uint64_t d3 = (uint64_t) _mm256_extract_epi64(diff, 3);

  /*  Compute (u - v) * omega_inv using scalar Montgomery multiply.  */
  uint64_t hi, lo;

  ntt_umul(&hi, &lo, d0, omega_inv0);
  uint64_t r0 = mont_redc(hi, lo, mont);

  ntt_umul(&hi, &lo, d1, omega_inv1);
  uint64_t r1 = mont_redc(hi, lo, mont);

  ntt_umul(&hi, &lo, d2, omega_inv2);
  uint64_t r2 = mont_redc(hi, lo, mont);

  ntt_umul(&hi, &lo, d3, omega_inv3);
  uint64_t r3 = mont_redc(hi, lo, mont);

  __m256i result_mul =
      _mm256_set_epi64x((long long) r3, (long long) r2, (long long) r1,
                        (long long) r0);

  /*  Store the multiplication result.  */
  _mm256_storeu_si256((__m256i *) &x[k + m_half], result_mul);

  /*  Compute omega_inv^4 for the next iteration.  */
  ntt_umul(&hi, &lo, omega_inv3, omega_m_inv);
  return mont_redc(hi, lo, mont);
}

/*  Prepare omega vector: [omega^0, omega^1, omega^2, omega^3] in Montgomery form.
    omega_start is the starting omega value (omega^j for some j).  */
static inline __m256i ntt_prepare_omega_vec(uint64_t omega_start,
                                            uint64_t omega_m,
                                            const arbint_mont_params_t * mont) {
  uint64_t hi, lo;
  uint64_t o0 = omega_start;

  ntt_umul(&hi, &lo, o0, omega_m);
  uint64_t o1 = mont_redc(hi, lo, mont);

  ntt_umul(&hi, &lo, o1, omega_m);
  uint64_t o2 = mont_redc(hi, lo, mont);

  ntt_umul(&hi, &lo, o2, omega_m);
  uint64_t o3 = mont_redc(hi, lo, mont);

  return _mm256_set_epi64x((long long) o3, (long long) o2, (long long) o1,
                           (long long) o0);
}

/*  Process 4 forward butterflies using precomputed omega table.
    Same as ntt_butterfly_forward_4 but loads omega values from table[j:j+4].  */
static inline void ntt_butterfly_forward_4_table(
    uint64_t * x, size_t k, size_t m_half, const uint64_t * omega_table,
    size_t j, const arbint_mont_params_t * mont) {
  uint64_t p = mont->p;
  __m256i p_vec = _mm256_set1_epi64x((long long) p);

  /*  Load u = x[k+j:k+j+4] and v = x[k+j+m_half:k+j+m_half+4]  */
  __m256i u = _mm256_loadu_si256((const __m256i *) &x[k + j]);
  __m256i v = _mm256_loadu_si256((const __m256i *) &x[k + j + m_half]);

  /*  Load omega values directly from precomputed table.  */
  __m256i omega_vec = _mm256_loadu_si256((const __m256i *) &omega_table[j]);

  /*  Extract omega values for scalar Montgomery multiplication.  */
  uint64_t omega0 = (uint64_t) _mm256_extract_epi64(omega_vec, 0);
  uint64_t omega1 = (uint64_t) _mm256_extract_epi64(omega_vec, 1);
  uint64_t omega2 = (uint64_t) _mm256_extract_epi64(omega_vec, 2);
  uint64_t omega3 = (uint64_t) _mm256_extract_epi64(omega_vec, 3);

  /*  Extract v values for scalar Montgomery multiplication.  */
  uint64_t v0 = (uint64_t) _mm256_extract_epi64(v, 0);
  uint64_t v1 = (uint64_t) _mm256_extract_epi64(v, 1);
  uint64_t v2 = (uint64_t) _mm256_extract_epi64(v, 2);
  uint64_t v3 = (uint64_t) _mm256_extract_epi64(v, 3);

  /*  Compute t[j] = omega[j] * v[j] using scalar Montgomery multiply.  */
  uint64_t hi, lo;

  ntt_umul(&hi, &lo, omega0, v0);
  uint64_t t0 = mont_redc(hi, lo, mont);

  ntt_umul(&hi, &lo, omega1, v1);
  uint64_t t1 = mont_redc(hi, lo, mont);

  ntt_umul(&hi, &lo, omega2, v2);
  uint64_t t2 = mont_redc(hi, lo, mont);

  ntt_umul(&hi, &lo, omega3, v3);
  uint64_t t3 = mont_redc(hi, lo, mont);

  __m256i t = _mm256_set_epi64x((long long) t3, (long long) t2, (long long) t1,
                                (long long) t0);

  /*  Vectorized mod_add and mod_sub.  */
  __m256i result_add = ntt_mod_add_avx2(u, t, p_vec);
  __m256i result_sub = ntt_mod_sub_avx2(u, t, p_vec);

  /*  Store results.  */
  _mm256_storeu_si256((__m256i *) &x[k + j], result_add);
  _mm256_storeu_si256((__m256i *) &x[k + j + m_half], result_sub);
}

/*  Process 4 inverse butterflies using precomputed omega_inv table.
    Same as ntt_butterfly_inverse_4 but loads omega_inv values from table[j:j+4].  */
static inline void ntt_butterfly_inverse_4_table(
    uint64_t * x, size_t k, size_t m_half, const uint64_t * omega_inv_table,
    size_t j, const arbint_mont_params_t * mont) {
  uint64_t p = mont->p;
  __m256i p_vec = _mm256_set1_epi64x((long long) p);

  /*  Load u = x[k+j:k+j+4] and v = x[k+j+m_half:k+j+m_half+4]  */
  __m256i u = _mm256_loadu_si256((const __m256i *) &x[k + j]);
  __m256i v = _mm256_loadu_si256((const __m256i *) &x[k + j + m_half]);

  /*  Vectorized mod_add for x[k + j] = (u + v) mod p.  */
  __m256i result_add = ntt_mod_add_avx2(u, v, p_vec);

  /*  Vectorized mod_sub for diff = (u - v) mod p.  */
  __m256i diff = ntt_mod_sub_avx2(u, v, p_vec);

  /*  Store the addition result immediately.  */
  _mm256_storeu_si256((__m256i *) &x[k + j], result_add);

  /*  Load omega_inv values directly from precomputed table.  */
  __m256i omega_inv_vec =
      _mm256_loadu_si256((const __m256i *) &omega_inv_table[j]);

  /*  Extract values for scalar Montgomery multiplication.  */
  uint64_t omega_inv0 = (uint64_t) _mm256_extract_epi64(omega_inv_vec, 0);
  uint64_t omega_inv1 = (uint64_t) _mm256_extract_epi64(omega_inv_vec, 1);
  uint64_t omega_inv2 = (uint64_t) _mm256_extract_epi64(omega_inv_vec, 2);
  uint64_t omega_inv3 = (uint64_t) _mm256_extract_epi64(omega_inv_vec, 3);

  uint64_t d0 = (uint64_t) _mm256_extract_epi64(diff, 0);
  uint64_t d1 = (uint64_t) _mm256_extract_epi64(diff, 1);
  uint64_t d2 = (uint64_t) _mm256_extract_epi64(diff, 2);
  uint64_t d3 = (uint64_t) _mm256_extract_epi64(diff, 3);

  /*  Compute (u - v) * omega_inv using scalar Montgomery multiply.  */
  uint64_t hi, lo;

  ntt_umul(&hi, &lo, d0, omega_inv0);
  uint64_t r0 = mont_redc(hi, lo, mont);

  ntt_umul(&hi, &lo, d1, omega_inv1);
  uint64_t r1 = mont_redc(hi, lo, mont);

  ntt_umul(&hi, &lo, d2, omega_inv2);
  uint64_t r2 = mont_redc(hi, lo, mont);

  ntt_umul(&hi, &lo, d3, omega_inv3);
  uint64_t r3 = mont_redc(hi, lo, mont);

  __m256i result_mul =
      _mm256_set_epi64x((long long) r3, (long long) r2, (long long) r1,
                        (long long) r0);

  /*  Store the multiplication result.  */
  _mm256_storeu_si256((__m256i *) &x[k + j + m_half], result_mul);
}

/* ========== AVX2 NTT Transforms ========== */

/*  Forward NTT with AVX2 optimization: Cooley-Tukey decimation-in-time.
    x must be in Montgomery form on entry; remains in Montgomery form on exit.

    When precomputed full omega tables are available (omega_full[s] != NULL),
    we use direct table lookups instead of computing omega values on the fly.
    This eliminates the omega accumulation chain in the inner loop.  */
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
    uint64_t omega_m = roots->omega[s];

    /*  Check if precomputed full omega table is available for this stage.  */
    int use_full_table = (roots->omega_full != NULL &&
                          s >= NTT_FULL_OMEGA_MIN_STAGE &&
                          s < roots->full_max_log2 &&
                          roots->omega_full[s] != NULL);

    if (m_half >= NTT_AVX2_BUTTERFLY_THRESHOLD) {
      /*  AVX2 path: process 4 butterflies at a time.  */
      for (size_t k = 0u; k < n; k += m) {
        size_t j = 0u;

        if (use_full_table) {
          /*  Use precomputed omega table for direct lookups.  */
          const uint64_t * omega_table = roots->omega_full[s];

          /*  Process groups of 4 using table lookups.  */
          for (; j + 4u <= m_half; j += 4u) {
            ntt_butterfly_forward_4_table(x, k, m_half, omega_table, j, mont);
          }

          /*  Handle remaining elements (0-3) with scalar code using table.  */
          for (; j < m_half; ++j) {
            uint64_t omega = omega_table[j];
            uint64_t t = mont_mul(omega, x[k + j + m_half], mont);
            uint64_t u = x[k + j];
            x[k + j] = mod_add(u, t, p);
            x[k + j + m_half] = mod_sub(u, t, p);
          }
        } else {
          /*  Compute omega values on the fly.  */
          uint64_t omega = to_mont(1u, mont);

          /*  Process groups of 4.  */
          for (; j + 4u <= m_half; j += 4u) {
            __m256i omega_vec = ntt_prepare_omega_vec(omega, omega_m, mont);
            omega = ntt_butterfly_forward_4(x, k + j, m_half, omega_vec, omega_m,
                                            mont);
          }

          /*  Handle remaining elements (0-3) with scalar code.  */
          for (; j < m_half; ++j) {
            uint64_t t = mont_mul(omega, x[k + j + m_half], mont);
            uint64_t u = x[k + j];
            x[k + j] = mod_add(u, t, p);
            x[k + j + m_half] = mod_sub(u, t, p);
            omega = mont_mul(omega, omega_m, mont);
          }
        }
      }
    } else {
      /*  Scalar path for small m_half.  */
      for (size_t k = 0u; k < n; k += m) {
        uint64_t omega = to_mont(1u, mont);
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
}

/*  Inverse NTT with AVX2 optimization: Gentleman-Sande decimation-in-frequency.
    x must be in Montgomery form on entry.
    On exit, x contains INTT(x) * n (not yet scaled by n^-1).

    When precomputed full omega_inv tables are available (omega_inv_full[s] != NULL),
    we use direct table lookups instead of computing omega_inv values on the fly.  */
static void ntt_inverse(uint64_t * x, size_t log2_n,
                        const arbint_ntt_roots_t * roots,
                        const arbint_mont_params_t * mont) {
  size_t n = (size_t) 1u << log2_n;
  uint64_t p = mont->p;

  /*  Gentleman-Sande butterfly iterations (reverse order of forward).  */
  for (size_t s = log2_n; s > 0u; --s) {
    size_t m = (size_t) 1u << s;
    size_t m_half = m >> 1;
    size_t stage_idx = s - 1u;  /* omega_inv index for this stage */
    uint64_t omega_m_inv = roots->omega_inv[stage_idx];

    /*  Check if precomputed full omega_inv table is available for this stage.  */
    int use_full_table = (roots->omega_inv_full != NULL &&
                          stage_idx >= NTT_FULL_OMEGA_MIN_STAGE &&
                          stage_idx < roots->full_max_log2 &&
                          roots->omega_inv_full[stage_idx] != NULL);

    if (m_half >= NTT_AVX2_BUTTERFLY_THRESHOLD) {
      /*  AVX2 path: process 4 butterflies at a time.  */
      for (size_t k = 0u; k < n; k += m) {
        size_t j = 0u;

        if (use_full_table) {
          /*  Use precomputed omega_inv table for direct lookups.  */
          const uint64_t * omega_inv_table = roots->omega_inv_full[stage_idx];

          /*  Process groups of 4 using table lookups.  */
          for (; j + 4u <= m_half; j += 4u) {
            ntt_butterfly_inverse_4_table(x, k, m_half, omega_inv_table, j, mont);
          }

          /*  Handle remaining elements (0-3) with scalar code using table.  */
          for (; j < m_half; ++j) {
            uint64_t omega_inv = omega_inv_table[j];
            uint64_t u = x[k + j];
            uint64_t v = x[k + j + m_half];
            x[k + j] = mod_add(u, v, p);
            x[k + j + m_half] = mont_mul(mod_sub(u, v, p), omega_inv, mont);
          }
        } else {
          /*  Compute omega_inv values on the fly.  */
          uint64_t omega_inv = to_mont(1u, mont);

          /*  Process groups of 4.  */
          for (; j + 4u <= m_half; j += 4u) {
            __m256i omega_inv_vec =
                ntt_prepare_omega_vec(omega_inv, omega_m_inv, mont);
            omega_inv = ntt_butterfly_inverse_4(x, k + j, m_half, omega_inv_vec,
                                                omega_m_inv, mont);
          }

          /*  Handle remaining elements (0-3) with scalar code.  */
          for (; j < m_half; ++j) {
            uint64_t u = x[k + j];
            uint64_t v = x[k + j + m_half];
            x[k + j] = mod_add(u, v, p);
            x[k + j + m_half] = mont_mul(mod_sub(u, v, p), omega_inv, mont);
            omega_inv = mont_mul(omega_inv, omega_m_inv, mont);
          }
        }
      }
    } else {
      /*  Scalar path for small m_half.  */
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
