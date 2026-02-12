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

#ifndef ARBINT_MUL_H
#define ARBINT_MUL_H

#include "arbint_addsub.h"
#include "arbint_base.h"

/*  Shared wide multiply helpers used by mul/div/gcd/ntt kernels.
    The BMI2 variants use _mulx_u64 only when the including translation unit
    defines ARBINT_USE_BMI2_INTRIN and is compiled with BMI2 support.  */

static inline void arbint_umul_limb_generic(arbint_limb_t * hi,
                                            arbint_limb_t * lo,
                                            arbint_limb_t a,
                                            arbint_limb_t b) {
#if ARBINT_LIMB_BITS == 32
  uint64_t p = (uint64_t) a * (uint64_t) b;
  *lo = (arbint_limb_t) p;
  *hi = (arbint_limb_t) (p >> 32);
#else
  arbint_limb_t a0 = a & ARBINT_HALF_MASK;
  arbint_limb_t a1 = a >> ARBINT_HALF_BITS;
  arbint_limb_t b0 = b & ARBINT_HALF_MASK;
  arbint_limb_t b1 = b >> ARBINT_HALF_BITS;
  arbint_limb_t w0 = a0 * b0;
  arbint_limb_t t = a1 * b0 + (w0 >> ARBINT_HALF_BITS);
  arbint_limb_t w1 = t & ARBINT_HALF_MASK;
  arbint_limb_t w2 = t >> ARBINT_HALF_BITS;
  w1 = a0 * b1 + w1;
  *hi = a1 * b1 + w2 + (w1 >> ARBINT_HALF_BITS);
  *lo = (w1 << ARBINT_HALF_BITS) | (w0 & ARBINT_HALF_MASK);
#endif
}

static inline void arbint_umul_u64_generic(uint64_t * hi, uint64_t * lo,
                                           uint64_t a, uint64_t b) {
  uint64_t a0 = a & UINT64_C(0xFFFFFFFF);
  uint64_t a1 = a >> 32;
  uint64_t b0 = b & UINT64_C(0xFFFFFFFF);
  uint64_t b1 = b >> 32;
  uint64_t w0 = a0 * b0;
  uint64_t t = a1 * b0 + (w0 >> 32);
  uint64_t w1 = t & UINT64_C(0xFFFFFFFF);
  uint64_t w2 = t >> 32;
  w1 = a0 * b1 + w1;
  *hi = a1 * b1 + w2 + (w1 >> 32);
  *lo = (w1 << 32) | (w0 & UINT64_C(0xFFFFFFFF));
}

#if defined(ARBINT_USE_BMI2_INTRIN) && HAS_BMI2 && ARBINT_LIMB_BITS == 64
  #include <immintrin.h>
static inline void arbint_umul_limb_bmi2(arbint_limb_t * hi,
                                         arbint_limb_t * lo,
                                         arbint_limb_t a,
                                         arbint_limb_t b) {
  unsigned long long hi64 = 0ull;
  unsigned long long lo64 =
      _mulx_u64((unsigned long long) a, (unsigned long long) b, &hi64);
  *lo = (arbint_limb_t) lo64;
  *hi = (arbint_limb_t) hi64;
}

static inline void arbint_umul_u64_bmi2(uint64_t * hi, uint64_t * lo,
                                        uint64_t a, uint64_t b) {
  unsigned long long hi64 = 0ull;
  unsigned long long lo64 = _mulx_u64((unsigned long long) a,
                                      (unsigned long long) b, &hi64);
  *lo = (uint64_t) lo64;
  *hi = (uint64_t) hi64;
}
#else
static inline void arbint_umul_limb_bmi2(arbint_limb_t * hi,
                                         arbint_limb_t * lo,
                                         arbint_limb_t a,
                                         arbint_limb_t b) {
  arbint_umul_limb_generic(hi, lo, a, b);
}

static inline void arbint_umul_u64_bmi2(uint64_t * hi, uint64_t * lo,
                                        uint64_t a, uint64_t b) {
  arbint_umul_u64_generic(hi, lo, a, b);
}
#endif

typedef struct arbint_mul_kernel_table {
  arbint_err_t (*mul_impl)(arbint_t rop, const arbint_t a, const arbint_t b);
  arbint_err_t (*sqr_impl)(arbint_t rop, const arbint_t a);
  size_t (*mul_limb_1)(arbint_limb_t * dst, const arbint_limb_t * a, size_t an,
                       arbint_limb_t b);
  arbint_err_t (*mul_mag)(arbint_limb_t * dst, size_t * out_used,
                          const arbint_limb_t * a, size_t an,
                          const arbint_limb_t * b, size_t bn,
                          const arbint_alloc_t * alloc);
  size_t (*mulacc)(arbint_limb_t * dst, size_t dst_n, size_t dst_cap,
                   const arbint_limb_t * a, size_t an,
                   const arbint_limb_t * c, size_t cn);
  size_t (*mulacc_1)(arbint_limb_t * dst, size_t dst_n, size_t dst_cap,
                     const arbint_limb_t * a, size_t an, arbint_limb_t b);
} arbint_mul_kernel_table_t;

const arbint_mul_kernel_table_t * arbint_mul_kernel_table_get(void);

int arbint_mul_cap(size_t an, size_t bn, size_t * out);

arbint_err_t arbint_mul_impl_generic(arbint_t rop, const arbint_t a,
                                     const arbint_t b);

arbint_err_t arbint_sqr_impl_generic(arbint_t rop, const arbint_t a);

size_t arbint_mul_limb_1_generic(arbint_limb_t * dst, const arbint_limb_t * a,
                                 size_t an, arbint_limb_t b);

/*  Recursive magnitude multiplication: dst = a * b (magnitudes only).
    Uses schoolbook, Karatsuba, or Toom-3 depending on operand size.
    Returns normalized result limb count via *out_used.
    dst must have capacity for at least an + bn + 1 limbs.  */
arbint_err_t arbint_mul_mag_generic(arbint_limb_t * dst, size_t * out_used,
                                    const arbint_limb_t * a, size_t an,
                                    const arbint_limb_t * b, size_t bn,
                                    const arbint_alloc_t * alloc);

/*  Fused multiply-accumulate for single limb multiplier.
    Computes dst[0..dst_n-1] += a[0..an-1] * b in place.
    dst must have capacity for at least max(dst_n, an) + 1 limbs.
    Returns new normalized limb count.  */
size_t arbint_mulacc_1_generic(arbint_limb_t * dst, size_t dst_n,
                               size_t dst_cap, const arbint_limb_t * a,
                               size_t an, arbint_limb_t b);

/*  Fused multiply-accumulate for multi-limb multiplier.
    Computes dst[0..] += a[0..an-1] * c[0..cn-1] in place.
    dst must have capacity for at least max(dst_n, an + cn) + 1 limbs.
    Returns new normalized limb count.  */
size_t arbint_mulacc_generic(arbint_limb_t * dst, size_t dst_n, size_t dst_cap,
                             const arbint_limb_t * a, size_t an,
                             const arbint_limb_t * c, size_t cn);

#if HAS_BMI2
arbint_err_t arbint_mul_impl_bmi2(arbint_t rop, const arbint_t a,
                                  const arbint_t b);

arbint_err_t arbint_sqr_impl_bmi2(arbint_t rop, const arbint_t a);

size_t arbint_mul_limb_1_bmi2(arbint_limb_t * dst, const arbint_limb_t * a,
                              size_t an, arbint_limb_t b);

arbint_err_t arbint_mul_mag_bmi2(arbint_limb_t * dst, size_t * out_used,
                                 const arbint_limb_t * a, size_t an,
                                 const arbint_limb_t * b, size_t bn,
                                 const arbint_alloc_t * alloc);

size_t arbint_mulacc_1_bmi2(arbint_limb_t * dst, size_t dst_n, size_t dst_cap,
                            const arbint_limb_t * a, size_t an,
                            arbint_limb_t b);

size_t arbint_mulacc_bmi2(arbint_limb_t * dst, size_t dst_n, size_t dst_cap,
                          const arbint_limb_t * a, size_t an,
                          const arbint_limb_t * c, size_t cn);
#endif /* HAS_BMI2 */

/*  Multiplication thresholds (limb counts).  */
#define ARBINT_KARATSUBA_THRESHOLD 15u
#define ARBINT_TOOM3_THRESHOLD 20u

/*  NTT multiplication threshold (limb counts).
    NTT provides O(n log n) complexity vs Toom-3's O(n^1.465).
    Below this threshold, Toom-3 is faster due to NTT overhead.
    Dispatch condition: bn >= THRESH && an <= 2*bn (bn is smaller operand).  */
#define ARBINT_NTT_THRESHOLD 1024u

/*  Squaring thresholds (limb counts).
    Squaring exploits symmetry, so Karatsuba may be beneficial at smaller sizes
    than for general multiplication. Toom-3 squaring is simpler than Toom-3
    multiplication (no sign tracking at point -1).  */
#define ARBINT_SQR_KARATSUBA_THRESHOLD 10u
#define ARBINT_SQR_TOOM3_THRESHOLD 18u

#endif /* ARBINT_MUL_H */
