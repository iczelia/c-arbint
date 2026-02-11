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
size_t arbint_mulacc_1_generic(arbint_limb_t * dst, size_t dst_n, size_t dst_cap,
                               const arbint_limb_t * a, size_t an,
                               arbint_limb_t b);

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

/*  Squaring thresholds (limb counts).
    Squaring exploits symmetry, so Karatsuba may be beneficial at smaller sizes
    than for general multiplication. Toom-3 squaring is simpler than Toom-3
    multiplication (no sign tracking at point -1).  */
#define ARBINT_SQR_KARATSUBA_THRESHOLD 10u
#define ARBINT_SQR_TOOM3_THRESHOLD 18u

#endif /* ARBINT_MUL_H */
