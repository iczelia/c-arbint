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
#include "arbint_ntt_segmented.h"
#include "arbint_cpu.h"
#include "arbint_internal_util.h"
#include "config.h"

/*  Function pointer type for NTT multiplication dispatch.  */
typedef arbint_err_t (*arbint_mul_mag_ntt_fn_t)(
    arbint_limb_t * dst, size_t * out_used, const arbint_limb_t * a, size_t an,
    const arbint_limb_t * b, size_t bn, const arbint_alloc_t * alloc);

/*  Function pointer type for NTT squaring dispatch.  */
typedef arbint_err_t (*arbint_sqr_mag_ntt_fn_t)(
    arbint_limb_t * dst, size_t * out_used, const arbint_limb_t * a, size_t an,
    const arbint_alloc_t * alloc);

/*  Select optimal NTT implementation based on CPU features.
    Three-tier dispatch: compile-time always, runtime detection, fallback.
    Priority: AVX2 > BMI2 > generic.  */
static arbint_mul_mag_ntt_fn_t arbint_select_mul_mag_ntt(void) {
  /*  AVX2 tier (highest priority - AVX2 implies BMI2).  */
#if HAS_AVX2_ALWAYS
  return arbint_mul_mag_ntt_avx2;
#elif HAS_AVX2
  if (arbint_cpu_has_feature(ARBINT_CPU_FEATURE_AVX2))
    return arbint_mul_mag_ntt_avx2;
#endif

  /*  BMI2 tier (intermediate - _mulx_u64 optimization).  */
#if HAS_BMI2_ALWAYS
  return arbint_mul_mag_ntt_bmi2;
#elif HAS_BMI2
  if (arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2))
    return arbint_mul_mag_ntt_bmi2;
#endif

  /*  Generic fallback (portable half-limb multiplication).  */
  return arbint_mul_mag_ntt_generic;
}

/*  Select optimal NTT squaring implementation based on CPU features.
    Three-tier dispatch: compile-time always, runtime detection, fallback.
    Priority: AVX2 > BMI2 > generic.  */
static arbint_sqr_mag_ntt_fn_t arbint_select_sqr_mag_ntt(void) {
  /*  AVX2 tier (highest priority - AVX2 implies BMI2).  */
#if HAS_AVX2_ALWAYS
  return arbint_sqr_mag_ntt_avx2;
#elif HAS_AVX2
  if (arbint_cpu_has_feature(ARBINT_CPU_FEATURE_AVX2))
    return arbint_sqr_mag_ntt_avx2;
#endif

  /*  BMI2 tier (intermediate - _mulx_u64 optimization).  */
#if HAS_BMI2_ALWAYS
  return arbint_sqr_mag_ntt_bmi2;
#elif HAS_BMI2
  if (arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2))
    return arbint_sqr_mag_ntt_bmi2;
#endif

  /*  Generic fallback (portable half-limb multiplication).  */
  return arbint_sqr_mag_ntt_generic;
}

/*  Cached function pointers for NTT multiplication and squaring.
    Lazily initialized on first use.  */
static arbint_mul_mag_ntt_fn_t g_mul_mag_ntt = NULL;
static arbint_sqr_mag_ntt_fn_t g_sqr_mag_ntt = NULL;

/*  Public dispatch function for NTT multiplication.
    Used by arbint_mul_mag_rec when operands exceed NTT threshold.

    If operands exceed ARBINT_NTT_MAX_SIZE, dispatches to segmented
    multiplication which recursively splits operands until sub-products
    fit within NTT limits.  */
arbint_err_t arbint_mul_mag_ntt(arbint_limb_t * dst, size_t * out_used,
                                const arbint_limb_t * a, size_t an,
                                const arbint_limb_t * b, size_t bn,
                                const arbint_alloc_t * alloc) {
  /*  Check if operands exceed single-NTT capacity.
      NTT requires transform size >= an + bn - 1, rounded up to power of 2.
      If this exceeds ARBINT_NTT_MAX_SIZE, use segmented multiplication.  */
  size_t conv_len = an + bn;
  if (conv_len > ARBINT_NTT_MAX_SIZE)
    return arbint_mul_mag_ntt_segmented(dst, out_used, a, an, b, bn, alloc);

  ARBINT_LAZY_INIT(g_mul_mag_ntt, arbint_select_mul_mag_ntt);

  return g_mul_mag_ntt(dst, out_used, a, an, b, bn, alloc);
}

/*  Public dispatch function for NTT squaring.
    Used by arbint_sqr_mag_rec when operands exceed NTT threshold.

    Squaring uses only one forward NTT (vs two for multiplication), providing
    ~33% performance improvement over multiplication for same-size operands.

    If operand exceeds ARBINT_NTT_MAX_SIZE, dispatches to segmented squaring
    which recursively splits operand until sub-products fit within NTT limits.  */
arbint_err_t arbint_sqr_mag_ntt(arbint_limb_t * dst, size_t * out_used,
                                const arbint_limb_t * a, size_t an,
                                const arbint_alloc_t * alloc) {
  /*  Check if operand exceeds single-NTT capacity.
      NTT requires transform size >= 2*an - 1, rounded up to power of 2.
      If this exceeds ARBINT_NTT_MAX_SIZE, use segmented squaring.  */
  size_t conv_len = 2u * an;
  if (conv_len > ARBINT_NTT_MAX_SIZE)
    return arbint_sqr_mag_ntt_segmented(dst, out_used, a, an, alloc);

  ARBINT_LAZY_INIT(g_sqr_mag_ntt, arbint_select_sqr_mag_ntt);

  return g_sqr_mag_ntt(dst, out_used, a, an, alloc);
}
