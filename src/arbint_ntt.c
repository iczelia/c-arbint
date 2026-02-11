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
#include "arbint_cpu.h"
#include "config.h"

/*  Function pointer type for NTT multiplication dispatch.  */
typedef arbint_err_t (*arbint_mul_mag_ntt_fn_t)(
    arbint_limb_t * dst, size_t * out_used, const arbint_limb_t * a, size_t an,
    const arbint_limb_t * b, size_t bn, const arbint_alloc_t * alloc);

/*  Select optimal NTT implementation based on CPU features.
    Three-tier dispatch: compile-time always, runtime detection, fallback.  */
static arbint_mul_mag_ntt_fn_t arbint_select_mul_mag_ntt(void) {
#if HAS_AVX2_ALWAYS
  return arbint_mul_mag_ntt_avx2;
#elif HAS_AVX2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_AVX2)
             ? arbint_mul_mag_ntt_avx2
             : arbint_mul_mag_ntt_generic;
#else
  return arbint_mul_mag_ntt_generic;
#endif
}

/*  Cached function pointer for NTT multiplication.
    Lazily initialized on first use.  */
static arbint_mul_mag_ntt_fn_t g_mul_mag_ntt = NULL;

/*  Public dispatch function for NTT multiplication.
    Used by arbint_mul_mag_rec when operands exceed NTT threshold.  */
arbint_err_t arbint_mul_mag_ntt(arbint_limb_t * dst, size_t * out_used,
                                const arbint_limb_t * a, size_t an,
                                const arbint_limb_t * b, size_t bn,
                                const arbint_alloc_t * alloc) {
  if (g_mul_mag_ntt == NULL)
    g_mul_mag_ntt = arbint_select_mul_mag_ntt();

  return g_mul_mag_ntt(dst, out_used, a, an, b, bn, alloc);
}
