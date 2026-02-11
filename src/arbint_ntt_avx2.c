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

    Current implementation: thin wrapper around generic.
    Future optimization: vectorize butterfly add/sub operations while keeping
    Montgomery multiplication scalar (since AVX2 lacks efficient 64x64->128
    multiply, unless...).

    Optimization strategy (for future work):
    1. Vectorize mod_add/mod_sub: 4 parallel 64-bit additions with conditional
       subtraction using _mm256_cmpgt_epi64 or unsigned comparison trick.
    2. Vectorize butterfly memory access: _mm256_loadu_si256 for loading
       4 elements at once.
    3. Keep Montgomery multiply scalar: use BMI2 _mulx_u64 or half-limb method.
    4. Consider cache-friendly 4-step FFT for very large transforms.  */

#include "arbint_ntt.h"
#include "config.h"

#include <immintrin.h>

/*  For now, delegate to generic implementation.
    The three-tier dispatch in arbint_ntt.c will select this function when
    AVX2 is available, but currently it provides no speedup over generic.

    TODO: Implement vectorized butterflies for actual performance benefit.  */
arbint_err_t arbint_mul_mag_ntt_avx2(arbint_limb_t * dst, size_t * out_used,
                                     const arbint_limb_t * a, size_t an,
                                     const arbint_limb_t * b, size_t bn,
                                     const arbint_alloc_t * alloc) {
  /*  Delegate to generic implementation for now.  */
  return arbint_mul_mag_ntt_generic(dst, out_used, a, an, b, bn, alloc);
}
