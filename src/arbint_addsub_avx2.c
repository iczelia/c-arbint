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

#include "arbint_addsub.h"

#include "config.h"

#include <immintrin.h>
#include <string.h>

/*  Helper function to align carry bits from AVX2 vector.

    Given carries = [c3, c2, c1, c0] (carry-out from each 64-bit limb),
    we need to align them as [c2, c1, c0, carry_scalar] where carry_scalar
    is the carry from the previous iteration.

    This is done by:
    1. Permuting the vector to shift carries right by one position
    2. Blending in the scalar carry at position 0  */
static inline __m256i arbint_avx2_align_carries(__m256i carries,
                                                uint64_t carry_scalar) {
  /*  Permute to shift right: [c3,c2,c1,c0] -> [c2,c1,c0,c3]
      The last element c3 will be overwritten by carry_scalar.  */
  __m256i shifted_carries =
      _mm256_permute4x64_epi64(carries, _MM_SHUFFLE(2, 1, 0, 3));

  /*  Create vector with carry_scalar in position 0, zeros elsewhere.  */
  __m256i scalar_carry_vec =
      _mm256_set_epi64x(0, 0, 0, (long long) carry_scalar);

  /*  Blend: keep shifted_carries[3:1], replace [0] with carry_scalar.
      Blend mask 0x03 = 0b00000011 selects lower 2 32-bit elements
      (which is the lowest 64-bit element in the 4x64-bit vector).  */
  __m256i aligned =
      _mm256_blend_epi32(shifted_carries, scalar_carry_vec, 0x03);

  return aligned;
}

/*  AVX2-optimized magnitude doubling implementation.

    Doubles a multi-limb magnitude (multiply by 2) using AVX2 SIMD intrinsics.
    Processes 4 limbs at a time using 256-bit vectors.

    Algorithm:
    - Extract carry-out bits (bit 63) from each limb
    - Shift all limbs left by 1
    - OR in the aligned carries from previous limbs
    - Handle remaining 0-3 limbs with scalar code

    Returns the number of limbs in the result (may be nx or nx+1 if final
    carry).  */
size_t arbint__dbl_mag_avx2(arbint_limb_t * dst, const arbint_limb_t * x,
                            size_t nx) {
#if ARBINT_LIMB_BITS == 64
  size_t i;
  uint64_t carry_scalar = 0u;
  size_t n_vec = nx & ~((size_t) 3); /*  Number of limbs to process with SIMD
                                        (multiple of 4).  */

  /*  Phase 1: Process complete 4-limb chunks with AVX2.  */
  for (i = 0u; i < n_vec; i += 4u) {
    /*  Load 4 limbs into AVX2 vector: [x[i+3], x[i+2], x[i+1], x[i+0]]  */
    __m256i vec = _mm256_loadu_si256((const __m256i *) &x[i]);

    /*  Extract carry-out bits (bit 63 of each limb).  */
    __m256i carries = _mm256_srli_epi64(vec, 63);

    /*  Shift each limb left by 1.  */
    __m256i shifted = _mm256_slli_epi64(vec, 1);

    /*  Align carries to match the shifted limbs.  */
    __m256i aligned_carries = arbint_avx2_align_carries(carries, carry_scalar);

    /*  Combine shifted limbs with aligned carries.  */
    __m256i result = _mm256_or_si256(shifted, aligned_carries);

    /*  Store result.  */
    _mm256_storeu_si256((__m256i *) &dst[i], result);

    /*  Extract the top carry for the next iteration (from limb i+3).  */
    carry_scalar = (uint64_t) _mm256_extract_epi64(carries, 3);
  }

  /*  Phase 2: Handle remaining 0-3 limbs with scalar code.  */
  for (; i < nx; ++i) {
    arbint_limb_t xi = x[i];
    arbint_limb_t d = (arbint_limb_t) (xi << 1);
    dst[i] = d + (arbint_limb_t) carry_scalar;
    carry_scalar = (uint64_t) (xi >> 63);
  }

  /*  Phase 3: Write final carry if non-zero.  */
  if (carry_scalar != 0u) {
    dst[nx] = (arbint_limb_t) carry_scalar;
    return nx + 1u;
  }

  return nx;
#else
  #error "AVX2 enabled but 64-bit mode not supported?"
#endif /* ARBINT_LIMB_BITS == 64 */
}
