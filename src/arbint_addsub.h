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

#ifndef ARBINT_ADDSUB_H
#define ARBINT_ADDSUB_H

#include "arbint_base.h"

#if ARBINT_TARGET_X86_FAMILY && ARBINT_COMPILER_GNU_CLANG
  #include <immintrin.h>
  #define ARBINT_HAVE_X86_INTRIN 1
#elif ARBINT_TARGET_X86_FAMILY && ARBINT_COMPILER_MSVC
  #include <intrin.h>
  #define ARBINT_HAVE_X86_INTRIN 1
#else
  #define ARBINT_HAVE_X86_INTRIN 0
#endif /* ARBINT_HAVE_X86_INTRIN */

#if ARBINT_HAVE_X86_INTRIN && ARBINT_LIMB_BITS == 64
  #define ARBINT_HAVE_X86_CARRY_KERNEL 1
typedef unsigned long long arbint_x86_carry_word_t;
  #define ARBINT_X86_ADDCARRY _addcarry_u64
  #define ARBINT_X86_SUBBORROW _subborrow_u64
#elif ARBINT_HAVE_X86_INTRIN && ARBINT_LIMB_BITS == 32
  #define ARBINT_HAVE_X86_CARRY_KERNEL 1
typedef unsigned int arbint_x86_carry_word_t;
  #define ARBINT_X86_ADDCARRY _addcarry_u32
  #define ARBINT_X86_SUBBORROW _subborrow_u32
#else
  #define ARBINT_HAVE_X86_CARRY_KERNEL 0
#endif /* ARBINT_HAVE_X86_CARRY_KERNEL */

/*  Magnitude arithmetic; callers ensure valid non-overflowing dimensions.
    - arbint__add_mag handles either operand order.
    - arbint__sub_mag requires |x| >= |y| and nx >= ny.  */
size_t arbint__add_mag(arbint_limb_t * dst, const arbint_limb_t * x, size_t nx,
                       const arbint_limb_t * y, size_t ny);
size_t arbint__sub_mag(arbint_limb_t * dst, const arbint_limb_t * x, size_t nx,
                       const arbint_limb_t * y, size_t ny);
size_t arbint__dbl_mag(arbint_limb_t * dst, const arbint_limb_t * x,
                       size_t nx);

/*  Threshold for AVX2 vectorization of doubling operation.
    Below this threshold, scalar implementation is used to avoid SIMD overhead.
 */
#define ARBINT_DBL_AVX2_THRESHOLD 16u

#if HAS_AVX2
/*  AVX2-optimized magnitude doubling implementation.  */
size_t arbint__dbl_mag_avx2(arbint_limb_t * dst, const arbint_limb_t * x,
                            size_t nx);
#endif /* HAS_AVX2 */

#endif /* ARBINT_ADDSUB_H */
