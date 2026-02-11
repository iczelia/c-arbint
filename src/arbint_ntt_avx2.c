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

    AVX2 implies BMI2 on all modern CPUs (Intel Haswell, AMD Zen and later),
    so we use _mulx_u64 for fast 64x64->128 multiplication.

    Future optimization opportunities:
    1. Vectorize mod_add/mod_sub: 4 parallel 64-bit additions with conditional
       subtraction using _mm256_cmpgt_epi64 or unsigned comparison trick.
    2. Vectorize butterfly memory access: _mm256_loadu_si256 for loading
       4 elements at once.
    3. Consider cache-friendly 4-step FFT for very large transforms.  */

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

#define ARBINT_NTT_MUL_MAG_FN arbint_mul_mag_ntt_avx2
#include "arbint_ntt_core.inc"
