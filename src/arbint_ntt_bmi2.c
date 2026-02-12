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

/*  BMI2-optimized NTT multiplication.

    Uses _mulx_u64 intrinsic for fast 64x64->128 bit multiplication,
    replacing the half-limb decomposition used in the generic version.
    This provides ~25-30% speedup in NTT operations on BMI2-capable CPUs
    (Intel Haswell/AMD Zen and later).  */

#define ARBINT_USE_BMI2_INTRIN 1

#include "arbint_mul.h"
#include "arbint_ntt.h"
#include "config.h"

#include <stdlib.h>
#include <string.h>

/* ========== Wide Multiplication ========== */

/*  Multiply two 64-bit integers producing full 128-bit result (hi:lo = a * b).

    Uses BMI2 _mulx_u64 intrinsic for single-instruction 64x64->128 multiply.
    This replaces 4 half-limb multiplications with one hardware multiply.  */
static inline void ntt_umul(uint64_t * hi, uint64_t * lo, uint64_t a,
                            uint64_t b) {
  arbint_umul_u64_bmi2(hi, lo, a, b);
}

#define ARBINT_NTT_CACHE_CLEAR_FN arbint_ntt_cache_clear_bmi2
#define ARBINT_NTT_MUL_MAG_FN arbint_mul_mag_ntt_bmi2
#include "arbint_ntt_core.inc"
