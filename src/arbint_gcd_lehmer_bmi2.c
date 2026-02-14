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

/*  Lehmer GCD - BMI2-optimized implementation.

    This file instantiates the Lehmer GCD template with BMI2-optimized
    multiply primitives. Uses _mulx_u64 for 64x64->128 multiplication.

    BMI2 optimization:
    -----------------
    Uses _mulx_u64 intrinsic for 64x64->128 multiplication:
      - Latency: 1 cycle (vs ~5 for half-limb decomposition)
      - Throughput: 1/cycle on modern Intel/AMD
      - No flag modification (unlike IMUL)

    Performance impact:
    - Matrix application: ~2.5x faster
    - Overall Lehmer GCD: ~1.5-2x faster on large operands

    For algorithm documentation, see arbint_gcd_lehmer_core.inc.  */

#define ARBINT_USE_BMI2_INTRIN 1

#include "arbint_gcd.h"

#include "arbint_addsub.h"
#include "arbint_base.h"
#include "arbint_div.h"
#include "arbint_mul.h"
#include "arbint_shift.h"

#include <string.h>

/*  Define template parameters for BMI2 implementation.  */
#define ARBINT_LEHMER_MUL_1_FN     arbint_lehmer_mul_1_bmi2
#define ARBINT_LEHMER_APPLY_FN    arbint_lehmer_apply_matrix_bmi2
#define ARBINT_LEHMER_FALLBACK_FN arbint_gcd_binary_fallback_bmi2
#define ARBINT_LEHMER_ENTRY_FN    arbint_gcd_lehmer_bmi2
#define ARBINT_LEHMER_UMUL_FN     arbint_umul_limb_bmi2

#include "arbint_gcd_lehmer_core.inc"
