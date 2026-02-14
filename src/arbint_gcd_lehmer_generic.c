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

/*  Lehmer GCD - Portable (generic) implementation.

    This file instantiates the Lehmer GCD template with generic (portable)
    multiply primitives. Uses half-limb decomposition for wide multiplication.

    For algorithm documentation, see arbint_gcd_lehmer_core.inc.  */

#include "arbint_gcd.h"

#include "arbint_addsub.h"
#include "arbint_base.h"
#include "arbint_div.h"
#include "arbint_mul.h"
#include "arbint_shift.h"

#include <string.h>

/*  Define template parameters for generic implementation.  */
#define ARBINT_LEHMER_MUL_1_FN     arbint_lehmer_mul_1_generic
#define ARBINT_LEHMER_APPLY_FN    arbint_lehmer_apply_matrix_generic
#define ARBINT_LEHMER_FALLBACK_FN arbint_gcd_binary_fallback
#define ARBINT_LEHMER_ENTRY_FN    arbint_gcd_lehmer_generic
#define ARBINT_LEHMER_UMUL_FN     arbint_umul_limb_generic

#include "arbint_gcd_lehmer_core.inc"
