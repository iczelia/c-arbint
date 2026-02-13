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

#include "arbint_mul.h"
#include "arbint_ntt.h"
#include "config.h"

#include <stdlib.h>
#include <string.h>

/* ========== Wide Multiplication ========== */

/*  Multiply two 64-bit integers producing full 128-bit result (hi:lo = a * b).

    Uses half-limb (32x32) multiplication to synthesize 128-bit product
    via Karatsuba-like decomposition (same pattern as seen in
    arbint_div_generic.c).  */
static inline void ntt_umul(uint64_t * hi, uint64_t * lo, uint64_t a,
                            uint64_t b) {
  arbint_umul_u64_generic(hi, lo, a, b);
}

#define ARBINT_NTT_CACHE_CLEAR_FN arbint_ntt_cache_clear_generic
#define ARBINT_NTT_MUL_MAG_FN arbint_mul_mag_ntt_generic
#define ARBINT_NTT_SQR_MAG_FN arbint_sqr_mag_ntt_generic
#include "arbint_ntt_core.inc"
