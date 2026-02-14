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

#ifndef ARBINT_DIV_BZ_H
#define ARBINT_DIV_BZ_H

#include "arbint_base.h"

/*  Burnikel-Ziegler Division Threshold.

    Minimum divisor size (in limbs) to use BZ divide-and-conquer.
    Below this, Knuth's Algorithm D is used as base case.

    BZ achieves O(n^1.58) complexity vs Knuth's O(n*m) by reducing
    division to multiplication via recursive 3n/2n structure.

    Threshold tuned based on:
    - BZ overhead: recursive calls, scratch allocation
    - Karatsuba threshold: 15 limbs
    - Typical crossover: ~2.5x Karatsuba threshold  */
#ifndef ARBINT_BZ_THRESHOLD
#define ARBINT_BZ_THRESHOLD 40u
#endif

/*  Burnikel-Ziegler division: compute q = floor(n/d) and r = n mod d.

    Preconditions:
      - np points to nn limbs of dividend (little-endian)
      - dp points to dn limbs of divisor (little-endian)
      - dp[dn-1] != 0 (divisor normalized to actual limb count)
      - dn >= ARBINT_BZ_THRESHOLD
      - qp has capacity for (nn - dn + 1) limbs if non-NULL
      - rp has capacity for dn limbs (required)

    Postconditions:
      - n = q * d + r with 0 <= r < d
      - If qp != NULL: quotient written to qp[0..nn-dn]
      - Remainder written to rp[0..dn-1]

    Complexity: O(M(n) * log(n/m)) where M(n) is multiplication cost.  */
arbint_err_t arbint_div_mag_bz(const arbint_limb_t * np, size_t nn,
                                const arbint_limb_t * dp, size_t dn,
                                arbint_limb_t * qp, arbint_limb_t * rp,
                                const arbint_alloc_t * alloc);

/*  Compute scratch size needed for arbint_div_mag_bz.  */
size_t arbint_div_bz_scratch_size(size_t nn, size_t dn);

#endif /*  ARBINT_DIV_BZ_H  */
