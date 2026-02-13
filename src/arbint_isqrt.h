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

#ifndef ARBINT_ISQRT_H
#define ARBINT_ISQRT_H

#include "arbint_base.h"

/*  Threshold for switching from Newton iteration to Karatsuba sqrt.
    Below this threshold, simple Newton is used.
    Above it, divide-and-conquer Karatsuba sqrt is used.

    The optimal threshold depends on the relative cost of division vs
    multiplication. A value of 4-8 limbs is typical.  */
#define ARBINT_ISQRT_KARATSUBA_THRESHOLD 8u

/*  Table for initial 1/sqrt approximation (8-bit precision).
    Entry i (for i in 0..127) approximates floor(256 / sqrt((128+i)/128)).
    This gives about 8 bits of precision for the top byte of the input.  */
extern const unsigned char arbint_invsqrt_tab[128];

/*  Compute floor(sqrt(x)) for a single limb x.
    x must be normalized (top bit set, i.e., x >= 2^(LIMB_BITS-1)).
    Returns the square root.  */
arbint_limb_t arbint_isqrt_1(arbint_limb_t x);

/*  Compute floor(sqrt({np, 2})) for a 2-limb input.
    Input {np, 2} must be normalized (np[1] >= 2^(LIMB_BITS-2)).
    Returns sqrt in *sp, remainder carry in return value.  */
arbint_limb_t arbint_isqrt_2(arbint_limb_t * sp, const arbint_limb_t * np);

/*  Shared Newton sqrt backend for small magnitudes.
    Computes rop = floor(sqrt(a)) for non-negative a.  */
arbint_err_t arbint_isqrt_newton_mag(arbint_t rop, const arbint_t a);

/*  Karatsuba divide-and-conquer square root.

    Computes floor(sqrt({np, 2n})) and writes result to {sp, n}.
    Returns remainder (0 or 1) indicating if input was a perfect square.

    Input {np, 2n} must be normalized: np[2n-1] >= 2^(LIMB_BITS-2).
    Output {sp, n} is the integer square root.

    scratch must have capacity for at least floor(n/2)+1 limbs.

    The algorithm is based on "Karatsuba Square Root" by Paul Zimmermann.
    Complexity: O(M(n)) where M(n) is multiplication time.  */
int arbint_isqrt_dc(arbint_limb_t * sp, arbint_limb_t * np, size_t n,
                    arbint_limb_t * scratch, const arbint_alloc_t * alloc);

#endif /*  ARBINT_ISQRT_H  */
