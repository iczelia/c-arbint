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

#ifndef ARBINT_BITOPS_H
#define ARBINT_BITOPS_H

#include "arbint_base.h"

/*  Threshold for POPCNT-accelerated popcount/hammingdist.  */
#define ARBINT_POPCOUNT_SSE_THRESHOLD 8u

#if HAS_POPCNT
size_t arbint__popcount_limbs_popcnt(const arbint_limb_t * x, size_t n);
size_t arbint__hamming_limbs_popcnt(const arbint_limb_t * a, size_t an,
                                    const arbint_limb_t * b, size_t bn);
#endif /* HAS_POPCNT */

#endif /* ARBINT_BITOPS_H */
