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

#include "arbint_bitops.h"

#include "config.h"

#include <nmmintrin.h>

/*  POPCNT-accelerated popcount over a limb array.
    Uses the _mm_popcnt_u64 (64-bit) or _mm_popcnt_u32 (32-bit) intrinsic
    which maps directly to the hardware POPCNT instruction.  */
size_t arbint__popcount_limbs_popcnt(const arbint_limb_t * x, size_t n) {
  size_t count = 0u;
  size_t i;
  for (i = 0u; i < n; ++i) {
#if ARBINT_LIMB_BITS == 64
    count += (size_t) _mm_popcnt_u64((unsigned long long) x[i]);
#else
    count += (size_t) _mm_popcnt_u32((unsigned int) x[i]);
#endif /* ARBINT_LIMB_BITS */
  }
  return count;
}

/*  POPCNT-accelerated hamming distance between two limb arrays.
    Computes popcount(a[i] ^ b[i]) for overlapping limbs, then
    popcount of the remaining tail from the longer array.  */
size_t arbint__hamming_limbs_popcnt(const arbint_limb_t * a, size_t an,
                                    const arbint_limb_t * b, size_t bn) {
  size_t count = 0u;
  size_t i;
  size_t min_n = (an < bn) ? an : bn;

  for (i = 0u; i < min_n; ++i) {
#if ARBINT_LIMB_BITS == 64
    count += (size_t) _mm_popcnt_u64((unsigned long long) (a[i] ^ b[i]));
#else
    count += (size_t) _mm_popcnt_u32((unsigned int) (a[i] ^ b[i]));
#endif /* ARBINT_LIMB_BITS */
  }
  for (i = min_n; i < an; ++i) {
#if ARBINT_LIMB_BITS == 64
    count += (size_t) _mm_popcnt_u64((unsigned long long) a[i]);
#else
    count += (size_t) _mm_popcnt_u32((unsigned int) a[i]);
#endif /* ARBINT_LIMB_BITS */
  }
  for (i = min_n; i < bn; ++i) {
#if ARBINT_LIMB_BITS == 64
    count += (size_t) _mm_popcnt_u64((unsigned long long) b[i]);
#else
    count += (size_t) _mm_popcnt_u32((unsigned int) b[i]);
#endif /* ARBINT_LIMB_BITS */
  }
  return count;
}
