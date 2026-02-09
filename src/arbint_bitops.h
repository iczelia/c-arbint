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

/*  Portable popcount of a single limb.  */
static inline unsigned arbint_popcount_limb(arbint_limb_t x) {
#if ARBINT_COMPILER_GNU_CLANG
  #if ARBINT_LIMB_BITS == 64
  return (unsigned) __builtin_popcountll((unsigned long long) x);
  #else
  return (unsigned) __builtin_popcount((unsigned) x);
  #endif /* ARBINT_LIMB_BITS */
#elif ARBINT_COMPILER_MSVC
  #if ARBINT_LIMB_BITS == 64
  return (unsigned) __popcnt64(x);
  #else
  return (unsigned) __popcnt(x);
  #endif /* ARBINT_LIMB_BITS */
#else
  {
    unsigned count = 0u;
    while (x != 0u) {
      x &= x - 1u;
      ++count;
    }
    return count;
  }
#endif /* ARBINT_COMPILER_GNU_CLANG */
}

/*  Count trailing zeros in a single limb.  Undefined if x == 0.  */
static inline unsigned arbint_ctz_limb(arbint_limb_t x) {
#if ARBINT_COMPILER_GNU_CLANG
  #if ARBINT_LIMB_BITS == 64
  return (unsigned) __builtin_ctzll((unsigned long long) x);
  #else
  return (unsigned) __builtin_ctz((unsigned) x);
  #endif /* ARBINT_LIMB_BITS */
#elif ARBINT_COMPILER_MSVC
  {
    unsigned long idx;
  #if ARBINT_LIMB_BITS == 64
    _BitScanForward64(&idx, x);
  #else
    _BitScanForward(&idx, x);
  #endif /* ARBINT_LIMB_BITS */
    return (unsigned) idx;
  }
#else
  {
    unsigned n = 0u;
    while ((x & 1u) == 0u) {
      ++n;
      x >>= 1u;
    }
    return n;
  }
#endif /* ARBINT_COMPILER_GNU_CLANG */
}

/*  Threshold for POPCNT-accelerated popcount/hammingdist.  */
#define ARBINT_POPCOUNT_SSE_THRESHOLD 8u

#if HAS_POPCNT
size_t arbint__popcount_limbs_popcnt(const arbint_limb_t * x, size_t n);
size_t arbint__hamming_limbs_popcnt(const arbint_limb_t * a, size_t an,
                                    const arbint_limb_t * b, size_t bn);
#endif /* HAS_POPCNT */

#endif /* ARBINT_BITOPS_H */
