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

#ifndef ARBINT_GCD_H
#define ARBINT_GCD_H

#include "arbint_base.h"

/*  Threshold: when min(an, bn) <= this, use Euclidean with tdiv_r.  */
#define ARBINT_GCD_EUCLID_THRESHOLD 4u

/*  Count trailing zeros across a limb array.
    Returns the total number of trailing zero bits.
    Precondition: n > 0 and the value is nonzero (at least one nonzero limb).
    For zero magnitude (n == 0 or all limbs zero), behavior is undefined.  */
static inline size_t arbint_gcd_mag_ctz(const arbint_limb_t * p, size_t n) {
  size_t i;
  size_t ctz = 0u;

  for (i = 0u; i < n && p[i] == 0u; ++i)
    ctz += ARBINT_LIMB_BITS;

  if (i < n)
    ctz += arbint_ctz_limb(p[i]);

  return ctz;
}

/*  Single-limb Euclidean GCD.  */
static inline arbint_limb_t arbint_gcd_limb(arbint_limb_t a, arbint_limb_t b) {
  arbint_limb_t t;
  while (b != 0u) {
    t = b;
    b = a % b;
    a = t;
  }
  return a;
}

/*  u32 Euclidean GCD helper.  */
static inline uint32_t arbint_gcd_u32u32(uint32_t a, uint32_t b) {
  uint32_t t;
  while (b != 0u) {
    t = b;
    b = a % b;
    a = t;
  }
  return a;
}

#endif /* ARBINT_GCD_H */
