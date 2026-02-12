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

#ifndef ARBINT_SHIFT_H
#define ARBINT_SHIFT_H

#include "arbint_base.h"

/*  In-place right shift of limb array by k bits.
    Returns new normalized size after shift.
    Safe for k == 0 (no-op) and k >= n * LIMB_BITS (returns 0).  */
static inline size_t arbint_rshift_limbs_inplace(arbint_limb_t * p, size_t n,
                                                 size_t k) {
  size_t limb_shift;
  unsigned bit_shift;
  size_t i;

  if (k == 0u || n == 0u)
    return n;

  limb_shift = k / ARBINT_LIMB_BITS;
  bit_shift = (unsigned) (k % ARBINT_LIMB_BITS);

  if (limb_shift >= n)
    return 0u;

  if (limb_shift > 0u) {
    for (i = 0u; i < n - limb_shift; ++i)
      p[i] = p[i + limb_shift];
    n -= limb_shift;
  }

  if (bit_shift > 0u) {
    unsigned lshift = ARBINT_LIMB_BITS - bit_shift;
    for (i = 0u; i < n - 1u; ++i)
      p[i] = (p[i] >> bit_shift) | (p[i + 1u] << lshift);
    p[n - 1u] >>= bit_shift;
  }

  return arbint_norm_used(p, n);
}

/*  In-place left shift of limb array by k bits.
    Caller must ensure p has capacity for at least (n + k/LIMB_BITS + 1) limbs.
    Returns new size after shift.  */
static inline size_t arbint_lshift_limbs_inplace(arbint_limb_t * p, size_t n,
                                                 size_t k, size_t cap) {
  size_t limb_shift;
  unsigned bit_shift;
  size_t new_n;
  size_t i;
  arbint_limb_t carry;

  if (k == 0u || n == 0u)
    return n;

  limb_shift = k / ARBINT_LIMB_BITS;
  bit_shift = (unsigned) (k % ARBINT_LIMB_BITS);

  new_n = n + limb_shift;
  if (bit_shift > 0u && (p[n - 1u] >> (ARBINT_LIMB_BITS - bit_shift)) != 0u)
    ++new_n;

  if (new_n > cap)
    return 0u;

  if (bit_shift > 0u) {
    carry = 0u;
    for (i = 0u; i < n; ++i) {
      arbint_limb_t new_carry = p[i] >> (ARBINT_LIMB_BITS - bit_shift);
      p[i] = (p[i] << bit_shift) | carry;
      carry = new_carry;
    }
    if (carry != 0u)
      p[n++] = carry;
  }

  if (limb_shift > 0u) {
    for (i = n; i != 0u; --i)
      p[i - 1u + limb_shift] = p[i - 1u];
    for (i = 0u; i < limb_shift; ++i)
      p[i] = 0u;
    n += limb_shift;
  }

  return n;
}

#endif /* ARBINT_SHIFT_H */
