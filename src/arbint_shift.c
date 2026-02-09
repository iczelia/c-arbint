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

#include "arbint_shift.h"

#include "arbint.h"

#include <string.h>

/*  Left-shift: rop = a * 2^k.
    Preserves sign.  Handles aliasing (rop == a).

    Algorithm:
      k = limb_shift * LIMB_BITS + bit_shift.
      The result has at most (an + limb_shift + 1) limbs.
      Low limb_shift limbs of result are zero.
      Remaining limbs are shifted copies of the source with carries
      propagated upward.  */
arbint_err_t arbint_shl(arbint_t rop, const arbint_t a, uint32_t k) {
  int sign;
  size_t an;
  size_t limb_shift;
  unsigned bit_shift;
  size_t need;
  size_t i;
  arbint_err_t rc;
  const arbint_limb_t * ap;
  arbint_limb_t * rp;
  arbint_limb_t carry;

  if (rop == NULL || a == NULL)
    return ARBINT_EINVAL;

  sign = (a[0]._sz > 0) - (a[0]._sz < 0);
  if (sign == 0 || k == 0u) {
    if (sign == 0) {
      arbint_zero(rop);
      return ARBINT_OK;
    }
    return arbint_set(rop, a);
  }

  an = arbint_abs_sz(a[0]._sz);
  limb_shift = (size_t) (k / ARBINT_LIMB_BITS);
  bit_shift = (unsigned) (k % ARBINT_LIMB_BITS);

  /*  need = an + limb_shift + (bit_shift ? 1 : 0).
      Check for overflow.  */
  need = an + limb_shift;
  if (need < an)
    return ARBINT_EOVERFLOW;
  if (bit_shift != 0u) {
    if (need == SIZE_MAX)
      return ARBINT_EOVERFLOW;
    ++need;
  }

  if (rop[0]._cap < need) {
    rc = arbint_resize(rop, need);
    if (rc != ARBINT_OK)
      return rc;
  }

  rp = ARBINT_LIMBS(rop);
  ap = ARBINT_CLIMBS(a);

  if (bit_shift == 0u) {
    /*  Whole-limb shift only.  Copy high-to-low so in-place aliasing
        is safe (destination indices are always >= source indices).  */
    for (i = an; i != 0u; --i)
      rp[i - 1u + limb_shift] = ap[i - 1u];
  } else {
    /*  Combined limb + bit shift.  Process high-to-low.  */
    unsigned rshift = ARBINT_LIMB_BITS - bit_shift;
    carry = 0u;

    /*  Top limb: carry from the highest source limb.  */
    carry = ap[an - 1u] >> rshift;
    rp[an + limb_shift] = carry;

    /*  Middle limbs: shift and combine adjacent source limbs.  */
    for (i = an - 1u; i != 0u; --i)
      rp[i + limb_shift] = (ap[i] << bit_shift) | (ap[i - 1u] >> rshift);

    /*  Lowest shifted limb.  */
    rp[limb_shift] = ap[0] << bit_shift;
  }

  /*  Zero the low limb_shift positions.  */
  if (limb_shift != 0u)
    memset(rp, 0, limb_shift * sizeof(arbint_limb_t));

  size_t used = arbint_norm_used(rp, need);
  if (!arbint_set_signed_sz(rop, used, sign))
    return ARBINT_EOVERFLOW;
  return ARBINT_OK;
}

/*  Right-shift: rop = trunc(a / 2^k) (truncation toward zero).
    Preserves sign.  Handles aliasing (rop == a).

    Algorithm:
      k = limb_shift * LIMB_BITS + bit_shift.
      If limb_shift >= an, the result is zero.
      Otherwise, discard low limb_shift limbs and shift remaining
      limbs rightward by bit_shift bits.  */
arbint_err_t arbint_shr(arbint_t rop, const arbint_t a, uint32_t k) {
  int sign;
  size_t an;
  size_t limb_shift;
  unsigned bit_shift;
  size_t remain;
  size_t i;
  arbint_err_t rc;
  const arbint_limb_t * ap;
  arbint_limb_t * rp;

  if (rop == NULL || a == NULL)
    return ARBINT_EINVAL;

  sign = (a[0]._sz > 0) - (a[0]._sz < 0);
  if (sign == 0 || k == 0u) {
    if (sign == 0) {
      arbint_zero(rop);
      return ARBINT_OK;
    }
    return arbint_set(rop, a);
  }

  an = arbint_abs_sz(a[0]._sz);
  limb_shift = (size_t) (k / ARBINT_LIMB_BITS);
  bit_shift = (unsigned) (k % ARBINT_LIMB_BITS);

  /*  If we shift away all limbs, result is zero.  */
  if (limb_shift >= an) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  remain = an - limb_shift;

  if (rop[0]._cap < remain) {
    rc = arbint_resize(rop, remain);
    if (rc != ARBINT_OK)
      return rc;
  }

  rp = ARBINT_LIMBS(rop);
  ap = ARBINT_CLIMBS(a);

  if (bit_shift == 0u) {
    /*  Whole-limb shift only.  Copy low-to-high so in-place aliasing
        is safe (destination indices are always <= source indices).  */
    for (i = 0u; i < remain; ++i)
      rp[i] = ap[i + limb_shift];
  } else {
    unsigned lshift = ARBINT_LIMB_BITS - bit_shift;

    /*  Shift and combine adjacent source limbs, low to high.  */
    for (i = 0u; i < remain - 1u; ++i)
      rp[i] = (ap[i + limb_shift] >> bit_shift) |
              (ap[i + limb_shift + 1u] << lshift);

    /*  Top limb: no higher limb to borrow from.  */
    rp[remain - 1u] = ap[remain - 1u + limb_shift] >> bit_shift;
  }

  size_t used = arbint_norm_used(rp, remain);
  if (used == 0u) {
    arbint_zero(rop);
    return ARBINT_OK;
  }
  if (!arbint_set_signed_sz(rop, used, sign))
    return ARBINT_EOVERFLOW;
  return ARBINT_OK;
}
