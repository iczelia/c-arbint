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

/*  Internal limb-level helpers for binary GCD.
    These are designed for tight inner loops with minimal overhead.  */

/*  Count trailing zeros across a limb array.
    Returns the total number of trailing zero bits.
    Precondition: n > 0 and the value is nonzero (at least one nonzero limb).
    For zero magnitude (n == 0 or all limbs zero), behavior is undefined. */
static inline size_t arbint_gcd_mag_ctz(const arbint_limb_t * p, size_t n) {
  size_t i;
  size_t ctz = 0u;

  /* Count full zero limbs. */
  for (i = 0u; i < n && p[i] == 0u; ++i)
    ctz += ARBINT_LIMB_BITS;

  /* Count trailing zeros in first nonzero limb. */
  if (i < n)
    ctz += arbint_ctz_limb(p[i]);

  return ctz;
}

/*  Strip leading zero limbs, return new normalized size.
    Returns 0 if all limbs are zero. */
static inline size_t arbint_gcd_normalize(const arbint_limb_t * p, size_t n) {
  while (n > 0u && p[n - 1u] == 0u)
    --n;
  return n;
}

/*  Compare magnitudes. Returns -1 if a < b, 0 if a == b, +1 if a > b.
    Both arrays must be normalized (no leading zeros). */
static inline int arbint_gcd_cmp_mag(const arbint_limb_t * a, size_t an,
                                     const arbint_limb_t * b, size_t bn) {
  size_t i;

  if (an != bn)
    return (an > bn) ? 1 : -1;

  /* Same length: compare from high to low. */
  for (i = an; i != 0u; --i) {
    if (a[i - 1u] != b[i - 1u])
      return (a[i - 1u] > b[i - 1u]) ? 1 : -1;
  }

  return 0;
}

/*  In-place right shift by k bits.
    Precondition: k < ARBINT_LIMB_BITS * n (shift doesn't exceed total bits).
    Returns new normalized size after shift. */
static inline size_t arbint_gcd_rshift_inplace(arbint_limb_t * p, size_t n,
                                               size_t k) {
  size_t limb_shift;
  unsigned bit_shift;
  size_t i;
  arbint_limb_t carry;

  if (k == 0u || n == 0u)
    return n;

  limb_shift = k / ARBINT_LIMB_BITS;
  bit_shift = (unsigned) (k % ARBINT_LIMB_BITS);

  /* Shift by whole limbs first. */
  if (limb_shift > 0u) {
    if (limb_shift >= n) {
      /* Shift exceeds size: result is zero. */
      p[0] = 0u;
      return 0u;
    }
    for (i = 0u; i < n - limb_shift; ++i)
      p[i] = p[i + limb_shift];
    n -= limb_shift;
  }

  /* Shift remaining bits within limbs. */
  if (bit_shift > 0u) {
    carry = 0u;
    for (i = n; i != 0u; --i) {
      arbint_limb_t new_carry = p[i - 1u] << (ARBINT_LIMB_BITS - bit_shift);
      p[i - 1u] = (p[i - 1u] >> bit_shift) | carry;
      carry = new_carry;
    }
  }

  /* Normalize: strip leading zeros. */
  return arbint_gcd_normalize(p, n);
}

/*  In-place left shift by k bits.
    Caller must ensure p has enough capacity for the expanded result.
    Returns new size after shift. */
static inline size_t arbint_gcd_lshift_inplace(arbint_limb_t * p, size_t n,
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

  /* Calculate new size. */
  new_n = n + limb_shift;
  if (bit_shift > 0u) {
    /* Check if we need an extra limb for overflow. */
    arbint_limb_t top = p[n - 1u];
    if ((top >> (ARBINT_LIMB_BITS - bit_shift)) != 0u)
      ++new_n;
  }

  if (new_n > cap)
    return 0u; /* Overflow: caller error. */

  /* Shift bits within limbs first (from high to low to avoid overwriting). */
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

  /* Shift by whole limbs (from high to low). */
  if (limb_shift > 0u) {
    for (i = n; i != 0u; --i)
      p[i - 1u + limb_shift] = p[i - 1u];
    for (i = 0u; i < limb_shift; ++i)
      p[i] = 0u;
    n += limb_shift;
  }

  return n;
}

/*  In-place subtraction: u -= v, assuming u >= v.
    Both u and v must be normalized.
    Returns new normalized size of u. */
static inline size_t arbint_gcd_sub_inplace(arbint_limb_t * u, size_t un,
                                            const arbint_limb_t * v,
                                            size_t vn) {
  size_t i;
  arbint_limb_t borrow = 0u;

  /* Subtract v from u. */
  for (i = 0u; i < vn; ++i) {
    arbint_limb_t ui = u[i];
    arbint_limb_t vi = v[i];
    arbint_limb_t diff = ui - vi - borrow;
    /* Borrow if ui < vi + borrow (with wrap-around check). */
    borrow = (ui < vi) || (borrow != 0u && ui == vi) ? 1u : 0u;
    u[i] = diff;
  }

  /* Propagate borrow through remaining limbs. */
  for (; i < un && borrow != 0u; ++i) {
    arbint_limb_t ui = u[i];
    u[i] = ui - borrow;
    borrow = (ui == 0u) ? 1u : 0u;
  }

  /* Normalize: strip leading zeros. */
  return arbint_gcd_normalize(u, un);
}

/*  Single-limb Euclidean GCD. */
static inline arbint_limb_t arbint_gcd_limb(arbint_limb_t a, arbint_limb_t b) {
  arbint_limb_t t;
  while (b != 0u) {
    t = b;
    b = a % b;
    a = t;
  }
  return a;
}

/*  u32 Euclidean GCD helper. */
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
