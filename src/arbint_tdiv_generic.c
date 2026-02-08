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

#include "arbint_tdiv.h"

#include "config.h"

#include <string.h>

#if ARBINT_LIMB_BITS == 64
static inline void arbint_umul(arbint_limb_t * hi, arbint_limb_t * lo,
                               arbint_limb_t a, arbint_limb_t b) {
  arbint_limb_t a0 = a & ARBINT_HALF_MASK;
  arbint_limb_t a1 = a >> ARBINT_HALF_BITS;
  arbint_limb_t b0 = b & ARBINT_HALF_MASK;
  arbint_limb_t b1 = b >> ARBINT_HALF_BITS;
  arbint_limb_t w0 = a0 * b0;
  arbint_limb_t t = a1 * b0 + (w0 >> ARBINT_HALF_BITS);
  arbint_limb_t w1 = t & ARBINT_HALF_MASK;
  arbint_limb_t w2 = t >> ARBINT_HALF_BITS;
  w1 = a0 * b1 + w1;
  *hi = a1 * b1 + w2 + (w1 >> ARBINT_HALF_BITS);
  *lo = (w1 << ARBINT_HALF_BITS) | (w0 & ARBINT_HALF_MASK);
}
#elif ARBINT_LIMB_BITS == 32
static inline void arbint_umul(arbint_limb_t * hi, arbint_limb_t * lo,
                               arbint_limb_t a, arbint_limb_t b) {
  uint64_t p = (uint64_t) a * (uint64_t) b;
  *lo = (arbint_limb_t) p;
  *hi = (arbint_limb_t) (p >> 32);
}
#else
  #error "Unsupported limb size"
#endif

/*  Reciprocal of a normalized limb d (MSB set).
    Returns v = floor((beta^2 - 1) / d) - beta.  */
static inline arbint_limb_t arbint_prepare_barrett(arbint_limb_t d) {
  arbint_limb_t d0;
  arbint_limb_t d1;
  arbint_limb_t v;
  arbint_limb_t p;
  arbint_limb_t r;
  arbint_limb_t t;
  arbint_limb_t ql;

  d1 = d >> ARBINT_HALF_BITS;
  d0 = d & ARBINT_HALF_MASK;

  /* Compute high half of reciprocal: qh = ~d / d1 (half-by-half). */
  v = (arbint_limb_t) ((~d) / d1);
  r = ((~d) - v * d1) << ARBINT_HALF_BITS;
  r |= ARBINT_HALF_MASK;

  /* Adjust for d0. */
  p = v * d0;
  if (r < p) {
    v--;
    r += d;
    if (r >= d && r < p) {
      v--;
      r += d;
    }
  }
  r -= p;

  /* Compute low half of reciprocal. */
  t = (r >> ARBINT_HALF_BITS) * v + r;
  ql = (t >> ARBINT_HALF_BITS) + (arbint_limb_t) 1u;

  r = (r << ARBINT_HALF_BITS) + ARBINT_HALF_MASK - ql * d;
  if (r >= (t << ARBINT_HALF_BITS)) {
    ql--;
    r += d;
  }

  v = (v << ARBINT_HALF_BITS) + ql;
  if (r >= d) {
    v++;
  }

  return v;
}

/*  Division step using precomputed reciprocal.
    Divides (nh : nl) by d, producing quotient *q and remainder *r.
    Requires: nh < d, d normalized (MSB set), di = invert_limb(d).  */
static inline void arbint_ubarrett(arbint_limb_t * q, arbint_limb_t * r,
                                   arbint_limb_t nh, arbint_limb_t nl,
                                   arbint_limb_t d, arbint_limb_t di) {
  arbint_limb_t qh;
  arbint_limb_t ql;
  arbint_limb_t _r;
  arbint_limb_t mask;

  arbint_umul(&qh, &ql, nh, di);

  arbint_limb_t lo = ql + nl;
  arbint_limb_t carry = (lo < ql) ? (arbint_limb_t) 1u : (arbint_limb_t) 0u;
  ql = lo;
  qh = qh + (nh + (arbint_limb_t) 1u) + carry;

  _r = nl - qh * d;

  /* First correction: if the estimate was 1 too high. */
  mask = -(arbint_limb_t) (_r > ql);
  qh += mask;
  _r += mask & d;

  /* Second correction: if the estimate was 1 too low. */
  if (_r >= d) {
    _r -= d;
    qh++;
  }

  *q = qh;
  *r = _r;
}

arbint_err_t arbint_tdiv_qr_u32_generic_impl(arbint_t q, arbint_t r,
                                             const arbint_t n, uint32_t dmag,
                                             int dsign) {
  const arbint_limb_t * np;
  size_t nn;
  int nsign;
  arbint_err_t rc;
  size_t i;
  arbint_limb_t rem;
  int qsign;
  int rsign;
  size_t q_used;
  arbint_limb_t d_norm;
  arbint_limb_t di;
  unsigned shift;

  if ((q == NULL && r == NULL) || n == NULL)
    return ARBINT_EINVAL;
  if (dmag == 0u)
    return ARBINT_EZERO;

  rc = arbint_get_mag_view(n, &np, &nn, &nsign);
  if (rc != ARBINT_OK)
    return rc;

  if (nn == 0u) {
    if (q != NULL)
      arbint_zero(q);
    if (r != NULL)
      arbint_zero(r);
    return ARBINT_OK;
  }

  if (q != NULL) {
    rc = arbint_resize(q, nn);
    if (rc != ARBINT_OK)
      return rc;
  }

  /*  Promote the 32-bit divisor to full limb width and normalize
      (shift left so the MSB is set). Precompute the reciprocal once.  */
  d_norm = (arbint_limb_t) dmag;
  shift = arbint_clz_limb(d_norm);
  d_norm <<= shift;
  di = arbint_prepare_barrett(d_norm);

  /*  Seed the remainder with the bits shifted out of the top limb.
      When the entire dividend is conceptually shifted left by 'shift'
      bits, the carry out of np[nn-1] is np[nn-1] >> (BITS - shift).
      This is always < d_norm, so it's a valid initial remainder.  */
  if (shift != 0u)
    rem = np[nn - 1u] >> (ARBINT_LIMB_BITS - shift);
  else
    rem = (arbint_limb_t) 0u;

  for (i = nn; i != 0u; --i) {
    arbint_limb_t nl;
    arbint_limb_t qi;

    if (shift != 0u) {
      nl = (np[i - 1u] << shift);
      if (i >= 2u)
        nl |= (np[i - 2u] >> (ARBINT_LIMB_BITS - shift));
    } else {
      nl = np[i - 1u];
    }

    arbint_ubarrett(&qi, &rem, rem, nl, d_norm, di);

    if (q != NULL)
      ARBINT_LIMBS(q)[i - 1u] = qi;
  }

  /* Un-normalize the remainder. */
  rem >>= shift;

  if (q != NULL) {
    q_used = arbint_norm_used(ARBINT_LIMBS(q), nn);
    qsign = (q_used == 0u || nsign == 0) ? 0 : ((nsign == dsign) ? 1 : -1);
    if (!arbint_set_signed_sz(q, q_used, qsign))
      return ARBINT_EOVERFLOW;
  }

  if (r != NULL) {
    rsign = (rem == 0u || nsign == 0) ? 0 : nsign;
    rc = arbint_resize(r, 1u);
    if (rc != ARBINT_OK)
      return rc;
    ARBINT_LIMBS(r)[0] = rem;
    if (!arbint_set_signed_sz(r, (rem != 0u) ? 1u : 0u, rsign))
      return ARBINT_EOVERFLOW;
  }

  return ARBINT_OK;
}
