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

#ifndef ARBINT_DIV_PREPARE_FN
  #error "ARBINT_DIV_PREPARE_FN must be defined"
#endif

#ifndef ARBINT_DIV_OMIT_STEP
  #ifndef ARBINT_DIV_STEP_FN
    #error "ARBINT_DIV_STEP_FN must be defined"
  #endif

  #ifndef ARBINT_DIV_UMUL_FN
    #error "ARBINT_DIV_UMUL_FN must be defined"
  #endif
#endif

static inline arbint_limb_t ARBINT_DIV_PREPARE_FN(arbint_limb_t d) {
  arbint_limb_t d0;
  arbint_limb_t d1;
  arbint_limb_t v;
  arbint_limb_t p;
  arbint_limb_t r;
  arbint_limb_t t;
  arbint_limb_t ql;

  d1 = d >> ARBINT_HALF_BITS;
  d0 = d & ARBINT_HALF_MASK;

  v = (arbint_limb_t) ((~d) / d1);
  r = ((~d) - v * d1) << ARBINT_HALF_BITS;
  r |= ARBINT_HALF_MASK;

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

  t = (r >> ARBINT_HALF_BITS) * v + r;
  ql = (t >> ARBINT_HALF_BITS) + (arbint_limb_t) 1u;

  r = (r << ARBINT_HALF_BITS) + ARBINT_HALF_MASK - ql * d;
  if (r >= (t << ARBINT_HALF_BITS)) {
    ql--;
    r += d;
  }

  v = (v << ARBINT_HALF_BITS) + ql;
  if (r >= d)
    v++;

  return v;
}

#ifndef ARBINT_DIV_OMIT_STEP
static inline void ARBINT_DIV_STEP_FN(arbint_limb_t * q, arbint_limb_t * r,
                                      arbint_limb_t nh, arbint_limb_t nl,
                                      arbint_limb_t d, arbint_limb_t di) {
  arbint_limb_t qh;
  arbint_limb_t ql;
  arbint_limb_t _r;
  arbint_limb_t mask;

  ARBINT_DIV_UMUL_FN(&qh, &ql, nh, di);

  {
    arbint_limb_t lo = ql + nl;
    arbint_limb_t carry = (lo < ql) ? (arbint_limb_t) 1u : (arbint_limb_t) 0u;
    ql = lo;
    qh = qh + (nh + (arbint_limb_t) 1u) + carry;
  }

  _r = nl - qh * d;

  mask = (arbint_limb_t) 0u - (arbint_limb_t) (_r > ql);
  qh += mask;
  _r += mask & d;

  if (_r >= d) {
    _r -= d;
    qh++;
  }

  *q = qh;
  *r = _r;
}
#endif

#undef ARBINT_DIV_PREPARE_FN
#undef ARBINT_DIV_OMIT_STEP
#undef ARBINT_DIV_STEP_FN
#undef ARBINT_DIV_UMUL_FN
