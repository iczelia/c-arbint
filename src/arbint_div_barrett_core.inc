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

/*  Two-limb reciprocal for 3-by-2 division.

    Preconditions:
      - d1 has MSB set (normalized divisor), i.e., d1 >= B/2.
      - d = d1*B + d0 is the full 2-limb divisor.

    Postcondition: *dinv = floor((B^3 - 1) / d) - B.

    Invariant: B <= floor((B^3-1)/d) <= 2*B - 1, so 0 <= dinv < B.  */
static inline void arbint_recip_2limb(arbint_limb_t * dinv, arbint_limb_t d1,
                                      arbint_limb_t d0) {
  arbint_limb_t v;
  arbint_limb_t p;
  arbint_limb_t t1;
  arbint_limb_t t0;
  arbint_limb_t mask;

  v = ARBINT_DIV_PREPARE_FN(d1);

  p = d1 * v;
  p += d0;
  if (p < d0) {
    v--;
    mask = (arbint_limb_t) 0u - (arbint_limb_t) (p >= d1);
    p -= d1;
    v += mask;
    p -= mask & d1;
  }

  ARBINT_DIV_UMUL_FN(&t1, &t0, d0, v);
  p += t1;
  if (p < t1) {
    v--;
    if (p >= d1) {
      if (p > d1 || t0 >= d0)
        v--;
    }
  }

  *dinv = v;
}

/*  Divide 3-limb by 2-limb, producing 1-limb quotient and 2-limb remainder.

    Preconditions:
      - d1 has MSB set (normalized divisor), i.e., d1 >= B/2.
      - d = d1*B + d0 is the 2-limb divisor.
      - n = n2*B^2 + n1*B + n0 is the 3-limb dividend.
      - CRITICAL: n < d*B, which means n2 < d1 (quotient fits in one limb).
      - dinv was computed by arbint_recip_2limb(d1, d0) for this divisor.

    Postconditions:
      - *q = floor(n / d), a single limb (guaranteed by precondition).
      - (*r1, *r0) = n mod d, a 2-limb remainder with 0 <= r < d.
      - Mathematical identity: n = q * d + r.  */
static inline void arbint_div_3by2(arbint_limb_t * q, arbint_limb_t * r1,
                                   arbint_limb_t * r0, arbint_limb_t n2,
                                   arbint_limb_t n1, arbint_limb_t n0,
                                   arbint_limb_t d1, arbint_limb_t d0,
                                   arbint_limb_t dinv) {
  arbint_limb_t _q;
  arbint_limb_t _q0;
  arbint_limb_t _r1;
  arbint_limb_t _r0;
  arbint_limb_t _t1;
  arbint_limb_t _t0;
  arbint_limb_t _mask;
#if ARBINT_HAVE_X86_CARRY_KERNEL
  arbint_x86_carry_word_t _lo;
  arbint_x86_carry_word_t _hi;
  unsigned char _cy;
#else
  arbint_limb_t _x;
#endif

  ARBINT_DIV_UMUL_FN(&_q, &_q0, n2, dinv);
#if ARBINT_HAVE_X86_CARRY_KERNEL
  _cy = ARBINT_X86_ADDCARRY(0, _q0, n1, &_lo);
  (void) ARBINT_X86_ADDCARRY(_cy, _q, n2, &_hi);
  _q0 = (arbint_limb_t) _lo;
  _q = (arbint_limb_t) _hi;
#else
  _x = _q0 + n1;
  _q = _q + n2 + (_x < _q0);
  _q0 = _x;
#endif

  _r1 = n1 - d1 * _q;

#if ARBINT_HAVE_X86_CARRY_KERNEL
  _cy = ARBINT_X86_SUBBORROW(0, n0, d0, &_lo);
  (void) ARBINT_X86_SUBBORROW(_cy, _r1, d1, &_hi);
  _r0 = (arbint_limb_t) _lo;
  _r1 = (arbint_limb_t) _hi;
#else
  _x = n0 - d0;
  _r1 = _r1 - d1 - (n0 < d0);
  _r0 = _x;
#endif

  ARBINT_DIV_UMUL_FN(&_t1, &_t0, d0, _q);
#if ARBINT_HAVE_X86_CARRY_KERNEL
  _cy = ARBINT_X86_SUBBORROW(0, _r0, _t0, &_lo);
  (void) ARBINT_X86_SUBBORROW(_cy, _r1, _t1, &_hi);
  _r0 = (arbint_limb_t) _lo;
  _r1 = (arbint_limb_t) _hi;
#else
  _x = _r0 - _t0;
  _r1 = _r1 - _t1 - (_r0 < _t0);
  _r0 = _x;
#endif

  _q++;

  _mask = (arbint_limb_t) 0u - (arbint_limb_t) (_r1 >= _q0);
  _q += _mask;
#if ARBINT_HAVE_X86_CARRY_KERNEL
  _cy = ARBINT_X86_ADDCARRY(0, _r0, _mask & d0, &_lo);
  (void) ARBINT_X86_ADDCARRY(_cy, _r1, _mask & d1, &_hi);
  _r0 = (arbint_limb_t) _lo;
  _r1 = (arbint_limb_t) _hi;
#else
  _x = _r0 + (_mask & d0);
  _r1 = _r1 + (_mask & d1) + (_x < _r0);
  _r0 = _x;
#endif

  if (_r1 >= d1) {
    if (_r1 > d1 || _r0 >= d0) {
      _q++;
#if ARBINT_HAVE_X86_CARRY_KERNEL
      _cy = ARBINT_X86_SUBBORROW(0, _r0, d0, &_lo);
      (void) ARBINT_X86_SUBBORROW(_cy, _r1, d1, &_hi);
      _r0 = (arbint_limb_t) _lo;
      _r1 = (arbint_limb_t) _hi;
#else
      _x = _r0 - d0;
      _r1 = _r1 - d1 - (_r0 < d0);
      _r0 = _x;
#endif
    }
  }

  *q = _q;
  *r1 = _r1;
  *r0 = _r0;
}
#endif /*  ARBINT_DIV_OMIT_STEP  */

#undef ARBINT_DIV_PREPARE_FN
#undef ARBINT_DIV_OMIT_STEP
#undef ARBINT_DIV_STEP_FN
#undef ARBINT_DIV_UMUL_FN
