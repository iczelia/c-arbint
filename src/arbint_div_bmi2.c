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

#define ARBINT_USE_BMI2_INTRIN 1
#include "arbint_div.h"
#include "arbint_mul.h"

#include "config.h"

#include <string.h>

/*  Multiply two limbs producing full double-width result (hi:lo = a * b).

    BMI2-optimized variant using _mulx_u64 intrinsic for 64-bit limbs, which
    provides faster wide multiplication on Haswell+ (latency 3-4 cycles vs
    4-5 for IMUL, throughput 0.5 vs 1.0). On 32-bit, falls back to uint64_t
    cast.

    Algorithm complexity: O(1) with 1 wide multiplication (BMI2) or hardware
    multiply (32-bit).

    Precondition: hi and lo must point to valid writable limbs.  */
static inline void arbint_umul(arbint_limb_t * hi, arbint_limb_t * lo,
                               arbint_limb_t a, arbint_limb_t b) {
  arbint_umul_limb_bmi2(hi, lo, a, b);
}

/*  Shared reciprocal, division-step, and 3-by-2 division core.  */
#define ARBINT_DIV_PREPARE_FN arbint_prepare_barrett
#define ARBINT_DIV_STEP_FN arbint_utdiv_barrett
#define ARBINT_DIV_UMUL_FN arbint_umul
#include "arbint_div_barrett_core.inc"

/*  BMI2-optimized single-limb division using reciprocal method.
    Uses _mulx_u64 for fast wide multiply in reciprocal-based division
    step.  */
arbint_err_t arbint_div_mag_single_limb_bmi2(const arbint_limb_t * np,
                                             size_t nn, arbint_limb_t d_limb,
                                             arbint_limb_t * qp,
                                             size_t * q_used,
                                             arbint_limb_t * rem_out) {
  arbint_limb_t d_norm;
  arbint_limb_t di;
  arbint_limb_t rem;
  unsigned shift;
  size_t i;

  if (np == NULL || d_limb == 0u || rem_out == NULL)
    return ARBINT_EINVAL;

  if (nn == 0u) {
    if (q_used != NULL)
      *q_used = 0u;
    *rem_out = 0u;
    return ARBINT_OK;
  }

  if (qp != NULL && q_used == NULL)
    return ARBINT_EINVAL;

  shift = arbint_clz_limb(d_limb);
  d_norm = d_limb << shift;
  di = arbint_prepare_barrett(d_norm);

  rem = (shift != 0u) ? (np[nn - 1u] >> (ARBINT_LIMB_BITS - shift)) : 0u;

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

    arbint_utdiv_barrett(&qi, &rem, rem, nl, d_norm, di);

    if (qp != NULL)
      qp[i - 1u] = qi;
  }

  rem >>= shift;

  if (qp != NULL && q_used != NULL)
    *q_used = arbint_norm_used(qp, nn);
  *rem_out = rem;

  return ARBINT_OK;
}

/*  Return normalized dividend limb u[idx] for u = np << shift.
    Valid for idx in [0, nn], where u[nn] is the top carry limb.  */
static inline arbint_limb_t
arbint_get_shifted_limb_bmi2(const arbint_limb_t * np, size_t nn,
                             unsigned shift, size_t idx) {
  if (shift == 0u)
    return (idx < nn) ? np[idx] : (arbint_limb_t) 0u;

  if (idx == nn)
    return np[nn - 1u] >> (ARBINT_LIMB_BITS - shift);
  if (idx == 0u)
    return np[0] << shift;
  return (np[idx] << shift) | (np[idx - 1u] >> (ARBINT_LIMB_BITS - shift));
}

/*  BMI2-optimized 3-by-2 long division for dn == 2.
    Uses fused arbint_div_3by2 with 2-limb reciprocal for ~3 wide muls
    per quotient digit instead of ~5-7 in the old implementation.  */
arbint_err_t arbint_div_mag_two_limb_bmi2(const arbint_limb_t * np, size_t nn,
                                          const arbint_limb_t * dp,
                                          arbint_limb_t * qp,
                                          arbint_limb_t * rp) {
  unsigned shift;
  arbint_limb_t d0;
  arbint_limb_t d1;
  arbint_limb_t dinv;
  arbint_limb_t r1;
  arbint_limb_t r0;
  size_t j;
  size_t qn;

  if (np == NULL || dp == NULL || rp == NULL)
    return ARBINT_EINVAL;
  if (nn < 2u || dp[1] == 0u)
    return ARBINT_EINVAL;

  /*  Normalize divisor: shift so MSB of d1 is set.  */
  shift = arbint_clz_limb(dp[1]);
  if (shift != 0u) {
    d0 = dp[0] << shift;
    d1 = (dp[1] << shift) | (dp[0] >> (ARBINT_LIMB_BITS - shift));
  } else {
    d0 = dp[0];
    d1 = dp[1];
  }

  /*  Compute 2-limb reciprocal ONCE (outside the loop).  */
  arbint_recip_2limb(&dinv, d1, d0);

  /*  Check for shift overflow: when shift > 0, the normalized dividend may
      have an extra leading limb (shifted[nn] = np[nn-1] >> (64 - shift)).
      If this is non-zero, the quotient may have nn limbs instead of nn-1.  */
  {
    arbint_limb_t overflow =
        (shift != 0u) ? (np[nn - 1u] >> (ARBINT_LIMB_BITS - shift)) : 0u;

    if (overflow != 0u) {
      /*  With overflow: quotient has up to nn limbs.
          Initialize (r1, r0) from (shifted[nn], shifted[nn-1]).  */
      r1 = overflow;
      r0 = arbint_get_shifted_limb_bmi2(np, nn, shift, nn - 1u);
      qn = nn;
    } else {
      /*  No overflow: quotient has up to nn-1 limbs.
          Initialize (r1, r0) from (shifted[nn-1], shifted[nn-2]).  */
      r1 = arbint_get_shifted_limb_bmi2(np, nn, shift, nn - 1u);
      r0 = arbint_get_shifted_limb_bmi2(np, nn, shift, nn - 2u);
      qn = nn - 1u;
    }
  }

  /*  Handle MSQ: ensure r1 < d1 before entering main loop.
      This satisfies arbint_div_3by2's precondition: n2 < d1.  */
  if (r1 > d1 || (r1 == d1 && r0 >= d0)) {
    arbint_x86_carry_word_t _lo;
    arbint_x86_carry_word_t _hi;
    unsigned char _cy;
    _cy = ARBINT_X86_SUBBORROW(0, r0, d0, &_lo);
    (void) ARBINT_X86_SUBBORROW(_cy, r1, d1, &_hi);
    r0 = (arbint_limb_t) _lo;
    r1 = (arbint_limb_t) _hi;
    if (qp != NULL)
      qp[qn - 1u] = 1u;
  } else {
    if (qp != NULL)
      qp[qn - 1u] = 0u;
  }

  /*  Main division loop: each iteration produces one quotient limb.
      After the MSQ handling above, r1 < d1 is guaranteed.

      Index mapping: for quotient position j, we bring in normalized dividend
      limb at index j (with/without overflow accounted for in qn).  */
  for (j = qn - 1u; j-- != 0u;) {
    arbint_limb_t u0 = arbint_get_shifted_limb_bmi2(np, nn, shift, j);
    arbint_limb_t q_digit;

    arbint_div_3by2(&q_digit, &r1, &r0, r1, r0, u0, d1, d0, dinv);

    if (qp != NULL)
      qp[j] = q_digit;
  }

  /*  Un-normalize remainder.  */
  if (shift != 0u) {
    rp[0] = (r0 >> shift) | (r1 << (ARBINT_LIMB_BITS - shift));
    rp[1] = r1 >> shift;
  } else {
    rp[0] = r0;
    rp[1] = r1;
  }

  return ARBINT_OK;
}

/*  BMI2-optimized truncated division by uint32_t using reciprocal algorithm.
    Uses _mulx_u64 for fast wide multiply in reciprocal-based division step.
    1.5-2x faster than generic implementation.  */
arbint_err_t arbint_div_qr_u32_bmi2_impl(arbint_t q, arbint_t r,
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

  /*  CRITICAL: Seed the remainder with bits shifted out of the top limb.  */
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

    arbint_utdiv_barrett(&qi, &rem, rem, nl, d_norm, di);

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

/*  BMI2-optimized Barrett reduction for modular exponentiation.

    Reduce x mod d in-place using precomputed Barrett parameters.
    Uses _mulx_u64 for fast wide multiply in reciprocal-based division step.

    Parameters:
      x      - Input/output: value to reduce in-place
      d_norm - Normalized divisor: d << shift (MSB set)
      di     - Precomputed reciprocal from arbint_prepare_barrett(d_norm)
      shift  - Normalization shift: arbint_clz_limb(d)

    Truncated division semantics: remainder sign matches dividend sign.  */
arbint_err_t arbint_mod_u32_barrett_bmi2(arbint_t x, arbint_limb_t d_norm,
                                         arbint_limb_t di, unsigned shift) {
  const arbint_limb_t * xp;
  size_t xn;
  int xsign;
  arbint_err_t rc;
  arbint_limb_t rem;
  size_t i;

  if (x == NULL)
    return ARBINT_EINVAL;

  rc = arbint_get_mag_view(x, &xp, &xn, &xsign);
  if (rc != ARBINT_OK)
    return rc;

  if (xn == 0u) {
    arbint_zero(x);
    return ARBINT_OK;
  }

  /* Seed remainder with carry from top limb (on-the-fly shifting). */
  rem = (shift != 0u) ? (xp[xn - 1u] >> (ARBINT_LIMB_BITS - shift)) : 0u;

  /* Process limbs high-to-low, discard quotient digits. */
  for (i = xn; i != 0u; --i) {
    arbint_limb_t nl;
    arbint_limb_t qi;

    if (shift != 0u) {
      nl = (xp[i - 1u] << shift);
      if (i >= 2u)
        nl |= (xp[i - 2u] >> (ARBINT_LIMB_BITS - shift));
    } else {
      nl = xp[i - 1u];
    }

    arbint_utdiv_barrett(&qi, &rem, rem, nl, d_norm, di);
    (void) qi; /* Quotient not needed. */
  }

  /* Un-normalize remainder. */
  rem >>= shift;

  /* Store result in x. */
  if (rem == 0u) {
    arbint_zero(x);
  } else {
    int rsign = (xsign == 0) ? 0 : xsign;
    rc = arbint_resize(x, 1u);
    if (rc != ARBINT_OK)
      return rc;
    ARBINT_LIMBS(x)[0] = rem;
    if (!arbint_set_signed_sz(x, 1u, rsign))
      return ARBINT_EOVERFLOW;
  }

  return ARBINT_OK;
}

/*  Fast truncated quotient by 3 using fixed Barrett constants.
    Computes q = trunc(n / 3) and discards the remainder.
    Supports aliasing (q may be n).  */
arbint_err_t arbint_tdiv_q_3_bmi2(arbint_t q, const arbint_t n) {
  const arbint_limb_t * np;
  size_t nn;
  size_t i;
  size_t q_used;
  int nsign;
  arbint_limb_t rem;
  arbint_err_t rc;

#if ARBINT_LIMB_BITS == 64
  static const arbint_limb_t d_norm = UINT64_C(0xC000000000000000);
  static const arbint_limb_t di = UINT64_C(0x5555555555555555);
  static const unsigned shift = 62u;
#elif ARBINT_LIMB_BITS == 32
  static const arbint_limb_t d_norm = UINT32_C(0xC0000000);
  static const arbint_limb_t di = UINT32_C(0x55555555);
  static const unsigned shift = 30u;
#else
  #error "Unsupported ARBINT_LIMB_BITS for div3"
#endif

  if (q == NULL || n == NULL)
    return ARBINT_EINVAL;

  nsign = (n[0]._sz > 0) - (n[0]._sz < 0);
  if (nsign == 0) {
    arbint_zero(q);
    return ARBINT_OK;
  }

  nn = arbint_abs_sz(n[0]._sz);
  rc = arbint_resize(q, nn);
  if (rc != ARBINT_OK)
    return rc;

  np = ARBINT_CLIMBS(n);
  rem = np[nn - 1u] >> (ARBINT_LIMB_BITS - shift);

  for (i = nn; i != 0u; --i) {
    arbint_limb_t nl = np[i - 1u] << shift;
    arbint_limb_t qi;
    if (i >= 2u)
      nl |= np[i - 2u] >> (ARBINT_LIMB_BITS - shift);

    arbint_utdiv_barrett(&qi, &rem, rem, nl, d_norm, di);
    ARBINT_LIMBS(q)[i - 1u] = qi;
  }

  q_used = arbint_norm_used(ARBINT_LIMBS(q), nn);
  if (!arbint_set_signed_sz(q, q_used, (q_used == 0u) ? 0 : nsign))
    return ARBINT_EOVERFLOW;
  return ARBINT_OK;
}
