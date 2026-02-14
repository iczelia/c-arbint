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

/*  Shared reciprocal and division-step core.  */
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
    Allocation-free rolling Knuth-D step with BMI2 reciprocal/multiply ops.  */
arbint_err_t arbint_div_mag_two_limb_bmi2(const arbint_limb_t * np, size_t nn,
                                          const arbint_limb_t * dp,
                                          arbint_limb_t * qp,
                                          arbint_limb_t * rp) {
  unsigned shift;
  arbint_limb_t d0;
  arbint_limb_t d1;
  arbint_limb_t di;
  arbint_limb_t r1;
  arbint_limb_t r0;
  size_t j;

  if (np == NULL || dp == NULL || rp == NULL)
    return ARBINT_EINVAL;
  if (nn < 2u || dp[1] == 0u)
    return ARBINT_EINVAL;

  shift = arbint_clz_limb(dp[1]);
  if (shift != 0u) {
    d0 = dp[0] << shift;
    d1 = (dp[1] << shift) | (dp[0] >> (ARBINT_LIMB_BITS - shift));
  } else {
    d0 = dp[0];
    d1 = dp[1];
  }
  di = arbint_prepare_barrett(d1);

  r1 = arbint_get_shifted_limb_bmi2(np, nn, shift, nn);
  r0 = arbint_get_shifted_limb_bmi2(np, nn, shift, nn - 1u);

  for (j = nn - 1u; j-- != 0u;) {
    arbint_limb_t u0 = arbint_get_shifted_limb_bmi2(np, nn, shift, j);
    arbint_limb_t qhat;

    if (r1 >= d1) {
      qhat = ~(arbint_limb_t) 0u;
    } else {
      arbint_limb_t rhat;
      arbint_limb_t prod_hi_test;
      arbint_limb_t prod_lo_test;

      arbint_utdiv_barrett(&qhat, &rhat, r1, r0, d1, di);
      while (1) {
        arbint_umul(&prod_hi_test, &prod_lo_test, qhat, d0);
        if (prod_hi_test > rhat ||
            (prod_hi_test == rhat && prod_lo_test > u0)) {
          qhat--;
          rhat += d1;
          if (rhat < d1)
            break;
        } else {
          break;
        }
      }
    }

    {
      arbint_limb_t v0;
      arbint_limb_t v1;
      arbint_limb_t prod0_hi;
      arbint_limb_t prod0_lo;
      arbint_limb_t prod1_hi;
      arbint_limb_t prod1_lo;
      arbint_limb_t p0;
      arbint_limb_t p1;
      arbint_limb_t p2;

      arbint_umul(&prod0_hi, &prod0_lo, qhat, d0);
      arbint_umul(&prod1_hi, &prod1_lo, qhat, d1);
      p0 = prod0_lo;
      p1 = prod1_lo + prod0_hi;
      p2 = prod1_hi + ((p1 < prod1_lo) ? 1u : 0u);

#if ARBINT_HAVE_X86_CARRY_KERNEL
      {
        unsigned char borrow;
        unsigned char carry;
        arbint_x86_carry_word_t out = (arbint_x86_carry_word_t) 0;

        borrow = ARBINT_X86_SUBBORROW(
            0u, (arbint_x86_carry_word_t) u0, (arbint_x86_carry_word_t) p0,
            &out);
        v0 = (arbint_limb_t) out;
        borrow = ARBINT_X86_SUBBORROW(
            borrow, (arbint_x86_carry_word_t) r0, (arbint_x86_carry_word_t) p1,
            &out);
        v1 = (arbint_limb_t) out;
        borrow = ARBINT_X86_SUBBORROW(
            borrow, (arbint_x86_carry_word_t) r1, (arbint_x86_carry_word_t) p2,
            &out);
        if (borrow != 0u) {
          carry = ARBINT_X86_ADDCARRY(0u, (arbint_x86_carry_word_t) v0,
                                      (arbint_x86_carry_word_t) d0, &out);
          v0 = (arbint_limb_t) out;
          carry = ARBINT_X86_ADDCARRY(carry, (arbint_x86_carry_word_t) v1,
                                      (arbint_x86_carry_word_t) d1, &out);
          v1 = (arbint_limb_t) out;
          qhat--;
        }
      }
#else
      {
        arbint_limb_t borrow = 0u;
        arbint_limb_t old_top;
        arbint_limb_t top_after_sub;

        {
          arbint_limb_t t = u0 - p0;
          arbint_limb_t borrow1 = (t > u0) ? 1u : 0u;
          arbint_limb_t t2 = t - borrow;
          arbint_limb_t borrow2 = (t2 > t) ? 1u : 0u;
          v0 = t2;
          borrow = borrow1 | borrow2;
        }
        {
          arbint_limb_t t = r0 - p1;
          arbint_limb_t borrow1 = (t > r0) ? 1u : 0u;
          arbint_limb_t t2 = t - borrow;
          arbint_limb_t borrow2 = (t2 > t) ? 1u : 0u;
          v1 = t2;
          borrow = borrow1 | borrow2;
        }

        old_top = r1;
        {
          arbint_limb_t borrow_in = borrow;
          arbint_limb_t t = old_top - p2;
          arbint_limb_t borrow1 = (t > old_top) ? 1u : 0u;
          top_after_sub = t - borrow_in;
          borrow = borrow1 | ((top_after_sub > t) ? 1u : 0u);
        }
        if (borrow != 0u) {
          arbint_limb_t carry = 0u;
          arbint_limb_t sum_lo = v0 + carry;
          arbint_limb_t carry1 = (sum_lo < carry) ? 1u : 0u;
          arbint_limb_t sum = sum_lo + d0;
          arbint_limb_t carry2 = (sum < sum_lo) ? 1u : 0u;
          v0 = sum;
          carry = carry1 + carry2;

          sum_lo = v1 + carry;
          carry1 = (sum_lo < carry) ? 1u : 0u;
          sum = sum_lo + d1;
          carry2 = (sum < sum_lo) ? 1u : 0u;
          v1 = sum;
          qhat--;
        }
      }
#endif

      r1 = v1;
      r0 = v0;
    }

    if (qp != NULL)
      qp[j] = qhat;
  }

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
