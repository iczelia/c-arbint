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

#include <stdlib.h>
#include <string.h>

/*  Multiply two limbs producing full double-width result (hi:lo = a * b).

    For 64-bit limbs without __uint128_t, uses half-limb (32x32) multiplication
    to synthesize a 128-bit product via Karatsuba-like decomposition:
      a = a1*2^32 + a0,  b = b1*2^32 + b0
      a*b = a1*b1*2^64 + ((a1*b0 + a0*b1)*2^32) + a0*b0

    For 32-bit limbs, casts to uint64_t and uses hardware multiply.

    Algorithm complexity: O(1) with 4 half-limb multiplications (64-bit case)

    Precondition: hi and lo must point to valid writable limbs.  */
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
#endif /* ARBINT_LIMB_BITS */

/*  Compute reciprocal of a normalized limb for Barrett reduction.

    Given a normalized divisor d (MSB set, d >= beta/2 where beta =
    2^LIMB_BITS), computes the reciprocal v = floor((beta^2 - 1) / d) - beta.

    This reciprocal allows replacing expensive hardware division with two
    multiplications in the division-by-reciprocal algorithm
    (Granlund-Montgomery method). The normalization requirement (MSB set)
    ensures sufficient precision in the reciprocal approximation.

    Algorithm: Computes the reciprocal using half-limb arithmetic via Newton
    iteration adapted from mini-gmp. First computes high half of reciprocal
    via ~d / d1 (where d1 is the high half of d), then refines to full
    precision by computing the low half with remainder tracking and correction
    steps.

    Parameters:
      d - Normalized divisor (MSB must be set, i.e., d >= 2^(LIMB_BITS-1))

    Returns:
      Reciprocal v = floor((2^(2*LIMB_BITS) - 1) / d) - 2^LIMB_BITS

    Precondition: d must be normalized (d >= 2^(LIMB_BITS-1)). If d is not
    normalized, the reciprocal will have insufficient precision and division
    results will be incorrect. Caller must normalize via left-shift before
    calling this function.

    Algorithm complexity: O(1) with ~10 half-limb operations.

    Note: The returned reciprocal v satisfies beta <= v < beta + (beta^2 -
    1)/d, where beta = 2^LIMB_BITS. This property is essential for the
    correctness of the division-by-reciprocal algorithm.  */
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

/*  Perform single division step using precomputed reciprocal (Barrett
   reduction).

    Divides the two-limb dividend (nh:nl) by the single-limb divisor d,
    producing quotient q and remainder r. Uses precomputed reciprocal di to
    replace the expensive hardware DIV instruction with two multiplications
    plus a few ALU operations (~5-10 cycles instead of ~40-80 cycles for
    hardware DIV on x86-64).

    This is the core primitive of the Granlund-Montgomery division algorithm.
    The reciprocal-based approach trades one expensive division (when computing
    the reciprocal) for many cheap multiplications (when dividing each limb),
    resulting in dramatic speedups for multi-limb division by a single limb.

    Algorithm: Estimates the quotient via qh = (nh * di + nh + 1 + carry) >>
    BITS, then applies two correction steps to handle the +/-1 error inherent
    in the reciprocal approximation. The first correction uses branchless
    arithmetic (mask trick) for better performance; the second uses a simple
    conditional.

    Parameters:
      q  - Output pointer for quotient (single limb)
      r  - Output pointer for remainder (single limb)
      nh - High limb of two-limb dividend (most significant)
      nl - Low limb of two-limb dividend (least significant)
      d  - Normalized divisor (MSB set, i.e., d >= 2^(LIMB_BITS-1))
      di - Precomputed reciprocal from arbint_prepare_barrett(d)

    Returns:
      *q = floor((nh * 2^LIMB_BITS + nl) / d)
      *r = (nh * 2^LIMB_BITS + nl) mod d

    Preconditions:
      - nh < d (ensures quotient fits in one limb)
      - d is normalized (MSB set, guaranteed by caller's shift)
      - di = arbint_prepare_barrett(d) (correct reciprocal)
      - q and r point to valid writable limbs

    Algorithm complexity: O(1) with 2 wide multiplications + ~5-10 ALU ops.

    Performance: On modern x86-64, this is ~5-8x faster than hardware DIV.
    On ARM, speedup is even more dramatic (~10-15x) due to slow DIV
    instruction.

    Correctness: The reciprocal approximation may be off by +/-1, hence the two
    correction steps. The first correction handles overestimation (quotient too
    high), the second handles underestimation (quotient too low). After both
    corrections, the result is always exact.  */
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
  mask = (arbint_limb_t) 0u - (arbint_limb_t) (_r > ql);
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

/*  Generic truncated division by uint32_t (magnitude) using Barrett
    reduction.  */
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

  /*  CRITICAL: Seed the remainder with bits shifted out of the top limb.

      The division algorithm normalizes the divisor by shifting it left so its
      MSB is set. To maintain mathematical correctness, the dividend must be
      shifted by the same amount. However, we perform this shift "on the fly"
      during the main loop to avoid allocating a temporary shifted copy.

      When conceptually shifting the entire dividend left by 'shift' bits, the
      carry out of the top limb np[nn-1] becomes the initial remainder:
        initial_rem = np[nn-1] >> (LIMB_BITS - shift)

      This value is always < d_norm (since shift > 0 means d_norm >= 2^(BITS-1)
      and the carry is at most 2^shift - 1), so it's a valid starting
      remainder.  */
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

/*  Single-limb division using reciprocal method.
    Divides n-limb dividend by single-limb divisor using precomputed
    reciprocal.  */
arbint_err_t
arbint_div_mag_single_limb_generic(const arbint_limb_t * np, size_t nn,
                                   arbint_limb_t d_limb, arbint_limb_t * qp,
                                   size_t * q_used, arbint_limb_t * rem_out) {
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

    arbint_ubarrett(&qi, &rem, rem, nl, d_norm, di);

    if (qp != NULL)
      qp[i - 1u] = qi;
  }

  rem >>= shift;

  if (qp != NULL && q_used != NULL)
    *q_used = arbint_norm_used(qp, nn);
  *rem_out = rem;

  return ARBINT_OK;
}

/*  Knuth Algorithm D: Multi-limb division by multi-limb divisor.

    Implements Knuth's Algorithm D from TAOCP Vol 2, Section 4.3.1.
    Divides n-limb dividend by m-limb divisor, producing quotient and
    remainder.

    Algorithm overview:
    1. Normalize: shift divisor left so MSB is set (improves quotient
       estimation)
    2. For each quotient position (high to low):
       a. Estimate quotient digit using top 2-3 limbs
       b. Multiply estimate by divisor
       c. Subtract from dividend
       d. Adjust if estimate was too large (rare)
    3. Un-normalize remainder

    Complexity: O(n * m) for n-limb dividend and m-limb divisor.
    This is dramatically faster than the binary algorithm's O(n^2 * bits).

    Parameters:
      np  - Dividend limbs (n limbs)
      nn  - Number of dividend limbs
      dp  - Divisor limbs (m limbs)
      dn  - Number of divisor limbs (must be >= 2)
      qp  - Output quotient buffer ((nn - dn + 1) limbs, can be NULL)
      rp  - Output remainder buffer (dn limbs, must not be NULL)

    Returns: ARBINT_OK on success, error code on failure.

    Preconditions:
    - nn >= dn >= 2
    - dp[dn-1] != 0 (divisor is normalized to actual used limbs)
    - qp has space for (nn - dn + 1) limbs if non-NULL
    - rp has space for dn limbs (always required)

    Postconditions:
    - If qp != NULL: quotient written to qp[0..nn-dn]
    - Remainder written to rp[0..dn-1]
    - Result is exact (no approximation)

    References:
    - Knuth TAOCP Vol 2, Algorithm 4.3.1.D
    - Hacker's Delight (2nd ed.), Chapter 9  */
arbint_err_t arbint_div_mag_knuth(const arbint_limb_t * np, size_t nn,
                                  const arbint_limb_t * dp, size_t dn,
                                  arbint_limb_t * qp, arbint_limb_t * rp) {
  arbint_limb_t di;
  unsigned shift;
  size_t i;
  size_t j;
  arbint_limb_t qhat;
  arbint_limb_t * u;
  arbint_limb_t * d_shifted;
  arbint_limb_t d_top;
  arbint_limb_t u_top2;
  arbint_limb_t u_top1;
  arbint_limb_t u_top0;
  arbint_limb_t carry;
  arbint_limb_t borrow;
  arbint_limb_t prod_hi;
  arbint_limb_t prod_lo;

  if (np == NULL || dp == NULL || rp == NULL)
    return ARBINT_EINVAL;
  if (dn < 2u || nn < dn)
    return ARBINT_EINVAL;
  if (dp[dn - 1u] == 0u)
    return ARBINT_EINVAL;

  /* Allocate working storage for normalized dividend and divisor.
     u needs nn+1 limbs (extra limb for carry from normalization).
     d_shifted needs dn limbs.  */
  u = (arbint_limb_t *) malloc((nn + 1u) * sizeof(arbint_limb_t));
  if (u == NULL)
    return ARBINT_ENOMEM;

  d_shifted = (arbint_limb_t *) malloc(dn * sizeof(arbint_limb_t));
  if (d_shifted == NULL) {
    free(u);
    return ARBINT_ENOMEM;
  }

  /* Step D1: Normalize.
     Shift divisor left so MSB is set. This ensures quotient estimation
     is accurate (off by at most 1). Same shift applied to dividend.  */
  shift = arbint_clz_limb(dp[dn - 1u]);

  /* Normalize divisor into d_shifted.  */
  if (shift != 0u) {
    carry = 0u;
    for (i = 0u; i < dn; ++i) {
      d_shifted[i] = (dp[i] << shift) | carry;
      carry = dp[i] >> (ARBINT_LIMB_BITS - shift);
    }
  } else {
    memcpy(d_shifted, dp, dn * sizeof(arbint_limb_t));
  }

  /* Normalize dividend into u.  */
  if (shift != 0u) {
    carry = 0u;
    for (i = 0u; i < nn; ++i) {
      u[i] = (np[i] << shift) | carry;
      carry = np[i] >> (ARBINT_LIMB_BITS - shift);
    }
    u[nn] = carry;
  } else {
    memcpy(u, np, nn * sizeof(arbint_limb_t));
    u[nn] = 0u;
  }

  /* Precompute reciprocal of divisor top limb for quotient estimation.  */
  d_top = d_shifted[dn - 1u];
  di = arbint_prepare_barrett(d_top);

  /* Step D2-D7: Main loop.
     Process quotient positions from high to low (j = nn - dn down to 0).  */
  for (j = nn - dn + 1u; j != 0u; --j) {
    size_t idx = j - 1u;

    /* Step D3: Estimate quotient digit.
       Use top 3 limbs of current dividend position: u[idx+dn], u[idx+dn-1],
       u[idx+dn-2]. Divide (u[idx+dn]:u[idx+dn-1]) by d_top using reciprocal.
     */
    u_top2 = u[idx + dn];
    u_top1 = u[idx + dn - 1u];
    u_top0 = u[idx + dn - 2u];

    if (u_top2 >= d_top) {
      /* Special case: dividend >= beta * divisor, quotient will be max limb
         value. Occurs when previous subtraction borrowed.  */
      qhat = (arbint_limb_t) (~(arbint_limb_t) 0u);
    } else {
      /* Normal case: use reciprocal to estimate quotient.  */
      arbint_limb_t rhat;
      arbint_ubarrett(&qhat, &rhat, u_top2, u_top1, d_top, di);

      /* Refine estimate using next limb (Knuth's correction step).
         Check if qhat * d[dn-2] > rhat * beta + u[idx+dn-2].
         If so, decrement qhat and update rhat. Loop until condition is false.
         This loop executes at most twice.  */
      arbint_limb_t d_second = d_shifted[dn - 2u];
      arbint_limb_t prod_hi_test;
      arbint_limb_t prod_lo_test;

      while (1) {
        arbint_umul(&prod_hi_test, &prod_lo_test, qhat, d_second);

        /* Test: qhat * d[dn-2] > rhat * beta + u[idx+dn-2]
           Equivalently: prod_hi_test > rhat or
                         (prod_hi_test == rhat and prod_lo_test > u_top0)  */
        if (prod_hi_test > rhat ||
            (prod_hi_test == rhat && prod_lo_test > u_top0)) {
          qhat--;
          rhat += d_top;
          /* If rhat overflowed (>= beta), stop: qhat is now definitely
           * correct.  */
          if (rhat < d_top)
            break;
          /* Otherwise, continue loop to check again.  */
        } else {
          break;
        }
      }
    }

    /* Step D4: Multiply and subtract.
       Compute u[idx..idx+dn] -= qhat * d_shifted[0..dn-1].
       We compute qhat * d_shifted and subtract it from u[idx..idx+dn].
       The product is at most (beta-1) * (beta-1) = beta^2 - 2*beta + 1,
       which requires dn+1 limbs. We subtract limb-by-limb with borrow
       propagation.  */
    borrow = 0u;
    for (i = 0u; i < dn; ++i) {
      arbint_limb_t ui = u[idx + i];
      arbint_umul(&prod_hi, &prod_lo, qhat, d_shifted[i]);

      /* Subtract (prod_hi:prod_lo) from (borrow:ui), giving new borrow and new
         ui. First: ui -= prod_lo  */
      arbint_limb_t new_ui = ui - prod_lo;
      arbint_limb_t borrow1 = (new_ui > ui) ? 1u : 0u;

      /* Second: new_ui -= borrow  */
      arbint_limb_t new_ui2 = new_ui - borrow;
      arbint_limb_t borrow2 = (new_ui2 > new_ui) ? 1u : 0u;

      /* New borrow is prod_hi + borrow1 + borrow2.  */
      borrow = prod_hi + borrow1 + borrow2;
      u[idx + i] = new_ui2;
    }

    /* Subtract final borrow from top limb.  */
    arbint_limb_t old_top = u[idx + dn];
    u[idx + dn] = old_top - borrow;

    /* Step D5: Test remainder.
       If subtraction caused borrow out (old_top < borrow), the quotient digit
       was too large by 1. This happens rarely (probability < 2/beta for
       normalized divisor).  */
    if (u[idx + dn] > old_top) {
      /* Step D6: Add back divisor and decrement quotient.  */
      qhat--;
      carry = 0u;
      for (i = 0u; i < dn; ++i) {
        arbint_limb_t sum_lo = u[idx + i] + carry;
        arbint_limb_t carry1 = (sum_lo < carry) ? 1u : 0u;
        arbint_limb_t sum = sum_lo + d_shifted[i];
        arbint_limb_t carry2 = (sum < sum_lo) ? 1u : 0u;
        u[idx + i] = sum;
        carry = carry1 + carry2;
      }
      u[idx + dn] += carry;
    }

    /* Store quotient digit.  */
    if (qp != NULL)
      qp[idx] = qhat;
  }

  /* Step D8: Un-normalize remainder.
     The remainder is in u[0..dn-1], but it's still shifted. Shift right to
     restore.  */
  if (shift != 0u) {
    for (i = 0u; i < dn - 1u; ++i)
      rp[i] = (u[i] >> shift) | (u[i + 1u] << (ARBINT_LIMB_BITS - shift));
    rp[dn - 1u] = u[dn - 1u] >> shift;
  } else {
    memcpy(rp, u, dn * sizeof(arbint_limb_t));
  }

  free(d_shifted);
  free(u);
  return ARBINT_OK;
}
