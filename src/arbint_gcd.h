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

/*  Threshold: when min(an, bn) >= this, use Lehmer GCD.
    Below this, Binary GCD overhead is lower than Lehmer setup cost.  */
#define ARBINT_LEHMER_THRESHOLD 32u

/*  Lehmer matrix: tracks quotient sequence transformations.
    After processing quotients q0, q1, ..., qk, represents:
      a_new = v0*a_old - v1*b_old  (if k is even)
      b_new = w1*b_old - w0*a_old  (alternating signs)
    Initialized to identity: v0=1, v1=0, w0=0, w1=1.  */
typedef struct {
  arbint_limb_t v0, v1;
  arbint_limb_t w0, w1;
} arbint_lehmer_matrix_t;

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

/*  Lehmer step: simulate Euclidean GCD on top limbs, building a matrix.

    Given the top 2 limbs of a and b (where a >= b), simulates quotient
    computation using single-limb arithmetic. Builds a 2x2 matrix M such
    that applying M to (a, b) produces reduced values.

    Algorithm overview:
    ------------------
    Standard Euclidean GCD computes a sequence of quotients and remainders:
      a = q0*b + r0,  b = q1*r0 + r1,  r0 = q2*r1 + r2, ...

    For multi-limb integers, each step requires expensive O(n) division.
    Lehmer's insight: the first ~LIMB_BITS quotients depend only on the top
    1-2 limbs, so we can simulate them with single-limb arithmetic.

    Instead of actually computing remainders, we track a 2x2 matrix that
    represents the cumulative transformation. After k quotients, we have:
      [u_new]   [v0  -v1] [u]       (if k even, signs flip if k odd)
      [v_new] = [-w0  w1] [v]

    The matrix uses MAGNITUDES only; signs alternate based on parity.

    Matrix coefficient recurrence (all positive):
    ---------------------------------------------
    Initialize: v0=1, v1=0, w0=0, w1=1  (identity matrix)
    For each quotient q:
      (v0, v1) <- (w0, w1)           // shift old w into v
      (w0, w1) <- (v0 + q*w0, v1 + q*w1)  // new w = old_v + q * old_w

    After k steps:
      - If k is even: u_new = v0*u - v1*v,  v_new = w1*v - w0*u
      - If k is odd:  u_new = v1*v - v0*u,  v_new = w0*u - w1*v

    Critical safety checks:
    ----------------------
    The simulation can diverge from true multi-limb quotients when:
    1. Matrix coefficients approach overflow (we stop at LIMB_MAX/2)
    2. The simulated quotient is uncertain (q from top limbs != true q)
    3. Simulated b becomes zero before the real b would

    When any safety check triggers, we return count=0 or the count so far,
    and the caller must perform a full-precision Euclidean step.

    Parameters:
      m     - output matrix (updated in place)
      a1:a0 - top 2 limbs of a (a1 is most significant)
      b1:b0 - top 2 limbs of b
      even  - output: 1 if even number of quotients processed, 0 if odd

    Returns: number of quotients incorporated (0 if simulation unsafe).

    Complexity: O(LIMB_BITS) iterations, each O(1), so O(LIMB_BITS) total.
    This replaces O(LIMB_BITS) expensive multi-limb divisions with one
    matrix application.  */
static inline unsigned arbint_lehmer_step(arbint_lehmer_matrix_t * m,
                                          arbint_limb_t a1, arbint_limb_t a0,
                                          arbint_limb_t b1, arbint_limb_t b0,
                                          int * even) {
  arbint_limb_t v0 = 1u, v1 = 0u;
  arbint_limb_t w0 = 0u, w1 = 1u;
  arbint_limb_t q;
  arbint_limb_t t0, t1;
  arbint_limb_t new_w0, new_w1;
  unsigned count = 0u;

  /*  Overflow threshold: stop if any coefficient exceeds this.  */
#if ARBINT_LIMB_BITS == 64
  const arbint_limb_t LIMIT = UINT64_C(0x7FFFFFFFFFFFFFFF);
#else
  const arbint_limb_t LIMIT = UINT32_C(0x7FFFFFFF);
#endif

  *even = 1;

  /*  Main simulation loop.  */
  while (b1 != 0u || b0 != 0u) {
    /*  Compute quotient q = floor(a / b) using top limbs.
        For safety, we use a conservative estimate.  */

    /*  If b1 == 0 but b0 != 0, the quotient could be very large.
        We can only safely continue if a1 == 0 too.  */
    if (b1 == 0u) {
      if (a1 != 0u)
        break; /*  Quotient too large to simulate safely.  */
      /*  Both a1 and b1 are zero; single-limb quotient.  */
      if (b0 == 0u)
        break;
      q = a0 / b0;
    } else {
      /*  Both have nonzero high limbs; quotient is at most a1/b1 + 1.  */
      q = a1 / b1;
    }

    if (q == 0u)
      break; /*  Would make no progress.  */

    /*  Check if matrix update would overflow.
        new_w = v + q*w, so check v + q*w <= LIMIT.  */
    if (w0 > 0u && q > (LIMIT - v0) / w0)
      break;
    if (w1 > 0u && q > (LIMIT - v1) / w1)
      break;

    /*  Compute new w values: w_new = v + q*w.  */
    new_w0 = v0 + q * w0;
    new_w1 = v1 + q * w1;

    /*  Update (a, b) = (b, a - q*b).  */
    t0 = b0;
    t1 = b1;

    /*  a - q*b with borrow.
        Compute q*b as a two-limb value. If q*b does not fit in two limbs,
        simulation is unsafe and we must stop (caller will do a full step).  */
    {
      arbint_limb_t prod_lo;
      arbint_limb_t prod_hi;

#if ARBINT_LIMB_BITS == 64
      arbint_limb_t q_lo = q & 0xFFFFFFFFu;
      arbint_limb_t q_hi = q >> 32;
      arbint_limb_t b_lo = b0 & 0xFFFFFFFFu;
      arbint_limb_t b_hi = b0 >> 32;
      arbint_limb_t x0 = q_lo * b_lo;
      arbint_limb_t x1 = q_lo * b_hi;
      arbint_limb_t x2 = q_hi * b_lo;
      arbint_limb_t x3 = q_hi * b_hi;
      arbint_limb_t mid = (x0 >> 32) + (x1 & 0xFFFFFFFFu) + (x2 & 0xFFFFFFFFu);
      arbint_limb_t qb1_lo;
      arbint_limb_t qb1_hi;
      arbint_limb_t y0, y1, y2, y3, ymid;
      prod_lo = (mid << 32) | (x0 & 0xFFFFFFFFu);
      prod_hi = x3 + (x1 >> 32) + (x2 >> 32) + (mid >> 32);

      /*  Compute q*b1 as a 128-bit product split into qb1_hi:qb1_lo.  */
      y0 = q_lo * (b1 & 0xFFFFFFFFu);
      y1 = q_lo * (b1 >> 32);
      y2 = q_hi * (b1 & 0xFFFFFFFFu);
      y3 = q_hi * (b1 >> 32);
      ymid = (y0 >> 32) + (y1 & 0xFFFFFFFFu) + (y2 & 0xFFFFFFFFu);
      qb1_lo = (ymid << 32) | (y0 & 0xFFFFFFFFu);
      qb1_hi = y3 + (y1 >> 32) + (y2 >> 32) + (ymid >> 32);

      /*  q*b = q*b0 + (q*b1)<<64 must fit in two limbs:
          - q*b1 must fit in one limb (qb1_hi == 0)
          - prod_hi + qb1_lo must not overflow limb.  */
      if (qb1_hi != 0u)
        break;
      if (prod_hi > ~(arbint_limb_t) 0 - qb1_lo)
        break;
      prod_hi += qb1_lo;
#else
      /*  32-bit limbs: use 64-bit arithmetic.  */
      uint64_t prod64 = (uint64_t) q * b0;
      uint64_t qb1_64;
      arbint_limb_t qb1_lo;
      prod_lo = (arbint_limb_t) prod64;
      prod_hi = (arbint_limb_t) (prod64 >> 32);

      /*  q*b = q*b0 + (q*b1)<<32 must fit in two limbs:
          - q*b1 must fit in one limb
          - prod_hi + q*b1 must not overflow limb.  */
      qb1_64 = (uint64_t) q * (uint64_t) b1;
      if ((qb1_64 >> 32) != 0u)
        break;
      qb1_lo = (arbint_limb_t) qb1_64;
      if (prod_hi > (arbint_limb_t) UINT32_MAX - qb1_lo)
        break;
      prod_hi += qb1_lo;
#endif

      /*  a - prod.  */
      if (a0 < prod_lo) {
        a0 = a0 - prod_lo; /*  Wraps, borrow into a1.  */
        if (a1 == 0u || a1 - 1u < prod_hi)
          break; /*  Underflow.  */
        a1 = a1 - 1u - prod_hi;
      } else {
        a0 = a0 - prod_lo;
        if (a1 < prod_hi)
          break; /*  Underflow.  */
        a1 = a1 - prod_hi;
      }
    }

    b0 = a0;
    b1 = a1;
    a0 = t0;
    a1 = t1;

    /*  Update matrix: (v0, v1) <- (w0, w1), (w0, w1) <- new_w.  */
    v0 = w0;
    v1 = w1;
    w0 = new_w0;
    w1 = new_w1;

    ++count;
    *even = !(*even);
  }

  m->v0 = v0;
  m->v1 = v1;
  m->w0 = w0;
  m->w1 = w1;

  return count;
}

/*  Lehmer GCD function declarations (implemented in platform-specific
    files).  */
arbint_err_t arbint_gcd_lehmer_generic(arbint_t g, const arbint_limb_t * ap,
                                       size_t an, const arbint_limb_t * bp,
                                       size_t bn,
                                       const arbint_alloc_t * alloc);

#if HAS_BMI2
arbint_err_t arbint_gcd_lehmer_bmi2(arbint_t g, const arbint_limb_t * ap,
                                    size_t an, const arbint_limb_t * bp,
                                    size_t bn, const arbint_alloc_t * alloc);
#endif

#endif /* ARBINT_GCD_H */
