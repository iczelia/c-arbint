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

#ifndef ARBINT_DIV_NEWTON_H
#define ARBINT_DIV_NEWTON_H

#include "arbint_base.h"

/*  Newton-Raphson Division Thresholds.

    These thresholds control when the Newton-Raphson division algorithm is
    used instead of Knuth's Algorithm D.  Tunable via compile-time defines.  */

/*  Minimum divisor size (in limbs) to use reciprocal-based division.
    Below this threshold, Knuth's Algorithm D is used.

    For divisors above this threshold, reciprocal-based division is used.
    The reciprocal is computed via Newton iteration (for large divisors)
    or direct division (for moderate divisors).  The crossover is tuned
    for where O(M(n)) complexity beats O(n^2) schoolbook.  */
#ifndef ARBINT_NEWTON_DIV_THRESHOLD
#define ARBINT_NEWTON_DIV_THRESHOLD 256u
#endif

/*  Minimum divisor size for divide-and-conquer inverse computation.
    Below this, direct division is used for the reciprocal.

    For balanced divisions, Newton inversion is preferred above this size.
    The implementation verifies Newton output and automatically falls back
    to direct inversion when precision requirements exceed the validated
    Newton range.  */
#ifndef ARBINT_INVERT_DC_THRESHOLD
#define ARBINT_INVERT_DC_THRESHOLD 256u
#endif

/*  Number of guard limbs for reciprocal computation.

    The reciprocal v = floor(2^p / d) is computed with p = (dn + g) * BITS
    where g = ARBINT_INVERT_GUARD_LIMBS.  With g = 2, the quotient
    approximation error is bounded such that at most 2 correction iterations
    are needed.

    Mathematical bound:
      Let v = floor(2^p / d), so 2^p = v*d + r where 0 <= r < d.
      q_approx = floor(n * v / 2^p)
      Error: |q - q_approx| <= ceil(n*r / (d*2^p)) + 1 < 2
      when p = (dn + 2) * BITS and n < 2^(nn * BITS).  */
#ifndef ARBINT_INVERT_GUARD_LIMBS
#define ARBINT_INVERT_GUARD_LIMBS 2u
#endif

/*  Lookup table for initial 1/x approximation.

    Entry i gives floor(2^16 / (2^8 + i)) - 2^8 for i in 0..127.
    The input d is normalized (MSB set), so top 8 bits are in [128..255].
    Index by (top_byte - 128), add 256 to get ~9-bit approximation of
    2^16 / top_byte.

    This table is analogous to arbint_invsqrt_tab used in isqrt.  */
extern const unsigned char arbint_inv_tab[128];

/*  Compute the reciprocal v = floor(2^p / d).

    Preconditions:
      - d is normalized: d[dn-1] has MSB set (i.e., d >= 2^((dn-1)*BITS + BITS-1))
      - dn >= 1
      - v has capacity for at least (dn + ARBINT_INVERT_GUARD_LIMBS + 1) limbs

    Postconditions:
      - v * d <= 2^p < (v+1) * d, where p = (dn + ARBINT_INVERT_GUARD_LIMBS) * BITS
      - *v_n contains the actual number of significant limbs in v

    Algorithm:
      For dn >= ARBINT_INVERT_DC_THRESHOLD, uses Newton iteration:
        x_{k+1} = x_k * (2 - d*x_k / 2^p).
      Starts from a single-limb inverse computed via lookup table refinement.
      Each Newton step doubles precision.  Complexity: O(M(n)).

      For smaller dn, uses direct division via Knuth Algorithm D.
      This has O(n^2) complexity but lower constant factors.  */
arbint_err_t arbint_invert_newton(arbint_limb_t * v, size_t * v_n,
                                   const arbint_limb_t * d, size_t dn,
                                   const arbint_alloc_t * alloc);

/*  Compute the reciprocal v = floor(2^p / d) using divide-and-conquer.

    Same interface as arbint_invert_newton, but uses recursive structure
    for O(M(n)) complexity (compared to O(M(n) log n) for iterative).

    Falls back to arbint_invert_newton for dn <= ARBINT_INVERT_DC_THRESHOLD.

    The scratch buffer must have at least arbint_invert_dc_scratch_size(dn)
    limbs of capacity.  */
arbint_err_t arbint_invert_dc(arbint_limb_t * v, size_t * v_n,
                               const arbint_limb_t * d, size_t dn,
                               arbint_limb_t * scratch, size_t scratch_size,
                               const arbint_alloc_t * alloc);

/*  Compute scratch size needed for arbint_invert_dc.  */
size_t arbint_invert_dc_scratch_size(size_t dn);

/*  Newton-Raphson division: compute q = floor(n/d) and r = n mod d.

    Preconditions:
      - np points to nn limbs of dividend (little-endian)
      - dp points to dn limbs of divisor (little-endian)
      - dp[dn-1] != 0 (divisor is normalized to actual limb count)
      - dn >= 2 (single-limb divisors should use Barrett reduction)
      - qp has capacity for (nn - dn + 1) limbs if non-NULL
      - rp has capacity for dn limbs (required)

    Postconditions:
      - n = q * d + r with 0 <= r < d
      - If qp != NULL: quotient written to qp[0..nn-dn]
      - Remainder written to rp[0..dn-1]

    Algorithm:
      1. Normalize divisor (shift left so MSB set)
      2. Compute reciprocal v = floor(2^p / d_norm)
      3. q_approx = floor((n_shifted * v) / 2^p)
      4. r = n - q_approx * d
      5. Correct q_approx (at most a few iterations with guard limbs)
      6. Un-shift remainder

    Complexity: O(n*m) for reciprocal + O(M(n)) for multiplication.
    For very large divisors where the O(M(n)) property matters, this
    amortizes well when the same divisor is used for multiple divisions.  */
arbint_err_t arbint_div_mag_newton(const arbint_limb_t * np, size_t nn,
                                    const arbint_limb_t * dp, size_t dn,
                                    arbint_limb_t * qp, arbint_limb_t * rp,
                                    const arbint_alloc_t * alloc);

/*  Compute scratch size needed for arbint_div_mag_newton.  */
size_t arbint_div_newton_scratch_size(size_t nn, size_t dn);

#endif /*  ARBINT_DIV_NEWTON_H  */
