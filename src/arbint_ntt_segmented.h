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

#ifndef ARBINT_NTT_SEGMENTED_H
#define ARBINT_NTT_SEGMENTED_H

#include "arbint_base.h"
#include "arbint_ntt.h"

/*  Segmented NTT multiplication for operands exceeding ARBINT_NTT_MAX_SIZE.

    When operand sizes exceed the NTT limit (2^24 limbs), we split them into
    segments and use Karatsuba-style decomposition:

      A = A_hi * B^k + A_lo
      B = B_hi * B^k + B_lo

      A * B = A_hi*B_hi * B^(2k)
            + (A_hi*B_lo + A_lo*B_hi) * B^k
            + A_lo*B_lo

    Each sub-product fits within NTT limits. For three-way splitting:

      A = A2 * B^(2k) + A1 * B^k + A0  (Toom-3 style)

    The segment size k is chosen as ARBINT_NTT_MAX_SIZE / 2 to ensure
    that each sub-product (up to 2k limbs) fits within NTT capacity.

    This approach allows multiplication of arbitrarily large numbers
    while leveraging O(n log n) NTT for the bulk of the computation.  */

/*  Segment size for splitting operands.
    Each sub-product is at most 2 * SEGMENT_SIZE limbs, which must fit
    within ARBINT_NTT_MAX_SIZE. We use MAX_SIZE/2 to be safe.  */
#define ARBINT_NTT_SEGMENT_SIZE (ARBINT_NTT_MAX_SIZE / 2u)

/*  Segmented NTT multiplication: dst = a * b for very large operands.

    Splits operands into segments of ARBINT_NTT_SEGMENT_SIZE limbs and
    uses Karatsuba decomposition with NTT for sub-products.

    Preconditions:
    - an >= bn (caller swaps if needed)
    - At least one of an, bn exceeds ARBINT_NTT_MAX_SIZE
    - dst must have capacity for an + bn limbs

    Returns result limb count via out_used.

    Complexity: O(n log n) for the dominant NTT operations, with
    O(n^1.58) overhead from Karatsuba decomposition levels.  */
arbint_err_t arbint_mul_mag_ntt_segmented(arbint_limb_t * dst, size_t * out_used,
                                          const arbint_limb_t * a, size_t an,
                                          const arbint_limb_t * b, size_t bn,
                                          const arbint_alloc_t * alloc);

/*  Segmented NTT squaring: dst = a^2 for very large operands.

    Splits operand into segments and uses Karatsuba decomposition:
      A = A_hi * B^k + A_lo
      A^2 = A_hi^2 * B^(2k) + 2*A_hi*A_lo * B^k + A_lo^2

    Using the identity for squaring:
      (A_lo + A_hi)^2 = A_lo^2 + A_hi^2 + 2*A_lo*A_hi
      => 2*A_lo*A_hi = (A_lo + A_hi)^2 - A_lo^2 - A_hi^2

    All three sub-products are squares (not general multiplications),
    so the squaring optimization applies recursively at every level.

    Preconditions:
    - 2*an > ARBINT_NTT_MAX_SIZE (requires segmentation)
    - dst must have capacity for 2*an + 1 limbs

    Returns result limb count via out_used.  */
arbint_err_t arbint_sqr_mag_ntt_segmented(arbint_limb_t * dst, size_t * out_used,
                                          const arbint_limb_t * a, size_t an,
                                          const arbint_alloc_t * alloc);

/*  Check if operands require segmented multiplication.
    Returns 1 if convolution length an + bn exceeds ARBINT_NTT_MAX_SIZE.  */
static inline int arbint_needs_segmented_ntt(size_t an, size_t bn) {
  return (an + bn > ARBINT_NTT_MAX_SIZE) ? 1 : 0;
}

/*  Check if operand requires segmented squaring.
    Returns 1 if convolution length 2*an exceeds ARBINT_NTT_MAX_SIZE.  */
static inline int arbint_needs_segmented_ntt_sqr(size_t an) {
  return (2u * an > ARBINT_NTT_MAX_SIZE) ? 1 : 0;
}

#endif /*  ARBINT_NTT_SEGMENTED_H  */
