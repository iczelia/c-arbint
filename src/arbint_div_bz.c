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

/*  Burnikel-Ziegler Divide-and-Conquer Division.

    This module implements the Burnikel-Ziegler algorithm for dividing
    an n-limb dividend by an m-limb divisor where m is in the range
    where Knuth's Algorithm D is suboptimal but Newton-Raphson is overkill.

    Algorithm complexity: O(M(n) * log(n/m)) where M(n) is multiplication cost.
    For Karatsuba multiplication, this achieves O(n^1.58) overall.

    References:
    - C. Burnikel and J. Ziegler, "Fast Recursive Division", MPI-I-98-1-022  */

#include "arbint_div_bz.h"

#include "arbint_addsub.h"
#include "arbint_base.h"
#include "arbint_cmp.h"
#include "arbint_div.h"
#include "arbint_mul.h"

#include <string.h>

/*  Subtraction with borrow detection: dst = a - b.
    Returns 1 if result is negative (borrow out), 0 otherwise.
    an must be >= bn. dst may alias a but not b.  */
static int arbint_sub_with_borrow(arbint_limb_t * dst,
                                   const arbint_limb_t * a, size_t an,
                                   const arbint_limb_t * b, size_t bn) {
  arbint_limb_t borrow = 0u;
  size_t i;

  for (i = 0u; i < bn; ++i) {
    arbint_limb_t ai = a[i];
    arbint_limb_t bi = b[i];
    arbint_limb_t diff = ai - bi - borrow;
    borrow = (ai < bi || (ai == bi && borrow)) ? 1u : 0u;
    dst[i] = diff;
  }
  for (; i < an; ++i) {
    arbint_limb_t ai = a[i];
    arbint_limb_t diff = ai - borrow;
    borrow = (ai < borrow) ? 1u : 0u;
    dst[i] = diff;
  }
  return (int) borrow;
}

/*  Addition: dst = a + b, returns carry.  */
static arbint_limb_t arbint_add_n(arbint_limb_t * dst,
                                   const arbint_limb_t * a,
                                   const arbint_limb_t * b, size_t n) {
  arbint_limb_t carry = 0u;
  size_t i;

  for (i = 0u; i < n; ++i) {
    arbint_limb_t sum = a[i] + b[i] + carry;
    carry = (sum < a[i] || (carry && sum == a[i])) ? 1u : 0u;
    dst[i] = sum;
  }
  return carry;
}

/*  Compute scratch size needed for BZ division.  */
size_t arbint_div_bz_scratch_size(size_t nn, size_t dn) {
  /*  Conservative estimate for recursive calls.
      Each level needs about 8*dn limbs, recursion depth is O(log(dn)).
      For unbalanced, we also need space for iterative processing.  */
  size_t n = (nn > dn) ? nn : dn;
  if (n > SIZE_MAX / 20u)
    return SIZE_MAX;
  return 20u * n + 64u;
}

/*  Forward declaration.  */
static arbint_err_t arbint_div_sb(arbint_limb_t * qp, arbint_limb_t * rp,
                                   const arbint_limb_t * ap, size_t an,
                                   const arbint_limb_t * bp, size_t bn);

/*  DIV_SB: Schoolbook division wrapper.
    Handles the case where the divisor's top limb might be zero due to padding.
    Finds the actual size and delegates to Knuth.  */
static arbint_err_t arbint_div_sb(arbint_limb_t * qp, arbint_limb_t * rp,
                                   const arbint_limb_t * ap, size_t an,
                                   const arbint_limb_t * bp, size_t bn) {
  size_t actual_bn = bn;
  size_t actual_an = an;
  size_t qn = (an >= bn) ? an - bn + 1u : 1u;
  size_t i;

  /*  Always clear full output capacities.
      Knuth writes only used quotient/remainder limbs. Callers expect the
      remaining higher limbs to be zero, especially in recursive BZ steps.  */
  if (qp != NULL)
    memset(qp, 0, qn * sizeof(arbint_limb_t));
  memset(rp, 0, bn * sizeof(arbint_limb_t));

  /*  Normalize: find actual size of divisor.  */
  while (actual_bn > 0u && bp[actual_bn - 1u] == 0u)
    actual_bn--;

  if (actual_bn == 0u)
    return ARBINT_EZERO;

  /*  Normalize: find actual size of dividend.  */
  while (actual_an > 0u && ap[actual_an - 1u] == 0u)
    actual_an--;

  if (actual_an == 0u) {
    /*  Dividend is zero: q = 0, r = 0.  */
    return ARBINT_OK;
  }

  if (actual_an < actual_bn) {
    /*  Dividend < divisor: q = 0, r = dividend.  */
    for (i = 0u; i < actual_an; ++i)
      rp[i] = ap[i];
    return ARBINT_OK;
  }

  /*  Call Knuth with normalized sizes.  */
  return arbint_div_mag_knuth(ap, actual_an, bp, actual_bn, qp, rp);
}

/*  DIV3N2N: Core 3n/2-by-n division subroutine.

    Divides a (3*half_n)-limb number by an n-limb number (where n = 2*half_n),
    producing a half_n-limb quotient and n-limb remainder.

    Preconditions:
      - n is even (n = 2 * half_n)
      - ap has 3*half_n limbs (A = A2:A1:A0, each half_n limbs)
      - bp has n limbs (B = B1:B0, each half_n limbs)
      - B is normalized (bp[n-1] != 0)
      - A < B * beta^(half_n), i.e., the quotient fits in half_n limbs

    Outputs:
      - qp receives half_n limbs of quotient
      - rp receives n limbs of remainder  */
static arbint_err_t arbint_div_3n2n(arbint_limb_t * qp, arbint_limb_t * rp,
                                     const arbint_limb_t * ap,
                                     const arbint_limb_t * bp, size_t n,
                                     arbint_limb_t * scratch,
                                     const arbint_alloc_t * alloc) {
  size_t half_n;
  const arbint_limb_t * b1;
  const arbint_limb_t * b0;
  const arbint_limb_t * a2;
  const arbint_limb_t * a1;
  const arbint_limb_t * a0;
  arbint_limb_t * q_tmp;
  arbint_limb_t * r_high;
  arbint_limb_t * d_prod;
  arbint_limb_t * r_full;
  size_t d_used;
  int borrow;
  arbint_err_t rc;
  size_t i;

  half_n = n / 2u;

  /*  B = B1:B0, each half_n limbs.  */
  b0 = bp;
  b1 = bp + half_n;

  /*  A = A2:A1:A0, each half_n limbs.  */
  a0 = ap;
  a1 = ap + half_n;
  a2 = ap + n;

  /*  Workspace layout:
      q_tmp:  half_n + 1 limbs (temporary quotient - may overflow)
      r_high: n limbs (R' from recursive division - may be larger in overflow case)
      d_prod: n + 2 limbs (Q * B0 product)
      r_full: n + half_n + 2 limbs (working remainder - for (R':A0) in overflow case)  */
  q_tmp = scratch;
  r_high = scratch + half_n + 1u;
  d_prod = r_high + n;
  r_full = d_prod + n + 2u;

  /*  Step 1: Check if A2 >= B1.
      If A2 >= B1, the quotient might overflow half_n limbs.
      In this case, set Q' = beta^half_n - 1 (maximum value).  */
  if (arbint_cmp_mag_limbs(a2, half_n, b1, half_n) >= 0) {
    /*  Q' = all ones.  */
    for (i = 0u; i < half_n; ++i)
      q_tmp[i] = ~(arbint_limb_t) 0u;

    /*  Compute R' = (A2:A1) - Q' * B1.
        Since Q' = beta^half_n - 1:
        Q' * B1 = B1 * beta^half_n - B1
        R' = A2*beta^half_n + A1 - B1*beta^half_n + B1
           = (A2 - B1)*beta^half_n + (A1 + B1)

        R' could be larger than half_n limbs in this case.
        Store the full result in r_high (which now has n limbs).  */
    rc = arbint_mul_mag_generic(d_prod, &d_used, q_tmp, half_n, b1, half_n,
                                 alloc);
    if (rc != ARBINT_OK)
      return rc;

    /*  r_high = A2:A1 (n limbs), then subtract d_prod.  */
    memcpy(r_high, a1, half_n * sizeof(arbint_limb_t));
    memcpy(r_high + half_n, a2, half_n * sizeof(arbint_limb_t));
    (void) arbint_sub_with_borrow(r_high, r_high, n, d_prod, d_used);
  } else {
    /*  Step 2: Normal case - (Q', R') = DIV(A2:A1, B1).
        This divides n limbs by half_n limbs, producing half_n-limb Q'
        and half_n-limb R'.

        Zero r_high first because arbint_div_sb may write fewer than half_n
        limbs if the divisor b1 has leading zeros (actual_bn < half_n).
        Without this, positions actual_bn..n-1 would contain garbage.  */
    memset(r_high, 0, n * sizeof(arbint_limb_t));
    rc = arbint_div_sb(q_tmp, r_high, a1, n, b1, half_n);
    if (rc != ARBINT_OK)
      return rc;
  }

  /*  Step 3: D = Q' * B0.  */
  rc = arbint_mul_mag_generic(d_prod, &d_used, q_tmp, half_n, b0, half_n,
                               alloc);
  if (rc != ARBINT_OK)
    return rc;

  /*  Step 4: R = (R':A0) - D.
      Construct R':A0 in r_full:
      - Low half_n limbs = A0
      - High part = R' (half_n limbs in normal case, n limbs in overflow case)

      For the overflow case, (R':A0) can be up to n + half_n limbs.
      r_full has n + half_n + 2 limbs allocated, which is sufficient.  */
  memcpy(r_full, a0, half_n * sizeof(arbint_limb_t));
  memcpy(r_full + half_n, r_high, n * sizeof(arbint_limb_t));
  /*  Clear any extra limbs to ensure clean state.  */
  r_full[n + half_n] = 0u;

  /*  Subtraction: R = (R':A0) - D.
      In the overflow case, (R':A0) can have up to n + half_n limbs.
      We need to subtract D and propagate borrow through all limbs.  */
  borrow = arbint_sub_with_borrow(r_full, r_full, n + half_n, d_prod, d_used);

  /*  Step 5: Correction loop - while R < 0: Q'--, R += B.
      The correction is needed at most twice.  */
  while (borrow) {
    /*  Q' -= 1.  */
    arbint_limb_t c = 1u;
    for (i = 0u; i < half_n && c != 0u; ++i) {
      arbint_limb_t old = q_tmp[i];
      q_tmp[i] = old - c;
      c = (q_tmp[i] > old) ? 1u : 0u;
    }

    /*  R += B.  Add B (n limbs) to r_full and propagate carry.  */
    arbint_limb_t carry = arbint_add_n(r_full, r_full, bp, n);
    /*  Propagate carry through remaining limbs.  */
    for (i = n; i < n + half_n && carry != 0u; ++i) {
      arbint_limb_t old = r_full[i];
      r_full[i] = old + carry;
      carry = (r_full[i] < old) ? 1u : 0u;
    }
    /*  If carry occurred (wrapped around), R became positive.  */
    borrow = (carry == 0u) ? 1 : 0;
    if (carry != 0u)
      borrow = 0;
  }

  /*  Output results.  */
  memcpy(qp, q_tmp, half_n * sizeof(arbint_limb_t));
  memcpy(rp, r_full, n * sizeof(arbint_limb_t));

  return ARBINT_OK;
}

/*  DIV2N1N: Balanced 2n-by-n division.

    Divides a 2n-limb number by an n-limb number using recursive 3n/2n division.

    Preconditions:
      - ap has 2n limbs
      - bp has n limbs, with bp[n-1] != 0
      - n >= 2

    Outputs:
      - qp receives n limbs of quotient
      - rp receives n limbs of remainder  */
static arbint_err_t arbint_div_2n1n(arbint_limb_t * qp, arbint_limb_t * rp,
                                     const arbint_limb_t * ap,
                                     const arbint_limb_t * bp, size_t n,
                                     arbint_limb_t * scratch,
                                     const arbint_alloc_t * alloc) {
  size_t half_n;
  arbint_limb_t * q_hi;
  arbint_limb_t * r_tmp;
  arbint_limb_t * a_concat;
  arbint_err_t rc;
  size_t i;

  /*  Base case: use Knuth for small n.  */
  if (n < ARBINT_BZ_THRESHOLD)
    return arbint_div_sb(qp, rp, ap, 2u * n, bp, n);

  /*  Handle odd n by rounding up.
      Use a limb-shift transform:
        A' = A * B, D' = D * B, with B = limb radix.
      Then A' / D' has the same quotient, and remainder is shifted by one limb.
      This converts odd n to even n+1 without losing normalization.  */
  if (n % 2u != 0u) {
    size_t ne = n + 1u;
    size_t ap2_n = 2u * ne;
    size_t bp2_n = ne;
    size_t q2_n = ne;
    size_t r2_n = ne;
    size_t work_n;
    size_t total_n;
    arbint_limb_t * buf = NULL;
    arbint_limb_t * ap2 = NULL;
    arbint_limb_t * bp2 = NULL;
    arbint_limb_t * q2 = NULL;
    arbint_limb_t * r2 = NULL;
    arbint_limb_t * work2 = NULL;
    arbint_err_t rc_odd;

    work_n = arbint_div_bz_scratch_size(ap2_n, bp2_n);
    if (work_n == SIZE_MAX)
      return ARBINT_EOVERFLOW;
    if (ap2_n > SIZE_MAX - bp2_n || ap2_n + bp2_n > SIZE_MAX - q2_n)
      return ARBINT_EOVERFLOW;
    total_n = ap2_n + bp2_n + q2_n;
    if (total_n > SIZE_MAX - r2_n || total_n + r2_n > SIZE_MAX - work_n)
      return ARBINT_EOVERFLOW;
    total_n += r2_n + work_n;

    buf = arbint_alloc_limbs(alloc, total_n);
    if (buf == NULL)
      return ARBINT_ENOMEM;

    ap2 = buf;
    bp2 = ap2 + ap2_n;
    q2 = bp2 + bp2_n;
    r2 = q2 + q2_n;
    work2 = r2 + r2_n;

    /*  A' = A * B, D' = D * B (little-endian limb shift).  */
    memset(ap2, 0, ap2_n * sizeof(arbint_limb_t));
    memcpy(ap2 + 1u, ap, (2u * n) * sizeof(arbint_limb_t));

    memset(bp2, 0, bp2_n * sizeof(arbint_limb_t));
    memcpy(bp2 + 1u, bp, n * sizeof(arbint_limb_t));

    rc_odd = arbint_div_2n1n(q2, r2, ap2, bp2, ne, work2, alloc);
    if (rc_odd != ARBINT_OK) {
      arbint_free_limbs(alloc, buf);
      return rc_odd;
    }

    memcpy(qp, q2, n * sizeof(arbint_limb_t));
    memcpy(rp, r2 + 1u, n * sizeof(arbint_limb_t));

    arbint_free_limbs(alloc, buf);
    return ARBINT_OK;
  }

  half_n = n / 2u;

  /*  Workspace layout:
      q_hi:     half_n limbs (high quotient digit)
      r_tmp:    n limbs (intermediate remainder)
      a_concat: 3*half_n limbs (for DIV3N2N input)
      rest:     workspace for recursive calls  */
  q_hi = scratch;
  r_tmp = scratch + half_n;
  a_concat = scratch + half_n + n;

  /*  Step 1: High quotient digit.
      (Q1, R1) = DIV3N2N(A[half_n..2n-1], B).
      A[half_n..2n-1] is 3*half_n limbs: A1:A2:A3 where A = A0:A1:A2:A3.
      This produces half_n-limb Q1 and n-limb R1.  */
  rc = arbint_div_3n2n(q_hi, r_tmp, ap + half_n, bp, n,
                        a_concat + 3u * half_n + 1u, alloc);
  if (rc != ARBINT_OK)
    return rc;

  /*  Step 2: Construct dividend for low quotient.
      We need a 3*half_n limb number: A0:R1[0..n-1].
      - a_concat[0..half_n-1] = A0 (lowest half_n limbs of A)
      - a_concat[half_n..3*half_n-1] = R1 (n limbs)  */
  memcpy(a_concat, ap, half_n * sizeof(arbint_limb_t));
  memcpy(a_concat + half_n, r_tmp, n * sizeof(arbint_limb_t));

  /*  Step 3: Low quotient digit.
      (Q0, R) = DIV3N2N(a_concat, B).
      This produces half_n-limb Q0 and n-limb R.  */
  rc = arbint_div_3n2n(qp, rp, a_concat, bp, n,
                        a_concat + 3u * half_n + 1u, alloc);
  if (rc != ARBINT_OK)
    return rc;

  /*  Step 4: Combine quotient.
      Q = Q1 * beta^half_n + Q0.
      Q0 is already in qp[0..half_n-1].
      Put Q1 in qp[half_n..n-1].  */
  for (i = 0u; i < half_n; ++i)
    qp[half_n + i] = q_hi[i];

  return ARBINT_OK;
}

/*  Main entry point: Burnikel-Ziegler division.  */
arbint_err_t arbint_div_mag_bz(const arbint_limb_t * np, size_t nn,
                                const arbint_limb_t * dp, size_t dn,
                                arbint_limb_t * qp, arbint_limb_t * rp,
                                const arbint_alloc_t * alloc) {
  arbint_limb_t * scratch = NULL;
  arbint_limb_t * d_norm = NULL;
  arbint_limb_t * n_norm = NULL;
  arbint_limb_t * q_tmp = NULL;
  arbint_limb_t * r_tmp = NULL;
  arbint_limb_t * work = NULL;
  unsigned shift;
  size_t scratch_size;
  size_t i;
  size_t qn;
  size_t n_effective;
  arbint_err_t rc = ARBINT_OK;

  /*  Validate inputs.  */
  if (np == NULL || dp == NULL || rp == NULL)
    return ARBINT_EINVAL;
  if (dn < 2u || nn < dn)
    return ARBINT_EINVAL;
  if (dp[dn - 1u] == 0u)
    return ARBINT_EINVAL;

  /*  Fall back to Knuth for small divisors.  */
  if (dn < ARBINT_BZ_THRESHOLD)
    return arbint_div_mag_knuth(np, nn, dp, dn, qp, rp);

  /*  Allocate scratch space.  */
  scratch_size = arbint_div_bz_scratch_size(nn, dn);
  scratch = arbint_alloc_limbs(alloc, scratch_size);
  if (scratch == NULL)
    return ARBINT_ENOMEM;

  /*  Workspace layout.  */
  d_norm = scratch;
  n_norm = d_norm + dn + 1u;
  q_tmp = n_norm + nn + 2u;
  r_tmp = q_tmp + nn + 2u;
  work = r_tmp + dn + 2u;

  /*  Normalize divisor (shift so MSB is set).  */
  shift = arbint_clz_limb(dp[dn - 1u]);
  if (shift != 0u) {
    arbint_limb_t carry = 0u;
    for (i = 0u; i < dn; ++i) {
      d_norm[i] = (dp[i] << shift) | carry;
      carry = dp[i] >> (ARBINT_LIMB_BITS - shift);
    }
  } else {
    memcpy(d_norm, dp, dn * sizeof(arbint_limb_t));
  }

  /*  Normalize dividend.  */
  if (shift != 0u) {
    arbint_limb_t carry = 0u;
    for (i = 0u; i < nn; ++i) {
      n_norm[i] = (np[i] << shift) | carry;
      carry = np[i] >> (ARBINT_LIMB_BITS - shift);
    }
    n_norm[nn] = carry;
    n_effective = (carry != 0u) ? nn + 1u : nn;
  } else {
    memcpy(n_norm, np, nn * sizeof(arbint_limb_t));
    n_norm[nn] = 0u;
    n_effective = nn;
  }

  qn = nn - dn + 1u;

  /*  Initialize outputs.  */
  for (i = 0u; i < qn; ++i)
    q_tmp[i] = 0u;
  for (i = 0u; i < dn; ++i)
    r_tmp[i] = 0u;

  /*  For balanced case (n_effective <= 2*dn), use single DIV2N1N.  */
  if (n_effective <= 2u * dn) {
    /*  Pad dividend to exactly 2*dn limbs.  */
    for (i = n_effective; i < 2u * dn; ++i)
      n_norm[i] = 0u;

    rc = arbint_div_2n1n(q_tmp, r_tmp, n_norm, d_norm, dn, work, alloc);
    if (rc != ARBINT_OK)
      goto cleanup;
  } else {
    /*  Unbalanced division: n_effective > 2*dn.
        Fall back to Knuth which handles this correctly.  */
    arbint_free_limbs(alloc, scratch);
    return arbint_div_mag_knuth(np, nn, dp, dn, qp, rp);
  }

  /*  Un-normalize remainder.  */
  if (shift != 0u) {
    for (i = 0u; i < dn - 1u; ++i)
      rp[i] = (r_tmp[i] >> shift) | (r_tmp[i + 1u] << (ARBINT_LIMB_BITS - shift));
    rp[dn - 1u] = r_tmp[dn - 1u] >> shift;
  } else {
    memcpy(rp, r_tmp, dn * sizeof(arbint_limb_t));
  }

  /*  Copy quotient.  */
  if (qp != NULL)
    memcpy(qp, q_tmp, qn * sizeof(arbint_limb_t));

cleanup:
  arbint_free_limbs(alloc, scratch);
  return rc;
}
