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

/*  Segmented NTT multiplication for operands beyond 2^24 limbs.

    Algorithm overview:
    ------------------
    When operands exceed NTT capacity, we split them into segments and
    recursively multiply. The decomposition follows Karatsuba's identity:

    For two-way split (k = segment size):
      A = A1 * B^k + A0
      B = B1 * B^k + B0

      A * B = A1*B1 * B^(2k) + (A1*B0 + A0*B1) * B^k + A0*B0

    We can optimize using Karatsuba's trick:
      A1*B0 + A0*B1 = (A0+A1)*(B0+B1) - A0*B0 - A1*B1

    This reduces 4 multiplications to 3.

    For operands requiring more than 2 segments, we recursively apply
    this decomposition until each sub-product fits within NTT limits.

    Memory layout:
    -------------
    Workspace allocation is carefully managed to minimize allocations.
    We preallocate space for all intermediate products and accumulate
    results directly into the output buffer where possible.

    Error handling:
    --------------
    All allocation failures are propagated. The function uses goto-based
    cleanup to ensure resources are freed on all exit paths.  */

#include "arbint_ntt_segmented.h"
#include "arbint_addsub.h"

#include <string.h>

/*  Forward declaration for recursive multiplication.
    This allows sub-products to use either NTT or segmented multiplication
    based on operand size.  */
static arbint_err_t arbint_mul_mag_segmented_rec(arbint_limb_t * dst,
                                                  size_t * out_used,
                                                  const arbint_limb_t * a,
                                                  size_t an,
                                                  const arbint_limb_t * b,
                                                  size_t bn,
                                                  const arbint_alloc_t * alloc);

/*  Add src to dst at offset, handling carries.
    dst[offset..] += src[0..src_n-1]
    dst must have capacity for offset + src_n + 1 limbs.
    Returns new used limb count.  */
static size_t add_at_offset(arbint_limb_t * dst, size_t dst_n, size_t dst_cap,
                            const arbint_limb_t * src, size_t src_n,
                            size_t offset) {
  size_t i;
  arbint_limb_t carry = 0u;
  size_t new_n;

  if (src_n == 0u)
    return dst_n;

  /*  Zero-extend dst if needed.  */
  for (i = dst_n; i < offset + src_n + 1u && i < dst_cap; ++i)
    dst[i] = 0u;

  /*  Add with carry propagation.  */
  for (i = 0u; i < src_n; ++i) {
    size_t pos = offset + i;
    arbint_limb_t d;
    arbint_limb_t s;
    arbint_limb_t sum;
    arbint_limb_t c1;
    arbint_limb_t sum2;
    arbint_limb_t c2;

    if (pos >= dst_cap)
      break;

    d = dst[pos];
    s = src[i];
    sum = d + s;
    c1 = (sum < d) ? 1u : 0u;
    sum2 = sum + carry;
    c2 = (sum2 < sum) ? 1u : 0u;
    dst[pos] = sum2;
    carry = c1 | c2;
  }

  /*  Propagate remaining carry.  */
  i = offset + src_n;
  while (carry != 0u && i < dst_cap) {
    arbint_limb_t sum = dst[i] + carry;
    carry = (sum < dst[i]) ? 1u : 0u;
    dst[i] = sum;
    ++i;
  }

  new_n = (offset + src_n > dst_n) ? offset + src_n : dst_n;
  if (i > new_n)
    new_n = i;
  if (new_n > dst_cap)
    new_n = dst_cap;

  return arbint_norm_used(dst, new_n);
}

/*  Karatsuba-style segmented multiplication.

    Splits operands at k = ARBINT_NTT_SEGMENT_SIZE:
      A = A1 * B^k + A0
      B = B1 * B^k + B0

    Computes:
      z0 = A0 * B0
      z2 = A1 * B1
      z1 = (A0 + A1) * (B0 + B1) - z0 - z2

    Result: z2 * B^(2k) + z1 * B^k + z0  */
static arbint_err_t karatsuba_segment(arbint_limb_t * dst, size_t * out_used,
                                       const arbint_limb_t * a, size_t an,
                                       const arbint_limb_t * b, size_t bn,
                                       size_t k,
                                       const arbint_alloc_t * alloc) {
  arbint_err_t rc = ARBINT_OK;
  arbint_limb_t * work = NULL;
  arbint_limb_t * z0;
  arbint_limb_t * z2;
  arbint_limb_t * z1;
  arbint_limb_t * sa;
  arbint_limb_t * sb;
  size_t z0_cap;
  size_t z2_cap;
  size_t z1_cap;
  size_t sa_cap;
  size_t sb_cap;
  size_t total;
  size_t a0n;
  size_t a1n;
  size_t b0n;
  size_t b1n;
  const arbint_limb_t * a0;
  const arbint_limb_t * a1;
  const arbint_limb_t * b0;
  const arbint_limb_t * b1;
  size_t z0n = 0u;
  size_t z2n = 0u;
  size_t z1n = 0u;
  size_t san;
  size_t sbn;
  size_t cap;
  size_t used;

  /*  Split operands.  */
  a0 = a;
  a0n = (k < an) ? k : an;
  a0n = arbint_norm_used(a0, a0n);

  if (an > k) {
    a1 = a + k;
    a1n = an - k;
    a1n = arbint_norm_used(a1, a1n);
  } else {
    a1 = a;
    a1n = 0u;
  }

  b0 = b;
  b0n = (k < bn) ? k : bn;
  b0n = arbint_norm_used(b0, b0n);

  if (bn > k) {
    b1 = b + k;
    b1n = bn - k;
    b1n = arbint_norm_used(b1, b1n);
  } else {
    b1 = b;
    b1n = 0u;
  }

  /*  Handle degenerate cases where one half is zero.  */
  if (a1n == 0u && b1n == 0u) {
    /*  Both are single-segment: just multiply.  */
    return arbint_mul_mag_segmented_rec(dst, out_used, a0, a0n, b0, b0n, alloc);
  }

  if (a1n == 0u) {
    /*  A is single segment: A * B = A0 * (B1 * B^k + B0)
                                   = A0 * B1 * B^k + A0 * B0  */
    cap = an + bn + 1u;
    memset(dst, 0u, cap * sizeof(arbint_limb_t));

    rc = arbint_mul_mag_segmented_rec(dst, &z0n, a0, a0n, b0, b0n, alloc);
    if (rc != ARBINT_OK)
      return rc;

    z2_cap = a0n + b1n + 1u;
    z2 = arbint_alloc_limbs(alloc, z2_cap);
    if (z2 == NULL)
      return ARBINT_ENOMEM;

    rc = arbint_mul_mag_segmented_rec(z2, &z2n, a0, a0n, b1, b1n, alloc);
    if (rc != ARBINT_OK) {
      arbint_free_limbs(alloc, z2);
      return rc;
    }

    used = add_at_offset(dst, z0n, cap, z2, z2n, k);
    arbint_free_limbs(alloc, z2);
    *out_used = arbint_norm_used(dst, used);
    return ARBINT_OK;
  }

  if (b1n == 0u) {
    /*  B is single segment: A * B = (A1 * B^k + A0) * B0
                                   = A1 * B0 * B^k + A0 * B0  */
    cap = an + bn + 1u;
    memset(dst, 0u, cap * sizeof(arbint_limb_t));

    rc = arbint_mul_mag_segmented_rec(dst, &z0n, a0, a0n, b0, b0n, alloc);
    if (rc != ARBINT_OK)
      return rc;

    z2_cap = a1n + b0n + 1u;
    z2 = arbint_alloc_limbs(alloc, z2_cap);
    if (z2 == NULL)
      return ARBINT_ENOMEM;

    rc = arbint_mul_mag_segmented_rec(z2, &z2n, a1, a1n, b0, b0n, alloc);
    if (rc != ARBINT_OK) {
      arbint_free_limbs(alloc, z2);
      return rc;
    }

    used = add_at_offset(dst, z0n, cap, z2, z2n, k);
    arbint_free_limbs(alloc, z2);
    *out_used = arbint_norm_used(dst, used);
    return ARBINT_OK;
  }

  /*  Full Karatsuba case.  */

  /*  Calculate workspace sizes.  */
  z0_cap = a0n + b0n + 1u;
  z2_cap = a1n + b1n + 1u;
  sa_cap = ((a0n > a1n) ? a0n : a1n) + 1u;
  sb_cap = ((b0n > b1n) ? b0n : b1n) + 1u;
  z1_cap = sa_cap + sb_cap + 1u;

  total = z0_cap + z2_cap + z1_cap + sa_cap + sb_cap;

  work = arbint_alloc_limbs(alloc, total);
  if (work == NULL)
    return ARBINT_ENOMEM;

  z0 = work;
  z2 = z0 + z0_cap;
  z1 = z2 + z2_cap;
  sa = z1 + z1_cap;
  sb = sa + sa_cap;

  /*  Compute z0 = A0 * B0.  */
  rc = arbint_mul_mag_segmented_rec(z0, &z0n, a0, a0n, b0, b0n, alloc);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Compute z2 = A1 * B1.  */
  rc = arbint_mul_mag_segmented_rec(z2, &z2n, a1, a1n, b1, b1n, alloc);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Compute sa = A0 + A1 and sb = B0 + B1.  */
  san = arbint__add_mag(sa, a0, a0n, a1, a1n);
  sbn = arbint__add_mag(sb, b0, b0n, b1, b1n);

  /*  Compute z1 = sa * sb.  */
  rc = arbint_mul_mag_segmented_rec(z1, &z1n, sa, san, sb, sbn, alloc);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  z1 = z1 - z0 - z2.  */
  z1n = arbint__sub_mag(z1, z1, z1n, z0, z0n);
  z1n = arbint__sub_mag(z1, z1, z1n, z2, z2n);

  /*  Assemble result: dst = z2 * B^(2k) + z1 * B^k + z0.  */
  cap = an + bn + 1u;
  memset(dst, 0u, cap * sizeof(arbint_limb_t));

  /*  Copy z0 to dst[0..].  */
  if (z0n > 0u)
    memcpy(dst, z0, z0n * sizeof(arbint_limb_t));

  /*  Add z2 at offset 2k.  */
  used = z0n;
  if (z2n > 0u) {
    if (2u * k + z2n > cap) {
      rc = ARBINT_EOVERFLOW;
      goto cleanup;
    }
    memcpy(dst + 2u * k, z2, z2n * sizeof(arbint_limb_t));
    if (2u * k + z2n > used)
      used = 2u * k + z2n;
  }

  /*  Add z1 at offset k.  */
  used = add_at_offset(dst, used, cap, z1, z1n, k);

  *out_used = arbint_norm_used(dst, used);
  rc = ARBINT_OK;

cleanup:
  arbint_free_limbs(alloc, work);
  return rc;
}

/*  Recursive segmented multiplication.
    Dispatches to NTT for operands that fit, or to Karatsuba segmentation
    for larger operands.  */
static arbint_err_t arbint_mul_mag_segmented_rec(arbint_limb_t * dst,
                                                  size_t * out_used,
                                                  const arbint_limb_t * a,
                                                  size_t an,
                                                  const arbint_limb_t * b,
                                                  size_t bn,
                                                  const arbint_alloc_t * alloc) {
  size_t k;
  size_t max_n;

  if (an == 0u || bn == 0u) {
    dst[0] = 0u;
    *out_used = 0u;
    return ARBINT_OK;
  }

  /*  Ensure an >= bn.  */
  if (an < bn) {
    const arbint_limb_t * t = a;
    size_t tn = an;
    a = b;
    an = bn;
    b = t;
    bn = tn;
  }

  /*  Check if operands fit within single NTT.  */
  max_n = (an > bn) ? an : bn;
  if (max_n <= ARBINT_NTT_MAX_SIZE && an + bn <= ARBINT_NTT_MAX_SIZE) {
    /*  Use standard NTT multiplication.  */
    return arbint_mul_mag_ntt(dst, out_used, a, an, b, bn, alloc);
  }

  /*  Choose segment size.
      We want each sub-product to fit within NTT limits.
      The largest sub-product is (a0 + a1) * (b0 + b1), which can be
      up to (k+1) + (k+1) = 2k+2 limbs for inputs of k limbs each.
      So we need 2k+2 <= ARBINT_NTT_MAX_SIZE, thus k <= (MAX_SIZE-2)/2.

      We use ARBINT_NTT_SEGMENT_SIZE = MAX_SIZE/2 which satisfies this.  */
  k = ARBINT_NTT_SEGMENT_SIZE;

  /*  If the smaller operand fits in one segment, use simplified path.  */
  if (bn <= k) {
    /*  B fits in one segment. Split A only.  */
    return karatsuba_segment(dst, out_used, a, an, b, bn, k, alloc);
  }

  /*  Both operands exceed one segment. Apply full Karatsuba.  */
  return karatsuba_segment(dst, out_used, a, an, b, bn, k, alloc);
}

/*  Public entry point for segmented NTT multiplication.  */
arbint_err_t arbint_mul_mag_ntt_segmented(arbint_limb_t * dst, size_t * out_used,
                                          const arbint_limb_t * a, size_t an,
                                          const arbint_limb_t * b, size_t bn,
                                          const arbint_alloc_t * alloc) {
  if (dst == NULL || out_used == NULL || a == NULL || b == NULL)
    return ARBINT_EINVAL;

  if (alloc == NULL || alloc->realloc == NULL)
    return ARBINT_EINVAL;

  return arbint_mul_mag_segmented_rec(dst, out_used, a, an, b, bn, alloc);
}

/*  Forward declaration for recursive squaring.  */
static arbint_err_t arbint_sqr_mag_segmented_rec(arbint_limb_t * dst,
                                                  size_t * out_used,
                                                  const arbint_limb_t * a,
                                                  size_t an,
                                                  const arbint_alloc_t * alloc);

/*  Karatsuba-style segmented squaring.

    Splits operand at k = ARBINT_NTT_SEGMENT_SIZE:
      A = A1 * B^k + A0

    Computes using Karatsuba identity for squaring:
      z0 = A0^2
      z2 = A1^2
      z1 = (A0 + A1)^2 - z0 - z2  (this equals 2*A0*A1)

    Result: z2 * B^(2k) + z1 * B^k + z0

    All three sub-products are squares, so the squaring optimization
    applies recursively.  */
static arbint_err_t karatsuba_sqr_segment(arbint_limb_t * dst, size_t * out_used,
                                           const arbint_limb_t * a, size_t an,
                                           size_t k,
                                           const arbint_alloc_t * alloc) {
  arbint_err_t rc = ARBINT_OK;
  arbint_limb_t * work = NULL;
  arbint_limb_t * z0;
  arbint_limb_t * z2;
  arbint_limb_t * z1;
  arbint_limb_t * sa;
  size_t z0_cap;
  size_t z2_cap;
  size_t z1_cap;
  size_t sa_cap;
  size_t total;
  size_t a0n;
  size_t a1n;
  const arbint_limb_t * a0;
  const arbint_limb_t * a1;
  size_t z0n = 0u;
  size_t z2n = 0u;
  size_t z1n = 0u;
  size_t san;
  size_t cap;
  size_t used;

  /*  Split operand.  */
  a0 = a;
  a0n = (k < an) ? k : an;
  a0n = arbint_norm_used(a0, a0n);

  if (an > k) {
    a1 = a + k;
    a1n = an - k;
    a1n = arbint_norm_used(a1, a1n);
  } else {
    a1 = a;
    a1n = 0u;
  }

  /*  Handle degenerate case where high half is zero.  */
  if (a1n == 0u) {
    /*  Single-segment: just square.  */
    return arbint_sqr_mag_segmented_rec(dst, out_used, a0, a0n, alloc);
  }

  /*  Full Karatsuba squaring case.  */

  /*  Calculate workspace sizes.
      z0 needs 2*a0n limbs, z2 needs 2*a1n limbs.
      sa = a0 + a1 needs max(a0n, a1n) + 1 limbs.
      z1 = sa^2 needs 2*sa_cap limbs.  */
  z0_cap = 2u * a0n + 1u;
  z2_cap = 2u * a1n + 1u;
  sa_cap = ((a0n > a1n) ? a0n : a1n) + 1u;
  z1_cap = 2u * sa_cap + 1u;

  total = z0_cap + z2_cap + z1_cap + sa_cap;

  work = arbint_alloc_limbs(alloc, total);
  if (work == NULL)
    return ARBINT_ENOMEM;

  z0 = work;
  z2 = z0 + z0_cap;
  z1 = z2 + z2_cap;
  sa = z1 + z1_cap;

  /*  Compute z0 = A0^2.  */
  rc = arbint_sqr_mag_segmented_rec(z0, &z0n, a0, a0n, alloc);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Compute z2 = A1^2.  */
  rc = arbint_sqr_mag_segmented_rec(z2, &z2n, a1, a1n, alloc);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Compute sa = A0 + A1.  */
  san = arbint__add_mag(sa, a0, a0n, a1, a1n);

  /*  Compute z1 = (A0 + A1)^2.  */
  rc = arbint_sqr_mag_segmented_rec(z1, &z1n, sa, san, alloc);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  z1 = z1 - z0 - z2 = 2*A0*A1.  */
  z1n = arbint__sub_mag(z1, z1, z1n, z0, z0n);
  z1n = arbint__sub_mag(z1, z1, z1n, z2, z2n);

  /*  Assemble result: dst = z2 * B^(2k) + z1 * B^k + z0.  */
  cap = 2u * an + 1u;
  memset(dst, 0u, cap * sizeof(arbint_limb_t));

  /*  Copy z0 to dst[0..].  */
  if (z0n > 0u)
    memcpy(dst, z0, z0n * sizeof(arbint_limb_t));

  /*  Add z2 at offset 2k.  */
  used = z0n;
  if (z2n > 0u) {
    if (2u * k + z2n > cap) {
      rc = ARBINT_EOVERFLOW;
      goto cleanup;
    }
    memcpy(dst + 2u * k, z2, z2n * sizeof(arbint_limb_t));
    if (2u * k + z2n > used)
      used = 2u * k + z2n;
  }

  /*  Add z1 at offset k.  */
  used = add_at_offset(dst, used, cap, z1, z1n, k);

  *out_used = arbint_norm_used(dst, used);
  rc = ARBINT_OK;

cleanup:
  arbint_free_limbs(alloc, work);
  return rc;
}

/*  Recursive segmented squaring.
    Dispatches to NTT for operands that fit, or to Karatsuba segmentation
    for larger operands.  */
static arbint_err_t arbint_sqr_mag_segmented_rec(arbint_limb_t * dst,
                                                  size_t * out_used,
                                                  const arbint_limb_t * a,
                                                  size_t an,
                                                  const arbint_alloc_t * alloc) {
  size_t k;

  if (an == 0u) {
    dst[0] = 0u;
    *out_used = 0u;
    return ARBINT_OK;
  }

  /*  Check if operand fits within single NTT.
      For squaring, conv_len = 2*an - 1, so we need 2*an <= MAX_SIZE.  */
  if (2u * an <= ARBINT_NTT_MAX_SIZE) {
    /*  Use standard NTT squaring.  */
    return arbint_sqr_mag_ntt(dst, out_used, a, an, alloc);
  }

  /*  Choose segment size.
      We want each sub-product to fit within NTT limits.
      The largest sub-product is (a0 + a1)^2, which can be
      up to 2*(k+1) = 2k+2 limbs for inputs of k limbs each.
      So we need 2k+2 <= ARBINT_NTT_MAX_SIZE, thus k <= (MAX_SIZE-2)/2.

      We use ARBINT_NTT_SEGMENT_SIZE = MAX_SIZE/2 which satisfies this.  */
  k = ARBINT_NTT_SEGMENT_SIZE;

  /*  Apply Karatsuba squaring decomposition.  */
  return karatsuba_sqr_segment(dst, out_used, a, an, k, alloc);
}

/*  Public entry point for segmented NTT squaring.  */
arbint_err_t arbint_sqr_mag_ntt_segmented(arbint_limb_t * dst, size_t * out_used,
                                          const arbint_limb_t * a, size_t an,
                                          const arbint_alloc_t * alloc) {
  if (dst == NULL || out_used == NULL || a == NULL)
    return ARBINT_EINVAL;

  if (alloc == NULL || alloc->realloc == NULL)
    return ARBINT_EINVAL;

  return arbint_sqr_mag_segmented_rec(dst, out_used, a, an, alloc);
}
