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

/*  Reciprocal-based Division with Newton Iteration.

    This module implements O(M(n)) division using Newton iteration for
    computing the reciprocal.

    Algorithm overview:
      1. Normalize divisor d so MSB is set
      2. Compute reciprocal v = floor(2^p / d) using Newton iteration
      3. Compute q_approx = floor((n * v) >> p)
      4. Correct q_approx (at most 2 iterations with guard limbs)
      5. Un-normalize remainder

    Newton iteration for 1/d (scaled):
      We compute v such that v * d is close to 2^p.
      Starting with ~32 bits of precision from table + one limb iteration,
      each Newton step doubles precision: v' = 2*v - v*v*d >> p

    Complexity: O(M(n)) where M(n) is the multiplication complexity.  */

#include "arbint_div_newton.h"
#include "arbint_div.h"
#include "arbint_mul.h"
#include "arbint_addsub.h"
#include "arbint_cmp.h"
#include "arbint_shift.h"

#include <string.h>

/*  Lookup table for initial 1/x approximation.

    For normalized input with top byte b in [128..255],
    entry i = arbint_inv_tab[b - 128] satisfies:
      (256 + i) * b <= 2^16 < (257 + i) * b

    So x0 = 256 + arbint_inv_tab[b - 128] approximates 2^16 / b.  */
const unsigned char arbint_inv_tab[128] = {
    0xff, 0xfc, 0xf8, 0xf4, 0xf0, 0xed, 0xe9, 0xe6,
    0xe2, 0xdf, 0xdc, 0xd9, 0xd6, 0xd3, 0xd0, 0xcd,
    0xca, 0xc7, 0xc4, 0xc2, 0xbf, 0xbc, 0xba, 0xb7,
    0xb5, 0xb2, 0xb0, 0xae, 0xab, 0xa9, 0xa7, 0xa5,
    0xa3, 0xa1, 0x9f, 0x9c, 0x9a, 0x99, 0x97, 0x95,
    0x93, 0x91, 0x8f, 0x8d, 0x8c, 0x8a, 0x88, 0x87,
    0x85, 0x83, 0x82, 0x80, 0x7f, 0x7d, 0x7c, 0x7a,
    0x79, 0x77, 0x76, 0x75, 0x73, 0x72, 0x70, 0x6f,
    0x6e, 0x6d, 0x6b, 0x6a, 0x69, 0x68, 0x66, 0x65,
    0x64, 0x63, 0x62, 0x61, 0x60, 0x5f, 0x5e, 0x5d,
    0x5c, 0x5b, 0x5a, 0x59, 0x58, 0x57, 0x56, 0x55,
    0x54, 0x53, 0x52, 0x51, 0x50, 0x50, 0x4f, 0x4e,
    0x4d, 0x4c, 0x4c, 0x4b, 0x4a, 0x49, 0x49, 0x48,
    0x47, 0x46, 0x46, 0x45, 0x44, 0x44, 0x43, 0x42,
    0x42, 0x41, 0x40, 0x40, 0x3f, 0x3e, 0x3e, 0x3d,
    0x3d, 0x3c, 0x3c, 0x3b, 0x3a, 0x3a, 0x39, 0x39
};

#ifndef ARBINT_NEWTON_MAX_CORRECTIONS
#define ARBINT_NEWTON_MAX_CORRECTIONS 64u
#endif

#ifndef ARBINT_MUL_HIGH_SCHOOLBOOK_THRESHOLD
#define ARBINT_MUL_HIGH_SCHOOLBOOK_THRESHOLD 128u
#endif

/*  Checked addition for size_t.  */
static inline int arbint_size_add(size_t a, size_t b, size_t * out) {
  size_t sum = a + b;
  if (sum < a)
    return 0;
  *out = sum;
  return 1;
}

/*  Checked multiplication for size_t.  */
static inline int arbint_size_mul(size_t a, size_t b, size_t * out) {
  if (a != 0 && b > SIZE_MAX / a)
    return 0;
  *out = a * b;
  return 1;
}

/*  Add one limb to a fixed-size little-endian accumulator.
    Returns 1 on success, 0 if the addition overflows acc_n limbs.  */
static inline int arbint_acc_add_limb(arbint_limb_t * acc, size_t acc_n,
                                      size_t idx, arbint_limb_t x) {
  while (x != 0) {
    arbint_limb_t prev;
    arbint_limb_t sum;

    if (idx >= acc_n)
      return 0;
    prev = acc[idx];
    sum = prev + x;
    acc[idx] = sum;
    x = (sum < prev) ? 1u : 0u;
    ++idx;
  }
  return 1;
}

/*  Fallback high product via full multiplication and extraction.
    Computes out = floor((a*b) / B^cut), truncated to out_cap limbs.  */
static arbint_err_t arbint_mul_high_full(arbint_limb_t * out, size_t * out_n,
                                         const arbint_limb_t * a, size_t an,
                                         const arbint_limb_t * b, size_t bn,
                                         size_t cut, size_t out_cap,
                                         const arbint_alloc_t * alloc) {
  arbint_limb_t * prod = NULL;
  size_t prod_cap;
  size_t prod_n;
  size_t used = 0;
  size_t i;
  arbint_err_t rc;

  if (!arbint_size_add(an, bn, &prod_cap))
    return ARBINT_EOVERFLOW;
  if (!arbint_size_add(prod_cap, 1, &prod_cap))
    return ARBINT_EOVERFLOW;

  prod = arbint_alloc_limbs(alloc, prod_cap);
  if (prod == NULL)
    return ARBINT_ENOMEM;

  rc = arbint_mul_mag_generic(prod, &prod_n, a, an, b, bn, alloc);
  if (rc != ARBINT_OK) {
    arbint_free_limbs(alloc, prod);
    return rc;
  }

  if (prod_n > cut) {
    used = prod_n - cut;
    if (used > out_cap)
      used = out_cap;
    memcpy(out, prod + cut, used * sizeof(arbint_limb_t));
  }

  for (i = used; i < out_cap; ++i)
    out[i] = 0;

  *out_n = arbint_norm_used(out, used);

  arbint_free_limbs(alloc, prod);
  return ARBINT_OK;
}

/*  Schoolbook high product without materializing low limbs.
    Exact, column-wise accumulation:
      - process columns [0, cut) to carry state only
      - emit columns [cut, cut + out_cap) into out[]
    The implementation is intended for small/medium operand sizes.  */
static arbint_err_t
arbint_mul_high_schoolbook(arbint_limb_t * out, size_t * out_n,
                           const arbint_limb_t * a, size_t an,
                           const arbint_limb_t * b, size_t bn, size_t cut,
                           size_t out_cap) {
  enum { ARBINT_MUL_HIGH_ACC_LIMBS = 5 };
  arbint_limb_t carry[ARBINT_MUL_HIGH_ACC_LIMBS];
  size_t full_cols;
  size_t upper;
  size_t k;
  size_t out_used = 0;
  size_t i;

  memset(carry, 0, sizeof(carry));

  if (out_cap == 0 || an == 0 || bn == 0) {
    *out_n = 0;
    return ARBINT_OK;
  }

  if (!arbint_size_add(an, bn, &full_cols))
    return ARBINT_EOVERFLOW;
  if (cut >= full_cols) {
    memset(out, 0, out_cap * sizeof(arbint_limb_t));
    *out_n = 0;
    return ARBINT_OK;
  }

  if (!arbint_size_add(cut, out_cap, &upper))
    upper = full_cols;
  if (upper > full_cols)
    upper = full_cols;

  for (k = 0; k < upper; ++k) {
    arbint_limb_t acc[ARBINT_MUL_HIGH_ACC_LIMBS];
    size_t i_lo;
    size_t i_hi;

    memcpy(acc, carry, sizeof(acc));

    if (k >= bn - 1)
      i_lo = k - (bn - 1);
    else
      i_lo = 0;
    i_hi = (k < an) ? k : (an - 1);

    if (i_lo <= i_hi) {
      size_t ii;
      for (ii = i_lo;; ++ii) {
        size_t jj = k - ii;
        arbint_limb_t hi;
        arbint_limb_t lo;

        arbint_umul_limb_generic(&hi, &lo, a[ii], b[jj]);
        if (!arbint_acc_add_limb(acc, ARBINT_MUL_HIGH_ACC_LIMBS, 0, lo) ||
            !arbint_acc_add_limb(acc, ARBINT_MUL_HIGH_ACC_LIMBS, 1, hi))
          return ARBINT_EOVERFLOW;
        if (ii == i_hi)
          break;
      }
    }

    if (k >= cut && out_used < out_cap)
      out[out_used++] = acc[0];

    for (i = 0; i + 1 < ARBINT_MUL_HIGH_ACC_LIMBS; ++i)
      carry[i] = acc[i + 1];
    carry[ARBINT_MUL_HIGH_ACC_LIMBS - 1] = 0;
  }

  for (i = out_used; i < out_cap; ++i)
    out[i] = 0;
  *out_n = arbint_norm_used(out, out_used);
  return ARBINT_OK;
}

/*  Compute out = floor((a*b) / B^cut), truncated to out_cap limbs.  */
static arbint_err_t arbint_mul_high(arbint_limb_t * out, size_t * out_n,
                                    const arbint_limb_t * a, size_t an,
                                    const arbint_limb_t * b, size_t bn,
                                    size_t cut, size_t out_cap,
                                    const arbint_alloc_t * alloc) {
  if (an <= ARBINT_MUL_HIGH_SCHOOLBOOK_THRESHOLD &&
      bn <= ARBINT_MUL_HIGH_SCHOOLBOOK_THRESHOLD) {
    arbint_err_t rc =
        arbint_mul_high_schoolbook(out, out_n, a, an, b, bn, cut, out_cap);
    if (rc == ARBINT_OK)
      return ARBINT_OK;
    /*  If fixed-width accumulator overflows on pathological sizes,
        fall back to full multiplication.  */
  }

  return arbint_mul_high_full(out, out_n, a, an, b, bn, cut, out_cap, alloc);
}

/*  Compare x against B^p_limbs.
    Returns -1 if x < B^p, 0 if x == B^p, 1 if x > B^p.  */
static int arbint_cmp_mag_pow_limb(const arbint_limb_t * x, size_t xn,
                                   size_t p_limbs) {
  size_t i;

  while (xn > 0 && x[xn - 1] == 0)
    --xn;

  if (xn < p_limbs + 1)
    return -1;
  if (xn > p_limbs + 1)
    return 1;

  if (x[p_limbs] < 1)
    return -1;
  if (x[p_limbs] > 1)
    return 1;

  for (i = 0; i < p_limbs; ++i) {
    if (x[i] != 0)
      return 1;
  }
  return 0;
}

/* ========================================================================= */
/*  Newton Iteration for Reciprocal                                          */
/* ========================================================================= */

/*  Compute floor(2^(2*BITS) / d1) where d1 is a single normalized limb.
    Uses the lookup table for initial approximation, then one Newton step.  */
static arbint_limb_t arbint_invert_limb(arbint_limb_t d1) {
  arbint_limb_t d0;
  arbint_limb_t v0, v1, v2;
  arbint_limb_t p, t0, t1;
  unsigned idx;

  /*  d1 is normalized, so MSB is set.  */
  idx = (unsigned)(d1 >> (ARBINT_LIMB_BITS - 8)) - 128;
  v0 = (arbint_limb_t)arbint_inv_tab[idx] + 256;

  /*  v0 approximates 2^16 / (d1 >> (BITS-16)) with ~8 bits precision.
      Refine to ~16 bits: v1 = v0 * (2^17 - v0 * (d1 >> (BITS-16))) >> 16.  */
#if ARBINT_LIMB_BITS == 32
  {
    uint64_t d0_64 = (uint64_t) (d1 >> 16);
    uint64_t v0_64 = (uint64_t) v0;
    uint64_t p64 = v0_64 * v0_64 * d0_64;
    uint64_t v1_64 = (v0_64 * ((uint64_t) 1 << 17) - p64) >> 16;
    v1 = (arbint_limb_t) v1_64;
  }
#else
  d0 = d1 >> (ARBINT_LIMB_BITS - 16);
  p = v0 * v0 * d0;
  v1 = (v0 * ((arbint_limb_t)1 << 17) - p) >> 16;
#endif

#if ARBINT_LIMB_BITS == 32
  /*  32-bit limbs: one wide Newton refinement reaches full-limb precision.  */
  {
    arbint_limb_t e, hi, lo;

    arbint_umul_limb_generic(&hi, &lo, v1, d1);
    e = ~hi;
    arbint_umul_limb_generic(&t1, &t0, v1, e);
    v1 += t1;
  }

  return v1;
#else
  /*  Refine to ~32 bits: v2 = v1 * (2^33 - v1 * (d1 >> (BITS-32))) >> 32.  */
  d0 = d1 >> (ARBINT_LIMB_BITS - 32);
  p = v1 * v1 * d0;
  v2 = (v1 * ((arbint_limb_t)1 << 33) - p) >> 32;

  /*  Final refinement to full limb precision using wide multiply.  */
  {
    arbint_limb_t e, hi, lo;

    /*  Compute v2 * d1, get high part.  */
    arbint_umul_limb_generic(&hi, &lo, v2, d1);
    e = ~hi;  /*  e = 2^BITS - 1 - hi, approximately the error.  */

    /*  v3 = v2 + (v2 * e >> BITS).  */
    arbint_umul_limb_generic(&t1, &t0, v2, e);
    v2 += t1;
  }

  return v2;
#endif
}

/*  Compute v = floor(2^(2*dn*BITS) / d) using Newton iteration.

    The algorithm works by:
    1. Computing an initial single-limb inverse of the top limb
    2. Iteratively doubling precision using v' = 2*v - v^2*d >> shift

    Parameters:
      v:       output (capacity: dn + 1 limbs)
      v_n:     output number of limbs
      d:       normalized divisor (dn limbs, MSB of d[dn-1] is set)
      dn:      number of limbs in d
      alloc:   allocator  */
static arbint_err_t arbint_invert_newton_core(arbint_limb_t * v, size_t * v_n,
                                               const arbint_limb_t * d,
                                               size_t dn,
                                               const arbint_alloc_t * alloc) {
  size_t k;
  size_t sizes[64];
  int num_sizes;
  int i;
  arbint_limb_t * scratch = NULL;
  arbint_limb_t * vv = NULL;
  arbint_limb_t * vvd = NULL;
  arbint_limb_t * two_v = NULL;
  size_t scratch_size;
  size_t vv_n, vvd_n, two_v_n;
  arbint_err_t rc = ARBINT_OK;

  if (dn == 0) {
    *v_n = 0;
    return ARBINT_EINVAL;
  }

  /*  Base case: single limb divisor.  */
  if (dn == 1) {
    v[0] = arbint_invert_limb(d[0]);
    v[1] = 0;
    *v_n = 1;
    return ARBINT_OK;
  }

  /*  Build list of sizes from dn down to 1.
      sizes[0] = dn, sizes[num_sizes-1] = 1.  */
  num_sizes = 0;
  k = dn;
  while (k > 1) {
    sizes[num_sizes++] = k;
    k = (k + 1) / 2;
  }
  sizes[num_sizes++] = 1;

  /*  Allocate scratch: v^2 (2*dn+2), v^2*d (3*dn+3), 2*v (dn+2).  */
  scratch_size = 6 * dn + 10;
  scratch = arbint_alloc_limbs(alloc, scratch_size);
  if (scratch == NULL)
    return ARBINT_ENOMEM;

  vv = scratch;
  vvd = vv + 2 * dn + 2;
  two_v = vvd + 3 * dn + 3;

  /*  Initialize with single-limb inverse.  */
  v[0] = arbint_invert_limb(d[dn - 1]);
  *v_n = 1;

  /*  Newton iterations, doubling precision each time.  */
  for (i = num_sizes - 2; i >= 0; --i) {
    size_t target_n = sizes[i];
    size_t d_offset = dn - target_n;
    const arbint_limb_t * d_top = d + d_offset;
    size_t shift_limbs;
    size_t new_v_n;
    size_t j;

    /*  v' = 2*v - v^2 * d_top >> (target_n * BITS).

        Here v approximates 2^((v_n)*BITS) / d_top[0..v_n-1].
        We want v' to approximate 2^((target_n)*BITS) / d_top[0..target_n-1].  */

    /*  Compute v^2.  */
    rc = arbint_mul_mag_generic(vv, &vv_n, v, *v_n, v, *v_n, alloc);
    if (rc != ARBINT_OK)
      goto cleanup;

    /*  Compute v^2 * d_top.  */
    rc = arbint_mul_mag_generic(vvd, &vvd_n, vv, vv_n, d_top, target_n, alloc);
    if (rc != ARBINT_OK)
      goto cleanup;

    /*  Shift right by (*v_n) limbs.  */
    shift_limbs = *v_n;
    if (vvd_n > shift_limbs) {
      memmove(vvd, vvd + shift_limbs,
              (vvd_n - shift_limbs) * sizeof(arbint_limb_t));
      vvd_n -= shift_limbs;
    } else {
      vvd_n = 0;
    }

    /*  Compute 2*v, shifted left by (target_n - *v_n) limbs.  */
    {
      size_t pad = target_n - *v_n;
      memset(two_v, 0, pad * sizeof(arbint_limb_t));
      two_v_n = arbint__dbl_mag(two_v + pad, v, *v_n);
      two_v_n += pad;
    }

    /*  v' = 2*v - v^2*d >> shift.  */
    if (arbint_cmp_mag_limbs(two_v, two_v_n, vvd, vvd_n) >= 0) {
      new_v_n = arbint__sub_mag(v, two_v, two_v_n, vvd, vvd_n);
    } else {
      /*  Underflow - shouldn't happen with correct iteration.  */
      memset(v, 0, target_n * sizeof(arbint_limb_t));
      new_v_n = 0;
    }

    /*  Trim to target_n limbs.  */
    if (new_v_n > target_n) {
      memmove(v, v + (new_v_n - target_n), target_n * sizeof(arbint_limb_t));
      new_v_n = target_n;
    }

    /*  Zero-pad if needed.  */
    for (j = new_v_n; j < target_n; ++j)
      v[j] = 0;

    *v_n = arbint_norm_used(v, target_n);
    if (*v_n == 0) {
      *v_n = 1;
      v[0] = 1;
    }
  }

  rc = ARBINT_OK;

cleanup:
  arbint_free_limbs(alloc, scratch);
  return rc;
}

/*  Compute v = floor(2^p / d) by direct division (fallback).  */
static arbint_err_t arbint_invert_direct(arbint_limb_t * v, size_t * v_n,
                                          const arbint_limb_t * d, size_t dn,
                                          size_t p_limbs,
                                          const arbint_alloc_t * alloc) {
  size_t num_limbs;
  size_t alloc_size;
  arbint_limb_t * num = NULL;
  arbint_limb_t * rem = NULL;
  size_t i;
  arbint_err_t rc = ARBINT_OK;

  if (v == NULL || v_n == NULL || d == NULL || dn == 0)
    return ARBINT_EINVAL;
  if (p_limbs < dn)
    return ARBINT_EINVAL;

  if (!arbint_size_add(p_limbs, 1, &num_limbs))
    return ARBINT_EOVERFLOW;
  if (!arbint_size_add(num_limbs, dn, &alloc_size))
    return ARBINT_EOVERFLOW;

  num = arbint_alloc_limbs(alloc, alloc_size);
  if (num == NULL)
    return ARBINT_ENOMEM;
  rem = num + num_limbs;

  memset(num, 0, num_limbs * sizeof(arbint_limb_t));
  num[p_limbs] = 1;

  rc = arbint_div_mag_knuth(num, num_limbs, d, dn, v, rem);
  if (rc != ARBINT_OK)
    goto cleanup;

  *v_n = arbint_norm_used(v, p_limbs - dn + 2);

  for (i = *v_n; i < p_limbs - dn + 2; ++i)
    v[i] = 0;

cleanup:
  arbint_free_limbs(alloc, num);
  return rc;
}

/*  Internal: compute reciprocal v = floor(2^(p_limbs * BITS) / d).  */
static arbint_err_t arbint_invert_newton_prec(arbint_limb_t * v, size_t * v_n,
                                              const arbint_limb_t * d,
                                              size_t dn, size_t p_limbs,
                                              const arbint_alloc_t * alloc) {
  if (p_limbs < dn)
    return ARBINT_EINVAL;

  /*  Use Newton for large divisors, direct for small.  */
  if (dn >= ARBINT_INVERT_DC_THRESHOLD) {
    arbint_err_t rc;
    size_t newton_v_n;
    size_t shift_limbs;
    size_t two_dn;
    size_t result_cap;
    size_t result_n;
    size_t adjust;
    int cmp_pow;
    arbint_limb_t * check_prod = NULL;
    size_t check_prod_n;

    if (!arbint_size_add(p_limbs - dn, 2, &result_cap))
      return ARBINT_EOVERFLOW;
    if (!arbint_size_mul(dn, 2, &two_dn))
      return ARBINT_EOVERFLOW;

    /*  Newton core builds B^(2*dn)/d. If target precision exceeds that,
        direct inversion is required for correctness.  */
    if (p_limbs > two_dn)
      return arbint_invert_direct(v, v_n, d, dn, p_limbs, alloc);

    rc = arbint_invert_newton_core(v, &newton_v_n, d, dn, alloc);
    if (rc != ARBINT_OK) {
      if (rc == ARBINT_ENOMEM)
        return rc;
      return arbint_invert_direct(v, v_n, d, dn, p_limbs, alloc);
    }

    /*  Newton core approximates floor(B^(2*dn) / d). Convert to floor(B^p / d)
        by right-shifting (2*dn - p) limbs when p <= 2*dn.
        For p > 2*dn the core output lacks required low precision.  */
    shift_limbs = two_dn - p_limbs;

    if (newton_v_n > shift_limbs) {
      result_n = newton_v_n - shift_limbs;
      memmove(v, v + shift_limbs, result_n * sizeof(arbint_limb_t));
    } else {
      return arbint_invert_direct(v, v_n, d, dn, p_limbs, alloc);
    }

    if (result_n == 0 || result_n > result_cap)
      return arbint_invert_direct(v, v_n, d, dn, p_limbs, alloc);

    check_prod = arbint_alloc_limbs(alloc, result_cap + dn + 4);
    if (check_prod == NULL)
      return ARBINT_ENOMEM;

    rc = arbint_mul_mag_generic(check_prod, &check_prod_n, v, result_n, d, dn,
                                 alloc);
    if (rc != ARBINT_OK)
      goto fallback_direct;

    /*  Enforce v*d <= B^p.  */
    adjust = 0;
    while ((cmp_pow = arbint_cmp_mag_pow_limb(check_prod, check_prod_n,
                                              p_limbs)) > 0) {
      if (adjust++ >= 64 || (result_n == 1 && v[0] == 0))
        goto fallback_direct;
      arbint_limb_sub_1(v, v, result_n, 1);
      result_n = arbint_norm_used(v, result_n);
      if (result_n == 0)
        goto fallback_direct;
      check_prod_n = arbint__sub_mag(check_prod, check_prod, check_prod_n, d,
                                     dn);
    }

    /*  Enforce (v+1)*d > B^p.  */
    adjust = 0;
    while (1) {
      size_t vp1_n;
      arbint_limb_t carry;

      if (adjust++ >= 64)
        goto fallback_direct;

      vp1_n = arbint__add_mag(check_prod, check_prod, check_prod_n, d, dn);
      if (arbint_cmp_mag_pow_limb(check_prod, vp1_n, p_limbs) > 0) {
        check_prod_n = arbint__sub_mag(check_prod, check_prod, vp1_n, d, dn);
        break;
      }
      carry = arbint_limb_add_1(v, v, result_n, 1);
      if (carry != 0) {
        if (result_n >= result_cap)
          goto fallback_direct;
        v[result_n++] = carry;
      }
      check_prod_n = vp1_n;
    }

    arbint_free_limbs(alloc, check_prod);

    *v_n = arbint_norm_used(v, result_n);
    if (*v_n == 0)
      *v_n = 1;
    for (adjust = *v_n; adjust < result_cap; ++adjust)
      v[adjust] = 0;
    return ARBINT_OK;

fallback_direct:
    arbint_free_limbs(alloc, check_prod);
    rc = arbint_invert_direct(v, v_n, d, dn, p_limbs, alloc);
    if (rc != ARBINT_OK)
      return rc;
    for (adjust = *v_n; adjust < result_cap; ++adjust)
      v[adjust] = 0;
    return ARBINT_OK;
  }

  return arbint_invert_direct(v, v_n, d, dn, p_limbs, alloc);
}

/*  Public: compute reciprocal v = floor(2^((dn + GUARD) * BITS) / d).  */
arbint_err_t arbint_invert_newton(arbint_limb_t * v, size_t * v_n,
                                   const arbint_limb_t * d, size_t dn,
                                   const arbint_alloc_t * alloc) {
  size_t p_limbs;

  if (!arbint_size_add(dn, ARBINT_INVERT_GUARD_LIMBS, &p_limbs))
    return ARBINT_EOVERFLOW;

  return arbint_invert_newton_prec(v, v_n, d, dn, p_limbs, alloc);
}

/*  Compute scratch size for D&C inverse.  */
size_t arbint_invert_dc_scratch_size(size_t dn) {
  return 8 * dn + 32;
}

/*  D&C inverse wrapper.  */
arbint_err_t arbint_invert_dc(arbint_limb_t * v, size_t * v_n,
                               const arbint_limb_t * d, size_t dn,
                               arbint_limb_t * scratch, size_t scratch_size,
                               const arbint_alloc_t * alloc) {
  (void)scratch;
  (void)scratch_size;
  return arbint_invert_newton(v, v_n, d, dn, alloc);
}

/*  Compute scratch size for Newton division.  */
size_t arbint_div_newton_scratch_size(size_t nn, size_t dn) {
  size_t n_used = nn + 1;
  size_t p_limbs = nn + ARBINT_INVERT_GUARD_LIMBS;
  size_t v_cap = p_limbs - dn + 3;
  size_t q_cap = nn - dn + 2;

  /*  Layout: d_norm[dn] | n_norm[n_used] | v[v_cap] | q_tmp[q_cap+1]
            | qd[q_cap+dn+1] | r_tmp[nn+2]  */
  return dn + n_used + v_cap + (q_cap + 1) + (q_cap + dn + 1) + (nn + 2);
}

/*  Newton-Raphson division.  */
arbint_err_t arbint_div_mag_newton(const arbint_limb_t * np, size_t nn,
                                    const arbint_limb_t * dp, size_t dn,
                                    arbint_limb_t * qp, arbint_limb_t * rp,
                                    const arbint_alloc_t * alloc) {
  arbint_limb_t * d_norm = NULL;
  arbint_limb_t * n_norm = NULL;
  arbint_limb_t * v = NULL;
  arbint_limb_t * qd = NULL;
  arbint_limb_t * q_tmp = NULL;
  arbint_limb_t * r_tmp = NULL;
  arbint_limb_t * scratch = NULL;
  size_t v_n;
  size_t qd_n;
  size_t q_n;
  size_t r_n;
  unsigned shift;
  size_t total_scratch;
  size_t p_limbs;
  size_t v_cap;
  size_t q_cap;
  size_t n_used;
  size_t i;
  int corrections;
  arbint_err_t rc = ARBINT_OK;

  if (np == NULL || dp == NULL || rp == NULL)
    return ARBINT_EINVAL;
  if (dn < 2 || nn < dn)
    return ARBINT_EINVAL;
  if (dp[dn - 1] == 0)
    return ARBINT_EINVAL;

  /*  p = (nn + GUARD_LIMBS) * BITS.  */
  if (!arbint_size_add(nn, ARBINT_INVERT_GUARD_LIMBS, &p_limbs))
    return ARBINT_EOVERFLOW;

  q_cap = nn - dn + 2;

  if (!arbint_size_add(nn, 1, &n_used))
    return ARBINT_EOVERFLOW;

  if (p_limbs < dn)
    return ARBINT_EINVAL;
  v_cap = p_limbs - dn + 3;

  /*  Compute scratch size.
      Layout: d_norm[dn] | n_norm[n_used] | v[v_cap] | q_tmp[q_cap+1]
            | qd[q_cap+dn+1] | r_tmp[nn+2]  */
  {
    size_t part1, part2, part3;

    if (!arbint_size_add(dn, n_used, &part1))
      return ARBINT_EOVERFLOW;
    if (!arbint_size_add(part1, v_cap, &part1))
      return ARBINT_EOVERFLOW;
    if (!arbint_size_add(q_cap, 1, &part2))
      return ARBINT_EOVERFLOW;
    if (!arbint_size_add(part1, part2, &part1))
      return ARBINT_EOVERFLOW;
    if (!arbint_size_add(q_cap, dn, &part2))
      return ARBINT_EOVERFLOW;
    if (!arbint_size_add(part2, 1, &part2))
      return ARBINT_EOVERFLOW;
    if (!arbint_size_add(part1, part2, &part1))
      return ARBINT_EOVERFLOW;
    if (!arbint_size_add(nn, 2, &part3))
      return ARBINT_EOVERFLOW;
    if (!arbint_size_add(part1, part3, &total_scratch))
      return ARBINT_EOVERFLOW;
  }

  scratch = arbint_alloc_limbs(alloc, total_scratch);
  if (scratch == NULL)
    return ARBINT_ENOMEM;

  d_norm = scratch;
  n_norm = d_norm + dn;
  v = n_norm + n_used;
  q_tmp = v + v_cap;
  qd = q_tmp + q_cap + 1;
  r_tmp = qd + q_cap + dn + 1;

  /*  Normalize divisor.  */
  shift = arbint_clz_limb(dp[dn - 1]);

  if (shift != 0) {
    arbint_limb_t carry = 0;
    for (i = 0; i < dn; ++i) {
      d_norm[i] = (dp[i] << shift) | carry;
      carry = dp[i] >> (ARBINT_LIMB_BITS - shift);
    }
  } else {
    memcpy(d_norm, dp, dn * sizeof(arbint_limb_t));
  }

  /*  Normalize dividend.  */
  if (shift != 0) {
    arbint_limb_t carry = 0;
    for (i = 0; i < nn; ++i) {
      n_norm[i] = (np[i] << shift) | carry;
      carry = np[i] >> (ARBINT_LIMB_BITS - shift);
    }
    n_norm[nn] = carry;
  } else {
    memcpy(n_norm, np, nn * sizeof(arbint_limb_t));
    n_norm[nn] = 0;
  }

  /*  Compute reciprocal at division precision p = (nn + guard) limbs.
      This uses Newton where available and falls back to direct inversion
      if Newton refinement or verification cannot guarantee correctness.  */
  rc = arbint_invert_newton_prec(v, &v_n, d_norm, dn, p_limbs, alloc);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Compute q_approx = (n_norm * v) >> p using high-product extraction.  */
  rc = arbint_mul_high(q_tmp, &q_n, n_norm, n_used, v, v_n, p_limbs, q_cap,
                       alloc);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Compute r = n_norm - q * d_norm.  */
  if (q_n > 0) {
    rc = arbint_mul_mag_generic(qd, &qd_n, q_tmp, q_n, d_norm, dn, alloc);
    if (rc != ARBINT_OK)
      goto cleanup;
  } else {
    qd_n = 0;
  }

  for (i = qd_n; i < n_used; ++i)
    qd[i] = 0;

  /*  Correction loop.  */
  corrections = 0;

  if (arbint_cmp_mag_limbs(n_norm, n_used, qd, n_used) >= 0) {
    r_n = arbint__sub_mag(r_tmp, n_norm, n_used, qd, n_used);

    while (arbint_cmp_mag_limbs(r_tmp, r_n, d_norm, dn) >= 0 &&
           corrections < (int) ARBINT_NEWTON_MAX_CORRECTIONS) {
      r_n = arbint__sub_mag(r_tmp, r_tmp, r_n, d_norm, dn);
      if (q_n < q_cap) {
        arbint_limb_t carry = arbint_limb_add_1(q_tmp, q_tmp, q_n, 1);
        if (carry != 0 && q_n < q_cap)
          q_tmp[q_n++] = carry;
      }
      ++corrections;
    }
  } else {
    while (arbint_cmp_mag_limbs(qd, qd_n > n_used ? qd_n : n_used,
                                 n_norm, n_used) > 0 &&
           corrections < (int) ARBINT_NEWTON_MAX_CORRECTIONS) {
      if (q_n > 0) {
        arbint_limb_sub_1(q_tmp, q_tmp, q_n, 1);
        q_n = arbint_norm_used(q_tmp, q_n);
      }

      if (qd_n >= dn) {
        qd_n = arbint__sub_mag(qd, qd, qd_n, d_norm, dn);
      } else {
        break;
      }
      ++corrections;
    }

    if (arbint_cmp_mag_limbs(n_norm, n_used, qd, qd_n) >= 0) {
      r_n = arbint__sub_mag(r_tmp, n_norm, n_used, qd, qd_n);
    } else {
      rc = ARBINT_EINVAL;
      goto cleanup;
    }

    while (arbint_cmp_mag_limbs(r_tmp, r_n, d_norm, dn) >= 0 &&
           corrections < (int) ARBINT_NEWTON_MAX_CORRECTIONS) {
      r_n = arbint__sub_mag(r_tmp, r_tmp, r_n, d_norm, dn);
      if (q_n < q_cap) {
        arbint_limb_t carry = arbint_limb_add_1(q_tmp, q_tmp, q_n, 1);
        if (carry != 0 && q_n < q_cap)
          q_tmp[q_n++] = carry;
      }
      ++corrections;
    }
  }

  if (corrections >= (int) ARBINT_NEWTON_MAX_CORRECTIONS) {
    rc = arbint_div_mag_knuth(np, nn, dp, dn, qp, rp);
    goto cleanup;
  }

  /*  Un-normalize remainder.  */
  if (shift != 0) {
    r_n = arbint_rshift_limbs_inplace(r_tmp, r_n, shift);
  }

  r_n = arbint_norm_used(r_tmp, r_n);

  if (r_n > dn)
    r_n = dn;
  memcpy(rp, r_tmp, r_n * sizeof(arbint_limb_t));
  for (i = r_n; i < dn; ++i)
    rp[i] = 0;

  if (qp != NULL) {
    q_n = arbint_norm_used(q_tmp, q_n);
    if (q_n > nn - dn + 1)
      q_n = nn - dn + 1;
    memcpy(qp, q_tmp, q_n * sizeof(arbint_limb_t));
    for (i = q_n; i < nn - dn + 1; ++i)
      qp[i] = 0;
  }

  rc = ARBINT_OK;

cleanup:
  arbint_free_limbs(alloc, scratch);
  return rc;
}
