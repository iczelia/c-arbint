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

#include "arbint_mul.h"

#include "config.h"

#include <string.h>

/*  Multiply two limbs producing full double-width result (hi:lo = x * y).
    Uses half-limb multiplication and combining when 2x-width types
    unavailable.  */
static void arbint_mul_wide_limb(arbint_limb_t x, arbint_limb_t y,
                                 arbint_limb_t * hi, arbint_limb_t * lo) {
#if ARBINT_LIMB_BITS == 32
  uint64_t p = (uint64_t) x * (uint64_t) y;
  *lo = (arbint_limb_t) p;
  *hi = (arbint_limb_t) (p >> 32);
#else
  const arbint_limb_t mask = ARBINT_HALF_MASK;
  arbint_limb_t x0 = x & mask;
  arbint_limb_t x1 = x >> ARBINT_HALF_BITS;
  arbint_limb_t y0 = y & mask;
  arbint_limb_t y1 = y >> ARBINT_HALF_BITS;

  arbint_limb_t w0 = x0 * y0;
  arbint_limb_t t = x1 * y0 + (w0 >> ARBINT_HALF_BITS);
  arbint_limb_t w1 = t & mask;
  arbint_limb_t w2 = t >> ARBINT_HALF_BITS;

  w1 = x0 * y1 + w1;

  *hi = x1 * y1 + w2 + (w1 >> ARBINT_HALF_BITS);
  *lo = (w1 << ARBINT_HALF_BITS) | (w0 & mask);
#endif /* ARBINT_LIMB_BITS */
}

/*  Compute x*y + acc + carry, returning high limb in result and low limb
    in *out. Essential primitive for schoolbook multiplication inner loop.  */
static arbint_limb_t arbint_muladd_limb(arbint_limb_t x, arbint_limb_t y,
                                        arbint_limb_t acc, arbint_limb_t carry,
                                        arbint_limb_t * out) {
  arbint_limb_t ph;
  arbint_limb_t pl;

  arbint_mul_wide_limb(x, y, &ph, &pl);

#if ARBINT_HAVE_X86_CARRY_KERNEL
  {
    arbint_x86_carry_word_t t = (arbint_x86_carry_word_t) 0;
    arbint_x86_carry_word_t o = (arbint_x86_carry_word_t) 0;
    unsigned char c1;
    unsigned char c2;

    c1 = ARBINT_X86_ADDCARRY((unsigned char) 0, (arbint_x86_carry_word_t) pl,
                             (arbint_x86_carry_word_t) acc, &t);
    c2 = ARBINT_X86_ADDCARRY((unsigned char) 0, t,
                             (arbint_x86_carry_word_t) carry, &o);
    *out = (arbint_limb_t) o;
    return ph + (arbint_limb_t) c1 + (arbint_limb_t) c2;
  }
#else
  {
    arbint_limb_t s;
    arbint_limb_t c1;
    arbint_limb_t s2;
    arbint_limb_t c2;

    s = pl + acc;
    c1 = (s < pl) ? 1u : 0u;
    s2 = s + carry;
    c2 = (s2 < s) ? 1u : 0u;

    *out = s2;
    return ph + c1 + c2;
  }
#endif /* ARBINT_HAVE_X86_CARRY_KERNEL */
}

/*  Classical O(n^2) schoolbook multiplication algorithm.
    Multiplies each limb of a by each limb of b, accumulating into result.
    Returns normalized result limb count (usually an+bn, may be an+bn-1 if no
    high carry).  */
static size_t arbint_mul_schoolbook(arbint_limb_t * dst,
                                    const arbint_limb_t * a, size_t an,
                                    const arbint_limb_t * b, size_t bn) {
  size_t i;
  size_t n;

  if (an == 0u || bn == 0u) {
    dst[0] = 0u;
    return 0u;
  }

  n = an + bn + 1u;
  memset(dst, 0, n * sizeof(arbint_limb_t));

  for (i = 0u; i < an; ++i) {
    size_t j;
    arbint_limb_t carry = 0u;

    for (j = 0u; j < bn; ++j) {
      arbint_limb_t out;
      carry = arbint_muladd_limb(a[i], b[j], dst[i + j], carry, &out);
      dst[i + j] = out;
    }

    {
      size_t k = i + bn;
#if ARBINT_HAVE_X86_CARRY_KERNEL
      arbint_x86_carry_word_t out = (arbint_x86_carry_word_t) 0;
      unsigned char c;

      c = ARBINT_X86_ADDCARRY((unsigned char) 0,
                              (arbint_x86_carry_word_t) dst[k],
                              (arbint_x86_carry_word_t) carry, &out);
      dst[k] = (arbint_limb_t) out;
      ++k;
      while (c != 0u) {
        c = ARBINT_X86_ADDCARRY(c, (arbint_x86_carry_word_t) dst[k],
                                (arbint_x86_carry_word_t) 0, &out);
        dst[k] = (arbint_limb_t) out;
        ++k;
      }
#else
      while (carry != 0u) {
        arbint_limb_t t = dst[k] + carry;
        dst[k] = t;
        carry = (t < carry) ? 1u : 0u;
        ++k;
      }
#endif /* ARBINT_HAVE_X86_CARRY_KERNEL */
    }
  }

  return arbint_norm_used(dst, n);
}

/*  Add src array to dst starting at position 'shift' (dst[shift:] += src).
    Used by Karatsuba to accumulate sub-products at different offsets.
    Returns 1 on success, 0 on overflow (carry extends beyond dst_n).  */
static int arbint_add_shifted(arbint_limb_t * dst, size_t dst_n,
                              const arbint_limb_t * src, size_t src_n,
                              size_t shift) {
  size_t i;

  if (src_n == 0u)
    return 1;
  if (shift > dst_n || src_n > dst_n - shift)
    return 0;

#if ARBINT_HAVE_X86_CARRY_KERNEL
  {
    unsigned char carry = 0u;
    arbint_x86_carry_word_t out = (arbint_x86_carry_word_t) 0;

    for (i = 0u; i < src_n; ++i) {
      size_t di = shift + i;
      carry = ARBINT_X86_ADDCARRY(carry, (arbint_x86_carry_word_t) dst[di],
                                  (arbint_x86_carry_word_t) src[i], &out);
      dst[di] = (arbint_limb_t) out;
    }

    i = shift + src_n;
    while (carry != 0u) {
      if (i >= dst_n)
        return 0;
      carry = ARBINT_X86_ADDCARRY(carry, (arbint_x86_carry_word_t) dst[i],
                                  (arbint_x86_carry_word_t) 0, &out);
      dst[i] = (arbint_limb_t) out;
      ++i;
    }
  }
#else
  {
    arbint_limb_t carry = 0u;

    for (i = 0u; i < src_n; ++i) {
      size_t di = shift + i;
      arbint_limb_t d = dst[di];
      arbint_limb_t s = src[i];
      arbint_limb_t t = d + s;
      arbint_limb_t c1 = (t < d) ? 1u : 0u;
      arbint_limb_t t2 = t + carry;
      arbint_limb_t c2 = (t2 < t) ? 1u : 0u;
      dst[di] = t2;
      carry = (c1 | c2);
    }

    i = shift + src_n;
    while (carry != 0u) {
      arbint_limb_t d;
      arbint_limb_t t;

      if (i >= dst_n)
        return 0;
      d = dst[i];
      t = d + carry;
      dst[i] = t;
      carry = (t < d) ? 1u : 0u;
      ++i;
    }
  }
#endif /* ARBINT_HAVE_X86_CARRY_KERNEL */

  return 1;
}

/*  Recursive multiplication using Karatsuba algorithm for large operands.
    Falls back to schoolbook multiplication below ARBINT_KARATSUBA_THRESHOLD.
    Splits operands in half and uses identity: (a1*B + a0) * (b1*B + b0) =
    a1*b1*B^2 + ((a1+a0)*(b1+b0) - a1*b1 - a0*b0)*B + a0*b0.  */
static arbint_err_t arbint_mul_mag_rec(arbint_limb_t * dst, size_t * out_used,
                                       const arbint_limb_t * a, size_t an,
                                       const arbint_limb_t * b, size_t bn,
                                       const arbint_alloc_t * alloc) {
  size_t cap;

  if (out_used == NULL || dst == NULL || a == NULL || b == NULL)
    return ARBINT_EINVAL;

  if (an == 0u || bn == 0u) {
    dst[0] = 0u;
    *out_used = 0u;
    return ARBINT_OK;
  }

  if (an < bn) {
    const arbint_limb_t * tp = a;
    size_t tn = an;
    a = b;
    an = bn;
    b = tp;
    bn = tn;
  }

  if (!arbint_mul_cap(an, bn, &cap))
    return ARBINT_EOVERFLOW;

  if (bn < ARBINT_KARATSUBA_THRESHOLD || an < ARBINT_KARATSUBA_THRESHOLD ||
      bn <= 1u || bn > SIZE_MAX - bn || an > bn + bn) {
    *out_used = arbint_mul_schoolbook(dst, a, an, b, bn);
    return ARBINT_OK;
  }

  {
    size_t k = an / 2u;
    size_t n0 = k;
    size_t n1 = an - k;
    size_t m0 = (bn < k) ? bn : k;
    size_t m1 = bn - m0;
    size_t z0_used = 0u;
    size_t z2_used = 0u;
    size_t sx_cap;
    size_t sy_cap;
    size_t p_cap;
    size_t tmp_cap;
    arbint_limb_t * tmp = NULL;
    arbint_limb_t * sx;
    arbint_limb_t * sy;
    arbint_limb_t * p;
    size_t sx_used;
    size_t sy_used;
    size_t p_used = 0u;
    arbint_err_t rc;

    if (m1 == 0u) {
      *out_used = arbint_mul_schoolbook(dst, a, an, b, bn);
      return ARBINT_OK;
    }

    sx_cap = ((n0 > n1) ? n0 : n1) + 1u;
    sy_cap = ((m0 > m1) ? m0 : m1) + 1u;
    if (!arbint_mul_cap(sx_cap, sy_cap, &p_cap))
      return ARBINT_EOVERFLOW;

    if (sx_cap > SIZE_MAX - sy_cap || sx_cap + sy_cap > SIZE_MAX - p_cap)
      return ARBINT_EOVERFLOW;
    tmp_cap = sx_cap + sy_cap + p_cap;

    tmp = arbint_alloc_limbs(alloc, tmp_cap);
    if (tmp == NULL)
      return ARBINT_ENOMEM;
    sx = tmp;
    sy = sx + sx_cap;
    p = sy + sy_cap;

    memset(dst, 0, cap * sizeof(arbint_limb_t));

    rc = arbint_mul_mag_rec(dst, &z0_used, a, n0, b, m0, alloc);
    if (rc != ARBINT_OK)
      goto cleanup;

    rc = arbint_mul_mag_rec(dst + (2u * k), &z2_used, a + k, n1, b + m0, m1,
                            alloc);
    if (rc != ARBINT_OK)
      goto cleanup;

    sx_used = arbint__add_mag(sx, a, n0, a + k, n1);
    sy_used = arbint__add_mag(sy, b, m0, b + m0, m1);

    rc = arbint_mul_mag_rec(p, &p_used, sx, sx_used, sy, sy_used, alloc);
    if (rc != ARBINT_OK)
      goto cleanup;

    if (arbint_cmp_mag_limbs(p, p_used, dst, z0_used) < 0 ||
        arbint_cmp_mag_limbs(p, p_used, dst + (2u * k), z2_used) < 0) {
      rc = ARBINT_EINVAL;
      goto cleanup;
    }

    p_used = arbint__sub_mag(p, p, p_used, dst, z0_used);
    p_used = arbint__sub_mag(p, p, p_used, dst + (2u * k), z2_used);

    if (!arbint_add_shifted(dst, cap, p, p_used, k)) {
      rc = ARBINT_EOVERFLOW;
      goto cleanup;
    }

    *out_used = arbint_norm_used(dst, cap);
    rc = ARBINT_OK;

  cleanup:
    arbint_free_limbs(alloc, tmp);
    return rc;
  }
}

/*  Multiply multi-limb integer by single limb (dst = a * b).
    Specialized version of multiplication for single-limb multiplier.
    Returns normalized result limb count.  */
size_t arbint_mul_limb_1_generic(arbint_limb_t * dst, const arbint_limb_t * a,
                                 size_t an, arbint_limb_t b) {
  size_t i;
  arbint_limb_t carry = 0u;

  for (i = 0u; i < an; ++i) {
    arbint_limb_t out;
    carry = arbint_muladd_limb(a[i], b, (arbint_limb_t) 0u, carry, &out);
    dst[i] = out;
  }

  if (carry != 0u) {
    dst[an] = carry;
    return an + 1u;
  }

  return arbint_norm_used(dst, an);
}

/*  Generic portable multiplication implementation (rop = a * b).
    Handles aliasing by copying inputs if needed. Uses Karatsuba algorithm
    for large operands (>= ARBINT_KARATSUBA_THRESHOLD limbs).  */
arbint_err_t arbint_mul_impl_generic(arbint_t rop, const arbint_t a,
                                     const arbint_t b) {
  int as;
  int bs;
  int sign;
  size_t an;
  size_t bn;
  size_t cap;
  size_t used;
  arbint_err_t rc;
  const arbint_limb_t * ap;
  const arbint_limb_t * bp;
  const arbint_limb_t * am;
  const arbint_limb_t * bm;
  arbint_limb_t * rp;
  const arbint_alloc_t * alloc;
  arbint_limb_t * a_copy = NULL;
  arbint_limb_t * b_copy = NULL;

  if (rop == NULL || a == NULL || b == NULL)
    return ARBINT_EINVAL;

  as = (a[0]._sz > 0) - (a[0]._sz < 0);
  bs = (b[0]._sz > 0) - (b[0]._sz < 0);
  if (as == 0 || bs == 0) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  an = arbint_abs_sz(a[0]._sz);
  bn = arbint_abs_sz(b[0]._sz);
  if (!arbint_mul_cap(an, bn, &cap))
    return ARBINT_EOVERFLOW;

  rc = arbint_resize(rop, cap);
  if (rc != ARBINT_OK)
    return rc;
  if (rop[0]._ctx == NULL || rop[0]._ctx->a.realloc == NULL)
    return ARBINT_EINVAL;
  alloc = &rop[0]._ctx->a;

  ap = ARBINT_CLIMBS(a);
  bp = ARBINT_CLIMBS(b);
  am = ap;
  bm = bp;

  if (rop == a) {
    a_copy = arbint_alloc_limbs(alloc, an);
    if (a_copy == NULL)
      return ARBINT_ENOMEM;
    memcpy(a_copy, ap, an * sizeof(arbint_limb_t));
    am = a_copy;
  }

  if (rop == b) {
    if (a == b) {
      bm = am;
    } else {
      b_copy = arbint_alloc_limbs(alloc, bn);
      if (b_copy == NULL) {
        arbint_free_limbs(alloc, a_copy);
        return ARBINT_ENOMEM;
      }
      memcpy(b_copy, bp, bn * sizeof(arbint_limb_t));
      bm = b_copy;
    }
  }

  rp = ARBINT_LIMBS(rop);
  rc = arbint_mul_mag_rec(rp, &used, am, an, bm, bn, alloc);

  arbint_free_limbs(alloc, a_copy);
  arbint_free_limbs(alloc, b_copy);
  if (rc != ARBINT_OK)
    return rc;

  sign = (as == bs) ? 1 : -1;
  if (!arbint_set_signed_sz(rop, used, sign))
    return ARBINT_EOVERFLOW;

  return ARBINT_OK;
}
