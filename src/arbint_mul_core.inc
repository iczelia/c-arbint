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

/*  Shared multiplication core included by arbint_mul_generic.c and
    arbint_mul_bmi2.c.

    Before including this file, the translation unit must define:
      - static arbint_mul_wide_limb(x, y, hi, lo)
      - static arbint_muladd_limb(x, y, acc, carry, out)
      - ARBINT_MUL_LIMB_1_FN   -- name for the single-limb multiply function
      - ARBINT_MUL_IMPL_FN     -- name for the multi-limb multiply function
      - ARBINT_SQR_IMPL_FN     -- name for the squaring function  */

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

/*  Propagate carry starting at dst[pos] upward.
    Used by schoolbook squaring after doubling cross-terms.  */
static void arbint_propagate_carry(arbint_limb_t * dst, size_t pos,
                                   arbint_limb_t carry) {
#if ARBINT_HAVE_X86_CARRY_KERNEL
  arbint_x86_carry_word_t out = (arbint_x86_carry_word_t) 0;
  unsigned char c;

  c = ARBINT_X86_ADDCARRY((unsigned char) 0,
                           (arbint_x86_carry_word_t) dst[pos],
                           (arbint_x86_carry_word_t) carry, &out);
  dst[pos] = (arbint_limb_t) out;
  ++pos;
  while (c != 0u) {
    c = ARBINT_X86_ADDCARRY(c, (arbint_x86_carry_word_t) dst[pos],
                            (arbint_x86_carry_word_t) 0, &out);
    dst[pos] = (arbint_limb_t) out;
    ++pos;
  }
#else
  while (carry != 0u) {
    arbint_limb_t t = dst[pos] + carry;
    dst[pos] = t;
    carry = (t < carry) ? 1u : 0u;
    ++pos;
  }
#endif /* ARBINT_HAVE_X86_CARRY_KERNEL */
}

/*  Classical schoolbook squaring exploiting symmetry.

    For a = sum(a[i] * B^i), a^2 = sum(a[i]*a[j] * B^(i+j)) where the
    off-diagonal terms (i != j) appear twice.  We compute:
      1. Off-diagonal cross products:  for i < j, dst[i+j] += a[i]*a[j]
      2. Double them:                  dst[] *= 2  (left-shift by 1)
      3. Add diagonal squares:         dst[2*i] += a[i]*a[i]

    This halves the number of multiplications vs general schoolbook:
    n*(n-1)/2 + n cross + diagonal instead of n^2.

    Returns normalized result limb count.  */
static size_t arbint_sqr_schoolbook(arbint_limb_t * dst,
                                    const arbint_limb_t * a, size_t an) {
  size_t i;
  size_t n;

  if (an == 0u) {
    dst[0] = 0u;
    return 0u;
  }

  if (an == 1u) {
    arbint_limb_t hi;
    arbint_limb_t lo;
    arbint_mul_wide_limb(a[0], a[0], &hi, &lo);
    dst[0] = lo;
    dst[1] = hi;
    if (hi != 0u)
      return 2u;
    if (lo != 0u)
      return 1u;
    return 0u;
  }

  n = 2u * an + 1u;
  memset(dst, 0, n * sizeof(arbint_limb_t));

  /*  Step 1: accumulate off-diagonal cross products (i < j).  */
  for (i = 0u; i < an - 1u; ++i) {
    size_t j;
    arbint_limb_t carry = 0u;

    for (j = i + 1u; j < an; ++j) {
      arbint_limb_t out;
      carry = arbint_muladd_limb(a[i], a[j], dst[i + j], carry, &out);
      dst[i + j] = out;
    }

    arbint_propagate_carry(dst, i + an, carry);
  }

  /*  Step 2: double the cross products (left-shift by 1 bit).  */
  {
    arbint_limb_t shift_carry = 0u;
    for (i = 0u; i < n; ++i) {
      arbint_limb_t v = dst[i];
      dst[i] = (v << 1u) | shift_carry;
      shift_carry = v >> (ARBINT_LIMB_BITS - 1u);
    }
  }

  /*  Step 3: add diagonal squares a[i]^2 at position 2*i.  */
  for (i = 0u; i < an; ++i) {
    arbint_limb_t hi;
    arbint_limb_t lo;
    arbint_limb_t s;
    arbint_limb_t c;

    arbint_mul_wide_limb(a[i], a[i], &hi, &lo);

    s = dst[2u * i] + lo;
    c = (s < lo) ? 1u : 0u;
    dst[2u * i] = s;

    s = dst[2u * i + 1u] + hi + c;
    c = (s < hi || (c != 0u && s == hi)) ? 1u : 0u;
    dst[2u * i + 1u] = s;

    if (c != 0u)
      arbint_propagate_carry(dst, 2u * i + 2u, c);
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

/*  Recursive squaring using Karatsuba identity for large operands.
    Falls back to schoolbook squaring below ARBINT_KARATSUBA_THRESHOLD.
    Splits a into halves a = a1*B^k + a0 and uses:
      a^2 = a1^2*B^(2k) + ((a0+a1)^2 - a0^2 - a1^2)*B^k + a0^2.
    All three sub-products are themselves squarings, so the optimization
    applies recursively at every level.  */
static arbint_err_t arbint_sqr_mag_rec(arbint_limb_t * dst, size_t * out_used,
                                       const arbint_limb_t * a, size_t an,
                                       const arbint_alloc_t * alloc) {
  size_t cap;

  if (out_used == NULL || dst == NULL || a == NULL)
    return ARBINT_EINVAL;

  if (an == 0u) {
    dst[0] = 0u;
    *out_used = 0u;
    return ARBINT_OK;
  }

  if (!arbint_mul_cap(an, an, &cap))
    return ARBINT_EOVERFLOW;

  if (an < ARBINT_KARATSUBA_THRESHOLD || an <= 1u) {
    *out_used = arbint_sqr_schoolbook(dst, a, an);
    return ARBINT_OK;
  }

  {
    size_t k = an / 2u;
    size_t n0 = k;
    size_t n1 = an - k;
    size_t z0_used = 0u;
    size_t z2_used = 0u;
    size_t sx_cap;
    size_t p_cap;
    size_t tmp_cap;
    arbint_limb_t * tmp = NULL;
    arbint_limb_t * sx;
    arbint_limb_t * p;
    size_t sx_used;
    size_t p_used = 0u;
    arbint_err_t rc;

    sx_cap = ((n0 > n1) ? n0 : n1) + 1u;
    if (!arbint_mul_cap(sx_cap, sx_cap, &p_cap))
      return ARBINT_EOVERFLOW;

    if (sx_cap > SIZE_MAX - p_cap)
      return ARBINT_EOVERFLOW;
    tmp_cap = sx_cap + p_cap;

    tmp = arbint_alloc_limbs(alloc, tmp_cap);
    if (tmp == NULL)
      return ARBINT_ENOMEM;
    sx = tmp;
    p = sx + sx_cap;

    memset(dst, 0, cap * sizeof(arbint_limb_t));

    /*  z0 = a0^2, stored at dst[0..].  */
    rc = arbint_sqr_mag_rec(dst, &z0_used, a, n0, alloc);
    if (rc != ARBINT_OK)
      goto cleanup;

    /*  z2 = a1^2, stored at dst[2k..].  */
    rc = arbint_sqr_mag_rec(dst + (2u * k), &z2_used, a + k, n1, alloc);
    if (rc != ARBINT_OK)
      goto cleanup;

    /*  sx = a0 + a1.  */
    sx_used = arbint__add_mag(sx, a, n0, a + k, n1);

    /*  p = (a0 + a1)^2.  */
    rc = arbint_sqr_mag_rec(p, &p_used, sx, sx_used, alloc);
    if (rc != ARBINT_OK)
      goto cleanup;

    /*  p -= z0.  */
    if (arbint_cmp_mag_limbs(p, p_used, dst, z0_used) < 0 ||
        arbint_cmp_mag_limbs(p, p_used, dst + (2u * k), z2_used) < 0) {
      rc = ARBINT_EINVAL;
      goto cleanup;
    }

    p_used = arbint__sub_mag(p, p, p_used, dst, z0_used);
    /*  p -= z2.  */
    p_used = arbint__sub_mag(p, p, p_used, dst + (2u * k), z2_used);

    /*  dst += p * B^k.  */
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
size_t ARBINT_MUL_LIMB_1_FN(arbint_limb_t * dst, const arbint_limb_t * a,
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

/*  Multi-limb multiplication (rop = a * b).
    Handles aliasing by copying inputs if needed. Uses Karatsuba algorithm
    for large operands (>= ARBINT_KARATSUBA_THRESHOLD limbs).  */
arbint_err_t ARBINT_MUL_IMPL_FN(arbint_t rop, const arbint_t a,
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

/*  Squaring (rop = a^2).
    Handles aliasing by copying input if needed. Uses specialized squaring
    kernels (schoolbook and Karatsuba) that exploit the symmetry of squaring
    for roughly 2x fewer multiplications than general multiplication.  */
arbint_err_t ARBINT_SQR_IMPL_FN(arbint_t rop, const arbint_t a) {
  int as;
  size_t an;
  size_t cap;
  size_t used;
  arbint_err_t rc;
  const arbint_limb_t * ap;
  const arbint_limb_t * am;
  arbint_limb_t * rp;
  const arbint_alloc_t * alloc;
  arbint_limb_t * a_copy = NULL;

  if (rop == NULL || a == NULL)
    return ARBINT_EINVAL;

  as = (a[0]._sz > 0) - (a[0]._sz < 0);
  if (as == 0) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  an = arbint_abs_sz(a[0]._sz);
  if (!arbint_mul_cap(an, an, &cap))
    return ARBINT_EOVERFLOW;

  rc = arbint_resize(rop, cap);
  if (rc != ARBINT_OK)
    return rc;
  if (rop[0]._ctx == NULL || rop[0]._ctx->a.realloc == NULL)
    return ARBINT_EINVAL;
  alloc = &rop[0]._ctx->a;

  ap = ARBINT_CLIMBS(a);
  am = ap;

  if (rop == a) {
    a_copy = arbint_alloc_limbs(alloc, an);
    if (a_copy == NULL)
      return ARBINT_ENOMEM;
    memcpy(a_copy, ap, an * sizeof(arbint_limb_t));
    am = a_copy;
  }

  rp = ARBINT_LIMBS(rop);
  rc = arbint_sqr_mag_rec(rp, &used, am, an, alloc);

  arbint_free_limbs(alloc, a_copy);
  if (rc != ARBINT_OK)
    return rc;

  /*  Result of squaring is always non-negative.  */
  if (!arbint_set_signed_sz(rop, used, 1))
    return ARBINT_EOVERFLOW;

  return ARBINT_OK;
}
