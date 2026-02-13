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

/*  Forward declaration for recursive multiplication (used by Toom-3).  */
static arbint_err_t arbint_mul_mag_rec(arbint_limb_t * dst, size_t * out_used,
                                       const arbint_limb_t * a, size_t an,
                                       const arbint_limb_t * b, size_t bn,
                                       const arbint_alloc_t * alloc);

/* ========== Toom-Cook-3 helpers ========== */

/*  Right-shift by 1 bit, in-place. Input must be even.
    Returns normalized result limb count.  */
static size_t arbint_rsh1(arbint_limb_t * x, size_t n) {
  size_t i;

  if (n == 0u)
    return 0u;

  for (i = 0u; i < n - 1u; ++i)
    x[i] = (x[i] >> 1u) | (x[i + 1u] << (ARBINT_LIMB_BITS - 1u));
  x[n - 1u] >>= 1u;

  return arbint_norm_used(x, n);
}

/*  dst = x * 3. Returns used count. dst must have room for n+1 limbs.  */
static size_t arbint_mul3(arbint_limb_t * dst, const arbint_limb_t * x,
                           size_t n) {
  size_t dn;

  if (n == 0u)
    return 0u;

  dn = arbint__dbl_mag(dst, x, n);
  return arbint__add_mag(dst, dst, dn, x, n);
}

/*  dst = x * 12. Returns used count. dst must have room for n+2 limbs.  */
static size_t arbint_mul12(arbint_limb_t * dst, const arbint_limb_t * x,
                            size_t n) {
  size_t dn;

  if (n == 0u)
    return 0u;

  dn = arbint_mul3(dst, x, n);
  dn = arbint__dbl_mag(dst, dst, dn);
  return arbint__dbl_mag(dst, dst, dn);
}

/*  Signed-magnitude addition: dst = a + b (algebraic).
    Signs: 0 = non-negative, 1 = negative. Zero is always non-negative.  */
static void toom3_signed_add(arbint_limb_t * dst, size_t * dst_used,
                              int * dst_neg,
                              const arbint_limb_t * a, size_t an, int a_neg,
                              const arbint_limb_t * b, size_t bn, int b_neg) {
  if (a_neg == b_neg) {
    /*  Same sign: add magnitudes, keep sign.  */
    *dst_used = arbint__add_mag(dst, a, an, b, bn);
    *dst_neg = a_neg;
  } else {
    /*  Different signs: subtract smaller from larger.  */
    int cmp = arbint_cmp_mag_limbs(a, an, b, bn);
    if (cmp >= 0) {
      *dst_used = arbint__sub_mag(dst, a, an, b, bn);
      *dst_neg = a_neg;
    } else {
      *dst_used = arbint__sub_mag(dst, b, bn, a, an);
      *dst_neg = b_neg;
    }
  }

  if (*dst_used == 0u)
    *dst_neg = 0;
}

/*  Signed-magnitude subtraction: dst = a - b (algebraic).
    Implemented as a + (-b).  */
static void toom3_signed_sub(arbint_limb_t * dst, size_t * dst_used,
                              int * dst_neg,
                              const arbint_limb_t * a, size_t an, int a_neg,
                              const arbint_limb_t * b, size_t bn, int b_neg) {
  int b_neg_flipped = (bn == 0u) ? 0 : !b_neg;
  toom3_signed_add(dst, dst_used, dst_neg, a, an, a_neg, b, bn,
                   b_neg_flipped);
}

/*  Toom-Cook-3 multiplication: dst = a * b for large balanced operands.

    Splits each operand into 3 pieces of k = ceil(an/3) limbs, evaluates
    at points {0, 1, -1, 2, inf}, performs 5 recursive multiplications,
    then interpolates and recomposes.

    Precondition: an >= bn, bn >= ARBINT_TOOM3_THRESHOLD, an/bn <= 3.
    dst must have capacity for an + bn + 1 limbs.  */
static arbint_err_t arbint_toom3_mul(arbint_limb_t * dst, size_t * out_used,
                                      const arbint_limb_t * a, size_t an,
                                      const arbint_limb_t * b, size_t bn,
                                      const arbint_alloc_t * alloc) {
  size_t k;
  size_t a0n;
  size_t a1n;
  size_t a2n;
  size_t b0n;
  size_t b1n;
  size_t b2n;
  const arbint_limb_t * a0;
  const arbint_limb_t * a1;
  const arbint_limb_t * a2;
  const arbint_limb_t * b0;
  const arbint_limb_t * b1;
  const arbint_limb_t * b2;
  size_t prod_cap;
  size_t r_cap;
  size_t scr_cap;
  size_t total;
  arbint_limb_t * work = NULL;
  arbint_limb_t * w0;
  arbint_limb_t * w1;
  arbint_limb_t * wm1;
  arbint_limb_t * w2;
  arbint_limb_t * winf;
  arbint_limb_t * r1;
  arbint_limb_t * r2;
  arbint_limb_t * r3;
  arbint_limb_t * scr;
  arbint_limb_t * scr_b;
  size_t w0n;
  size_t w1n;
  size_t wm1n;
  size_t w2n;
  size_t winfn;
  size_t r1n;
  size_t r2n;
  size_t r3n;
  int wm1_neg;
  size_t va1n;
  size_t vb1n;
  size_t va_neg1n;
  size_t vb_neg1n;
  int a_neg1_sign;
  int b_neg1_sign;
  size_t va2n;
  size_t vb2n;
  size_t tn;
  int t_neg;
  size_t scrn;
  size_t cap;
  arbint_err_t rc;

  /*  Splitting: k = ceil(an / 3).  */
  k = (an + 2u) / 3u;

  /*  a = a2 * B^(2k) + a1 * B^k + a0.  */
  a0 = a;
  a0n = (k < an) ? k : an;
  a0n = arbint_norm_used(a0, a0n);

  if (an > k) {
    a1 = a + k;
    a1n = ((2u * k) < an) ? k : (an - k);
    a1n = arbint_norm_used(a1, a1n);
  } else {
    a1 = a;
    a1n = 0u;
  }

  if (an > 2u * k) {
    a2 = a + 2u * k;
    a2n = an - 2u * k;
    a2n = arbint_norm_used(a2, a2n);
  } else {
    a2 = a;
    a2n = 0u;
  }

  /*  b = b2 * B^(2k) + b1 * B^k + b0.  */
  b0 = b;
  b0n = (k < bn) ? k : bn;
  b0n = arbint_norm_used(b0, b0n);

  if (bn > k) {
    b1 = b + k;
    b1n = ((2u * k) < bn) ? k : (bn - k);
    b1n = arbint_norm_used(b1, b1n);
  } else {
    b1 = b;
    b1n = 0u;
  }

  if (bn > 2u * k) {
    b2 = b + 2u * k;
    b2n = bn - 2u * k;
    b2n = arbint_norm_used(b2, b2n);
  } else {
    b2 = b;
    b2n = 0u;
  }

  /*  Degenerate split: if either high piece is empty, the three-way split
      is wasteful. Fall back to Karatsuba via arbint_mul_mag_rec which will
      use schoolbook or Karatsuba as appropriate.  */
  if (a2n == 0u || b2n == 0u) {
    *out_used = arbint_mul_schoolbook(dst, a, an, b, bn);
    return ARBINT_OK;
  }

  /*  Allocate workspace.
      product buffers: 2k+7 limbs (product of two (k+3)-limb values).
      w0/winf: 2k+2 limbs (products of k-limb pieces).
      r1/r2/r3: 2k+3 limbs each.
      scr/scr_b: 2k+7 limbs each (scratch for scalar multiplies).  */
  prod_cap = 2u * k + 7u;
  r_cap = 2u * k + 3u;
  scr_cap = 2u * k + 7u;

  /*  Layout: w0 + w1 + wm1 + w2 + winf + r1 + r2 + r3 + scr + scr_b.
      w0 = 2k+2, w1 = 2k+7, wm1 = 2k+7, w2 = 2k+7, winf = 2k+2,
      r1 = 2k+3, r2 = 2k+3, r3 = 2k+3, scr = 2k+7, scr_b = 2k+7.
      Total = 20k + 48.  */
  total = 20u * k + 48u;

  work = arbint_alloc_limbs(alloc, total);
  if (work == NULL)
    return ARBINT_ENOMEM;

  w0   = work;
  w1   = w0   + (2u * k + 2u);
  wm1  = w1   + prod_cap;
  w2   = wm1  + prod_cap;
  winf = w2   + prod_cap;
  r1   = winf + (2u * k + 2u);
  r2   = r1   + r_cap;
  r3   = r2   + r_cap;
  scr  = r3   + r_cap;
  scr_b = scr  + scr_cap;

  /*  No memset: all buffers are written before being read.  */

  /* ---- Pointwise multiplications (cache-friendly order) ---- */

  /*  w0 = a0 * b0 and winf = a2 * b2 first (smallest inputs).  */
  rc = arbint_mul_mag_rec(w0, &w0n, a0, a0n, b0, b0n, alloc);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_mul_mag_rec(winf, &winfn, a2, a2n, b2, b2n, alloc);
  if (rc != ARBINT_OK)
    goto cleanup;

  /* ---- Evaluation and remaining products ---- */

  /*  Point 1: a(1) = a0 + a1 + a2, b(1) = b0 + b1 + b2.  */
  va1n = arbint__add_mag(scr, a0, a0n, a2, a2n);
  va1n = arbint__add_mag(scr, scr, va1n, a1, a1n);

  vb1n = arbint__add_mag(scr_b, b0, b0n, b2, b2n);
  vb1n = arbint__add_mag(scr_b, scr_b, vb1n, b1, b1n);

  /*  w1 = a(1) * b(1). Skip if either eval is zero.  */
  if (va1n == 0u || vb1n == 0u) {
    w1n = 0u;
  } else {
    rc = arbint_mul_mag_rec(w1, &w1n, scr, va1n, scr_b, vb1n, alloc);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

  /*  Point -1: a(-1) = (a0 + a2) - a1, b(-1) = (b0 + b2) - b1.
      Track signs explicitly.  */
  tn = arbint__add_mag(scr, a0, a0n, a2, a2n);
  if (arbint_cmp_mag_limbs(scr, tn, a1, a1n) >= 0) {
    va_neg1n = arbint__sub_mag(scr, scr, tn, a1, a1n);
    a_neg1_sign = 0;
  } else {
    va_neg1n = arbint__sub_mag(scr, a1, a1n, scr, tn);
    a_neg1_sign = 1;
  }
  if (va_neg1n == 0u)
    a_neg1_sign = 0;

  tn = arbint__add_mag(scr_b, b0, b0n, b2, b2n);
  if (arbint_cmp_mag_limbs(scr_b, tn, b1, b1n) >= 0) {
    vb_neg1n = arbint__sub_mag(scr_b, scr_b, tn, b1, b1n);
    b_neg1_sign = 0;
  } else {
    vb_neg1n = arbint__sub_mag(scr_b, b1, b1n, scr_b, tn);
    b_neg1_sign = 1;
  }
  if (vb_neg1n == 0u)
    b_neg1_sign = 0;

  /*  wm1 = |a(-1)| * |b(-1)|, sign tracked separately.  */
  wm1_neg = a_neg1_sign ^ b_neg1_sign;
  if (va_neg1n == 0u || vb_neg1n == 0u) {
    wm1n = 0u;
    wm1_neg = 0;
  } else {
    rc = arbint_mul_mag_rec(wm1, &wm1n, scr, va_neg1n, scr_b, vb_neg1n,
                            alloc);
    if (rc != ARBINT_OK)
      goto cleanup;
    if (wm1n == 0u)
      wm1_neg = 0;
  }

  /*  Point 2: a(2) = a0 + 2*(a1 + 2*a2), b(2) = b0 + 2*(b1 + 2*b2).  */
  tn = arbint__dbl_mag(scr, a2, a2n);
  tn = arbint__add_mag(scr, scr, tn, a1, a1n);
  tn = arbint__dbl_mag(scr, scr, tn);
  va2n = arbint__add_mag(scr, scr, tn, a0, a0n);

  tn = arbint__dbl_mag(scr_b, b2, b2n);
  tn = arbint__add_mag(scr_b, scr_b, tn, b1, b1n);
  tn = arbint__dbl_mag(scr_b, scr_b, tn);
  vb2n = arbint__add_mag(scr_b, scr_b, tn, b0, b0n);

  /*  w2 = a(2) * b(2). Skip if either eval is zero.  */
  if (va2n == 0u || vb2n == 0u) {
    w2n = 0u;
  } else {
    rc = arbint_mul_mag_rec(w2, &w2n, scr, va2n, scr_b, vb2n, alloc);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

  /* ---- Interpolation ---- */

  /*  r0 = w0 (already in w0).
      r4 = winf (already in winf).  */

  /*  Compute r2 = (w1 + wm1_algebraic - 2*w0 - 2*winf) / 2.
      Accumulate directly into r2 buffer, use scr as temp.  */

  /*  r2 = w1 + wm1_algebraic (signed).  */
  toom3_signed_add(r2, &r2n, &t_neg, w1, w1n, 0, wm1, wm1n, wm1_neg);

  /*  r2 -= 2*w0.  */
  scrn = arbint__dbl_mag(scr, w0, w0n);
  toom3_signed_sub(r2, &r2n, &t_neg, r2, r2n, t_neg, scr, scrn, 0);

  /*  r2 -= 2*winf.  */
  scrn = arbint__dbl_mag(scr, winf, winfn);
  toom3_signed_sub(r2, &r2n, &t_neg, r2, r2n, t_neg, scr, scrn, 0);

  /*  r2 numerator must be non-negative and even.  */
  assert(r2n == 0u || t_neg == 0);
  assert(r2n == 0u || (r2[0] & 1u) == 0u);

  /*  r2 = t / 2 (exact).  */
  r2n = arbint_rsh1(r2, r2n);

  /*  Compute r3 = (w2 - 3*w1 - wm1_algebraic + 3*w0 - 12*winf) / 6.
      Accumulate directly into r3 buffer, use scr as temp.  */

  /*  r3 = w2 + 3*w0.  */
  scrn = arbint_mul3(scr, w0, w0n);
  toom3_signed_add(r3, &r3n, &t_neg, w2, w2n, 0, scr, scrn, 0);

  /*  r3 -= 3*w1.  */
  scrn = arbint_mul3(scr, w1, w1n);
  toom3_signed_sub(r3, &r3n, &t_neg, r3, r3n, t_neg, scr, scrn, 0);

  /*  r3 -= wm1_algebraic.  */
  toom3_signed_sub(r3, &r3n, &t_neg, r3, r3n, t_neg,
                   wm1, wm1n, wm1_neg);

  /*  r3 -= 12*winf.  */
  scrn = arbint_mul12(scr, winf, winfn);
  toom3_signed_sub(r3, &r3n, &t_neg, r3, r3n, t_neg, scr, scrn, 0);

  /*  r3 numerator must be non-negative and even.  */
  assert(r3n == 0u || t_neg == 0);
  assert(r3n == 0u || (r3[0] & 1u) == 0u);

  /*  r3 = t / 6 = (t / 2) / 3.  */
  r3n = arbint_rsh1(r3, r3n);
  r3n = ARBINT_DIVEXACT3_FN(r3, r3n);

  /*  Compute r1 = w1 - w0 - r2 - r3 - winf.
      All subtractions are provably non-negative after correct r2, r3.  */
  assert(arbint_cmp_mag_limbs(w1, w1n, w0, w0n) >= 0);
  r1n = arbint__sub_mag(r1, w1, w1n, w0, w0n);

  assert(arbint_cmp_mag_limbs(r1, r1n, r2, r2n) >= 0);
  r1n = arbint__sub_mag(r1, r1, r1n, r2, r2n);

  assert(arbint_cmp_mag_limbs(r1, r1n, r3, r3n) >= 0);
  r1n = arbint__sub_mag(r1, r1, r1n, r3, r3n);

  assert(arbint_cmp_mag_limbs(r1, r1n, winf, winfn) >= 0);
  r1n = arbint__sub_mag(r1, r1, r1n, winf, winfn);

  /* ---- Recomposition ---- */

  /*  result = r0 + r1*B^k + r2*B^(2k) + r3*B^(3k) + r4*B^(4k).  */
  cap = an + bn + 1u;
  memset(dst, 0, cap * sizeof(arbint_limb_t));

  /*  r0 = w0.  */
  if (w0n > 0u)
    memcpy(dst, w0, w0n * sizeof(arbint_limb_t));

  /*  r4 = winf at offset 4k.  */
  if (winfn > 0u)
    memcpy(dst + 4u * k, winf, winfn * sizeof(arbint_limb_t));

  if (!arbint_add_shifted(dst, cap, r1, r1n, k)) {
    rc = ARBINT_EOVERFLOW;
    goto cleanup;
  }

  if (!arbint_add_shifted(dst, cap, r2, r2n, 2u * k)) {
    rc = ARBINT_EOVERFLOW;
    goto cleanup;
  }

  if (!arbint_add_shifted(dst, cap, r3, r3n, 3u * k)) {
    rc = ARBINT_EOVERFLOW;
    goto cleanup;
  }

  *out_used = arbint_norm_used(dst, cap);
  rc = ARBINT_OK;

cleanup:
  arbint_free_limbs(alloc, work);
  return rc;
}

/*  Forward declaration for recursive squaring (used by Toom-3 squaring).  */
static arbint_err_t arbint_sqr_mag_rec(arbint_limb_t * dst, size_t * out_used,
                                       const arbint_limb_t * a, size_t an,
                                       const arbint_alloc_t * alloc);

/*  Toom-Cook-3 squaring: dst = a^2 for large operands.

    Splits a into 3 pieces of k = ceil(an/3) limbs, evaluates at points
    {0, 1, -1, 2, inf}, performs 5 recursive squarings, then interpolates
    and recomposes.

    Key simplification vs multiplication: at point -1, we compute
    (a0 - a1 + a2)^2 = |a0 - a1 + a2|^2, so sign tracking is unnecessary.

    Precondition: an >= ARBINT_SQR_TOOM3_THRESHOLD.
    dst must have capacity for 2*an + 1 limbs.  */
static arbint_err_t arbint_toom3_sqr(arbint_limb_t * dst, size_t * out_used,
                                      const arbint_limb_t * a, size_t an,
                                      const arbint_alloc_t * alloc) {
  size_t k;
  size_t a0n;
  size_t a1n;
  size_t a2n;
  const arbint_limb_t * a0;
  const arbint_limb_t * a1;
  const arbint_limb_t * a2;
  size_t prod_cap;
  size_t r_cap;
  size_t scr_cap;
  size_t total;
  arbint_limb_t * work = NULL;
  arbint_limb_t * w0;
  arbint_limb_t * w1;
  arbint_limb_t * wm1;
  arbint_limb_t * w2;
  arbint_limb_t * winf;
  arbint_limb_t * r1;
  arbint_limb_t * r2;
  arbint_limb_t * r3;
  arbint_limb_t * scr;
  size_t w0n;
  size_t w1n;
  size_t wm1n;
  size_t w2n;
  size_t winfn;
  size_t r1n;
  size_t r2n;
  size_t r3n;
  size_t va1n;
  size_t va_neg1n;
  size_t va2n;
  size_t tn;
  int t_neg;
  size_t scrn;
  size_t cap;
  arbint_err_t rc;

  /*  Splitting: k = ceil(an / 3).  */
  k = (an + 2u) / 3u;

  /*  a = a2 * B^(2k) + a1 * B^k + a0.  */
  a0 = a;
  a0n = (k < an) ? k : an;
  a0n = arbint_norm_used(a0, a0n);

  if (an > k) {
    a1 = a + k;
    a1n = ((2u * k) < an) ? k : (an - k);
    a1n = arbint_norm_used(a1, a1n);
  } else {
    a1 = a;
    a1n = 0u;
  }

  if (an > 2u * k) {
    a2 = a + 2u * k;
    a2n = an - 2u * k;
    a2n = arbint_norm_used(a2, a2n);
  } else {
    a2 = a;
    a2n = 0u;
  }

  /*  Degenerate split: if high piece is empty, fall back to Karatsuba.  */
  if (a2n == 0u) {
    *out_used = arbint_sqr_schoolbook(dst, a, an);
    return ARBINT_OK;
  }

  /*  Allocate workspace.
      product buffers: 2k+7 limbs (square of (k+3)-limb value).
      w0/winf: 2k+2 limbs (squares of k-limb pieces).
      r1/r2/r3: 2k+3 limbs each.
      scr: 2k+7 limbs (scratch for evaluation).

      Layout: w0 + w1 + wm1 + w2 + winf + r1 + r2 + r3 + scr.
      w0 = 2k+2, w1 = 2k+7, wm1 = 2k+7, w2 = 2k+7, winf = 2k+2,
      r1 = 2k+3, r2 = 2k+3, r3 = 2k+3, scr = 2k+7.
      Total = 18k + 41.  */
  prod_cap = 2u * k + 7u;
  r_cap = 2u * k + 3u;
  scr_cap = 2u * k + 7u;
  total = 18u * k + 41u;

  work = arbint_alloc_limbs(alloc, total);
  if (work == NULL)
    return ARBINT_ENOMEM;

  w0   = work;
  w1   = w0   + (2u * k + 2u);
  wm1  = w1   + prod_cap;
  w2   = wm1  + prod_cap;
  winf = w2   + prod_cap;
  r1   = winf + (2u * k + 2u);
  r2   = r1   + r_cap;
  r3   = r2   + r_cap;
  scr  = r3   + r_cap;

  /* ---- Pointwise squarings ---- */

  /*  w0 = a0^2 and winf = a2^2 first (smallest inputs).  */
  rc = arbint_sqr_mag_rec(w0, &w0n, a0, a0n, alloc);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_sqr_mag_rec(winf, &winfn, a2, a2n, alloc);
  if (rc != ARBINT_OK)
    goto cleanup;

  /* ---- Evaluation and remaining squarings ---- */

  /*  Point 1: a(1) = a0 + a1 + a2.  */
  va1n = arbint__add_mag(scr, a0, a0n, a2, a2n);
  va1n = arbint__add_mag(scr, scr, va1n, a1, a1n);

  /*  w1 = a(1)^2.  */
  if (va1n == 0u) {
    w1n = 0u;
  } else {
    rc = arbint_sqr_mag_rec(w1, &w1n, scr, va1n, alloc);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

  /*  Point -1: a(-1) = a0 - a1 + a2 = (a0 + a2) - a1.
      For squaring, sign doesn't matter: (a(-1))^2 = |a(-1)|^2.  */
  tn = arbint__add_mag(scr, a0, a0n, a2, a2n);
  if (arbint_cmp_mag_limbs(scr, tn, a1, a1n) >= 0) {
    va_neg1n = arbint__sub_mag(scr, scr, tn, a1, a1n);
  } else {
    va_neg1n = arbint__sub_mag(scr, a1, a1n, scr, tn);
  }

  /*  wm1 = |a(-1)|^2.  */
  if (va_neg1n == 0u) {
    wm1n = 0u;
  } else {
    rc = arbint_sqr_mag_rec(wm1, &wm1n, scr, va_neg1n, alloc);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

  /*  Point 2: a(2) = a0 + 2*a1 + 4*a2 = a0 + 2*(a1 + 2*a2).  */
  tn = arbint__dbl_mag(scr, a2, a2n);
  tn = arbint__add_mag(scr, scr, tn, a1, a1n);
  tn = arbint__dbl_mag(scr, scr, tn);
  va2n = arbint__add_mag(scr, scr, tn, a0, a0n);

  /*  w2 = a(2)^2.  */
  if (va2n == 0u) {
    w2n = 0u;
  } else {
    rc = arbint_sqr_mag_rec(w2, &w2n, scr, va2n, alloc);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

  /* ---- Interpolation ---- */

  /*  r0 = w0 (already in w0).
      r4 = winf (already in winf).

      The interpolation formulas are the same as for multiplication.
      For squaring, wm1 is always non-negative (it's a square).  */

  /*  Compute r2 = (w1 + wm1 - 2*w0 - 2*winf) / 2.  */
  r2n = arbint__add_mag(r2, w1, w1n, wm1, wm1n);
  t_neg = 0;

  scrn = arbint__dbl_mag(scr, w0, w0n);
  toom3_signed_sub(r2, &r2n, &t_neg, r2, r2n, t_neg, scr, scrn, 0);

  scrn = arbint__dbl_mag(scr, winf, winfn);
  toom3_signed_sub(r2, &r2n, &t_neg, r2, r2n, t_neg, scr, scrn, 0);

  assert(r2n == 0u || t_neg == 0);
  assert(r2n == 0u || (r2[0] & 1u) == 0u);

  r2n = arbint_rsh1(r2, r2n);

  /*  Compute r3 = (w2 - 3*w1 - wm1 + 3*w0 - 12*winf) / 6.
      Note: wm1 is subtracted (not added with sign) because for squaring
      the product at -1 is always positive.  */
  scrn = arbint_mul3(scr, w0, w0n);
  toom3_signed_add(r3, &r3n, &t_neg, w2, w2n, 0, scr, scrn, 0);

  scrn = arbint_mul3(scr, w1, w1n);
  toom3_signed_sub(r3, &r3n, &t_neg, r3, r3n, t_neg, scr, scrn, 0);

  toom3_signed_sub(r3, &r3n, &t_neg, r3, r3n, t_neg, wm1, wm1n, 0);

  scrn = arbint_mul12(scr, winf, winfn);
  toom3_signed_sub(r3, &r3n, &t_neg, r3, r3n, t_neg, scr, scrn, 0);

  assert(r3n == 0u || t_neg == 0);
  assert(r3n == 0u || (r3[0] & 1u) == 0u);

  r3n = arbint_rsh1(r3, r3n);
  r3n = ARBINT_DIVEXACT3_FN(r3, r3n);

  /*  Compute r1 = w1 - w0 - r2 - r3 - winf.  */
  assert(arbint_cmp_mag_limbs(w1, w1n, w0, w0n) >= 0);
  r1n = arbint__sub_mag(r1, w1, w1n, w0, w0n);

  assert(arbint_cmp_mag_limbs(r1, r1n, r2, r2n) >= 0);
  r1n = arbint__sub_mag(r1, r1, r1n, r2, r2n);

  assert(arbint_cmp_mag_limbs(r1, r1n, r3, r3n) >= 0);
  r1n = arbint__sub_mag(r1, r1, r1n, r3, r3n);

  assert(arbint_cmp_mag_limbs(r1, r1n, winf, winfn) >= 0);
  r1n = arbint__sub_mag(r1, r1, r1n, winf, winfn);

  /* ---- Recomposition ---- */

  cap = 2u * an + 1u;
  memset(dst, 0, cap * sizeof(arbint_limb_t));

  if (w0n > 0u)
    memcpy(dst, w0, w0n * sizeof(arbint_limb_t));

  if (winfn > 0u)
    memcpy(dst + 4u * k, winf, winfn * sizeof(arbint_limb_t));

  if (!arbint_add_shifted(dst, cap, r1, r1n, k)) {
    rc = ARBINT_EOVERFLOW;
    goto cleanup;
  }

  if (!arbint_add_shifted(dst, cap, r2, r2n, 2u * k)) {
    rc = ARBINT_EOVERFLOW;
    goto cleanup;
  }

  if (!arbint_add_shifted(dst, cap, r3, r3n, 3u * k)) {
    rc = ARBINT_EOVERFLOW;
    goto cleanup;
  }

  *out_used = arbint_norm_used(dst, cap);
  rc = ARBINT_OK;

cleanup:
  arbint_free_limbs(alloc, work);
  return rc;
}

/* ========== End Toom-Cook-3 helpers ========== */

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

  /*  NTT for very large balanced operands.
      bn is the smaller operand (swap above guarantees an >= bn).
      Ratio check: an <= 2*bn ensures reasonably balanced operands.
      For highly unbalanced products, Toom-3 with unbalanced splitting is better.  */
  if (bn >= ARBINT_NTT_THRESHOLD && an <= 2u * bn)
    return arbint_mul_mag_ntt(dst, out_used, a, an, b, bn, alloc);

  /*  Toom-Cook-3 for large balanced operands (an <= 3*bn).
      Use subtraction form to avoid both division and multiplication overflow:
      an >= bn is guaranteed by the swap above, so an - bn >= 0.
      bn >= 96 so bn + bn cannot overflow size_t.  */
  if (bn >= ARBINT_TOOM3_THRESHOLD && an - bn <= bn + bn)
    return arbint_toom3_mul(dst, out_used, a, an, b, bn, alloc);

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
    Falls back to schoolbook squaring below ARBINT_SQR_KARATSUBA_THRESHOLD.
    Uses Toom-Cook-3 for operands >= ARBINT_SQR_TOOM3_THRESHOLD.
    Uses NTT for very large operands >= ARBINT_NTT_THRESHOLD.
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

  if (an < ARBINT_SQR_KARATSUBA_THRESHOLD || an <= 1u) {
    *out_used = arbint_sqr_schoolbook(dst, a, an);
    return ARBINT_OK;
  }

  /*  NTT for very large operands.
      NTT provides O(n log n) complexity vs Toom-3's O(n^1.465).
      Dispatch when an >= NTT threshold.  */
  if (an >= ARBINT_NTT_THRESHOLD)
    return arbint_sqr_mag_ntt(dst, out_used, a, an, alloc);

  /*  Toom-Cook-3 for large operands.  */
  if (an >= ARBINT_SQR_TOOM3_THRESHOLD)
    return arbint_toom3_sqr(dst, out_used, a, an, alloc);

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

/*  Fused multiply-accumulate for single limb multiplier.
    Computes dst[0..dst_n-1] += a[0..an-1] * b in place.
    dst must have capacity for at least max(dst_n, an) + 1 limbs.
    Returns new normalized limb count.

    The function accumulates the product directly into the destination,
    avoiding a temporary allocation for the intermediate product.  */
static size_t arbint_mulacc_1(arbint_limb_t * dst, size_t dst_n, size_t dst_cap,
                               const arbint_limb_t * a, size_t an,
                               arbint_limb_t b) {
  size_t i;
  arbint_limb_t carry = 0u;
  size_t new_n;

  if (an == 0u || b == 0u)
    return dst_n;

  /*  Zero-extend dst if the product may extend beyond current size.  */
  for (i = dst_n; i < an + 1u && i < dst_cap; ++i)
    dst[i] = 0u;

  /*  Core loop: for each limb of a, accumulate a[i] * b into dst[i].  */
  for (i = 0u; i < an; ++i) {
    arbint_limb_t out;
    arbint_limb_t d = (i < dst_cap) ? dst[i] : 0u;
    carry = arbint_muladd_limb(a[i], b, d, carry, &out);
    if (i < dst_cap)
      dst[i] = out;
  }

  /*  Propagate remaining carry through higher limbs.  */
  while (carry != 0u && i < dst_cap) {
    arbint_limb_t sum = dst[i] + carry;
    carry = (sum < dst[i]) ? 1u : 0u;
    dst[i] = sum;
    ++i;
  }

  new_n = (i > dst_n) ? i : dst_n;
  return arbint_norm_used(dst, new_n);
}

/*  Fused multiply-accumulate for multi-limb multiplier.
    Computes dst[0..] += a[0..an-1] * c[0..cn-1] in place.
    dst must have capacity for at least max(dst_n, an + cn) + 1 limbs.
    Returns new normalized limb count.

    Uses schoolbook algorithm: for each limb c[j], multiply a by c[j]
    and add at offset j. This is more efficient than computing the full
    product and adding it afterward because we avoid an intermediate
    allocation.  */
static size_t arbint_mulacc(arbint_limb_t * dst, size_t dst_n, size_t dst_cap,
                             const arbint_limb_t * a, size_t an,
                             const arbint_limb_t * c, size_t cn) {
  size_t j;
  size_t new_n;

  if (an == 0u || cn == 0u)
    return dst_n;

  /*  Zero-extend dst to accommodate the full product.  */
  for (j = dst_n; j < dst_cap; ++j)
    dst[j] = 0u;

  /*  For each limb of c, accumulate a * c[j] at offset j.  */
  for (j = 0u; j < cn; ++j) {
    size_t i;
    arbint_limb_t carry = 0u;

    for (i = 0u; i < an; ++i) {
      arbint_limb_t out;
      size_t pos = i + j;
      arbint_limb_t d = (pos < dst_cap) ? dst[pos] : 0u;
      carry = arbint_muladd_limb(a[i], c[j], d, carry, &out);
      if (pos < dst_cap)
        dst[pos] = out;
    }

    /*  Propagate carry through higher limbs.  */
    for (i = j + an; carry != 0u && i < dst_cap; ++i) {
      arbint_limb_t sum = dst[i] + carry;
      carry = (sum < dst[i]) ? 1u : 0u;
      dst[i] = sum;
    }
  }

  new_n = an + cn;
  if (new_n < dst_n)
    new_n = dst_n;
  if (new_n > dst_cap)
    new_n = dst_cap;

  return arbint_norm_used(dst, new_n);
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
