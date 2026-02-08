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

#include "arbint_addsub.h"

#include "config.h"

#include "arbint_cpu.h"

#include <assert.h>
#include <limits.h>
#include <string.h>

/*  Add magnitudes of two multi-limb integers.
    Handles operands in any order (swaps if nx < ny for efficiency).
    Returns result limb count (max(nx,ny) or max(nx,ny)+1 if final carry).  */
size_t arbint__add_mag(arbint_limb_t * dst, const arbint_limb_t * x, size_t nx,
                       const arbint_limb_t * y, size_t ny) {
  size_t i;
#if ARBINT_HAVE_X86_CARRY_KERNEL
  unsigned char carry = 0u;
  arbint_x86_carry_word_t out = (arbint_x86_carry_word_t) 0;
#else
  arbint_limb_t carry = 0u;
#endif /* ARBINT_HAVE_X86_CARRY_KERNEL */

  if (nx < ny) {
    const arbint_limb_t * tp = x;
    size_t tn = nx;
    x = y;
    nx = ny;
    y = tp;
    ny = tn;
  }

#if ARBINT_HAVE_X86_CARRY_KERNEL
  for (i = 0u; i < ny; ++i) {
    carry = ARBINT_X86_ADDCARRY(carry, (arbint_x86_carry_word_t) x[i],
                                (arbint_x86_carry_word_t) y[i], &out);
    dst[i] = (arbint_limb_t) out;
  }
  for (; i < nx; ++i) {
    carry = ARBINT_X86_ADDCARRY(carry, (arbint_x86_carry_word_t) x[i],
                                (arbint_x86_carry_word_t) 0, &out);
    dst[i] = (arbint_limb_t) out;
  }
  if (carry != 0u) {
    dst[nx] = 1u;
    return nx + 1u;
  }
  return nx;
#else
  for (i = 0u; i < ny; ++i) {
    arbint_limb_t xi = x[i];
    arbint_limb_t yi = y[i];
    arbint_limb_t s = xi + yi;
    arbint_limb_t c1 = (s < xi) ? 1u : 0u;
    arbint_limb_t s2 = s + carry;
    arbint_limb_t c2 = (s2 < s) ? 1u : 0u;
    dst[i] = s2;
    carry = (c1 | c2);
  }

  for (; i < nx; ++i) {
    arbint_limb_t xi = x[i];
    arbint_limb_t s = xi + carry;
    dst[i] = s;
    carry = (s < xi) ? 1u : 0u;
  }

  if (carry != 0u) {
    dst[nx] = carry;
    return nx + 1u;
  }
  return nx;
#endif /* ARBINT_HAVE_X86_CARRY_KERNEL */
}

/*  Forward declaration of scalar implementation.  */
static size_t arbint__dbl_mag_scalar(arbint_limb_t * dst,
                                     const arbint_limb_t * x, size_t nx);

/*  Function pointer type for doubling magnitude implementation.  */
typedef size_t (*arbint_dbl_mag_fn_t)(arbint_limb_t * dst,
                                      const arbint_limb_t * x, size_t nx);

/*  Select the optimal doubling implementation based on CPU features.  */
static arbint_dbl_mag_fn_t arbint_select_dbl_mag(void) {
#if HAS_AVX2_ALWAYS
  return arbint__dbl_mag_avx2;
#elif HAS_AVX2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_AVX2)
             ? arbint__dbl_mag_avx2
             : arbint__dbl_mag_scalar;
#else
  return arbint__dbl_mag_scalar;
#endif /* HAS_AVX2_ALWAYS */
}

/*  Scalar implementation of magnitude doubling (original implementation).  */
static size_t arbint__dbl_mag_scalar(arbint_limb_t * dst,
                                     const arbint_limb_t * x, size_t nx) {
#if ARBINT_HAVE_X86_CARRY_KERNEL
  size_t i;
  unsigned char carry = 0u;
  arbint_x86_carry_word_t out = (arbint_x86_carry_word_t) 0;

  for (i = 0u; i < nx; ++i) {
    carry = ARBINT_X86_ADDCARRY(carry, (arbint_x86_carry_word_t) x[i],
                                (arbint_x86_carry_word_t) x[i], &out);
    dst[i] = (arbint_limb_t) out;
  }
  if (carry != 0u) {
    dst[nx] = 1u;
    return nx + 1u;
  }
  return nx;
#else
  size_t i;
  arbint_limb_t carry = 0u;

  for (i = 0u; i < nx; ++i) {
    arbint_limb_t xi = x[i];
    arbint_limb_t d = (arbint_limb_t) (xi << 1);
    dst[i] = d + carry;
    carry = (arbint_limb_t) (xi >> (ARBINT_LIMB_BITS - 1));
  }
  if (carry != 0u) {
    dst[nx] = carry;
    return nx + 1u;
  }
  return nx;
#endif /* ARBINT_HAVE_X86_CARRY_KERNEL */
}

/*  Public doubling function with runtime dispatch and threshold check.  */
size_t arbint__dbl_mag(arbint_limb_t * dst, const arbint_limb_t * x,
                       size_t nx) {
  static arbint_dbl_mag_fn_t impl = NULL;

  /*  One-time initialization: select implementation based on CPU features.  */
  if (impl == NULL)
    impl = arbint_select_dbl_mag();

#if HAS_AVX2
  /*  Threshold check: avoid SIMD overhead for small operands.  */
  if (impl == arbint__dbl_mag_avx2 && nx < ARBINT_DBL_AVX2_THRESHOLD)
    return arbint__dbl_mag_scalar(dst, x, nx);
#endif /* HAS_AVX2 */

  return impl(dst, x, nx);
}

/*  Subtract magnitudes of two multi-limb integers.
    Requires |x| >= |y| and nx >= ny (caller ensures this).
    Returns normalized result limb count (may be less than nx if leading zeros
    created).  */
size_t arbint__sub_mag(arbint_limb_t * dst, const arbint_limb_t * x, size_t nx,
                       const arbint_limb_t * y, size_t ny) {
  assert(nx >= ny);
#if ARBINT_HAVE_X86_CARRY_KERNEL
  size_t i;
  unsigned char borrow = 0u;
  arbint_x86_carry_word_t out = (arbint_x86_carry_word_t) 0;

  for (i = 0u; i < ny; ++i) {
    borrow = ARBINT_X86_SUBBORROW(borrow, (arbint_x86_carry_word_t) x[i],
                                  (arbint_x86_carry_word_t) y[i], &out);
    dst[i] = (arbint_limb_t) out;
  }
  for (; i < nx; ++i) {
    borrow = ARBINT_X86_SUBBORROW(borrow, (arbint_x86_carry_word_t) x[i],
                                  (arbint_x86_carry_word_t) 0, &out);
    dst[i] = (arbint_limb_t) out;
  }
  (void) borrow;
  return arbint_norm_used(dst, nx);
#else
  size_t i;
  arbint_limb_t borrow = 0u;

  for (i = 0u; i < ny; ++i) {
    arbint_limb_t xi = x[i];
    arbint_limb_t yi = y[i];
    arbint_limb_t d = xi - yi;
    arbint_limb_t b1 = (xi < yi) ? 1u : 0u;
    arbint_limb_t d2 = d - borrow;
    arbint_limb_t b2 = (d < borrow) ? 1u : 0u;
    dst[i] = d2;
    borrow = (b1 | b2);
  }

  for (; i < nx; ++i) {
    arbint_limb_t xi = x[i];
    arbint_limb_t d = xi - borrow;
    dst[i] = d;
    borrow = (xi < borrow) ? 1u : 0u;
  }

  return arbint_norm_used(dst, nx);
#endif /* ARBINT_HAVE_X86_CARRY_KERNEL */
}

/*  Compare magnitude of multi-limb integer with uint32_t.
    Returns -1 if |x| < y, 0 if |x| == y, +1 if |x| > y.  */
static int arbint_cmp_mag_u32(const arbint_limb_t * x, size_t nx, uint32_t y) {
  if (nx == 0u)
    return (y == 0u) ? 0 : -1;
  if (nx > 1u)
    return 1;
  if (x[0] < (arbint_limb_t) y)
    return -1;
  if (x[0] > (arbint_limb_t) y)
    return 1;
  return 0;
}

/*  Add uint32_t to magnitude of multi-limb integer.
    Propagates carry and terminates early when carry becomes zero.
    Returns result limb count (nx or nx+1 if final carry).  */
static size_t arbint_add_mag_u32(arbint_limb_t * dst, const arbint_limb_t * x,
                                 size_t nx, uint32_t y) {
  size_t i;
  arbint_limb_t carry = (arbint_limb_t) y;

  for (i = 0u; i < nx; ++i) {
    arbint_limb_t xi = x[i];
    arbint_limb_t s = xi + carry;
    dst[i] = s;
    if (s >= xi) {
      carry = 0u;
      ++i;
      break;
    }
    carry = 1u;
  }

  if (i < nx) {
    if (dst != x)
      memcpy(dst + i, x + i, (nx - i) * sizeof(arbint_limb_t));
    return nx;
  }

  if (carry != 0u) {
    dst[nx] = carry;
    return nx + 1u;
  }
  return nx;
}

/*  Subtract uint32_t from magnitude of multi-limb integer.
    Assumes |x| >= y. Propagates borrow and terminates early when borrow
    becomes zero. Returns normalized result limb count (may be less than nx if
    leading zeros created).  */
static size_t arbint_sub_mag_u32(arbint_limb_t * dst, const arbint_limb_t * x,
                                 size_t nx, uint32_t y) {
  size_t i;
  arbint_limb_t borrow = (arbint_limb_t) y;

  for (i = 0u; i < nx; ++i) {
    arbint_limb_t xi = x[i];
    arbint_limb_t d = xi - borrow;
    dst[i] = d;
    if (xi >= borrow) {
      borrow = 0u;
      ++i;
      break;
    }
    borrow = 1u;
  }

  if (i < nx && dst != x)
    memcpy(dst + i, x + i, (nx - i) * sizeof(arbint_limb_t));

  return arbint_norm_used(dst, nx);
}

arbint_err_t arbint_add(arbint_t rop, const arbint_t a, const arbint_t b) {
  int as, bs, sign;
  size_t an, bn, need, used;
  arbint_err_t rc;
  const arbint_limb_t * ap;
  const arbint_limb_t * bp;
  arbint_limb_t * rp;
  int cmp;

  if (rop == NULL || a == NULL || b == NULL)
    return ARBINT_EINVAL;

  as = (a[0]._sz > 0) - (a[0]._sz < 0);
  bs = (b[0]._sz > 0) - (b[0]._sz < 0);
  if (as == 0)
    return arbint_set(rop, b);
  if (bs == 0)
    return arbint_set(rop, a);
  if (a == b) {
    an = arbint_abs_sz(a[0]._sz);
    if (an == SIZE_MAX)
      return ARBINT_EOVERFLOW;
    need = an + 1u;
    if (rop[0]._cap < need) {
      rc = arbint_resize(rop, need);
      if (rc != ARBINT_OK)
        return rc;
    }
    rp = ARBINT_LIMBS(rop);
    ap = ARBINT_CLIMBS(a);
    used = arbint__dbl_mag(rp, ap, an);
    if (!arbint_set_signed_sz(rop, used, as))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  an = arbint_abs_sz(a[0]._sz);
  bn = arbint_abs_sz(b[0]._sz);
  ap = ARBINT_CLIMBS(a);
  bp = ARBINT_CLIMBS(b);

  if (as == bs) {
    need = (an >= bn) ? an : bn;
    if (need == SIZE_MAX)
      return ARBINT_EOVERFLOW;
    ++need;
    if (rop[0]._cap < need) {
      rc = arbint_resize(rop, need);
      if (rc != ARBINT_OK)
        return rc;
    }
    rp = ARBINT_LIMBS(rop);
    ap = ARBINT_CLIMBS(a);
    bp = ARBINT_CLIMBS(b);
    if (an >= bn)
      used = arbint__add_mag(rp, ap, an, bp, bn);
    else
      used = arbint__add_mag(rp, bp, bn, ap, an);
    if (!arbint_set_signed_sz(rop, used, as))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  cmp = arbint_cmp_mag_limbs(ap, an, bp, bn);
  if (cmp == 0) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  if (cmp > 0) {
    need = an;
    sign = as;
  } else {
    need = bn;
    sign = bs;
  }
  if (rop[0]._cap < need) {
    rc = arbint_resize(rop, need);
    if (rc != ARBINT_OK)
      return rc;
  }

  rp = ARBINT_LIMBS(rop);
  ap = ARBINT_CLIMBS(a);
  bp = ARBINT_CLIMBS(b);
  if (cmp > 0)
    used = arbint__sub_mag(rp, ap, an, bp, bn);
  else
    used = arbint__sub_mag(rp, bp, bn, ap, an);

  if (!arbint_set_signed_sz(rop, used, sign))
    return ARBINT_EOVERFLOW;
  return ARBINT_OK;
}

arbint_err_t arbint_sub(arbint_t rop, const arbint_t a, const arbint_t b) {
  int as, bs, sign;
  size_t an, bn, need, used;
  arbint_err_t rc;
  const arbint_limb_t * ap;
  const arbint_limb_t * bp;
  arbint_limb_t * rp;
  int cmp;

  if (rop == NULL || a == NULL || b == NULL)
    return ARBINT_EINVAL;

  as = (a[0]._sz > 0) - (a[0]._sz < 0);
  bs = (b[0]._sz > 0) - (b[0]._sz < 0);
  if (a == b) {
    arbint_zero(rop);
    return ARBINT_OK;
  }
  if (bs == 0)
    return arbint_set(rop, a);
  if (as == 0) {
    rc = arbint_set(rop, b);
    if (rc != ARBINT_OK)
      return rc;
    rop[0]._sz = -rop[0]._sz;
    return ARBINT_OK;
  }

  an = arbint_abs_sz(a[0]._sz);
  bn = arbint_abs_sz(b[0]._sz);
  ap = ARBINT_CLIMBS(a);
  bp = ARBINT_CLIMBS(b);

  if (as != bs) {
    need = (an >= bn) ? an : bn;
    if (need == SIZE_MAX)
      return ARBINT_EOVERFLOW;
    ++need;
    if (rop[0]._cap < need) {
      rc = arbint_resize(rop, need);
      if (rc != ARBINT_OK)
        return rc;
    }
    rp = ARBINT_LIMBS(rop);
    ap = ARBINT_CLIMBS(a);
    bp = ARBINT_CLIMBS(b);
    if (an >= bn)
      used = arbint__add_mag(rp, ap, an, bp, bn);
    else
      used = arbint__add_mag(rp, bp, bn, ap, an);
    if (!arbint_set_signed_sz(rop, used, as))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  cmp = arbint_cmp_mag_limbs(ap, an, bp, bn);
  if (cmp == 0) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  if (cmp > 0) {
    need = an;
    sign = as;
  } else {
    need = bn;
    sign = -as;
  }
  if (rop[0]._cap < need) {
    rc = arbint_resize(rop, need);
    if (rc != ARBINT_OK)
      return rc;
  }

  rp = ARBINT_LIMBS(rop);
  ap = ARBINT_CLIMBS(a);
  bp = ARBINT_CLIMBS(b);
  if (cmp > 0)
    used = arbint__sub_mag(rp, ap, an, bp, bn);
  else
    used = arbint__sub_mag(rp, bp, bn, ap, an);
  if (!arbint_set_signed_sz(rop, used, sign))
    return ARBINT_EOVERFLOW;
  return ARBINT_OK;
}

arbint_err_t arbint_add_u32(arbint_t rop, const arbint_t a, uint32_t b) {
  int as;
  size_t an, need, used;
  arbint_err_t rc;
  const arbint_limb_t * ap;
  arbint_limb_t * rp;
  int cmp;

  if (rop == NULL || a == NULL)
    return ARBINT_EINVAL;
  if (b == 0u)
    return arbint_set(rop, a);

  as = (a[0]._sz > 0) - (a[0]._sz < 0);
  if (as == 0)
    return arbint_set_u32(rop, b);

  an = arbint_abs_sz(a[0]._sz);
  ap = ARBINT_CLIMBS(a);

  if (as > 0) {
    if (an == SIZE_MAX)
      return ARBINT_EOVERFLOW;
    need = an + 1u;
    if (rop[0]._cap < need) {
      rc = arbint_resize(rop, need);
      if (rc != ARBINT_OK)
        return rc;
    }
    rp = ARBINT_LIMBS(rop);
    ap = ARBINT_CLIMBS(a);
    used = arbint_add_mag_u32(rp, ap, an, b);
    if (!arbint_set_signed_sz(rop, used, 1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  cmp = arbint_cmp_mag_u32(ap, an, b);
  if (cmp == 0) {
    arbint_zero(rop);
    return ARBINT_OK;
  }
  if (cmp > 0) {
    need = an;
    if (rop[0]._cap < need) {
      rc = arbint_resize(rop, need);
      if (rc != ARBINT_OK)
        return rc;
    }
    rp = ARBINT_LIMBS(rop);
    ap = ARBINT_CLIMBS(a);
    used = arbint_sub_mag_u32(rp, ap, an, b);
    if (!arbint_set_signed_sz(rop, used, -1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  return arbint_set_u32(rop, b - (uint32_t) ap[0]);
}

arbint_err_t arbint_sub_u32(arbint_t rop, const arbint_t a, uint32_t b) {
  int as;
  size_t an, need, used;
  arbint_err_t rc;
  const arbint_limb_t * ap;
  arbint_limb_t * rp;
  int cmp;
  uint32_t diff;

  if (rop == NULL || a == NULL)
    return ARBINT_EINVAL;
  if (b == 0u)
    return arbint_set(rop, a);

  as = (a[0]._sz > 0) - (a[0]._sz < 0);
  if (as == 0) {
    rc = arbint_set_u32(rop, b);
    if (rc != ARBINT_OK)
      return rc;
    rop[0]._sz = -rop[0]._sz;
    return ARBINT_OK;
  }

  an = arbint_abs_sz(a[0]._sz);
  ap = ARBINT_CLIMBS(a);

  if (as < 0) {
    if (an == SIZE_MAX)
      return ARBINT_EOVERFLOW;
    need = an + 1u;
    if (rop[0]._cap < need) {
      rc = arbint_resize(rop, need);
      if (rc != ARBINT_OK)
        return rc;
    }
    rp = ARBINT_LIMBS(rop);
    ap = ARBINT_CLIMBS(a);
    used = arbint_add_mag_u32(rp, ap, an, b);
    if (!arbint_set_signed_sz(rop, used, -1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  cmp = arbint_cmp_mag_u32(ap, an, b);
  if (cmp == 0) {
    arbint_zero(rop);
    return ARBINT_OK;
  }
  if (cmp > 0) {
    need = an;
    if (rop[0]._cap < need) {
      rc = arbint_resize(rop, need);
      if (rc != ARBINT_OK)
        return rc;
    }
    rp = ARBINT_LIMBS(rop);
    ap = ARBINT_CLIMBS(a);
    used = arbint_sub_mag_u32(rp, ap, an, b);
    if (!arbint_set_signed_sz(rop, used, 1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  diff = b - (uint32_t) ap[0];
  rc = arbint_set_u32(rop, diff);
  if (rc != ARBINT_OK)
    return rc;
  rop[0]._sz = -rop[0]._sz;
  return ARBINT_OK;
}

arbint_err_t arbint_add_i32(arbint_t rop, const arbint_t a, int32_t b) {
  uint32_t mag;
  if (b >= 0)
    return arbint_add_u32(rop, a, (uint32_t) b);
  mag = (uint32_t) (-(b + 1)) + 1u;
  return arbint_sub_u32(rop, a, mag);
}

arbint_err_t arbint_sub_i32(arbint_t rop, const arbint_t a, int32_t b) {
  uint32_t mag;
  if (b >= 0)
    return arbint_sub_u32(rop, a, (uint32_t) b);
  mag = (uint32_t) (-(b + 1)) + 1u;
  return arbint_add_u32(rop, a, mag);
}
