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

#include "arbint_cpu.h"

typedef arbint_err_t (*arbint_mul_impl_fn_t)(arbint_t rop, const arbint_t a,
                                             const arbint_t b);

typedef arbint_err_t (*arbint_sqr_impl_fn_t)(arbint_t rop, const arbint_t a);

typedef size_t (*arbint_mul_limb_1_fn_t)(arbint_limb_t * dst,
                                         const arbint_limb_t * a, size_t an,
                                         arbint_limb_t b);

/*  Select optimal single-limb multiplication implementation.
    Prefers BMI2 when available for faster wide multiply.  */
static arbint_mul_limb_1_fn_t arbint_select_mul_limb_1(void) {
#if HAS_BMI2_ALWAYS
  return arbint_mul_limb_1_bmi2;
#elif HAS_BMI2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2)
             ? arbint_mul_limb_1_bmi2
             : arbint_mul_limb_1_generic;
#else
  return arbint_mul_limb_1_generic;
#endif /* HAS_BMI2_ALWAYS */
}

/*  Select optimal multi-limb multiplication implementation.
    Prefers BMI2 when available (uses _mulx_u64 for faster multiply).  */
static arbint_mul_impl_fn_t arbint_select_mul_impl(void) {
#if HAS_BMI2_ALWAYS
  return arbint_mul_impl_bmi2;
#elif HAS_BMI2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2)
             ? arbint_mul_impl_bmi2
             : arbint_mul_impl_generic;
#else
  return arbint_mul_impl_generic;
#endif /* HAS_BMI2_ALWAYS */
}

/*  Select optimal squaring implementation.
    Prefers BMI2 when available for faster wide multiply.  */
static arbint_sqr_impl_fn_t arbint_select_sqr_impl(void) {
#if HAS_BMI2_ALWAYS
  return arbint_sqr_impl_bmi2;
#elif HAS_BMI2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2)
             ? arbint_sqr_impl_bmi2
             : arbint_sqr_impl_generic;
#else
  return arbint_sqr_impl_generic;
#endif /* HAS_BMI2_ALWAYS */
}

/*  Compute required capacity for multiplication result with overflow check.
    Returns 1 on success (*out = an + bn + 1), 0 on overflow.  */
int arbint_mul_cap(size_t an, size_t bn, size_t * out) {
  if (out == NULL)
    return 0;
  if (bn > SIZE_MAX - 1u)
    return 0;
  if (an > SIZE_MAX - bn - 1u)
    return 0;
  *out = an + bn + 1u;
  return 1;
}

arbint_err_t arbint_mul(arbint_t rop, const arbint_t a, const arbint_t b) {
  static arbint_mul_impl_fn_t impl = NULL;

  if (impl == NULL)
    impl = arbint_select_mul_impl();

  return impl(rop, a, b);
}

arbint_err_t arbint_sqr(arbint_t rop, const arbint_t a) {
  static arbint_sqr_impl_fn_t impl = NULL;

  if (impl == NULL)
    impl = arbint_select_sqr_impl();

  return impl(rop, a);
}

/*  Count trailing zeros in a 32-bit unsigned integer.
    Precondition: v != 0.  */
static unsigned arbint_ctz32(uint32_t v) {
#if ARBINT_COMPILER_GNU_CLANG
  return (unsigned) __builtin_ctz(v);
#elif ARBINT_COMPILER_MSVC
  {
    unsigned long idx;
    _BitScanForward(&idx, (unsigned long) v);
    return (unsigned) idx;
  }
#else
  {
    unsigned n = 0u;
    while ((v & 1u) == 0u) {
      v >>= 1u;
      ++n;
    }
    return n;
  }
#endif /* ARBINT_COMPILER_GNU_CLANG */
}

/*  Multiply arbint by uint32_t (rop = a * b).
    Optimized path for single-limb multiplier with early exit for 0, 1,
    and powers of two (delegated to shift).  */
arbint_err_t arbint_mul_u32(arbint_t rop, const arbint_t a, uint32_t b) {
  static arbint_mul_limb_1_fn_t impl = NULL;
  int as;
  int sign;
  size_t an;
  size_t cap;
  size_t used;
  arbint_err_t rc;
  const arbint_limb_t * ap;
  arbint_limb_t * rp;

  if (rop == NULL || a == NULL)
    return ARBINT_EINVAL;

  as = (a[0]._sz > 0) - (a[0]._sz < 0);
  if (as == 0 || b == 0u) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  if (b == 1u)
    return arbint_set(rop, a);

  /*  Power-of-two fast path: delegate to left shift.  */
  if ((b & (b - 1u)) == 0u)
    return arbint_shl(rop, a, arbint_ctz32(b));

  an = arbint_abs_sz(a[0]._sz);
  cap = an + 1u;
  if (cap < an)
    return ARBINT_EOVERFLOW;

  rc = arbint_resize(rop, cap);
  if (rc != ARBINT_OK)
    return rc;

  ap = ARBINT_CLIMBS(a);
  rp = ARBINT_LIMBS(rop);

  if (impl == NULL)
    impl = arbint_select_mul_limb_1();

  used = impl(rp, ap, an, (arbint_limb_t) b);

  sign = as;
  if (!arbint_set_signed_sz(rop, used, sign))
    return ARBINT_EOVERFLOW;

  return ARBINT_OK;
}

/*  Multiply arbint by int32_t (rop = a * b).
    Handles sign extraction and delegates to mul_u32 for magnitude.  */
arbint_err_t arbint_mul_i32(arbint_t rop, const arbint_t a, int32_t b) {
  uint32_t mag;
  arbint_err_t rc;

  if (b == 0) {
    if (rop == NULL || a == NULL)
      return ARBINT_EINVAL;
    arbint_zero(rop);
    return ARBINT_OK;
  }

  if (b == 1)
    return arbint_set(rop, a);

  if (b == -1)
    return arbint_neg(rop, a);

  if (b > 0) {
    return arbint_mul_u32(rop, a, (uint32_t) b);
  }

  mag = (uint32_t) (-(b + 1)) + 1u;
  rc = arbint_mul_u32(rop, a, mag);
  if (rc != ARBINT_OK)
    return rc;

  rop[0]._sz = -rop[0]._sz;
  return ARBINT_OK;
}
