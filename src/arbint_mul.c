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
#include "arbint_div.h"

typedef arbint_err_t (*arbint_mul_impl_fn_t)(arbint_t rop, const arbint_t a,
                                             const arbint_t b);

typedef arbint_err_t (*arbint_sqr_impl_fn_t)(arbint_t rop, const arbint_t a);

typedef size_t (*arbint_mul_limb_1_fn_t)(arbint_limb_t * dst,
                                         const arbint_limb_t * a, size_t an,
                                         arbint_limb_t b);

typedef arbint_err_t (*arbint_mul_mag_fn_t)(arbint_limb_t * dst,
                                            size_t * out_used,
                                            const arbint_limb_t * a, size_t an,
                                            const arbint_limb_t * b, size_t bn,
                                            const arbint_alloc_t * alloc);

typedef size_t (*arbint_mulacc_fn_t)(arbint_limb_t * dst, size_t dst_n,
                                     size_t dst_cap, const arbint_limb_t * a,
                                     size_t an, const arbint_limb_t * c,
                                     size_t cn);

typedef size_t (*arbint_mulacc_1_fn_t)(arbint_limb_t * dst, size_t dst_n,
                                       size_t dst_cap, const arbint_limb_t * a,
                                       size_t an, arbint_limb_t b);

typedef arbint_err_t (*arbint_tdiv_q_3_fn_t)(arbint_t q, const arbint_t n);

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

/*  Select optimal magnitude multiplication (schoolbook/Karatsuba/Toom-3).  */
static arbint_mul_mag_fn_t arbint_select_mul_mag(void) {
#if HAS_BMI2_ALWAYS
  return arbint_mul_mag_bmi2;
#elif HAS_BMI2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2)
             ? arbint_mul_mag_bmi2
             : arbint_mul_mag_generic;
#else
  return arbint_mul_mag_generic;
#endif /* HAS_BMI2_ALWAYS */
}

/*  Select optimal fused multiply-accumulate (multi-limb).  */
static arbint_mulacc_fn_t arbint_select_mulacc(void) {
#if HAS_BMI2_ALWAYS
  return arbint_mulacc_bmi2;
#elif HAS_BMI2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2)
             ? arbint_mulacc_bmi2
             : arbint_mulacc_generic;
#else
  return arbint_mulacc_generic;
#endif /* HAS_BMI2_ALWAYS */
}

/*  Select optimal fused multiply-accumulate (single limb).  */
static arbint_mulacc_1_fn_t arbint_select_mulacc_1(void) {
#if HAS_BMI2_ALWAYS
  return arbint_mulacc_1_bmi2;
#elif HAS_BMI2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2)
             ? arbint_mulacc_1_bmi2
             : arbint_mulacc_1_generic;
#else
  return arbint_mulacc_1_generic;
#endif /* HAS_BMI2_ALWAYS */
}

/*  Select optimal /3 truncated quotient implementation.
    Prefers BMI2 when available for faster wide multiply in Barrett step.  */
static arbint_tdiv_q_3_fn_t arbint_select_tdiv_q_3(void) {
#if HAS_BMI2_ALWAYS
  return arbint_tdiv_q_3_bmi2;
#elif HAS_BMI2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2)
             ? arbint_tdiv_q_3_bmi2
             : arbint_tdiv_q_3_generic;
#else
  return arbint_tdiv_q_3_generic;
#endif /* HAS_BMI2_ALWAYS */
}

/*  Cached function pointers for addmul/submul dispatch.
    Lazily initialized on first use.  */
static arbint_mul_mag_fn_t g_mul_mag = NULL;
static arbint_mulacc_fn_t g_mulacc = NULL;
static arbint_mulacc_1_fn_t g_mulacc_1 = NULL;
static arbint_mul_limb_1_fn_t g_mul_limb_1 = NULL;

/*  Ensure all addmul dispatch function pointers are initialized.  */
static void arbint_init_addmul_dispatch(void) {
  if (g_mul_mag == NULL)
    g_mul_mag = arbint_select_mul_mag();
  if (g_mulacc == NULL)
    g_mulacc = arbint_select_mulacc();
  if (g_mulacc_1 == NULL)
    g_mulacc_1 = arbint_select_mulacc_1();
  if (g_mul_limb_1 == NULL)
    g_mul_limb_1 = arbint_select_mul_limb_1();
}

/*  Compute required capacity for multiplication result with overflow check.

    Calculates the capacity needed to store the product of n-limb and m-limb
    numbers, which is at most n + m limbs (actually n + m or n + m - 1, but
    we allocate +1 for safety and normalize after the operation).

    The function validates that (an + bn + 1) fits in size_t before computing
    it, preventing silent integer overflow. The checks are performed in an
    order-dependent sequence to avoid underflow in intermediate calculations.

    CRITICAL: Check ordering matters! The condition (an > SIZE_MAX - bn - 1u)
    relies on the previous check (bn > SIZE_MAX - 1u) having succeeded to avoid
    underflow. If bn == SIZE_MAX, then SIZE_MAX - bn == 0, and SIZE_MAX - bn -
    1u would underflow to SIZE_MAX, making the check pass incorrectly. The
    prior check catches bn >= SIZE_MAX, so this is safe.

    Parameters:
      an  - Number of limbs in first operand
      bn  - Number of limbs in second operand
      out - Output pointer for computed capacity (receives an + bn + 1)

    Returns:
      1 on success (*out = an + bn + 1), 0 on overflow or NULL out.

    Precondition: out must be non-NULL (checked).

    Overflow conditions detected:
      - out == NULL
      - bn > SIZE_MAX - 1u
      - an > SIZE_MAX - bn - 1 (i.e., an + bn + 1 would overflow)  */
int arbint_mul_cap(size_t an, size_t bn, size_t * out) {
  if (out == NULL)
    return 0;
  /*  Order-dependent overflow checks. Check bn first to prevent underflow
      in the second condition. See detailed comment above.  */
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

/*  Dispatch /3 truncated quotient to platform-specific implementation.  */
static arbint_err_t arbint_tdiv_q_3(arbint_t q, const arbint_t n) {
  static arbint_tdiv_q_3_fn_t impl = NULL;

  if (impl == NULL)
    impl = arbint_select_tdiv_q_3();

  return impl(q, n);
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

/*  Compute rop = base^exp via binary exponentiation (repeated squaring).

    Uses right-to-left scanning of exponent bits:
      result = 1
      while exp > 0:
        if exp is odd: result *= base
        base *= base
        exp >>= 1

    Edge cases:
      - base^0 = 1 (including 0^0 = 1)
      - base^1 = base
      - 0^exp = 0 for exp > 0

    Aliasing: rop may alias base. We work entirely in temporaries and
    swap the result into rop at the end.  */
arbint_err_t arbint_pow_u32(arbint_t rop, const arbint_t base, uint32_t exp) {
  arbint_ctx_t * ctx;
  arbint_t acc;
  arbint_t b;
  arbint_err_t rc;

  if (rop == NULL || base == NULL)
    return ARBINT_EINVAL;

  /*  base^0 = 1 (including 0^0 = 1).  */
  if (exp == 0u)
    return arbint_set_i32(rop, 1);

  /*  base^1 = base.  */
  if (exp == 1u)
    return arbint_set(rop, base);

  /*  0^exp = 0 for exp > 0.  */
  if (base[0]._sz == 0) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  ctx = rop[0]._ctx;
  if (ctx == NULL)
    ctx = base[0]._ctx;

  rc = arbint_init(acc, ctx);
  if (rc != ARBINT_OK)
    return rc;

  rc = arbint_init(b, ctx);
  if (rc != ARBINT_OK) {
    arbint_clear(acc);
    return rc;
  }

  /*  acc = 1, b = base.  */
  rc = arbint_set_i32(acc, 1);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_set(b, base);
  if (rc != ARBINT_OK)
    goto cleanup;

  while (exp > 1u) {
    if ((exp & 1u) != 0u) {
      rc = arbint_mul(acc, acc, b);
      if (rc != ARBINT_OK)
        goto cleanup;
    }
    rc = arbint_sqr(b, b);
    if (rc != ARBINT_OK)
      goto cleanup;
    exp >>= 1u;
  }

  /*  Final multiply for the remaining bit (exp == 1 here).  */
  rc = arbint_mul(acc, acc, b);
  if (rc != ARBINT_OK)
    goto cleanup;

  arbint_swap(rop, acc);
  rc = ARBINT_OK;

cleanup:
  arbint_clear(b);
  arbint_clear(acc);
  return rc;
}

/*  Local copy of Barrett reciprocal computation for modular exponentiation.

    Given a normalized divisor d (MSB set, d >= beta/2), computes the
    reciprocal v = floor((beta^2 - 1) / d) - beta. Used to precompute
    the reciprocal once and reuse for all reductions in the exponentiation
    loop.

    See arbint_tdiv_generic.c arbint_prepare_barrett() for full algorithm
    documentation.

    Precondition: d must be normalized (d >= 2^(LIMB_BITS-1)).  */
static inline arbint_limb_t arbint_pow_prepare_barrett(arbint_limb_t d) {
  arbint_limb_t d0;
  arbint_limb_t d1;
  arbint_limb_t v;
  arbint_limb_t p;
  arbint_limb_t r;
  arbint_limb_t t;
  arbint_limb_t ql;

  d1 = d >> ARBINT_HALF_BITS;
  d0 = d & ARBINT_HALF_MASK;

  /* Compute high half of reciprocal: qh = ~d / d1 (half-by-half). */
  v = (arbint_limb_t) ((~d) / d1);
  r = ((~d) - v * d1) << ARBINT_HALF_BITS;
  r |= ARBINT_HALF_MASK;

  /* Adjust for d0. */
  p = v * d0;
  if (r < p) {
    v--;
    r += d;
    if (r >= d && r < p) {
      v--;
      r += d;
    }
  }
  r -= p;

  /* Compute low half of reciprocal. */
  t = (r >> ARBINT_HALF_BITS) * v + r;
  ql = (t >> ARBINT_HALF_BITS) + (arbint_limb_t) 1u;

  r = (r << ARBINT_HALF_BITS) + ARBINT_HALF_MASK - ql * d;
  if (r >= (t << ARBINT_HALF_BITS)) {
    ql--;
    r += d;
  }

  v = (v << ARBINT_HALF_BITS) + ql;
  if (r >= d)
    v++;

  return v;
}

/*  Modular exponentiation with u32 modulus: rop = base^exp mod mod.

    Uses precomputed Barrett reduction - computes reciprocal once and
    reuses for all ~2*log2(exp) reduction operations. This is significantly
    faster than calling arbint_tdiv_r_u32() for each reduction, which would
    recompute the reciprocal each time.

    Edge cases:
      - mod == 0: ARBINT_EZERO
      - mod == 1: result is 0 (any integer mod 1 is 0)
      - exp == 0: result is 1 (including 0^0 = 1)
      - base == 0, exp > 0: result is 0

    Truncated division semantics: remainder sign matches base sign for odd
    exponents.

    Aliasing: rop may alias base. We work in temporaries and swap at end.  */
arbint_err_t arbint_pow_u32u32_tmod(arbint_t rop, const arbint_t base,
                                    uint32_t exp, uint32_t mod) {
  arbint_ctx_t * ctx;
  arbint_t acc;
  arbint_t b;
  arbint_err_t rc;
  arbint_limb_t d_norm;
  arbint_limb_t di;
  unsigned shift;

  if (rop == NULL || base == NULL)
    return ARBINT_EINVAL;
  if (mod == 0u)
    return ARBINT_EZERO;

  /* mod == 1: any integer mod 1 is 0. */
  if (mod == 1u) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  /* exp == 0: result is 1 (including 0^0). */
  if (exp == 0u)
    return arbint_set_i32(rop, 1);

  /* base == 0: 0^exp = 0 for exp > 0. */
  if (base[0]._sz == 0) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  /* Precompute Barrett parameters ONCE. This is the key optimization:
     we avoid recomputing shift/d_norm/di for each of the ~2*log2(exp)
     modular reductions. */
  d_norm = (arbint_limb_t) mod;
  shift = arbint_clz_limb(d_norm);
  d_norm <<= shift;
  di = arbint_pow_prepare_barrett(d_norm);

  ctx = rop[0]._ctx;
  if (ctx == NULL)
    ctx = base[0]._ctx;

  rc = arbint_init(acc, ctx);
  if (rc != ARBINT_OK)
    return rc;

  rc = arbint_init(b, ctx);
  if (rc != ARBINT_OK) {
    arbint_clear(acc);
    return rc;
  }

  /* acc = 1, b = base mod mod. */
  rc = arbint_set_i32(acc, 1);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_set(b, base);
  if (rc != ARBINT_OK)
    goto cleanup;

  /* Reduce base mod m first. */
  rc = arbint_mod_u32_barrett(b, d_norm, di, shift);
  if (rc != ARBINT_OK)
    goto cleanup;

  /* Binary exponentiation with Barrett reduction. */
  while (exp > 1u) {
    if ((exp & 1u) != 0u) {
      rc = arbint_mul(acc, acc, b);
      if (rc != ARBINT_OK)
        goto cleanup;
      rc = arbint_mod_u32_barrett(acc, d_norm, di, shift);
      if (rc != ARBINT_OK)
        goto cleanup;
    }
    rc = arbint_sqr(b, b);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_mod_u32_barrett(b, d_norm, di, shift);
    if (rc != ARBINT_OK)
      goto cleanup;
    exp >>= 1u;
  }

  /* Final multiply for remaining bit (exp == 1). */
  rc = arbint_mul(acc, acc, b);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_mod_u32_barrett(acc, d_norm, di, shift);
  if (rc != ARBINT_OK)
    goto cleanup;

  arbint_swap(rop, acc);
  rc = ARBINT_OK;

cleanup:
  arbint_clear(b);
  arbint_clear(acc);
  return rc;
}

/*  Compute rop = floor(sqrt(a)) via Newton's method (integer Heron).

    The iteration x_{n+1} = floor((x_n + floor(a / x_n)) / 2) converges
    quadratically from above. Termination: when x_{n+1} >= x_n, x_n is
    the answer.

    Initial guess: x_0 = 1 << ((nbits(a) + 1) / 2), which is always
    >= floor(sqrt(a)) and at most 2x too large.

    Returns ARBINT_EDOM if a < 0.  */
arbint_err_t arbint_isqrt(arbint_t rop, const arbint_t a) {
  arbint_ctx_t * ctx;
  arbint_t x;
  arbint_t t;
  arbint_err_t rc;
  size_t nbits;

  if (rop == NULL || a == NULL)
    return ARBINT_EINVAL;

  if (a[0]._sz < 0)
    return ARBINT_EDOM;

  /*  isqrt(0) = 0.  */
  if (a[0]._sz == 0) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  ctx = rop[0]._ctx;
  if (ctx == NULL)
    ctx = a[0]._ctx;

  rc = arbint_init(x, ctx);
  if (rc != ARBINT_OK)
    return rc;

  rc = arbint_init(t, ctx);
  if (rc != ARBINT_OK) {
    arbint_clear(x);
    return rc;
  }

  /*  Initial guess: x = 1 << ((nbits(a) + 1) / 2).  */
  nbits = arbint_nbits(a);
  rc = arbint_set_i32(x, 1);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_shl(x, x, (uint32_t) ((nbits + 1u) / 2u));
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Newton iteration: t = (x + a/x) / 2.
      Terminate when t >= x (sequence is monotonically decreasing).  */
  for (;;) {
    rc = arbint_tdiv_q(t, a, x);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_add(t, t, x);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_shr(t, t, 1u);
    if (rc != ARBINT_OK)
      goto cleanup;

    if (arbint_cmp(t, x) >= 0)
      break;

    rc = arbint_set(x, t);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

  rc = arbint_set(rop, x);

cleanup:
  arbint_clear(t);
  arbint_clear(x);
  return rc;
}

/*  Compute rop = floor(a^(1/k)) via integer Newton iteration.

    General Newton iteration for k-th root of a positive number:
      x_{n+1} = floor(((k-1)*x_n + floor(a / x_n^(k-1))) / k)

    Converges quadratically from above when starting with x_0 >= a^(1/k).
    Terminates when x_{n+1} >= x_n, at which point x_n is the answer.

    Initial guess: x_0 = 1 << ceil(nbits(a) / k), which is always
    >= floor(a^(1/k)).

    Special cases:
      - k=0: EDOM (undefined)
      - k=1: rop = a (trivial)
      - k=2: delegates to arbint_isqrt for efficiency
      - a=0: rop = 0
      - a=1: rop = 1
      - a<0: EDOM (API specifies a >= 0)

    Aliasing: rop may alias a. We work in temporaries and assign at the end.  */
arbint_err_t arbint_root(arbint_t rop, const arbint_t a, uint32_t k) {
  arbint_ctx_t * ctx;
  arbint_t x, t, xk1;
  arbint_err_t rc;
  size_t nbits;
  size_t shift;

  /*  Input validation.  */
  if (rop == NULL || a == NULL)
    return ARBINT_EINVAL;

  if (k == 0u)
    return ARBINT_EDOM;

  if (a[0]._sz < 0)
    return ARBINT_EDOM;

  /*  Trivial case: k=1 means rop = a.  */
  if (k == 1u)
    return arbint_set(rop, a);

  /*  Delegate k=2 to the optimized square root.  */
  if (k == 2u)
    return arbint_isqrt(rop, a);

  /*  root_k(0) = 0 for any k > 0.  */
  if (a[0]._sz == 0) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  /*  root_k(1) = 1 for any k > 0.  */
  if (arbint_is_one(a))
    return arbint_set_i32(rop, 1);

  /*  Select context for temporaries.  */
  ctx = rop[0]._ctx;
  if (ctx == NULL)
    ctx = a[0]._ctx;

  /*  Initialize temporaries.  */
  rc = arbint_init(x, ctx);
  if (rc != ARBINT_OK)
    return rc;

  rc = arbint_init(t, ctx);
  if (rc != ARBINT_OK) {
    arbint_clear(x);
    return rc;
  }

  rc = arbint_init(xk1, ctx);
  if (rc != ARBINT_OK) {
    arbint_clear(t);
    arbint_clear(x);
    return rc;
  }

  /*  Initial guess: x = 1 << ceil(nbits(a) / k).
      Overflow-safe computation of ceil(nbits / k).  */
  nbits = arbint_nbits(a);
  if (nbits > SIZE_MAX - (size_t)(k - 1u)) {
    /*  nbits + k - 1 would overflow; compute manually.  */
    shift = nbits / (size_t) k;
    if (nbits % (size_t) k != 0u)
      shift += 1u;
  } else {
    shift = (nbits + (size_t)(k - 1u)) / (size_t) k;
  }

  if (shift > (size_t) UINT32_MAX) {
    rc = ARBINT_EOVERFLOW;
    goto cleanup;
  }

  rc = arbint_set_i32(x, 1);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_shl(x, x, (uint32_t) shift);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Newton iteration loop.
      - k=3 specialized path:
          t = floor((2*x + floor(a / x^2)) / 3)
        Uses squaring and u32-division by 3 fast path.
      - k=4 specialized path:
          t = floor((3*x + floor(a / x^3)) / 4)
        Uses final division by 4 via right shift by 2.
      - generic path:
          t = floor(((k-1)*x + floor(a / x^(k-1))) / k)
      Terminate when t >= x.  */
  for (;;) {
    if (k == 3u) {
      /*  Compute x^2.  */
      rc = arbint_sqr(xk1, x);
      if (rc != ARBINT_OK)
        goto cleanup;

      /*  Compute floor(a / x^2).  */
      rc = arbint_tdiv_q(t, a, xk1);
      if (rc != ARBINT_OK)
        goto cleanup;

      /*  Compute 2*x and reuse xk1.  */
      rc = arbint_shl(xk1, x, 1u);
      if (rc != ARBINT_OK)
        goto cleanup;

      /*  t = 2*x + floor(a / x^2).  */
      rc = arbint_add(t, t, xk1);
      if (rc != ARBINT_OK)
        goto cleanup;

      /*  t = floor(t / 3), specialized fixed-divisor Barrett path.  */
      rc = arbint_tdiv_q_3(t, t);
      if (rc != ARBINT_OK)
        goto cleanup;
    } else {
      /*  Compute x^(k-1).  */
      rc = arbint_pow_u32(xk1, x, k - 1u);
      if (rc != ARBINT_OK)
        goto cleanup;

      /*  Compute floor(a / x^(k-1)).  */
      rc = arbint_tdiv_q(t, a, xk1);
      if (rc != ARBINT_OK)
        goto cleanup;

      /*  Compute (k-1) * x and reuse xk1.  */
      rc = arbint_mul_u32(xk1, x, k - 1u);
      if (rc != ARBINT_OK)
        goto cleanup;

      /*  t = (k-1)*x + floor(a / x^(k-1)).  */
      rc = arbint_add(t, t, xk1);
      if (rc != ARBINT_OK)
        goto cleanup;

      /*  t = floor(t / k).  */
      rc = arbint_tdiv_q_u32(t, t, k);
      if (rc != ARBINT_OK)
        goto cleanup;
    }

    /*  Check termination: if t >= x, we have converged.  */
    if (arbint_cmp(t, x) >= 0)
      break;

    /*  Update x = t for next iteration.  */
    rc = arbint_set(x, t);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

  /*  Store the result.  */
  rc = arbint_set(rop, x);

cleanup:
  arbint_clear(xk1);
  arbint_clear(t);
  arbint_clear(x);
  return rc;
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

  mag = arbint_i32_mag(b);
  rc = arbint_mul_u32(rop, a, mag);
  if (rc != ARBINT_OK)
    return rc;

  rop[0]._sz = -rop[0]._sz;
  return ARBINT_OK;
}

/*  Internal fused multiply-accumulate/subtract implementation.
    sign_flip: 0 for addmul (a += b*c), 1 for submul (a -= b*c).

    Uses platform-dispatched multiplication (BMI2 when available) and
    selects Karatsuba/Toom-3 for large operands.

    For same-sign operands, the schoolbook fused mulacc path avoids a
    temporary allocation. Above the Karatsuba threshold, we multiply
    via the recursive algorithm (which includes Karatsuba and Toom-3)
    into a temporary and then add, since the sub-quadratic algorithm
    wins over the fused O(n^2) accumulate.

    For different-sign operands, we always compute the product into a
    temporary (using the recursive multiplier) and subtract.  */
static arbint_err_t arbint_addsubmul_impl(arbint_t a, const arbint_t b,
                                          const arbint_t c, int sign_flip) {
  int as;
  int bs;
  int cs;
  int ps;
  size_t an;
  size_t bn;
  size_t cn;
  size_t min_bc;
  size_t cap;
  size_t prod_n;
  arbint_err_t rc;
  const arbint_alloc_t * alloc;
  arbint_limb_t * bp;
  arbint_limb_t * cp;
  arbint_limb_t * ap;
  arbint_limb_t * b_copy = NULL;
  arbint_limb_t * c_copy = NULL;
  arbint_limb_t * prod = NULL;

  if (a == NULL || b == NULL || c == NULL)
    return ARBINT_EINVAL;

  arbint_init_addmul_dispatch();

  bs = (b[0]._sz > 0) - (b[0]._sz < 0);
  cs = (c[0]._sz > 0) - (c[0]._sz < 0);
  ps = bs * cs;

  /*  b*c = 0, nothing to add/subtract.  */
  if (ps == 0)
    return ARBINT_OK;

  /*  For subtraction, flip the product sign.  */
  if (sign_flip)
    ps = -ps;

  as = (a[0]._sz > 0) - (a[0]._sz < 0);

  /*  a = 0: result is b*c (or -b*c for submul).  */
  if (as == 0) {
    rc = arbint_mul(a, b, c);
    if (rc != ARBINT_OK)
      return rc;
    if (sign_flip)
      a[0]._sz = -a[0]._sz;
    return ARBINT_OK;
  }

  an = arbint_abs_sz(a[0]._sz);
  bn = arbint_abs_sz(b[0]._sz);
  cn = arbint_abs_sz(c[0]._sz);

  /*  Calculate capacity needed: max(an, bn+cn) + 1 for possible carry.  */
  if (!arbint_mul_cap(bn, cn, &prod_n))
    return ARBINT_EOVERFLOW;
  cap = (prod_n > an) ? prod_n : an;
  if (cap > SIZE_MAX - 1u)
    return ARBINT_EOVERFLOW;
  cap += 1u;

  rc = arbint_resize(a, cap);
  if (rc != ARBINT_OK)
    return rc;
  if (a[0]._ctx == NULL || a[0]._ctx->a.realloc == NULL)
    return ARBINT_EINVAL;
  alloc = &a[0]._ctx->a;

  bp = (arbint_limb_t *) ARBINT_CLIMBS(b);
  cp = (arbint_limb_t *) ARBINT_CLIMBS(c);

  /*  Handle aliasing: copy b if a == b.  */
  if (a == b) {
    b_copy = arbint_alloc_limbs(alloc, bn);
    if (b_copy == NULL)
      return ARBINT_ENOMEM;
    memcpy(b_copy, bp, bn * sizeof(arbint_limb_t));
    bp = b_copy;
  }

  /*  Handle aliasing: copy c if a == c.  */
  if (a == c) {
    if (a == b) {
      cp = bp;
    } else {
      c_copy = arbint_alloc_limbs(alloc, cn);
      if (c_copy == NULL) {
        arbint_free_limbs(alloc, b_copy);
        return ARBINT_ENOMEM;
      }
      memcpy(c_copy, cp, cn * sizeof(arbint_limb_t));
      cp = c_copy;
    }
  }

  ap = ARBINT_LIMBS(a);
  min_bc = (bn < cn) ? bn : cn;

  if (as == ps) {
    /*  Same sign: add magnitudes.
        For small operands, fused schoolbook mulacc avoids a temp alloc.
        For large operands, Karatsuba/Toom-3 into temp + add is faster.  */
    if (min_bc < ARBINT_KARATSUBA_THRESHOLD) {
      size_t used = g_mulacc(ap, an, cap, bp, bn, cp, cn);
      if (!arbint_set_signed_sz(a, used, as)) {
        rc = ARBINT_EOVERFLOW;
        goto cleanup;
      }
    } else {
      size_t used;

      prod = arbint_alloc_limbs(alloc, prod_n);
      if (prod == NULL) {
        rc = ARBINT_ENOMEM;
        goto cleanup;
      }

      rc = g_mul_mag(prod, &prod_n, bp, bn, cp, cn, alloc);
      if (rc != ARBINT_OK)
        goto cleanup;

      used = arbint__add_mag(ap, ap, an, prod, prod_n);
      if (!arbint_set_signed_sz(a, used, as)) {
        rc = ARBINT_EOVERFLOW;
        goto cleanup;
      }
    }
  } else {
    /*  Different signs: compute product to temporary and subtract.  */
    size_t used;
    int cmp;

    prod = arbint_alloc_limbs(alloc, prod_n);
    if (prod == NULL) {
      rc = ARBINT_ENOMEM;
      goto cleanup;
    }

    rc = g_mul_mag(prod, &prod_n, bp, bn, cp, cn, alloc);
    if (rc != ARBINT_OK)
      goto cleanup;

    cmp = arbint_cmp_mag_limbs(ap, an, prod, prod_n);

    if (cmp >= 0) {
      /*  |a| >= |product|: result = |a| - |product|, keep sign of a.  */
      used = arbint__sub_mag(ap, ap, an, prod, prod_n);
      if (!arbint_set_signed_sz(a, used, as)) {
        rc = ARBINT_EOVERFLOW;
        goto cleanup;
      }
    } else {
      /*  |a| < |product|: result = |product| - |a|, use sign of product.  */
      used = arbint__sub_mag(ap, prod, prod_n, ap, an);
      if (!arbint_set_signed_sz(a, used, ps)) {
        rc = ARBINT_EOVERFLOW;
        goto cleanup;
      }
    }
  }

  rc = ARBINT_OK;

cleanup:
  arbint_free_limbs(alloc, prod);
  arbint_free_limbs(alloc, c_copy);
  arbint_free_limbs(alloc, b_copy);
  return rc;
}

/*  Fused multiply-accumulate: a += b * c.  */
arbint_err_t arbint_addmul(arbint_t a, const arbint_t b, const arbint_t c) {
  return arbint_addsubmul_impl(a, b, c, 0);
}

/*  Fused multiply-subtract: a -= b * c.  */
arbint_err_t arbint_submul(arbint_t a, const arbint_t b, const arbint_t c) {
  return arbint_addsubmul_impl(a, b, c, 1);
}

/*  Internal fused multiply-accumulate/subtract with u32 multiplier.
    sign_flip: 0 for addmul (a += b*c), 1 for submul (a -= b*c).
    Uses platform-dispatched single-limb multiply.  */
static arbint_err_t arbint_addsubmul_u32_impl(arbint_t a, const arbint_t b,
                                              uint32_t c, int sign_flip) {
  int as;
  int bs;
  int ps;
  size_t an;
  size_t bn;
  size_t cap;
  arbint_err_t rc;
  const arbint_alloc_t * alloc;
  const arbint_limb_t * bp;
  arbint_limb_t * ap;
  arbint_limb_t * b_copy = NULL;
  arbint_limb_t * prod = NULL;

  if (a == NULL || b == NULL)
    return ARBINT_EINVAL;

  if (c == 0u)
    return ARBINT_OK;

  if (c == 1u)
    return sign_flip ? arbint_sub(a, a, b) : arbint_add(a, a, b);

  bs = (b[0]._sz > 0) - (b[0]._sz < 0);

  if (bs == 0)
    return ARBINT_OK;

  /*  For subtraction, flip the product sign.  */
  ps = sign_flip ? -bs : bs;

  as = (a[0]._sz > 0) - (a[0]._sz < 0);

  /*  a = 0: result is b*c (or -b*c for submul).  */
  if (as == 0) {
    rc = arbint_mul_u32(a, b, c);
    if (rc != ARBINT_OK)
      return rc;
    if (sign_flip)
      a[0]._sz = -a[0]._sz;
    return ARBINT_OK;
  }

  arbint_init_addmul_dispatch();

  an = arbint_abs_sz(a[0]._sz);
  bn = arbint_abs_sz(b[0]._sz);

  /*  Product of bn limbs by single limb has at most bn+1 limbs.  */
  cap = (bn + 1u > an) ? bn + 2u : an + 2u;

  rc = arbint_resize(a, cap);
  if (rc != ARBINT_OK)
    return rc;
  if (a[0]._ctx == NULL || a[0]._ctx->a.realloc == NULL)
    return ARBINT_EINVAL;
  alloc = &a[0]._ctx->a;

  bp = ARBINT_CLIMBS(b);

  /*  Handle aliasing.  */
  if (a == b) {
    b_copy = arbint_alloc_limbs(alloc, bn);
    if (b_copy == NULL)
      return ARBINT_ENOMEM;
    memcpy(b_copy, bp, bn * sizeof(arbint_limb_t));
    bp = b_copy;
  }

  ap = ARBINT_LIMBS(a);

  if (as == ps) {
    /*  Same sign: use fused mulacc.  */
    size_t used = g_mulacc_1(ap, an, cap, bp, bn, (arbint_limb_t) c);
    if (!arbint_set_signed_sz(a, used, as)) {
      rc = ARBINT_EOVERFLOW;
      goto cleanup;
    }
  } else {
    /*  Different signs: compute product and subtract.  */
    size_t prod_n;
    size_t used;
    int cmp;

    prod = arbint_alloc_limbs(alloc, bn + 1u);
    if (prod == NULL) {
      rc = ARBINT_ENOMEM;
      goto cleanup;
    }

    prod_n = g_mul_limb_1(prod, bp, bn, (arbint_limb_t) c);

    cmp = arbint_cmp_mag_limbs(ap, an, prod, prod_n);

    if (cmp >= 0) {
      used = arbint__sub_mag(ap, ap, an, prod, prod_n);
      if (!arbint_set_signed_sz(a, used, as)) {
        rc = ARBINT_EOVERFLOW;
        goto cleanup;
      }
    } else {
      used = arbint__sub_mag(ap, prod, prod_n, ap, an);
      if (!arbint_set_signed_sz(a, used, ps)) {
        rc = ARBINT_EOVERFLOW;
        goto cleanup;
      }
    }
  }

  rc = ARBINT_OK;

cleanup:
  arbint_free_limbs(alloc, prod);
  arbint_free_limbs(alloc, b_copy);
  return rc;
}

/*  Fused multiply-accumulate with u32 multiplier: a += b * c.  */
arbint_err_t arbint_addmul_u32(arbint_t a, const arbint_t b, uint32_t c) {
  return arbint_addsubmul_u32_impl(a, b, c, 0);
}

/*  Fused multiply-subtract with u32 multiplier: a -= b * c.  */
arbint_err_t arbint_submul_u32(arbint_t a, const arbint_t b, uint32_t c) {
  return arbint_addsubmul_u32_impl(a, b, c, 1);
}

/*  Fused multiply-accumulate with i32 multiplier: a += b * c.  */
arbint_err_t arbint_addmul_i32(arbint_t a, const arbint_t b, int32_t c) {
  if (c == 0)
    return ARBINT_OK;

  if (c > 0)
    return arbint_addmul_u32(a, b, (uint32_t) c);

  /*  c < 0: a += b*c = a - b*|c|.  */
  return arbint_submul_u32(a, b, arbint_i32_mag(c));
}

/*  Fused multiply-subtract with i32 multiplier: a -= b * c.  */
arbint_err_t arbint_submul_i32(arbint_t a, const arbint_t b, int32_t c) {
  if (c == 0)
    return ARBINT_OK;

  if (c > 0)
    return arbint_submul_u32(a, b, (uint32_t) c);

  /*  c < 0: a -= b*c = a - b*c = a + b*|c|.  */
  return arbint_addmul_u32(a, b, arbint_i32_mag(c));
}
