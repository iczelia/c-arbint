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

#include "arbint_gcd.h"

#include "arbint_addsub.h"
#include "arbint_cpu.h"
#include "arbint_div.h"
#include "arbint_internal_util.h"
#include "arbint_mul.h"
#include "arbint_shift.h"
#include "config.h"

#include <string.h>

/*  Lehmer GCD dispatch function pointer type.

    Lehmer GCD has platform-specific implementations that use different
    multiply primitives. The dispatch selects the optimal implementation
    based on CPU features detected at runtime.  */
typedef arbint_err_t (*arbint_gcd_lehmer_fn_t)(
    arbint_t g, const arbint_limb_t * ap, size_t an, const arbint_limb_t * bp,
    size_t bn, const arbint_alloc_t * alloc);

/*  Select optimal Lehmer GCD implementation based on CPU features.

    Three-tier dispatch pattern (same as multiplication):
    1. HAS_BMI2_ALWAYS: library compiled with -mbmi2, always use BMI2
    2. HAS_BMI2: BMI2 code compiled separately, runtime CPUID check
    3. Neither: only generic implementation available

    BMI2 provides _mulx_u64 for fast 64x64->128 multiply (~3x faster
    than half-limb decomposition used in generic path).  */
static arbint_gcd_lehmer_fn_t arbint_select_gcd_lehmer(void) {
#if HAS_BMI2_ALWAYS
  return arbint_gcd_lehmer_bmi2;
#elif HAS_BMI2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2)
             ? arbint_gcd_lehmer_bmi2
             : arbint_gcd_lehmer_generic;
#else
  return arbint_gcd_lehmer_generic;
#endif
}

/*  Cached Lehmer GCD function pointer (lazily initialized).

    First call to arbint_gcd_lehmer_dispatch triggers CPU feature
    detection and caches the result. Subsequent calls use the cached
    pointer directly (no synchronization needed for pointer reads).  */
static arbint_gcd_lehmer_fn_t g_gcd_lehmer = NULL;

/*  Dispatch wrapper for Lehmer GCD.

    Called from arbint_gcd when min(an, bn) >= ARBINT_LEHMER_THRESHOLD.
    Performs lazy initialization of the function pointer on first call.  */
static arbint_err_t
arbint_gcd_lehmer_dispatch(arbint_t g, const arbint_limb_t * ap, size_t an,
                           const arbint_limb_t * bp, size_t bn,
                           const arbint_alloc_t * alloc) {
  ARBINT_LAZY_INIT(g_gcd_lehmer, arbint_select_gcd_lehmer);
  return g_gcd_lehmer(g, ap, an, bp, bn, alloc);
}

/*  Binary GCD (Stein's algorithm) for medium-size operands.

    Used when ARBINT_GCD_EUCLID_THRESHOLD < min_n < ARBINT_LEHMER_THRESHOLD
    (i.e., 5-31 limbs on 64-bit). For this size range, Stein's algorithm
    outperforms Euclidean GCD (avoids division) but doesn't justify
    Lehmer's matrix overhead.

    Algorithm: repeatedly subtract smaller from larger, stripping
    trailing zeros after each subtraction. Exploits the facts that:
    1. gcd(2a, 2b) = 2 * gcd(a, b)
    2. gcd(2a, b) = gcd(a, b) when b is odd
    3. gcd(a, b) = gcd(a-b, b) when a > b and both odd

    Uses internal limb-level operations for performance.
    Single combined allocation for u and v workspace.  */
static arbint_err_t arbint_gcd_binary(arbint_t g, const arbint_limb_t * ap,
                                      size_t an, const arbint_limb_t * bp,
                                      size_t bn,
                                      const arbint_alloc_t * alloc) {
  arbint_limb_t * scratch = NULL;
  arbint_limb_t * up;
  arbint_limb_t * vp;
  arbint_limb_t * tmp;
  size_t un;
  size_t vn;
  size_t max_n;
  size_t ctz_u;
  size_t ctz_v;
  size_t common;
  size_t cap;
  arbint_err_t rc = ARBINT_OK;

  max_n = (an > bn) ? an : bn;

  /*  Allocate combined workspace for u and v.
      Each needs cap = max_n + 1 limbs for final left shift.  */
  cap = max_n + 1u;
  scratch = arbint_alloc_limbs(alloc, cap * 2u);
  if (scratch == NULL)
    return ARBINT_ENOMEM;

  up = scratch;
  vp = scratch + cap;

  /*  Copy magnitudes to workspace.  */
  memcpy(up, ap, an * sizeof(arbint_limb_t));
  memcpy(vp, bp, bn * sizeof(arbint_limb_t));
  un = an;
  vn = bn;

  /*  Extract common power of 2.  */
  ctz_u = arbint_gcd_mag_ctz(up, un);
  ctz_v = arbint_gcd_mag_ctz(vp, vn);
  common = (ctz_u < ctz_v) ? ctz_u : ctz_v;

  /*  Divide both by their trailing zeros (make both odd).  */
  un = arbint_rshift_limbs_inplace(up, un, ctz_u);
  vn = arbint_rshift_limbs_inplace(vp, vn, ctz_v);

  /*  Main loop: both u and v are odd.  */
  while (vn != 0u) {
    int cmp;
    size_t shift;

    /*  Ensure u >= v (swap if needed).  */
    cmp = arbint_cmp_mag_limbs(up, un, vp, vn);
    if (cmp == 0)
      break;

    if (cmp < 0) {
      tmp = up;
      up = vp;
      vp = tmp;
      {
        size_t t = un;
        un = vn;
        vn = t;
      }
    }

    /*  u = u - v (u > v, both odd, so result is even and nonzero).  */
    un = arbint__sub_mag(up, up, un, vp, vn);

    /*  Strip trailing zeros from u (guaranteed at least 1).  */
    if (un > 0u) {
      shift = arbint_gcd_mag_ctz(up, un);
      un = arbint_rshift_limbs_inplace(up, un, shift);
    }
  }

  /*  Result = u << common.  */
  if (common > 0u)
    un = arbint_lshift_limbs_inplace(up, un, common, cap);

  /*  Copy result to g.  */
  rc = arbint_resize(g, un);
  if (rc != ARBINT_OK)
    goto cleanup;

  if (un > 0u)
    memcpy(ARBINT_LIMBS(g), up, un * sizeof(arbint_limb_t));

  g[0]._sz = (ptrdiff_t) un;

cleanup:
  arbint_free_limbs(alloc, scratch);
  return rc;
}

/*  Euclidean GCD for small operands (uses tdiv_r).  */
static arbint_err_t arbint_gcd_euclid(arbint_t g, const arbint_t a,
                                      const arbint_t b) {
  arbint_t u, v, r;
  arbint_err_t rc;

  rc = arbint_init_all(g[0]._ctx, u, v, r, (arbint_t *) NULL);
  if (rc != ARBINT_OK)
    return rc;

  /*  Set u = |a|, v = |b|.  */
  rc = arbint_abs(u, a);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_abs(v, b);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Euclidean algorithm: gcd(u, v) = gcd(v, u mod v).  */
  while (!arbint_is_zero(v)) {
    rc = arbint_tdiv_r(r, u, v);
    if (rc != ARBINT_OK)
      goto cleanup;
    arbint_swap(u, v);
    arbint_swap(v, r);
  }

  /*  Result is in u.  */
  rc = arbint_set(g, u);

cleanup:
  arbint_clear_all(u, v, r, (arbint_t *) NULL);
  return rc;
}

/*  arbint_gcd: compute greatest common divisor.

    Algorithm selection based on operand size (min_n = min limb count):
    - min_n == 1: single-limb fast path using arbint_gcd_limb
    - min_n <= 4: Euclidean algorithm with tdiv_r (simple, low overhead)
    - min_n < 32: Binary GCD / Stein's algorithm (avoids division)
    - min_n >= 32: Lehmer GCD (amortizes matrix overhead with O(n*M(n)))

    Handles aliasing: if g aliases a or b, makes a copy before computing.
    Result is always non-negative: gcd(-12, 8) = gcd(12, -8) = 4.  */
ARBINT_API arbint_err_t arbint_gcd(arbint_t g, const arbint_t a,
                                   const arbint_t b) {
  size_t an;
  size_t bn;
  size_t min_n;
  const arbint_limb_t * ap;
  const arbint_limb_t * bp;
  const arbint_alloc_t * alloc;
  arbint_limb_t * a_copy = NULL;
  arbint_limb_t * b_copy = NULL;
  arbint_err_t rc;

  if (g == NULL || a == NULL || b == NULL)
    return ARBINT_EINVAL;

  an = arbint_abs_sz(a[0]._sz);
  bn = arbint_abs_sz(b[0]._sz);

  /*  Edge cases: gcd(0, x) = |x|, gcd(x, 0) = |x|.  */
  if (an == 0u)
    return arbint_abs(g, b);
  if (bn == 0u)
    return arbint_abs(g, a);

  ap = ARBINT_CLIMBS(a);
  bp = ARBINT_CLIMBS(b);

  /*  Single-limb fast path.  */
  if (an == 1u && bn == 1u) {
    arbint_limb_t result = arbint_gcd_limb(ap[0], bp[0]);
    rc = arbint_resize(g, 1u);
    if (rc != ARBINT_OK)
      return rc;
    if (result == 0u) {
      g[0]._sz = 0;
    } else {
      ARBINT_LIMBS(g)[0] = result;
      g[0]._sz = 1;
    }
    return ARBINT_OK;
  }

  min_n = (an < bn) ? an : bn;

  /*  Small operand path: use Euclidean with tdiv_r.  */
  if (min_n <= ARBINT_GCD_EUCLID_THRESHOLD)
    return arbint_gcd_euclid(g, a, b);

  /*  Get allocator for large operand paths.  */
  alloc = arbint_pick_alloc3(g, a, b);
  if (alloc == NULL)
    return ARBINT_EINVAL;

  /*  Handle aliasing: copy operands if they alias g.  */
  if (g == a) {
    a_copy = arbint_alloc_limbs(alloc, an);
    if (a_copy == NULL)
      return ARBINT_ENOMEM;
    memcpy(a_copy, ap, an * sizeof(arbint_limb_t));
    ap = a_copy;
  }
  if (g == b && g != a) {
    b_copy = arbint_alloc_limbs(alloc, bn);
    if (b_copy == NULL) {
      arbint_free_limbs(alloc, a_copy);
      return ARBINT_ENOMEM;
    }
    memcpy(b_copy, bp, bn * sizeof(arbint_limb_t));
    bp = b_copy;
  }

  /*  Very large operand path: use Lehmer GCD for O(n * M(n)) complexity.  */
  if (min_n >= ARBINT_LEHMER_THRESHOLD)
    rc = arbint_gcd_lehmer_dispatch(g, ap, an, bp, bn, alloc);
  else
    /*  Medium operand path: binary GCD with internal limb ops.  */
    rc = arbint_gcd_binary(g, ap, an, bp, bn, alloc);

  arbint_free_limbs(alloc, b_copy);
  arbint_free_limbs(alloc, a_copy);

  return rc;
}

/*  arbint_gcd_u32: GCD with a uint32_t operand.
    Allocation-free: uses limb-level division to get remainder directly.  */
ARBINT_API arbint_err_t arbint_gcd_u32(arbint_t g, const arbint_t a,
                                       uint32_t b) {
  size_t an;
  arbint_limb_t rem_limb;
  arbint_err_t rc;

  if (g == NULL || a == NULL)
    return ARBINT_EINVAL;

  an = arbint_abs_sz(a[0]._sz);

  /*  Edge cases.  */
  if (b == 0u)
    return arbint_abs(g, a);
  if (an == 0u)
    return arbint_set_u32(g, b);
  if (b == 1u)
    return arbint_set_u32(g, 1u);

  /*  Compute |a| mod b directly using dispatched limb-level division.  */
  rc = arbint_mod_mag_single_limb(ARBINT_CLIMBS(a), an, (arbint_limb_t) b,
                                  &rem_limb);
  if (rc != ARBINT_OK)
    return rc;

  /*  Now compute gcd(rem, b) with hardware division.  */
  return arbint_set_u32(g, arbint_gcd_u32u32((uint32_t) rem_limb, b));
}

/*  arbint_lcm: compute least common multiple.
    lcm(a, b) = (|a| / gcd(a, b)) * |b|
    Result is always non-negative.

    Fast shortcuts:
    - gcd == 1: lcm = |a| * |b| (skip division)
    - one divides other: lcm = max(|a|, |b|) (skip multiplication)  */
ARBINT_API arbint_err_t arbint_lcm(arbint_t l, const arbint_t a,
                                   const arbint_t b) {
  arbint_t gcd_val, quotient;
  arbint_err_t rc;
  arbint_ctx_t * ctx;
  size_t an;
  size_t bn;
  size_t gn;

  if (l == NULL || a == NULL || b == NULL)
    return ARBINT_EINVAL;

  an = arbint_abs_sz(a[0]._sz);
  bn = arbint_abs_sz(b[0]._sz);

  /*  lcm(0, x) = lcm(x, 0) = 0.  */
  if (an == 0u || bn == 0u) {
    arbint_zero(l);
    return ARBINT_OK;
  }

  ctx = l[0]._ctx;

  rc = arbint_init_all(ctx, gcd_val, quotient, (arbint_t *) NULL);
  if (rc != ARBINT_OK)
    return rc;

  /*  Compute gcd(a, b).  */
  rc = arbint_gcd(gcd_val, a, b);
  if (rc != ARBINT_OK)
    goto cleanup;

  gn = arbint_abs_sz(gcd_val[0]._sz);

  /*  Fast path: gcd == 1 means lcm = |a| * |b|.  */
  if (gn == 1u && ARBINT_CLIMBS(gcd_val)[0] == 1u) {
    rc = arbint_abs(l, a);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_abs(quotient, b);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_mul(l, l, quotient);
    goto cleanup;
  }

  /*  Fast path: gcd == |a| means a divides b, lcm = |b|.  */
  if (gn == an && arbint_cmpabs(gcd_val, a) == 0) {
    rc = arbint_abs(l, b);
    goto cleanup;
  }

  /*  Fast path: gcd == |b| means b divides a, lcm = |a|.  */
  if (gn == bn && arbint_cmpabs(gcd_val, b) == 0) {
    rc = arbint_abs(l, a);
    goto cleanup;
  }

  /*  quotient = |a| / gcd (exact division).  */
  rc = arbint_abs(quotient, a);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_tdiv_q(quotient, quotient, gcd_val);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  l = quotient * |b|.  */
  rc = arbint_abs(l, b);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_mul(l, quotient, l);

cleanup:
  arbint_clear_all(gcd_val, quotient, (arbint_t *) NULL);
  return rc;
}

/*  arbint_lcm_u32: LCM with a uint32_t operand.
    Optimized: uses allocation-free mod for gcd computation.  */
ARBINT_API arbint_err_t arbint_lcm_u32(arbint_t l, const arbint_t a,
                                       uint32_t b) {
  arbint_t tmp;
  arbint_err_t rc;
  arbint_limb_t rem_limb;
  uint32_t g;
  size_t an;

  if (l == NULL || a == NULL)
    return ARBINT_EINVAL;

  an = arbint_abs_sz(a[0]._sz);

  /*  lcm(0, b) = lcm(a, 0) = 0.  */
  if (an == 0u || b == 0u) {
    arbint_zero(l);
    return ARBINT_OK;
  }

  /*  Fast path: lcm(a, 1) = |a|.  */
  if (b == 1u)
    return arbint_abs(l, a);

  /*  Compute |a| mod b directly using dispatched limb-level division.  */
  rc = arbint_mod_mag_single_limb(ARBINT_CLIMBS(a), an, (arbint_limb_t) b,
                                  &rem_limb);
  if (rc != ARBINT_OK)
    return rc;

  /*  Compute gcd(rem, b).  */
  if (rem_limb == 0u) {
    /*  b divides |a|, so lcm = |a|.  */
    return arbint_abs(l, a);
  }
  g = arbint_gcd_u32u32((uint32_t) rem_limb, b);

  /*  Fast path: if gcd == 1, lcm = |a| * b.  */
  if (g == 1u) {
    rc = arbint_abs(l, a);
    if (rc != ARBINT_OK)
      return rc;
    return arbint_mul_u32(l, l, b);
  }

  /*  l = (|a| / g) * b.  */
  rc = arbint_init(tmp, l[0]._ctx);
  if (rc != ARBINT_OK)
    return rc;

  rc = arbint_abs(tmp, a);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_tdiv_q_u32(tmp, tmp, g);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_mul_u32(l, tmp, b);

cleanup:
  arbint_clear(tmp);
  return rc;
}

/*  Multiply signed arbint by a single limb (non-negative multiplier).  */
static arbint_err_t arbint_xgcd_mul_limb(arbint_t rop, const arbint_t a,
                                         arbint_limb_t b) {
  const arbint_mul_kernel_table_t * k;
  int as;
  size_t an;
  size_t cap;
  size_t used;
  arbint_err_t rc;

  as = (a[0]._sz > 0) - (a[0]._sz < 0);
  if (as == 0 || b == 0u) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  if (b == 1u)
    return arbint_set(rop, a);

  an = arbint_abs_sz(a[0]._sz);
  cap = an + 1u;
  if (cap < an)
    return ARBINT_EOVERFLOW;

  rc = arbint_resize(rop, cap);
  if (rc != ARBINT_OK)
    return rc;

  k = arbint_mul_kernel_table_get();
  used = k->mul_limb_1(ARBINT_LIMBS(rop), ARBINT_CLIMBS(a), an, b);

  if (!arbint_set_signed_sz(rop, used, as))
    return ARBINT_EOVERFLOW;
  return ARBINT_OK;
}

/*  out = lhs_mul*lhs - rhs_mul*rhs, with non-negative limb multipliers.  */
static arbint_err_t arbint_xgcd_lincomb_sub(arbint_t out, const arbint_t lhs,
                                            arbint_limb_t lhs_mul,
                                            const arbint_t rhs,
                                            arbint_limb_t rhs_mul,
                                            arbint_t lhs_term,
                                            arbint_t rhs_term) {
  arbint_err_t rc;

  rc = arbint_xgcd_mul_limb(lhs_term, lhs, lhs_mul);
  if (rc != ARBINT_OK)
    return rc;

  rc = arbint_xgcd_mul_limb(rhs_term, rhs, rhs_mul);
  if (rc != ARBINT_OK)
    return rc;

  return arbint_sub(out, lhs_term, rhs_term);
}

/*  Apply one Lehmer block transform to pair (x0, x1).  */
static arbint_err_t arbint_xgcd_apply_lehmer_pair(
    arbint_t x0, arbint_t x1, const arbint_lehmer_matrix_t * m, int even,
    arbint_t out0, arbint_t out1, arbint_t tmp0, arbint_t tmp1) {
  arbint_err_t rc;

  if (even) {
    /*  (x0, x1) <- (v0*x0 - v1*x1, w1*x1 - w0*x0).  */
    rc = arbint_xgcd_lincomb_sub(out0, x0, m->v0, x1, m->v1, tmp0, tmp1);
    if (rc != ARBINT_OK)
      return rc;
    return arbint_xgcd_lincomb_sub(out1, x1, m->w1, x0, m->w0, tmp0, tmp1);
  }

  /*  (x0, x1) <- (v1*x1 - v0*x0, w0*x0 - w1*x1).  */
  rc = arbint_xgcd_lincomb_sub(out0, x1, m->v1, x0, m->v0, tmp0, tmp1);
  if (rc != ARBINT_OK)
    return rc;
  return arbint_xgcd_lincomb_sub(out1, x0, m->w0, x1, m->w1, tmp0, tmp1);
}

/*  Extended GCD: compute g = gcd(a, b) and Bezout coefficients x, y
    such that a*x + b*y = g.

    Uses extended Euclid with Lehmer block acceleration on large inputs:
    - Track (r0, r1) = remainders, starting with (|a|, |b|)
    - Track (s0, s1) = coefficients for first input
    - Track (t0, t1) = coefficients for second input
    - Large close-sized operands: apply one Lehmer matrix block
    - Otherwise: one exact Euclidean step q = floor(r0/r1)

    Sign handling: we compute with |a|, |b|, then adjust signs at the end.
    If a < 0: x = -s0. If b < 0: y = -t0.

    Aliasing: per arbint.h lines 213-218, outputs are written in order
    g, x, y. "The final value is whichever output is written last."
    So g==x gets x's value, g==y gets y's value, x==y gets y's value.

    Performance note: large inputs use the same Lehmer step simulation
    as arbint_gcd to reduce multi-limb divisions.  */
ARBINT_API arbint_err_t arbint_xgcd(arbint_t g, arbint_t x, arbint_t y,
                                    const arbint_t a, const arbint_t b) {
  arbint_t r0, r1, s0, s1, t0, t1, q, tmp;
  arbint_t pair0, pair1, prod0, prod1;
  arbint_t res_g, res_x, res_y;
  arbint_err_t rc;
  arbint_ctx_t * ctx;
  int a_neg, b_neg;
  int init_work = 0;

  if (g == NULL || x == NULL || y == NULL || a == NULL || b == NULL)
    return ARBINT_EINVAL;

  a_neg = (a[0]._sz < 0);
  b_neg = (b[0]._sz < 0);

  /*  Find context from any input.  */
  ctx = g[0]._ctx;
  if (ctx == NULL)
    ctx = x[0]._ctx;
  if (ctx == NULL)
    ctx = y[0]._ctx;
  if (ctx == NULL)
    ctx = a[0]._ctx;
  if (ctx == NULL)
    ctx = b[0]._ctx;

  /*  Handle gcd(0, 0) = 0 specially.  */
  if (a[0]._sz == 0 && b[0]._sz == 0) {
    arbint_zero(g);
    arbint_zero(x);
    arbint_zero(y);
    return ARBINT_OK;
  }

  /*  Allocate result temporaries for aliasing safety.  */
  rc = arbint_init_all(ctx, res_g, res_x, res_y, (arbint_t *) NULL);
  if (rc != ARBINT_OK)
    return rc;

  /*  Handle gcd(0, b) = |b| with x=0, y=sign(b).  */
  if (a[0]._sz == 0) {
    rc = arbint_abs(res_g, b);
    if (rc != ARBINT_OK)
      goto cleanup;
    arbint_zero(res_x);
    rc = arbint_set_i32(res_y, b_neg ? -1 : 1);
    if (rc != ARBINT_OK)
      goto cleanup;
    goto copy_results;
  }

  /*  Handle gcd(a, 0) = |a| with x=sign(a), y=0.  */
  if (b[0]._sz == 0) {
    rc = arbint_abs(res_g, a);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_set_i32(res_x, a_neg ? -1 : 1);
    if (rc != ARBINT_OK)
      goto cleanup;
    arbint_zero(res_y);
    goto copy_results;
  }

  /*  Initialize working temporaries.  */
  rc = arbint_init_all(ctx, r0, r1, s0, s1, t0, t1, q, tmp, pair0, pair1,
                       prod0, prod1, (arbint_t *) NULL);
  if (rc != ARBINT_OK)
    goto cleanup;
  init_work = 1;

  /*  Set initial values: r0 = |a|, r1 = |b|, s0 = 1, s1 = 0, t0 = 0,
      t1 = 1.  */
  rc = arbint_abs(r0, a);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_abs(r1, b);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_set_i32(s0, 1);
  if (rc != ARBINT_OK)
    goto cleanup;
  arbint_zero(s1);
  arbint_zero(t0);
  rc = arbint_set_i32(t1, 1);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Keep r0 >= r1 for Lehmer simulation precondition.  */
  if (arbint_cmpabs(r0, r1) < 0) {
    arbint_swap(r0, r1);
    arbint_swap(s0, s1);
    arbint_swap(t0, t1);
  }

  /*  Extended Euclidean loop.  */
  while (!arbint_is_zero(r1)) {
    size_t r0_n = arbint_abs_sz(r0[0]._sz);
    size_t r1_n = arbint_abs_sz(r1[0]._sz);
    unsigned count = 0u;
    int even = 1;
    arbint_lehmer_matrix_t m;

    /*  Lehmer simulation is reliable only when top lengths are close.  */
    if (r1_n >= ARBINT_LEHMER_THRESHOLD && r0_n <= r1_n + 1u) {
      const arbint_limb_t * r0p = ARBINT_CLIMBS(r0);
      const arbint_limb_t * r1p = ARBINT_CLIMBS(r1);
      arbint_limb_t a1 = r0p[r0_n - 1u];
      arbint_limb_t a0 = (r0_n >= 2u) ? r0p[r0_n - 2u] : 0u;
      arbint_limb_t b1 = r1p[r1_n - 1u];
      arbint_limb_t b0 = (r1_n >= 2u) ? r1p[r1_n - 2u] : 0u;

      count = arbint_lehmer_step(&m, a1, a0, b1, b0, &even);
    }

    if (count != 0u) {
      /*  Apply same Lehmer block to remainders and both coefficient pairs.  */
      rc = arbint_xgcd_apply_lehmer_pair(r0, r1, &m, even, pair0, pair1, prod0,
                                         prod1);
      if (rc != ARBINT_OK)
        goto cleanup;
      arbint_swap(r0, pair0);
      arbint_swap(r1, pair1);

      rc = arbint_xgcd_apply_lehmer_pair(s0, s1, &m, even, pair0, pair1, prod0,
                                         prod1);
      if (rc != ARBINT_OK)
        goto cleanup;
      arbint_swap(s0, pair0);
      arbint_swap(s1, pair1);

      rc = arbint_xgcd_apply_lehmer_pair(t0, t1, &m, even, pair0, pair1, prod0,
                                         prod1);
      if (rc != ARBINT_OK)
        goto cleanup;
      arbint_swap(t0, pair0);
      arbint_swap(t1, pair1);
      continue;
    }

    /*  q = floor(r0 / r1), r0 = r0 mod r1.  */
    rc = arbint_tdiv_qr(q, r0, r0, r1);
    if (rc != ARBINT_OK)
      goto cleanup;
    arbint_swap(r0, r1);

    /*  Update s: tmp = s0 - q*s1, (s0, s1) = (s1, tmp).  */
    rc = arbint_mul(tmp, q, s1);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_sub(tmp, s0, tmp);
    if (rc != ARBINT_OK)
      goto cleanup;
    arbint_swap(s0, s1);
    arbint_swap(s1, tmp);

    /*  Update t: tmp = t0 - q*t1, (t0, t1) = (t1, tmp).  */
    rc = arbint_mul(tmp, q, t1);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_sub(tmp, t0, tmp);
    if (rc != ARBINT_OK)
      goto cleanup;
    arbint_swap(t0, t1);
    arbint_swap(t1, tmp);
  }

  /*  Copy results with sign adjustment.
      We computed with |a|, |b|, so s0*|a| + t0*|b| = g.
      For signed a, b: a*x + b*y = g where x = s0*(a<0?-1:1), y =
      t0*(b<0?-1:1).  */
  rc = arbint_set(res_g, r0);
  if (rc != ARBINT_OK)
    goto cleanup;

  if (a_neg)
    rc = arbint_neg(res_x, s0);
  else
    rc = arbint_set(res_x, s0);
  if (rc != ARBINT_OK)
    goto cleanup;

  if (b_neg)
    rc = arbint_neg(res_y, t0);
  else
    rc = arbint_set(res_y, t0);
  if (rc != ARBINT_OK)
    goto cleanup;

copy_results:
  /*  Write outputs in order: g, x, y.
      Per aliasing contract, last write wins for overlapping outputs.  */
  rc = arbint_set(g, res_g);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_set(x, res_x);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_set(y, res_y);
  /*  Fall through to cleanup.  */

cleanup:
  if (init_work)
    arbint_clear_all(r0, r1, s0, s1, t0, t1, q, tmp, pair0, pair1, prod0,
                     prod1, (arbint_t *) NULL);
  arbint_clear_all(res_g, res_x, res_y, (arbint_t *) NULL);
  return rc;
}
