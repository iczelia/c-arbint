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
#include "arbint_shift.h"
#include "config.h"

#include <string.h>

/*  Lehmer GCD dispatch function pointer type.

    Lehmer GCD has platform-specific implementations that use different
    multiply primitives. The dispatch selects the optimal implementation
    based on CPU features detected at runtime.  */
typedef arbint_err_t (*arbint_gcd_lehmer_fn_t)(arbint_t g,
                                               const arbint_limb_t * ap,
                                               size_t an,
                                               const arbint_limb_t * bp,
                                               size_t bn,
                                               const arbint_alloc_t * alloc);

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
static arbint_err_t arbint_gcd_lehmer_dispatch(arbint_t g,
                                               const arbint_limb_t * ap,
                                               size_t an,
                                               const arbint_limb_t * bp,
                                               size_t bn,
                                               const arbint_alloc_t * alloc) {
  if (g_gcd_lehmer == NULL)
    g_gcd_lehmer = arbint_select_gcd_lehmer();
  return g_gcd_lehmer(g, ap, an, bp, bn, alloc);
}

/*  Allocator selection helpers.  */

static const arbint_alloc_t * arbint_gcd_get_alloc_from(const arbint_t x) {
  if (x == NULL || x[0]._ctx == NULL || x[0]._ctx->a.realloc == NULL)
    return NULL;
  return &x[0]._ctx->a;
}

static const arbint_alloc_t * arbint_gcd_pick_alloc(const arbint_t a,
                                                    const arbint_t b,
                                                    const arbint_t c) {
  const arbint_alloc_t * alloc;

  alloc = arbint_gcd_get_alloc_from(a);
  if (alloc != NULL)
    return alloc;
  alloc = arbint_gcd_get_alloc_from(b);
  if (alloc != NULL)
    return alloc;
  return arbint_gcd_get_alloc_from(c);
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
  arbint_ctx_t * ctx;
  int init_u = 0;
  int init_v = 0;
  int init_r = 0;

  ctx = g[0]._ctx;

  rc = arbint_init(u, ctx);
  if (rc != ARBINT_OK)
    return rc;
  init_u = 1;

  rc = arbint_init(v, ctx);
  if (rc != ARBINT_OK)
    goto cleanup;
  init_v = 1;

  rc = arbint_init(r, ctx);
  if (rc != ARBINT_OK)
    goto cleanup;
  init_r = 1;

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
  if (init_r)
    arbint_clear(r);
  if (init_v)
    arbint_clear(v);
  if (init_u)
    arbint_clear(u);
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
  alloc = arbint_gcd_pick_alloc(g, a, b);
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
  int init_gcd = 0;
  int init_quot = 0;
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

  rc = arbint_init(gcd_val, ctx);
  if (rc != ARBINT_OK)
    return rc;
  init_gcd = 1;

  /*  Compute gcd(a, b).  */
  rc = arbint_gcd(gcd_val, a, b);
  if (rc != ARBINT_OK)
    goto cleanup;

  gn = arbint_abs_sz(gcd_val[0]._sz);

  /*  Fast path: gcd == 1 means lcm = |a| * |b|.  */
  if (gn == 1u && ARBINT_CLIMBS(gcd_val)[0] == 1u) {
    arbint_clear(gcd_val);
    init_gcd = 0;

    rc = arbint_abs(l, a);
    if (rc != ARBINT_OK)
      return rc;
    rc = arbint_init(quotient, ctx);
    if (rc != ARBINT_OK)
      return rc;
    rc = arbint_abs(quotient, b);
    if (rc != ARBINT_OK) {
      arbint_clear(quotient);
      return rc;
    }
    rc = arbint_mul(l, l, quotient);
    arbint_clear(quotient);
    return rc;
  }

  /*  Fast path: gcd == |a| means a divides b, lcm = |b|.  */
  if (gn == an && arbint_cmpabs(gcd_val, a) == 0) {
    arbint_clear(gcd_val);
    return arbint_abs(l, b);
  }

  /*  Fast path: gcd == |b| means b divides a, lcm = |a|.  */
  if (gn == bn && arbint_cmpabs(gcd_val, b) == 0) {
    arbint_clear(gcd_val);
    return arbint_abs(l, a);
  }

  rc = arbint_init(quotient, ctx);
  if (rc != ARBINT_OK)
    goto cleanup;
  init_quot = 1;

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
  if (init_quot)
    arbint_clear(quotient);
  if (init_gcd)
    arbint_clear(gcd_val);
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
  int init_tmp = 0;

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
  init_tmp = 1;

  rc = arbint_abs(tmp, a);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_tdiv_q_u32(tmp, tmp, g);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_mul_u32(l, tmp, b);

cleanup:
  if (init_tmp)
    arbint_clear(tmp);
  return rc;
}
