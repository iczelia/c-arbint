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

/*  Lehmer GCD - BMI2-optimized implementation.

    This is the x86-64 optimized version of Lehmer GCD, selected at runtime
    when BMI2 instructions are available. The algorithm is identical to the
    generic implementation; only the low-level multiply primitives differ.

    BMI2 optimization:
    -----------------
    Uses _mulx_u64 intrinsic for 64x64->128 multiplication:
      - Latency: 1 cycle (vs ~5 for half-limb decomposition)
      - Throughput: 1/cycle on modern Intel/AMD
      - No flag modification (unlike IMUL)

    Performance impact:
    - Matrix application (arbint_lehmer_mul_1_bmi2): ~2.5x faster
    - Overall Lehmer GCD: ~1.5-2x faster on large operands

    The algorithm overview, fallback strategy, and memory layout are
    identical to arbint_gcd_lehmer_generic.c. See that file for detailed
    documentation of the Lehmer algorithm itself.

    Dispatch:
    --------
    Selected by arbint_gcd.c via three-tier dispatch:
    1. HAS_BMI2_ALWAYS: compiled with -mbmi2, always use BMI2
    2. HAS_BMI2: runtime check via CPUID, dispatch accordingly
    3. Otherwise: use generic implementation  */

#include "arbint_gcd.h"

#include "arbint_addsub.h"
#include "arbint_base.h"
#include "arbint_div.h"
#include "arbint_shift.h"
#include "config.h"

#include <string.h>

#if HAS_BMI2 && ARBINT_LIMB_BITS == 64
#include <immintrin.h>
#endif

/*  BMI2-optimized wide multiply: hi:lo = a * b.

    Uses _mulx_u64 intrinsic which computes a 128-bit product in a single
    cycle without modifying flags. This file is only compiled when BMI2
    is available, so no fallback paths are needed.  */
static inline void arbint_lehmer_umul_bmi2(arbint_limb_t * hi,
                                           arbint_limb_t * lo,
                                           arbint_limb_t a,
                                           arbint_limb_t b) {
  *lo = _mulx_u64(a, b, (unsigned long long *) hi);
}

/*  Multiply limb array by single limb: dst = src * mult.
    Returns result size (may be n or n+1).  */
static size_t arbint_lehmer_mul_1_bmi2(arbint_limb_t * dst,
                                       const arbint_limb_t * src, size_t n,
                                       arbint_limb_t mult) {
  size_t i;
  arbint_limb_t carry = 0u;
  arbint_limb_t hi, lo;

  for (i = 0u; i < n; ++i) {
    arbint_lehmer_umul_bmi2(&hi, &lo, src[i], mult);
    lo += carry;
    if (lo < carry)
      ++hi;
    dst[i] = lo;
    carry = hi;
  }

  if (carry != 0u) {
    dst[n] = carry;
    return n + 1u;
  }
  return n;
}

/*  Subtract limb arrays: dst = a - b (assumes a >= b).
    Returns result size.  */
static size_t arbint_lehmer_sub_bmi2(arbint_limb_t * dst,
                                     const arbint_limb_t * a, size_t an,
                                     const arbint_limb_t * b, size_t bn) {
  return arbint__sub_mag(dst, a, an, b, bn);
}

/*  Apply Lehmer matrix to operands.
    scratch layout: [v0*u | v1*v | w0*u | w1*v] each cap/4 limbs.
    scratch must have capacity for 4 * (max_n + 2) limbs.  */
static void arbint_lehmer_apply_matrix_bmi2(arbint_limb_t * up, size_t * un,
                                            arbint_limb_t * vp, size_t * vn,
                                            const arbint_lehmer_matrix_t * m,
                                            int even,
                                            arbint_limb_t * scratch,
                                            size_t scratch_cap) {
  size_t slot = scratch_cap / 4u;
  arbint_limb_t * v0u = scratch;
  arbint_limb_t * v1v = scratch + slot;
  arbint_limb_t * w0u = scratch + slot * 2u;
  arbint_limb_t * w1v = scratch + slot * 3u;
  size_t v0u_n, v1v_n, w0u_n, w1v_n;
  size_t new_un, new_vn;

  /*  Compute all four products using original u and v.  */
  if (m->v0 == 0u) {
    v0u_n = 0u;
  } else if (m->v0 == 1u) {
    memcpy(v0u, up, *un * sizeof(arbint_limb_t));
    v0u_n = *un;
  } else {
    v0u_n = arbint_lehmer_mul_1_bmi2(v0u, up, *un, m->v0);
  }

  if (m->v1 == 0u) {
    v1v_n = 0u;
  } else if (m->v1 == 1u) {
    memcpy(v1v, vp, *vn * sizeof(arbint_limb_t));
    v1v_n = *vn;
  } else {
    v1v_n = arbint_lehmer_mul_1_bmi2(v1v, vp, *vn, m->v1);
  }

  if (m->w0 == 0u) {
    w0u_n = 0u;
  } else if (m->w0 == 1u) {
    memcpy(w0u, up, *un * sizeof(arbint_limb_t));
    w0u_n = *un;
  } else {
    w0u_n = arbint_lehmer_mul_1_bmi2(w0u, up, *un, m->w0);
  }

  if (m->w1 == 0u) {
    w1v_n = 0u;
  } else if (m->w1 == 1u) {
    memcpy(w1v, vp, *vn * sizeof(arbint_limb_t));
    w1v_n = *vn;
  } else {
    w1v_n = arbint_lehmer_mul_1_bmi2(w1v, vp, *vn, m->w1);
  }

  /*  Compute new_u = v0*u - v1*v (or v1*v - v0*u if odd).  */
  if (even) {
    if (v0u_n > v1v_n ||
        (v0u_n == v1v_n && arbint_cmp_mag_limbs(v0u, v0u_n, v1v, v1v_n) >= 0)) {
      new_un = arbint_lehmer_sub_bmi2(up, v0u, v0u_n, v1v, v1v_n);
    } else {
      new_un = arbint_lehmer_sub_bmi2(up, v1v, v1v_n, v0u, v0u_n);
    }
  } else {
    if (v1v_n > v0u_n ||
        (v1v_n == v0u_n && arbint_cmp_mag_limbs(v1v, v1v_n, v0u, v0u_n) >= 0)) {
      new_un = arbint_lehmer_sub_bmi2(up, v1v, v1v_n, v0u, v0u_n);
    } else {
      new_un = arbint_lehmer_sub_bmi2(up, v0u, v0u_n, v1v, v1v_n);
    }
  }

  /*  Compute new_v = w1*v - w0*u (or w0*u - w1*v if odd).  */
  if (even) {
    if (w1v_n > w0u_n ||
        (w1v_n == w0u_n && arbint_cmp_mag_limbs(w1v, w1v_n, w0u, w0u_n) >= 0)) {
      new_vn = arbint_lehmer_sub_bmi2(vp, w1v, w1v_n, w0u, w0u_n);
    } else {
      new_vn = arbint_lehmer_sub_bmi2(vp, w0u, w0u_n, w1v, w1v_n);
    }
  } else {
    if (w0u_n > w1v_n ||
        (w0u_n == w1v_n && arbint_cmp_mag_limbs(w0u, w0u_n, w1v, w1v_n) >= 0)) {
      new_vn = arbint_lehmer_sub_bmi2(vp, w0u, w0u_n, w1v, w1v_n);
    } else {
      new_vn = arbint_lehmer_sub_bmi2(vp, w1v, w1v_n, w0u, w0u_n);
    }
  }

  *un = new_un;
  *vn = new_vn;
}

/*  Binary GCD (Stein's algorithm) - fallback for small operands.  */
static arbint_err_t arbint_gcd_binary_fallback_bmi2(arbint_t g,
                                                    arbint_limb_t * up,
                                                    size_t un,
                                                    arbint_limb_t * vp,
                                                    size_t vn,
                                                    const arbint_alloc_t * alloc) {
  arbint_limb_t * tmp;
  size_t ctz_u, ctz_v, common;
  size_t cap;
  arbint_err_t rc;

  (void) alloc;

  /*  Handle zero cases: gcd(u, 0) = u, gcd(0, v) = v.  */
  if (vn == 0u) {
    rc = arbint_resize(g, un);
    if (rc != ARBINT_OK)
      return rc;
    if (un > 0u)
      memcpy(ARBINT_LIMBS(g), up, un * sizeof(arbint_limb_t));
    g[0]._sz = (ptrdiff_t) un;
    return ARBINT_OK;
  }
  if (un == 0u) {
    rc = arbint_resize(g, vn);
    if (rc != ARBINT_OK)
      return rc;
    if (vn > 0u)
      memcpy(ARBINT_LIMBS(g), vp, vn * sizeof(arbint_limb_t));
    g[0]._sz = (ptrdiff_t) vn;
    return ARBINT_OK;
  }

  ctz_u = arbint_gcd_mag_ctz(up, un);
  ctz_v = arbint_gcd_mag_ctz(vp, vn);
  common = (ctz_u < ctz_v) ? ctz_u : ctz_v;

  un = arbint_rshift_limbs_inplace(up, un, ctz_u);
  vn = arbint_rshift_limbs_inplace(vp, vn, ctz_v);

  while (vn != 0u) {
    int cmp;
    size_t shift;

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

    un = arbint__sub_mag(up, up, un, vp, vn);

    if (un > 0u) {
      shift = arbint_gcd_mag_ctz(up, un);
      un = arbint_rshift_limbs_inplace(up, un, shift);
    }
  }

  cap = un + (common + ARBINT_LIMB_BITS - 1u) / ARBINT_LIMB_BITS + 1u;
  if (common > 0u)
    un = arbint_lshift_limbs_inplace(up, un, common, cap);

  rc = arbint_resize(g, un);
  if (rc != ARBINT_OK)
    return rc;

  if (un > 0u)
    memcpy(ARBINT_LIMBS(g), up, un * sizeof(arbint_limb_t));

  g[0]._sz = (ptrdiff_t) un;
  return ARBINT_OK;
}

/*  Lehmer GCD - main entry point (BMI2 implementation).  */
arbint_err_t arbint_gcd_lehmer_bmi2(arbint_t g, const arbint_limb_t * ap,
                                    size_t an, const arbint_limb_t * bp,
                                    size_t bn, const arbint_alloc_t * alloc) {
  arbint_limb_t * scratch = NULL;
  arbint_limb_t * up;
  arbint_limb_t * vp;
  arbint_limb_t * work;
  size_t un, vn;
  size_t max_n;
  size_t cap;
  arbint_err_t rc = ARBINT_OK;

  max_n = (an > bn) ? an : bn;

  cap = max_n + 2u;
  scratch = arbint_alloc_limbs(alloc, cap * 6u);
  if (scratch == NULL)
    return ARBINT_ENOMEM;

  up = scratch;
  vp = scratch + cap;
  work = scratch + cap * 2u;

  memcpy(up, ap, an * sizeof(arbint_limb_t));
  memcpy(vp, bp, bn * sizeof(arbint_limb_t));
  un = an;
  vn = bn;

  if (arbint_cmp_mag_limbs(up, un, vp, vn) < 0) {
    arbint_limb_t * tmp = up;
    up = vp;
    vp = tmp;
    size_t t = un;
    un = vn;
    vn = t;
  }

  while (vn >= ARBINT_GCD_EUCLID_THRESHOLD && un > 0u && vn > 0u) {
    arbint_lehmer_matrix_t m;
    arbint_limb_t a1, a0, b1, b0;
    unsigned count;
    int even;

    /*  Lehmer simulation is only reliable when operand lengths are close.
        If u has at least two more limbs than v, the true quotient may be
        very large while top-limb simulation suggests a tiny q; force an
        exact division step in that case.  */
    if (un > vn + 1u) {
      count = 0u;
    } else {
      a1 = up[un - 1u];
      a0 = (un >= 2u) ? up[un - 2u] : 0u;
      b1 = vp[vn - 1u];
      b0 = (vn >= 2u) ? vp[vn - 2u] : 0u;

      count = arbint_lehmer_step(&m, a1, a0, b1, b0, &even);
    }

    if (count == 0u) {
      /*  Simulation unsafe; do one Euclidean step (u = u mod v).
          Use division rather than repeated subtraction to avoid O(q)
          performance degradation for large quotients.  */
      if (vn == 1u) {
        /*  Single-limb divisor: use fast single-limb remainder.  */
        arbint_limb_t rem;
        (void) arbint_mod_mag_single_limb(up, un, vp[0], &rem);
        if (rem == 0u) {
          un = 0u;
        } else {
          up[0] = rem;
          un = 1u;
        }
      } else if (vn >= 2u && un >= vn) {
        /*  Multi-limb divisor: use Knuth division to compute remainder.
            work area provides scratch space; quotient discarded.  */
        arbint_limb_t * qp_scratch = work;
        arbint_limb_t * rp_scratch = work + (un - vn + 1u);
        arbint_err_t div_rc;

        div_rc = arbint_div_mag_knuth(up, un, vp, vn, qp_scratch, rp_scratch);
        if (div_rc == ARBINT_OK) {
          /*  Copy remainder back to up; normalize.  */
          memcpy(up, rp_scratch, vn * sizeof(arbint_limb_t));
          un = arbint_norm_used(up, vn);
        } else {
          /*  Division failed; fall back to subtraction-based.  */
          while (arbint_cmp_mag_limbs(up, un, vp, vn) >= 0) {
            un = arbint__sub_mag(up, up, un, vp, vn);
            if (un == 0u)
              break;
          }
        }
      } else {
        /*  un < vn: no division needed, u is already reduced.  */
      }
      /*  Swap u and v.  */
      {
        arbint_limb_t * tmp = up;
        up = vp;
        vp = tmp;
        size_t t = un;
        un = vn;
        vn = t;
      }
    } else {
      arbint_lehmer_apply_matrix_bmi2(up, &un, vp, &vn, &m, even,
                                      work, cap * 4u);

      un = arbint_norm_used(up, un);
      vn = arbint_norm_used(vp, vn);

      if (arbint_cmp_mag_limbs(up, un, vp, vn) < 0) {
        arbint_limb_t * tmp = up;
        up = vp;
        vp = tmp;
        size_t t = un;
        un = vn;
        vn = t;
      }
    }
  }

  rc = arbint_gcd_binary_fallback_bmi2(g, up, un, vp, vn, alloc);

  arbint_free_limbs(alloc, scratch);
  return rc;
}
