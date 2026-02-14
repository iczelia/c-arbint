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

/*  Lehmer GCD - shared template (included by generic and BMI2 variants).

    This file contains the implementation of Lehmer GCD that is shared between
    the generic and BMI2-optimized versions. Before including this file, define:

      ARBINT_LEHMER_MUL_1_FN       - limb-by-scalar multiply function name
      ARBINT_LEHMER_APPLY_FN      - matrix application function name
      ARBINT_LEHMER_FALLBACK_FN   - binary GCD fallback function name
      ARBINT_LEHMER_ENTRY_FN      - main entry point function name
      ARBINT_LEHMER_UMUL_FN       - wide multiply (hi:lo = a*b) function

    Algorithm overview:
    ------------------
    Lehmer's algorithm accelerates Euclidean GCD for large integers by
    simulating multiple quotient steps using only the top 1-2 limbs. This
    reduces the number of expensive multi-limb divisions from O(n) per step
    to O(1) matrix application, with O(log n) Lehmer steps total.

    See arbint_gcd_lehmer_generic.c for full algorithm documentation.  */

#ifndef ARBINT_LEHMER_MUL_1_FN
#error "ARBINT_LEHMER_MUL_1_FN must be defined before including this file"
#endif

#ifndef ARBINT_LEHMER_APPLY_FN
#error "ARBINT_LEHMER_APPLY_FN must be defined before including this file"
#endif

#ifndef ARBINT_LEHMER_FALLBACK_FN
#error "ARBINT_LEHMER_FALLBACK_FN must be defined before including this file"
#endif

#ifndef ARBINT_LEHMER_ENTRY_FN
#error "ARBINT_LEHMER_ENTRY_FN must be defined before including this file"
#endif

#ifndef ARBINT_LEHMER_UMUL_FN
#error "ARBINT_LEHMER_UMUL_FN must be defined before including this file"
#endif

/*  Multiply limb array by single limb: dst = src * mult.

    Used to compute the four products v0*u, v1*v, w0*u, w1*v when applying
    the Lehmer matrix. Each matrix coefficient is a single limb (limited
    to LIMB_MAX/2 by the overflow check in arbint_lehmer_step).

    Returns result size (may be n or n+1 if carry propagates).  */
static size_t ARBINT_LEHMER_MUL_1_FN(arbint_limb_t * dst,
                                     const arbint_limb_t * src, size_t n,
                                     arbint_limb_t mult) {
  size_t i;
  arbint_limb_t carry = 0u;
  arbint_limb_t hi, lo;

  for (i = 0u; i < n; ++i) {
    ARBINT_LEHMER_UMUL_FN(&hi, &lo, src[i], mult);
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

/*  Apply Lehmer matrix to operands.

    Given the matrix M = [[v0, v1], [w0, w1]] built by arbint_lehmer_step,
    computes the transformed operands:

      If even (k quotients where k is even):
        u_new = v0*u - v1*v
        v_new = w1*v - w0*u

      If odd (k quotients where k is odd):
        u_new = v1*v - v0*u
        v_new = w0*u - w1*v

    Optimization for small coefficients:
    - If coefficient == 0: product is zero (skip computation)
    - If coefficient == 1: product equals operand (memcpy, no multiply)
    - Otherwise: full limb-by-single-limb multiply

    Scratch layout: [v0*u | v1*v | w0*u | w1*v], each slot = cap/4 limbs.  */
static void ARBINT_LEHMER_APPLY_FN(arbint_limb_t * up, size_t * un,
                                   arbint_limb_t * vp, size_t * vn,
                                   const arbint_lehmer_matrix_t * m, int even,
                                   arbint_limb_t * scratch, size_t scratch_cap) {
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
    v0u_n = ARBINT_LEHMER_MUL_1_FN(v0u, up, *un, m->v0);
  }

  if (m->v1 == 0u) {
    v1v_n = 0u;
  } else if (m->v1 == 1u) {
    memcpy(v1v, vp, *vn * sizeof(arbint_limb_t));
    v1v_n = *vn;
  } else {
    v1v_n = ARBINT_LEHMER_MUL_1_FN(v1v, vp, *vn, m->v1);
  }

  if (m->w0 == 0u) {
    w0u_n = 0u;
  } else if (m->w0 == 1u) {
    memcpy(w0u, up, *un * sizeof(arbint_limb_t));
    w0u_n = *un;
  } else {
    w0u_n = ARBINT_LEHMER_MUL_1_FN(w0u, up, *un, m->w0);
  }

  if (m->w1 == 0u) {
    w1v_n = 0u;
  } else if (m->w1 == 1u) {
    memcpy(w1v, vp, *vn * sizeof(arbint_limb_t));
    w1v_n = *vn;
  } else {
    w1v_n = ARBINT_LEHMER_MUL_1_FN(w1v, vp, *vn, m->w1);
  }

  /*  Compute new_u = v0*u - v1*v (or v1*v - v0*u if odd).  */
  if (even) {
    if (v0u_n > v1v_n ||
        (v0u_n == v1v_n && arbint_cmp_mag_limbs(v0u, v0u_n, v1v, v1v_n) >= 0)) {
      new_un = arbint__sub_mag(up, v0u, v0u_n, v1v, v1v_n);
    } else {
      /*  Result would be negative; shouldn't happen with valid matrix.  */
      new_un = arbint__sub_mag(up, v1v, v1v_n, v0u, v0u_n);
    }
  } else {
    if (v1v_n > v0u_n ||
        (v1v_n == v0u_n && arbint_cmp_mag_limbs(v1v, v1v_n, v0u, v0u_n) >= 0)) {
      new_un = arbint__sub_mag(up, v1v, v1v_n, v0u, v0u_n);
    } else {
      new_un = arbint__sub_mag(up, v0u, v0u_n, v1v, v1v_n);
    }
  }

  /*  Compute new_v = w1*v - w0*u (or w0*u - w1*v if odd).  */
  if (even) {
    if (w1v_n > w0u_n ||
        (w1v_n == w0u_n && arbint_cmp_mag_limbs(w1v, w1v_n, w0u, w0u_n) >= 0)) {
      new_vn = arbint__sub_mag(vp, w1v, w1v_n, w0u, w0u_n);
    } else {
      new_vn = arbint__sub_mag(vp, w0u, w0u_n, w1v, w1v_n);
    }
  } else {
    if (w0u_n > w1v_n ||
        (w0u_n == w1v_n && arbint_cmp_mag_limbs(w0u, w0u_n, w1v, w1v_n) >= 0)) {
      new_vn = arbint__sub_mag(vp, w0u, w0u_n, w1v, w1v_n);
    } else {
      new_vn = arbint__sub_mag(vp, w1v, w1v_n, w0u, w0u_n);
    }
  }

  *un = new_un;
  *vn = new_vn;
}

/*  Binary GCD (Stein's algorithm) - fallback for small operands.

    After Lehmer reduction brings operands below ARBINT_GCD_EUCLID_THRESHOLD,
    we finish with Binary GCD which is efficient for small operands.  */
static arbint_err_t ARBINT_LEHMER_FALLBACK_FN(arbint_t g, arbint_limb_t * up,
                                              size_t un, arbint_limb_t * vp,
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

  /*  Result = u << common.  */
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

/*  Lehmer GCD - main entry point.

    Entry point called from arbint_gcd() when operands exceed
    ARBINT_LEHMER_THRESHOLD (32 limbs).

    Parameters:
      g     - output: receives gcd(|a|, |b|)
      ap/an - input a as limb array with size
      bp/bn - input b as limb array with size
      alloc - allocator for temporary workspace

    Returns: ARBINT_OK on success, ARBINT_ENOMEM if allocation fails.  */
arbint_err_t ARBINT_LEHMER_ENTRY_FN(arbint_t g, const arbint_limb_t * ap,
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

  /*  Allocate workspace: u, v, and scratch for matrix application.
      Each of u and v needs max_n + 2 limbs.
      Scratch needs 4 * (max_n + 2) for four temporary products.  */
  cap = max_n + 2u;
  scratch = arbint_alloc_limbs(alloc, cap * 6u);
  if (scratch == NULL)
    return ARBINT_ENOMEM;

  up = scratch;
  vp = scratch + cap;
  work = scratch + cap * 2u;

  /*  Copy inputs to workspace.  */
  memcpy(up, ap, an * sizeof(arbint_limb_t));
  memcpy(vp, bp, bn * sizeof(arbint_limb_t));
  un = an;
  vn = bn;

  /*  Ensure u >= v.  */
  if (arbint_cmp_mag_limbs(up, un, vp, vn) < 0) {
    arbint_limb_t * tmp = up;
    up = vp;
    vp = tmp;
    {
      size_t t = un;
      un = vn;
      vn = t;
    }
  }

  /*  Main Lehmer loop.  */
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
      /*  Extract top 2 limbs of u and v.  */
      a1 = up[un - 1u];
      a0 = (un >= 2u) ? up[un - 2u] : 0u;
      b1 = vp[vn - 1u];
      b0 = (vn >= 2u) ? vp[vn - 2u] : 0u;

      /*  Try Lehmer simulation.  */
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
        {
          size_t t = un;
          un = vn;
          vn = t;
        }
      }
    } else {
      /*  Apply matrix transformation.  */
      ARBINT_LEHMER_APPLY_FN(up, &un, vp, &vn, &m, even, work, cap * 4u);

      /*  Normalize.  */
      un = arbint_norm_used(up, un);
      vn = arbint_norm_used(vp, vn);

      /*  Ensure u >= v.  */
      if (arbint_cmp_mag_limbs(up, un, vp, vn) < 0) {
        arbint_limb_t * tmp = up;
        up = vp;
        vp = tmp;
        {
          size_t t = un;
          un = vn;
          vn = t;
        }
      }
    }
  }

  /*  Finish with binary GCD on reduced operands.  */
  rc = ARBINT_LEHMER_FALLBACK_FN(g, up, un, vp, vn, alloc);

  arbint_free_limbs(alloc, scratch);
  return rc;
}

/*  Clean up macro definitions.  */
#undef ARBINT_LEHMER_MUL_1_FN
#undef ARBINT_LEHMER_APPLY_FN
#undef ARBINT_LEHMER_FALLBACK_FN
#undef ARBINT_LEHMER_ENTRY_FN
#undef ARBINT_LEHMER_UMUL_FN
