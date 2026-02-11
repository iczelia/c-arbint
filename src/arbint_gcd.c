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
#include "arbint_shift.h"
#include "config.h"

#include <string.h>

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

/*  Binary GCD (Stein's algorithm) for large operands.
    Uses internal limb-level operations for performance.  */
static arbint_err_t arbint_gcd_binary(arbint_t g, const arbint_limb_t * ap,
                                      size_t an, const arbint_limb_t * bp,
                                      size_t bn,
                                      const arbint_alloc_t * alloc) {
  arbint_limb_t * up = NULL;
  arbint_limb_t * vp = NULL;
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

  /*  Allocate workspace for u and v.
      Need extra space for the final left shift.  */
  cap = max_n + 1u;
  up = arbint_alloc_limbs(alloc, cap);
  if (up == NULL) {
    rc = ARBINT_ENOMEM;
    goto cleanup;
  }
  vp = arbint_alloc_limbs(alloc, cap);
  if (vp == NULL) {
    rc = ARBINT_ENOMEM;
    goto cleanup;
  }

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
  arbint_free_limbs(alloc, vp);
  arbint_free_limbs(alloc, up);
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
    Result is always non-negative.  */
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

  /*  Large operand path: binary GCD with internal limb ops.  */
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

  rc = arbint_gcd_binary(g, ap, an, bp, bn, alloc);

  arbint_free_limbs(alloc, b_copy);
  arbint_free_limbs(alloc, a_copy);

  return rc;
}

/*  arbint_gcd_u32: GCD with a uint32_t operand.  */
ARBINT_API arbint_err_t arbint_gcd_u32(arbint_t g, const arbint_t a,
                                       uint32_t b) {
  arbint_t rem;
  size_t an;
  arbint_err_t rc;
  uint32_t r;
  int init_rem = 0;

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

  /*  Reduce |a| mod b to get a uint32_t remainder.  */
  rc = arbint_init(rem, g[0]._ctx);
  if (rc != ARBINT_OK)
    return rc;
  init_rem = 1;

  rc = arbint_tdiv_r_u32(rem, a, b);
  if (rc != ARBINT_OK)
    goto cleanup;

  if (arbint_is_zero(rem)) {
    arbint_clear(rem);
    return arbint_set_u32(g, b);
  }

  /*  Get absolute value of remainder as u32.  */
  rc = arbint_abs(rem, rem);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_get_u32(rem, &r);
  if (rc != ARBINT_OK)
    goto cleanup;

  arbint_clear(rem);
  init_rem = 0;

  /*  Now compute gcd(r, b) with hardware division.  */
  return arbint_set_u32(g, arbint_gcd_u32u32(r, b));

cleanup:
  if (init_rem)
    arbint_clear(rem);
  return rc;
}

/*  arbint_lcm: compute least common multiple.
    lcm(a, b) = (|a| / gcd(a, b)) * |b|
    Result is always non-negative.  */
ARBINT_API arbint_err_t arbint_lcm(arbint_t l, const arbint_t a,
                                   const arbint_t b) {
  arbint_t gcd_val, quotient;
  arbint_err_t rc;
  arbint_ctx_t * ctx;
  int init_gcd = 0;
  int init_quot = 0;
  size_t an;
  size_t bn;

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

  rc = arbint_init(quotient, ctx);
  if (rc != ARBINT_OK)
    goto cleanup;
  init_quot = 1;

  /*  Compute gcd(a, b).  */
  rc = arbint_gcd(gcd_val, a, b);
  if (rc != ARBINT_OK)
    goto cleanup;

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

/*  arbint_lcm_u32: LCM with a uint32_t operand.  */
ARBINT_API arbint_err_t arbint_lcm_u32(arbint_t l, const arbint_t a,
                                       uint32_t b) {
  arbint_t tmp, rem;
  arbint_err_t rc;
  arbint_ctx_t * ctx;
  uint32_t g;
  uint32_t r;
  size_t an;
  int init_tmp = 0;
  int init_rem = 0;

  if (l == NULL || a == NULL)
    return ARBINT_EINVAL;

  an = arbint_abs_sz(a[0]._sz);

  /*  lcm(0, b) = lcm(a, 0) = 0.  */
  if (an == 0u || b == 0u) {
    arbint_zero(l);
    return ARBINT_OK;
  }

  ctx = l[0]._ctx;

  /*  Compute gcd(|a| mod b, b) to get a u32 gcd.  */
  rc = arbint_init(rem, ctx);
  if (rc != ARBINT_OK)
    return rc;
  init_rem = 1;

  rc = arbint_tdiv_r_u32(rem, a, b);
  if (rc != ARBINT_OK)
    goto cleanup;

  if (arbint_is_zero(rem)) {
    g = b;
  } else {
    /*  Get absolute value of remainder as u32.  */
    rc = arbint_abs(rem, rem);
    if (rc != ARBINT_OK)
      goto cleanup;

    rc = arbint_get_u32(rem, &r);
    if (rc != ARBINT_OK)
      goto cleanup;
    g = arbint_gcd_u32u32(r, b);
  }

  arbint_clear(rem);
  init_rem = 0;

  /*  l = (|a| / g) * b.  */
  rc = arbint_init(tmp, ctx);
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
  if (init_rem)
    arbint_clear(rem);
  return rc;
}
