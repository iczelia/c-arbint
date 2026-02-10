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

#include "arbint_div.h"

#include "config.h"

#include "arbint_cpu.h"

#include <string.h>

/* ------------------------------------------------------------------ */
/*  Function pointer types for runtime dispatch.                      */
/* ------------------------------------------------------------------ */

typedef arbint_err_t (*arbint_div_qr_u32_impl_fn_t)(arbint_t q, arbint_t r,
                                                    const arbint_t n,
                                                    uint32_t dmag, int dsign);

typedef arbint_err_t (*arbint_div_mag_single_limb_fn_t)(
    const arbint_limb_t * np, size_t nn, arbint_limb_t d_limb,
    arbint_limb_t * qp, size_t * q_used, arbint_limb_t * rem_out);

typedef arbint_err_t (*arbint_mod_u32_barrett_fn_t)(arbint_t x,
                                                    arbint_limb_t d_norm,
                                                    arbint_limb_t di,
                                                    unsigned shift);

/* ------------------------------------------------------------------ */
/*  Runtime dispatch selectors.                                       */
/* ------------------------------------------------------------------ */

static arbint_div_qr_u32_impl_fn_t arbint_select_div_qr_u32_impl(void) {
#if HAS_BMI2_ALWAYS
  return arbint_div_qr_u32_bmi2_impl;
#elif HAS_BMI2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2)
             ? arbint_div_qr_u32_bmi2_impl
             : arbint_div_qr_u32_generic_impl;
#else
  return arbint_div_qr_u32_generic_impl;
#endif
}

static arbint_div_mag_single_limb_fn_t arbint_select_div_mag_single_limb(void) {
#if HAS_BMI2_ALWAYS
  return arbint_div_mag_single_limb_bmi2;
#elif HAS_BMI2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2)
             ? arbint_div_mag_single_limb_bmi2
             : arbint_div_mag_single_limb_generic;
#else
  return arbint_div_mag_single_limb_generic;
#endif
}

static arbint_mod_u32_barrett_fn_t arbint_select_mod_u32_barrett(void) {
#if HAS_BMI2_ALWAYS
  return arbint_mod_u32_barrett_bmi2;
#elif HAS_BMI2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2)
             ? arbint_mod_u32_barrett_bmi2
             : arbint_mod_u32_barrett_generic;
#else
  return arbint_mod_u32_barrett_generic;
#endif
}

/* ------------------------------------------------------------------ */
/*  Dispatched Barrett reduction.                                     */
/* ------------------------------------------------------------------ */

arbint_err_t arbint_mod_u32_barrett(arbint_t x, arbint_limb_t d_norm,
                                    arbint_limb_t di, unsigned shift) {
  static arbint_mod_u32_barrett_fn_t impl = NULL;

  if (impl == NULL)
    impl = arbint_select_mod_u32_barrett();

  return impl(x, d_norm, di, shift);
}

/* ------------------------------------------------------------------ */
/*  Dispatched u32 division.                                          */
/* ------------------------------------------------------------------ */

arbint_err_t arbint_div_qr_u32_dispatch(arbint_t q, arbint_t r,
                                        const arbint_t n, uint32_t dmag,
                                        int dsign) {
  static arbint_div_qr_u32_impl_fn_t impl = NULL;

  if (impl == NULL)
    impl = arbint_select_div_qr_u32_impl();

  return impl(q, r, n, dmag, dsign);
}

/* ------------------------------------------------------------------ */
/*  Allocator helpers.                                                */
/* ------------------------------------------------------------------ */

static const arbint_alloc_t * arbint_get_alloc_from(const arbint_t x) {
  if (x == NULL || x[0]._ctx == NULL || x[0]._ctx->a.realloc == NULL)
    return NULL;
  return &x[0]._ctx->a;
}

static const arbint_alloc_t * arbint_pick_alloc(const arbint_t a,
                                                const arbint_t b,
                                                const arbint_t c,
                                                const arbint_t d) {
  const arbint_alloc_t * alloc;

  alloc = arbint_get_alloc_from(a);
  if (alloc != NULL)
    return alloc;
  alloc = arbint_get_alloc_from(b);
  if (alloc != NULL)
    return alloc;
  alloc = arbint_get_alloc_from(c);
  if (alloc != NULL)
    return alloc;
  return arbint_get_alloc_from(d);
}

/* ------------------------------------------------------------------ */
/*  Magnitude view and signed assignment.                             */
/* ------------------------------------------------------------------ */

arbint_err_t arbint_get_mag_view(const arbint_t x, const arbint_limb_t ** xp,
                                 size_t * xn, int * sign) {
  size_t used;

  if (x == NULL || xp == NULL || xn == NULL || sign == NULL)
    return ARBINT_EINVAL;

  used = arbint_abs_sz(x[0]._sz);
  if (used == 0u) {
    *xp = NULL;
    *xn = 0u;
    *sign = 0;
    return ARBINT_OK;
  }

  if (x[0]._ptr == NULL || used > x[0]._cap)
    return ARBINT_EINVAL;

  *xp = ARBINT_CLIMBS(x);
  *xn = arbint_norm_used(*xp, used);
  *sign = (x[0]._sz < 0) ? -1 : 1;
  if (*xn == 0u)
    *sign = 0;
  return ARBINT_OK;
}

static arbint_err_t arbint_set_mag_signed(arbint_t x,
                                          const arbint_limb_t * mag,
                                          size_t used, int sign) {
  arbint_err_t rc;

  if (x == NULL)
    return ARBINT_OK;

  if (used == 0u) {
    arbint_zero(x);
    return ARBINT_OK;
  }

  rc = arbint_resize(x, used);
  if (rc != ARBINT_OK)
    return rc;

  memcpy(ARBINT_LIMBS(x), mag, used * sizeof(arbint_limb_t));
  if (!arbint_set_signed_sz(x, used, sign))
    return ARBINT_EOVERFLOW;

  return ARBINT_OK;
}

/* ------------------------------------------------------------------ */
/*  Core truncated magnitude division.                                */
/* ------------------------------------------------------------------ */

static arbint_err_t
arbint_tdiv_qr_mag_impl(arbint_t q, arbint_t r, const arbint_limb_t * np,
                        size_t nn, int nsign, const arbint_limb_t * dp,
                        size_t dn, int dsign, const arbint_alloc_t * alloc) {
  size_t qcap;
  size_t rcap;
  arbint_limb_t * qmag;
  arbint_limb_t * rmag;
  size_t q_used = 0u;
  size_t r_used;
  int qsign;
  int rsign;
  arbint_err_t rc;

  if (dn == 0u)
    return ARBINT_EZERO;
  if ((q == NULL && r == NULL) || alloc == NULL)
    return ARBINT_EINVAL;

  if (nn == 0u) {
    if (q != NULL)
      arbint_zero(q);
    if (r != NULL)
      arbint_zero(r);
    return ARBINT_OK;
  }

  if (nn < dn) {
    if (q != NULL)
      arbint_zero(q);
    if (r != NULL)
      return arbint_set_mag_signed(r, np, nn, nsign);
    return ARBINT_OK;
  }

  qcap = nn;
  if (dn > SIZE_MAX - 1u)
    return ARBINT_EOVERFLOW;
  rcap = dn + 1u;

  qmag = NULL;
  if (q != NULL) {
    qmag = arbint_alloc_limbs(alloc, qcap);
    if (qmag == NULL)
      return ARBINT_ENOMEM;
  }

  rmag = arbint_alloc_limbs(alloc, rcap);
  if (rmag == NULL) {
    arbint_free_limbs(alloc, qmag);
    return ARBINT_ENOMEM;
  }

  if (dn >= 2u) {
    rc = arbint_div_mag_knuth(np, nn, dp, dn, qmag, rmag);
    if (rc != ARBINT_OK) {
      arbint_free_limbs(alloc, rmag);
      arbint_free_limbs(alloc, qmag);
      return rc;
    }
    if (qmag != NULL)
      q_used = arbint_norm_used(qmag, nn - dn + 1u);
    else
      q_used = 0u;
    r_used = arbint_norm_used(rmag, dn);
  } else {
    static arbint_div_mag_single_limb_fn_t impl = NULL;
    arbint_limb_t rem_limb = 0u;

    if (impl == NULL)
      impl = arbint_select_div_mag_single_limb();

    rc = impl(np, nn, dp[0], qmag, &q_used, &rem_limb);
    if (rc != ARBINT_OK) {
      arbint_free_limbs(alloc, rmag);
      arbint_free_limbs(alloc, qmag);
      return rc;
    }
    rmag[0] = rem_limb;
    r_used = (rem_limb != 0u) ? 1u : 0u;
  }

  qsign = (q == NULL || q_used == 0u || nsign == 0)
              ? 0
              : ((nsign == dsign) ? 1 : -1);
  rsign = (r_used == 0u || nsign == 0) ? 0 : nsign;

  if (q != NULL) {
    rc = arbint_set_mag_signed(q, qmag, q_used, qsign);
    if (rc != ARBINT_OK) {
      arbint_free_limbs(alloc, rmag);
      arbint_free_limbs(alloc, qmag);
      return rc;
    }
  }

  if (r != NULL) {
    rc = arbint_set_mag_signed(r, rmag, r_used, rsign);
    if (rc != ARBINT_OK) {
      arbint_free_limbs(alloc, rmag);
      arbint_free_limbs(alloc, qmag);
      return rc;
    }
  }

  arbint_free_limbs(alloc, rmag);
  arbint_free_limbs(alloc, qmag);
  return ARBINT_OK;
}

/* ------------------------------------------------------------------ */
/*  Truncated division implementation (arbint / arbint).              */
/* ------------------------------------------------------------------ */

arbint_err_t arbint_tdiv_qr_impl(arbint_t q, arbint_t r, const arbint_t n,
                                 const arbint_t d) {
  const arbint_limb_t * np;
  const arbint_limb_t * dp;
  size_t nn;
  size_t dn;
  int nsign;
  int dsign;
  const arbint_alloc_t * alloc;
  arbint_err_t rc;

  if ((q == NULL && r == NULL) || n == NULL || d == NULL)
    return ARBINT_EINVAL;

  rc = arbint_get_mag_view(n, &np, &nn, &nsign);
  if (rc != ARBINT_OK)
    return rc;

  rc = arbint_get_mag_view(d, &dp, &dn, &dsign);
  if (rc != ARBINT_OK)
    return rc;
  if (dn == 0u)
    return ARBINT_EZERO;

  alloc = arbint_pick_alloc(q, r, n, d);
  if (alloc == NULL)
    return ARBINT_EINVAL;

  return arbint_tdiv_qr_mag_impl(q, r, np, nn, nsign, dp, dn, dsign, alloc);
}

/* ------------------------------------------------------------------ */
/*  Generalized division with rounding mode.                          */
/* ------------------------------------------------------------------ */

arbint_err_t arbint_div_qr_mode_impl(arbint_t q, arbint_t r, const arbint_t n,
                                     const arbint_t d, arbint_div_mode_t mode) {
  arbint_t q_tmp;
  arbint_t r_tmp;
  arbint_err_t rc;
  int nsign;
  int dsign;
  int need_adjust;

  if (mode == ARBINT_DIV_TRUNC)
    return arbint_tdiv_qr_impl(q, r, n, d);

  nsign = arbint_signum(n);
  dsign = arbint_signum(d);

  if (dsign == 0)
    return ARBINT_EZERO;

  if (nsign == 0) {
    if (q != NULL)
      arbint_zero(q);
    if (r != NULL)
      arbint_zero(r);
    return ARBINT_OK;
  }

  if (q == NULL) {
    rc = arbint_init(q_tmp, n[0]._ctx);
    if (rc != ARBINT_OK)
      return rc;
  } else {
    q_tmp[0] = q[0];
  }

  if (r == NULL) {
    rc = arbint_init(r_tmp, n[0]._ctx);
    if (rc != ARBINT_OK) {
      if (q == NULL)
        arbint_clear(q_tmp);
      return rc;
    }
  } else {
    r_tmp[0] = r[0];
  }

  rc = arbint_tdiv_qr_impl(q_tmp, r_tmp, n, d);
  if (rc != ARBINT_OK) {
    if (q == NULL)
      arbint_clear(q_tmp);
    if (r == NULL)
      arbint_clear(r_tmp);
    return rc;
  }

  if (arbint_is_zero(r_tmp)) {
    need_adjust = 0;
  } else if (mode == ARBINT_DIV_FLOOR) {
    need_adjust = (nsign != dsign);
  } else {
    need_adjust = (nsign == dsign);
  }

  if (need_adjust) {
    if (mode == ARBINT_DIV_FLOOR) {
      rc = arbint_sub_i32(q_tmp, q_tmp, 1);
      if (rc == ARBINT_OK)
        rc = arbint_add(r_tmp, r_tmp, d);
    } else {
      rc = arbint_add_i32(q_tmp, q_tmp, 1);
      if (rc == ARBINT_OK)
        rc = arbint_sub(r_tmp, r_tmp, d);
    }
    if (rc != ARBINT_OK) {
      if (q == NULL)
        arbint_clear(q_tmp);
      if (r == NULL)
        arbint_clear(r_tmp);
      return rc;
    }
  }

  if (q == NULL) {
    arbint_clear(q_tmp);
  } else {
    q[0] = q_tmp[0];
  }

  if (r == NULL) {
    arbint_clear(r_tmp);
  } else {
    r[0] = r_tmp[0];
  }

  return ARBINT_OK;
}
