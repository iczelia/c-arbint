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

/*  Floor division public API.
    Rounds quotient toward -infinity; remainder has same sign as divisor.  */

#include "arbint_div.h"

/* ------------------------------------------------------------------ */
/*  arbint / arbint floor division.                                   */
/* ------------------------------------------------------------------ */

arbint_err_t arbint_fdiv_qr(arbint_t q, arbint_t r, const arbint_t n,
                            const arbint_t d) {
  if (q == NULL || r == NULL)
    return ARBINT_EINVAL;
  return arbint_div_qr_mode_impl(q, r, n, d, ARBINT_DIV_FLOOR);
}

arbint_err_t arbint_fdiv_q(arbint_t q, const arbint_t n, const arbint_t d) {
  return arbint_div_qr_mode_impl(q, NULL, n, d, ARBINT_DIV_FLOOR);
}

arbint_err_t arbint_fdiv_r(arbint_t r, const arbint_t n, const arbint_t d) {
  return arbint_div_qr_mode_impl(NULL, r, n, d, ARBINT_DIV_FLOOR);
}

/* ------------------------------------------------------------------ */
/*  arbint / u32 floor division.                                      */
/*                                                                    */
/*  For positive divisor, fdiv adjustment needed when:                */
/*  r != 0 && n < 0: q -= 1, r += d                                   */
/* ------------------------------------------------------------------ */

static arbint_err_t arbint_fdiv_qr_u32_impl(arbint_t q, arbint_t r,
                                            const arbint_t n, uint32_t d) {
  arbint_t q_tmp;
  arbint_t r_tmp;
  arbint_err_t rc;
  int nsign;

  if (d == 0u)
    return ARBINT_EZERO;

  nsign = arbint_signum(n);

  /*  If n >= 0, fdiv == tdiv for positive divisor.  */
  if (nsign >= 0)
    return arbint_div_qr_u32_dispatch(q, r, n, d, 1);

  /*  n < 0, d > 0: need to compute both q and r for adjustment check.  */
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

  rc = arbint_div_qr_u32_dispatch(q_tmp, r_tmp, n, d, 1);
  if (rc != ARBINT_OK) {
    if (q == NULL)
      arbint_clear(q_tmp);
    if (r == NULL)
      arbint_clear(r_tmp);
    return rc;
  }

  /*  Adjust if r != 0: q -= 1, r += d.  */
  if (!arbint_is_zero(r_tmp)) {
    rc = arbint_sub_i32(q_tmp, q_tmp, 1);
    if (rc == ARBINT_OK)
      rc = arbint_add_u32(r_tmp, r_tmp, d);
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

arbint_err_t arbint_fdiv_qr_u32(arbint_t q, arbint_t r, const arbint_t n,
                                uint32_t d) {
  if (q == NULL || r == NULL)
    return ARBINT_EINVAL;
  return arbint_fdiv_qr_u32_impl(q, r, n, d);
}

arbint_err_t arbint_fdiv_q_u32(arbint_t q, const arbint_t n, uint32_t d) {
  return arbint_fdiv_qr_u32_impl(q, NULL, n, d);
}

arbint_err_t arbint_fdiv_r_u32(arbint_t r, const arbint_t n, uint32_t d) {
  return arbint_fdiv_qr_u32_impl(NULL, r, n, d);
}
