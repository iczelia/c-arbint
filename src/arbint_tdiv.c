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

/*  Truncated division public API.
    Rounds quotient toward zero; remainder has same sign as dividend.  */

#include "arbint_div.h"

/* ------------------------------------------------------------------ */
/*  arbint / arbint truncated division.                               */
/* ------------------------------------------------------------------ */

arbint_err_t arbint_tdiv_qr(arbint_t q, arbint_t r, const arbint_t n,
                            const arbint_t d) {
  if (q == NULL || r == NULL)
    return ARBINT_EINVAL;
  return arbint_tdiv_qr_impl(q, r, n, d);
}

arbint_err_t arbint_tdiv_q(arbint_t q, const arbint_t n, const arbint_t d) {
  return arbint_tdiv_qr_impl(q, NULL, n, d);
}

arbint_err_t arbint_tdiv_r(arbint_t r, const arbint_t n, const arbint_t d) {
  return arbint_tdiv_qr_impl(NULL, r, n, d);
}

/* ------------------------------------------------------------------ */
/*  arbint / u32 truncated division.                                  */
/* ------------------------------------------------------------------ */

arbint_err_t arbint_tdiv_qr_u32(arbint_t q, arbint_t r, const arbint_t n,
                                uint32_t d) {
  if (q == NULL || r == NULL)
    return ARBINT_EINVAL;
  if (d == 0u)
    return ARBINT_EZERO;
  return arbint_div_qr_u32_dispatch(q, r, n, d, 1);
}

arbint_err_t arbint_tdiv_q_u32(arbint_t q, const arbint_t n, uint32_t d) {
  if (d == 0u)
    return ARBINT_EZERO;
  return arbint_div_qr_u32_dispatch(q, NULL, n, d, 1);
}

arbint_err_t arbint_tdiv_r_u32(arbint_t r, const arbint_t n, uint32_t d) {
  if (d == 0u)
    return ARBINT_EZERO;
  return arbint_div_qr_u32_dispatch(NULL, r, n, d, 1);
}

/* ------------------------------------------------------------------ */
/*  arbint / i32 truncated division.                                  */
/* ------------------------------------------------------------------ */

arbint_err_t arbint_tdiv_qr_i32(arbint_t q, arbint_t r, const arbint_t n,
                                int32_t d) {
  if (q == NULL || r == NULL)
    return ARBINT_EINVAL;
  if (d == 0)
    return ARBINT_EZERO;
  return arbint_div_qr_u32_dispatch(q, r, n, arbint_i32_mag(d),
                                    (d < 0) ? -1 : 1);
}

arbint_err_t arbint_tdiv_q_i32(arbint_t q, const arbint_t n, int32_t d) {
  if (d == 0)
    return ARBINT_EZERO;
  return arbint_div_qr_u32_dispatch(q, NULL, n, arbint_i32_mag(d),
                                    (d < 0) ? -1 : 1);
}

arbint_err_t arbint_tdiv_r_i32(arbint_t r, const arbint_t n, int32_t d) {
  if (d == 0)
    return ARBINT_EZERO;
  return arbint_div_qr_u32_dispatch(NULL, r, n, arbint_i32_mag(d),
                                    (d < 0) ? -1 : 1);
}
