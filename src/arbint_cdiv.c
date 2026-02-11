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

/*  Ceiling division public API.
    Rounds quotient toward +Inf; remainder has opposite sign to divisor.  */

#include "arbint_div.h"

/*  arbint / arbint ceiling division.  */

arbint_err_t arbint_cdiv_qr(arbint_t q, arbint_t r, const arbint_t n,
                            const arbint_t d) {
  if (q == NULL || r == NULL)
    return ARBINT_EINVAL;
  return arbint_div_qr_mode_impl(q, r, n, d, ARBINT_DIV_CEIL);
}

arbint_err_t arbint_cdiv_q(arbint_t q, const arbint_t n, const arbint_t d) {
  return arbint_div_qr_mode_impl(q, NULL, n, d, ARBINT_DIV_CEIL);
}

arbint_err_t arbint_cdiv_r(arbint_t r, const arbint_t n, const arbint_t d) {
  return arbint_div_qr_mode_impl(NULL, r, n, d, ARBINT_DIV_CEIL);
}
