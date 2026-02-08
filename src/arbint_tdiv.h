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

#ifndef ARBINT_TDIV_H
#define ARBINT_TDIV_H

#include "arbint_addsub.h"

arbint_err_t arbint_get_mag_view(const arbint_t x, const arbint_limb_t ** xp,
                                 size_t * xn, int * sign);

arbint_err_t arbint_tdiv_qr_u32_generic_impl(arbint_t q, arbint_t r,
                                             const arbint_t n, uint32_t dmag,
                                             int dsign);

#if HAS_BMI2
arbint_err_t arbint_tdiv_qr_u32_bmi2_impl(arbint_t q, arbint_t r,
                                          const arbint_t n, uint32_t dmag,
                                          int dsign);
#endif

#endif
