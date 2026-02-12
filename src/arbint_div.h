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

#ifndef ARBINT_DIV_H
#define ARBINT_DIV_H

#include "arbint_addsub.h"

/*  Division rounding modes.  */
typedef enum {
  ARBINT_DIV_TRUNC, /*  Truncate toward zero (tdiv).  */
  ARBINT_DIV_FLOOR, /*  Floor toward -infinity (fdiv).  */
  ARBINT_DIV_CEIL   /*  Ceiling toward +infinity (cdiv).  */
} arbint_div_mode_t;

/*  Get read-only view of arbint magnitude and sign.  */
arbint_err_t arbint_get_mag_view(const arbint_t x, const arbint_limb_t ** xp,
                                 size_t * xn, int * sign);

/*  Core truncated division (used by all rounding modes).  */
arbint_err_t arbint_tdiv_qr_impl(arbint_t q, arbint_t r, const arbint_t n,
                                 const arbint_t d);

/*  Generalized division with rounding mode adjustment.  */
arbint_err_t arbint_div_qr_mode_impl(arbint_t q, arbint_t r, const arbint_t n,
                                     const arbint_t d, arbint_div_mode_t mode);

/*  Platform-specific u32 division implementations.  */
arbint_err_t arbint_div_qr_u32_generic_impl(arbint_t q, arbint_t r,
                                            const arbint_t n, uint32_t dmag,
                                            int dsign);

arbint_err_t
arbint_div_mag_single_limb_generic(const arbint_limb_t * np, size_t nn,
                                   arbint_limb_t d_limb, arbint_limb_t * qp,
                                   size_t * q_used, arbint_limb_t * rem_out);

arbint_err_t arbint_div_mag_knuth(const arbint_limb_t * np, size_t nn,
                                  const arbint_limb_t * dp, size_t dn,
                                  arbint_limb_t * qp, arbint_limb_t * rp);

/*  Barrett reduction for modular exponentiation.  */
arbint_err_t arbint_mod_u32_barrett_generic(arbint_t x, arbint_limb_t d_norm,
                                            arbint_limb_t di, unsigned shift);

arbint_limb_t arbint_div_prepare_barrett_limb(arbint_limb_t d_norm);

arbint_err_t arbint_mod_u32_barrett(arbint_t x, arbint_limb_t d_norm,
                                    arbint_limb_t di, unsigned shift);

arbint_err_t arbint_tdiv_q_3_generic(arbint_t q, const arbint_t n);
arbint_err_t arbint_tdiv_q_3_dispatch(arbint_t q, const arbint_t n);

/*  Dispatched u32 division (selects optimal implementation).  */
arbint_err_t arbint_div_qr_u32_dispatch(arbint_t q, arbint_t r,
                                        const arbint_t n, uint32_t dmag,
                                        int dsign);

/*  Dispatched single-limb remainder (remainder only, no quotient).
    Uses BMI2 path when available.  */
arbint_err_t arbint_mod_mag_single_limb(const arbint_limb_t * np, size_t nn,
                                        arbint_limb_t d_limb,
                                        arbint_limb_t * rem_out);

#if HAS_BMI2
arbint_err_t arbint_div_qr_u32_bmi2_impl(arbint_t q, arbint_t r,
                                         const arbint_t n, uint32_t dmag,
                                         int dsign);

arbint_err_t arbint_div_mag_single_limb_bmi2(const arbint_limb_t * np,
                                             size_t nn, arbint_limb_t d_limb,
                                             arbint_limb_t * qp,
                                             size_t * q_used,
                                             arbint_limb_t * rem_out);

arbint_err_t arbint_mod_u32_barrett_bmi2(arbint_t x, arbint_limb_t d_norm,
                                         arbint_limb_t di, unsigned shift);
arbint_err_t arbint_tdiv_q_3_bmi2(arbint_t q, const arbint_t n);
#endif /* HAS_BMI2 */

#endif /* ARBINT_DIV_H */
