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

arbint_err_t
arbint_div_mag_single_limb_generic(const arbint_limb_t * np, size_t nn,
                                   arbint_limb_t d_limb, arbint_limb_t * qp,
                                   size_t * q_used, arbint_limb_t * rem_out);

arbint_err_t arbint_div_mag_knuth(const arbint_limb_t * np, size_t nn,
                                  const arbint_limb_t * dp, size_t dn,
                                  arbint_limb_t * qp, arbint_limb_t * rp);

/*  Reduce x mod d in-place using precomputed Barrett parameters.

    This is a remainder-only reduction optimized for modular exponentiation,
    where the reciprocal is precomputed once and reused for many reductions.

    Parameters:
      x      - Input/output: value to reduce in-place
      d_norm - Normalized divisor: d << shift (MSB set)
      di     - Precomputed reciprocal from arbint_prepare_barrett(d_norm)
      shift  - Normalization shift: arbint_clz_limb(d)

    Truncated division semantics: remainder sign matches dividend sign.  */
arbint_err_t arbint_mod_u32_barrett_generic(arbint_t x, arbint_limb_t d_norm,
                                             arbint_limb_t di, unsigned shift);

/*  Dispatched wrapper for Barrett reduction (selects generic or BMI2).  */
arbint_err_t arbint_mod_u32_barrett(arbint_t x, arbint_limb_t d_norm,
                                     arbint_limb_t di, unsigned shift);

#if HAS_BMI2
arbint_err_t arbint_tdiv_qr_u32_bmi2_impl(arbint_t q, arbint_t r,
                                          const arbint_t n, uint32_t dmag,
                                          int dsign);

arbint_err_t arbint_div_mag_single_limb_bmi2(const arbint_limb_t * np,
                                             size_t nn, arbint_limb_t d_limb,
                                             arbint_limb_t * qp,
                                             size_t * q_used,
                                             arbint_limb_t * rem_out);

arbint_err_t arbint_mod_u32_barrett_bmi2(arbint_t x, arbint_limb_t d_norm,
                                          arbint_limb_t di, unsigned shift);
#endif /* HAS_BMI2 */

#endif /* ARBINT_TDIV_H */
