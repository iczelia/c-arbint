/* arbint - portable arbitrary-precision computation library
 *
 * Copyright (C) 2026 Kamila Szewczyk (k@iczelia.net)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef ARBINT_MUL_H
#define ARBINT_MUL_H

#include "arbint_addsub.h"
#include "arbint_base.h"

int arbint_mul_cap(size_t an, size_t bn, size_t * out);

arbint_err_t arbint_mul_impl_generic(arbint_t rop, const arbint_t a,
                                     const arbint_t b);

#if HAS_BMI2
arbint_err_t arbint_mul_impl_bmi2(arbint_t rop, const arbint_t a,
                                  const arbint_t b);
#endif

#define ARBINT_KARATSUBA_THRESHOLD 32u
#define ARBINT_HALF_BITS (ARBINT_LIMB_BITS / 2u)
#define ARBINT_HALF_MASK ((((arbint_limb_t) 1u) << ARBINT_HALF_BITS) - 1u)

#endif
