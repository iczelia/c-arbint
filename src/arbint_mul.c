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

#include "arbint_mul.h"

#include "config.h"

#include "arbint_cpu.h"

typedef arbint_err_t (*arbint_mul_impl_fn_t)(arbint_t rop, const arbint_t a,
                                             const arbint_t b);

static arbint_mul_impl_fn_t arbint_select_mul_impl(void) {
#if HAS_BMI2_ALWAYS
  return arbint_mul_impl_bmi2;
#elif HAS_BMI2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2)
             ? arbint_mul_impl_bmi2
             : arbint_mul_impl_generic;
#else
  return arbint_mul_impl_generic;
#endif
}

int arbint_mul_cap(size_t an, size_t bn, size_t * out) {
  if (out == NULL)
    return 0;
  if (bn > SIZE_MAX - 1u)
    return 0;
  if (an > SIZE_MAX - bn - 1u)
    return 0;
  *out = an + bn + 1u;
  return 1;
}

arbint_err_t arbint_mul(arbint_t rop, const arbint_t a, const arbint_t b) {
  static arbint_mul_impl_fn_t impl = NULL;

  if (impl == NULL)
    impl = arbint_select_mul_impl();

  return impl(rop, a, b);
}

arbint_err_t arbint_sqr(arbint_t rop, const arbint_t a) {
  if (rop == NULL || a == NULL)
    return ARBINT_EINVAL;
  return arbint_mul(rop, a, a);
}
