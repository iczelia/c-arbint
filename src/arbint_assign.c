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

#include "arbint_assign.h"

#include <limits.h>
#include <string.h>

static int arbint_get_abs_u64(const arbint_t x, uint64_t * out) {
  size_t used;
  const arbint_limb_t * xp;

  if (x == NULL || out == NULL)
    return 0;

  used = arbint_abs_sz(x[0]._sz);
  if (used == 0) {
    *out = 0u;
    return 1;
  }

  if (x[0]._ptr == NULL || used > x[0]._cap)
    return 0;

  xp = ARBINT_CLIMBS(x);
#if ARBINT_LIMB_BITS == 64
  if (used > 1u)
    return 0;
  *out = (uint64_t) xp[0];
#elif ARBINT_LIMB_BITS == 32
  if (used > 2u)
    return 0;
  *out = (uint64_t) xp[0];
  if (used == 2u)
    *out |= (uint64_t) xp[1] << 32;
#else
  #error "Unsupported ARBINT_LIMB_BITS"
#endif
  return 1;
}

arbint_err_t arbint_set(arbint_t rop, const arbint_t op) {
  size_t used;

  if (rop == NULL || op == NULL)
    return ARBINT_EINVAL;

  if (rop == op)
    return ARBINT_OK;

  used = arbint_abs_sz(op[0]._sz);
  if (used == 0u) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  if (op[0]._ptr == NULL || used > op[0]._cap)
    return ARBINT_EINVAL;

  if (rop[0]._cap < used) {
    arbint_err_t rc = arbint_resize(rop, used);
    if (rc != ARBINT_OK)
      return rc;
  }

  memcpy(ARBINT_LIMBS(rop), ARBINT_CLIMBS(op), used * sizeof(arbint_limb_t));
  rop[0]._sz = op[0]._sz;
  return ARBINT_OK;
}

arbint_err_t arbint_set_i32(arbint_t rop, int32_t v) {
  uint32_t mag;

  if (rop == NULL)
    return ARBINT_EINVAL;

  if (v == 0) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  if (rop[0]._cap < 1u) {
    arbint_err_t rc = arbint_resize(rop, 1u);
    if (rc != ARBINT_OK)
      return rc;
  }

  mag = (v < 0) ? (uint32_t) (-(v + 1)) + 1u : (uint32_t) v;
  ARBINT_LIMBS(rop)[0] = (arbint_limb_t) mag;
  rop[0]._sz = (v < 0) ? -1 : 1;
  return ARBINT_OK;
}

arbint_err_t arbint_set_u32(arbint_t rop, uint32_t v) {
  if (rop == NULL)
    return ARBINT_EINVAL;

  if (v == 0u) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  if (rop[0]._cap < 1u) {
    arbint_err_t rc = arbint_resize(rop, 1u);
    if (rc != ARBINT_OK)
      return rc;
  }

  ARBINT_LIMBS(rop)[0] = (arbint_limb_t) v;
  rop[0]._sz = 1;
  return ARBINT_OK;
}

arbint_err_t arbint_get_i32(const arbint_t op, int32_t * out) {
  size_t used;
  uint64_t mag;

  if (op == NULL || out == NULL)
    return ARBINT_EINVAL;

  used = arbint_abs_sz(op[0]._sz);
  if (used == 0u) {
    *out = 0;
    return ARBINT_OK;
  }

  if (!arbint_get_abs_u64(op, &mag))
    return ARBINT_EOVERFLOW;

  if (op[0]._sz > 0) {
    if (mag > (uint64_t) INT32_MAX)
      return ARBINT_EOVERFLOW;
    *out = (int32_t) mag;
    return ARBINT_OK;
  }

  if (mag > (uint64_t) INT32_MAX + 1u)
    return ARBINT_EOVERFLOW;

  if (mag == (uint64_t) INT32_MAX + 1u) {
    *out = INT32_MIN;
    return ARBINT_OK;
  }

  *out = -(int32_t) mag;
  return ARBINT_OK;
}

arbint_err_t arbint_get_u32(const arbint_t op, uint32_t * out) {
  uint64_t mag;

  if (op == NULL || out == NULL)
    return ARBINT_EINVAL;

  if (op[0]._sz < 0)
    return ARBINT_ESIGN;

  if (!arbint_get_abs_u64(op, &mag))
    return ARBINT_EOVERFLOW;

  if (mag > (uint64_t) UINT32_MAX)
    return ARBINT_EOVERFLOW;

  *out = (uint32_t) mag;
  return ARBINT_OK;
}

static int arbint_fits_u_impl(const arbint_t x, uint64_t maxv) {
  uint64_t mag;

  if (x == NULL || x[0]._sz < 0)
    return 0;

  if (!arbint_get_abs_u64(x, &mag))
    return 0;

  return mag <= maxv;
}

static int arbint_fits_i_impl(const arbint_t x, uint64_t pos_max,
                              uint64_t neg_abs_max) {
  uint64_t mag;

  if (x == NULL)
    return 0;

  if (!arbint_get_abs_u64(x, &mag))
    return 0;

  if (x[0]._sz < 0)
    return mag <= neg_abs_max;
  return mag <= pos_max;
}

int arbint_fits_u8(const arbint_t x) {
  return arbint_fits_u_impl(x, (uint64_t) UINT8_MAX);
}

int arbint_fits_u16(const arbint_t x) {
  return arbint_fits_u_impl(x, (uint64_t) UINT16_MAX);
}

int arbint_fits_u32(const arbint_t x) {
  return arbint_fits_u_impl(x, (uint64_t) UINT32_MAX);
}

int arbint_fits_u64(const arbint_t x) {
  return arbint_fits_u_impl(x, UINT64_MAX);
}

int arbint_fits_i8(const arbint_t x) {
  return arbint_fits_i_impl(x, (uint64_t) INT8_MAX, (uint64_t) INT8_MAX + 1u);
}

int arbint_fits_i16(const arbint_t x) {
  return arbint_fits_i_impl(x, (uint64_t) INT16_MAX,
                            (uint64_t) INT16_MAX + 1u);
}

int arbint_fits_i32(const arbint_t x) {
  return arbint_fits_i_impl(x, (uint64_t) INT32_MAX,
                            (uint64_t) INT32_MAX + 1u);
}

int arbint_fits_i64(const arbint_t x) {
  return arbint_fits_i_impl(x, (uint64_t) INT64_MAX,
                            (uint64_t) INT64_MAX + 1u);
}
