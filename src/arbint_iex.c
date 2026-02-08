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

#include "arbint_iex.h"
#include "arbint_addsub.h"

#include <string.h>

#define ARBINT_LIMB_BYTES (ARBINT_LIMB_BITS / 8u)

static size_t arbint_limb_used_bytes(arbint_limb_t x) {
  size_t n = 0u;
  while (x != (arbint_limb_t) 0u) {
    ++n;
    x >>= 8u;
  }
  return n;
}

arbint_err_t arbint_import(arbint_t rop, const void * buf, size_t nbytes) {
  const uint8_t * bp;
  size_t used_bytes;
  size_t limbs_needed;
  size_t i;

  if (rop == NULL)
    return ARBINT_EINVAL;
  if (buf == NULL && nbytes != 0u)
    return ARBINT_EINVAL;

  bp = (const uint8_t *) buf;
  used_bytes = nbytes;
  while (used_bytes != 0u && bp[used_bytes - 1u] == 0u)
    --used_bytes;

  if (used_bytes == 0u) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  if (used_bytes > SIZE_MAX - (ARBINT_LIMB_BYTES - 1u))
    return ARBINT_EOVERFLOW;
  limbs_needed = (used_bytes + ARBINT_LIMB_BYTES - 1u) / ARBINT_LIMB_BYTES;
  if (limbs_needed > PTRDIFF_MAX)
    return ARBINT_EOVERFLOW;

  {
    arbint_err_t rc = arbint_resize(rop, limbs_needed);
    if (rc != ARBINT_OK)
      return rc;
  }

  memset(ARBINT_LIMBS(rop), 0, limbs_needed * sizeof(arbint_limb_t));

  for (i = 0u; i < used_bytes; ++i) {
    size_t li = i / ARBINT_LIMB_BYTES;
    size_t bi = i % ARBINT_LIMB_BYTES;
    ARBINT_LIMBS(rop)[li] |= ((arbint_limb_t) bp[i]) << (bi * 8u);
  }

  limbs_needed = arbint_norm_used(ARBINT_LIMBS(rop), limbs_needed);
  if (limbs_needed > PTRDIFF_MAX)
    return ARBINT_EOVERFLOW;

  rop[0]._sz = (ptrdiff_t) limbs_needed;
  return ARBINT_OK;
}

arbint_err_t arbint_export(const arbint_t op, void ** out_buf,
                           size_t * out_nbytes) {
  size_t used;
  size_t top_nbytes;
  size_t nbytes;
  uint8_t * dst;
  const arbint_limb_t * p;
  size_t i;
  size_t j;
  size_t o;

  if (op == NULL || out_buf == NULL || out_nbytes == NULL)
    return ARBINT_EINVAL;

  *out_buf = NULL;
  *out_nbytes = 0u;

  used = arbint_abs_sz(op[0]._sz);
  if (used == 0u)
    return ARBINT_OK;
  if (op[0]._ptr == NULL || used > op[0]._cap)
    return ARBINT_EINVAL;
  if (op[0]._ctx == NULL || op[0]._ctx->a.realloc == NULL)
    return ARBINT_EINVAL;

  p = ARBINT_CLIMBS(op);
  top_nbytes = arbint_limb_used_bytes(p[used - 1u]);
  if (top_nbytes == 0u)
    return ARBINT_EINVAL;

  if (used - 1u > (SIZE_MAX - top_nbytes) / ARBINT_LIMB_BYTES)
    return ARBINT_EOVERFLOW;
  nbytes = (used - 1u) * ARBINT_LIMB_BYTES + top_nbytes;

  dst = (uint8_t *) op[0]._ctx->a.realloc(op[0]._ctx->a.ud, NULL, nbytes);
  if (dst == NULL)
    return ARBINT_ENOMEM;

  o = 0u;
  for (i = 0u; i < used; ++i) {
    arbint_limb_t w = p[i];
    size_t lim_nbytes = (i + 1u == used) ? top_nbytes : ARBINT_LIMB_BYTES;
    for (j = 0u; j < lim_nbytes; ++j) {
      dst[o++] = (uint8_t) (w & (arbint_limb_t) 0xffu);
      w >>= 8u;
    }
  }

  *out_buf = dst;
  *out_nbytes = nbytes;
  return ARBINT_OK;
}
