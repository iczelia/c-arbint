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

#ifndef ARBINT_INTERNAL_UTIL_H
#define ARBINT_INTERNAL_UTIL_H

#include "arbint_base.h"

#include <stddef.h>
#include <stdint.h>

/*  Error propagation with goto cleanup.
    Requires: 'rc' variable and 'cleanup' label in scope.  */
#define ARBINT_TRY(expr)                                                       \
  do {                                                                         \
    rc = (expr);                                                               \
    if (rc != ARBINT_OK)                                                       \
      goto cleanup;                                                            \
  } while (0)

/*  NULL pointer validation macros for public API entry points.  */
#define ARBINT_CHECK_NONNULL_1(a)                                              \
  do {                                                                         \
    if ((a) == NULL)                                                           \
      return ARBINT_EINVAL;                                                    \
  } while (0)

#define ARBINT_CHECK_NONNULL_2(a, b)                                           \
  do {                                                                         \
    if ((a) == NULL || (b) == NULL)                                            \
      return ARBINT_EINVAL;                                                    \
  } while (0)

#define ARBINT_CHECK_NONNULL_3(a, b, c)                                        \
  do {                                                                         \
    if ((a) == NULL || (b) == NULL || (c) == NULL)                             \
      return ARBINT_EINVAL;                                                    \
  } while (0)

#define ARBINT_CHECK_NONNULL_4(a, b, c, d)                                     \
  do {                                                                         \
    if ((a) == NULL || (b) == NULL || (c) == NULL || (d) == NULL)              \
      return ARBINT_EINVAL;                                                    \
  } while (0)

/*  Shared lazy-init helper for static dispatch pointers.  */
#define ARBINT_LAZY_INIT(impl_var, selector_fn)                               \
  do {                                                                        \
    if ((impl_var) == NULL)                                                   \
      (impl_var) = (selector_fn) ();                                          \
  } while (0)

/*  Extract allocator from a single arbint_t.  */
static inline const arbint_alloc_t *
arbint_get_alloc_from_obj(const arbint_t x) {
  if (x == NULL || x[0]._ctx == NULL || x[0]._ctx->a.realloc == NULL)
    return NULL;
  return &x[0]._ctx->a;
}

/*  Extract context from a single arbint_t.  */
static inline arbint_ctx_t * arbint_get_ctx_from_obj(const arbint_t x) {
  if (x == NULL)
    return NULL;
  return x[0]._ctx;
}

/*  Select first non-NULL context from 2 arbint_t values.  */
static inline arbint_ctx_t * arbint_pick_ctx2(const arbint_t a,
                                               const arbint_t b) {
  arbint_ctx_t * ctx = arbint_get_ctx_from_obj(a);
  if (ctx != NULL)
    return ctx;
  return arbint_get_ctx_from_obj(b);
}

/*  Select first non-NULL context from 3 arbint_t values.  */
static inline arbint_ctx_t * arbint_pick_ctx3(const arbint_t a,
                                               const arbint_t b,
                                               const arbint_t c) {
  arbint_ctx_t * ctx = arbint_pick_ctx2(a, b);
  if (ctx != NULL)
    return ctx;
  return arbint_get_ctx_from_obj(c);
}

static inline const arbint_alloc_t *
arbint_pick_alloc3(const arbint_t a, const arbint_t b, const arbint_t c) {
  const arbint_alloc_t * alloc;

  alloc = arbint_get_alloc_from_obj(a);
  if (alloc != NULL)
    return alloc;
  alloc = arbint_get_alloc_from_obj(b);
  if (alloc != NULL)
    return alloc;
  return arbint_get_alloc_from_obj(c);
}

static inline const arbint_alloc_t * arbint_pick_alloc4(const arbint_t a,
                                                        const arbint_t b,
                                                        const arbint_t c,
                                                        const arbint_t d) {
  const arbint_alloc_t * alloc;

  alloc = arbint_get_alloc_from_obj(a);
  if (alloc != NULL)
    return alloc;
  alloc = arbint_get_alloc_from_obj(b);
  if (alloc != NULL)
    return alloc;
  alloc = arbint_get_alloc_from_obj(c);
  if (alloc != NULL)
    return alloc;
  return arbint_get_alloc_from_obj(d);
}

/*  Convert bit index to limb index.  */
static inline size_t arbint_bit_to_limb_idx(size_t bit_idx) {
  return bit_idx / ARBINT_LIMB_BITS;
}

/*  Get bit position within a limb (0 to ARBINT_LIMB_BITS-1).  */
static inline unsigned arbint_bit_in_limb(size_t bit_idx) {
  return (unsigned) (bit_idx % ARBINT_LIMB_BITS);
}

/*  Number of limbs needed for k bits (ceiling division).  */
static inline size_t arbint_limbs_for_bits(size_t k) {
  return (k + ARBINT_LIMB_BITS - 1u) / ARBINT_LIMB_BITS;
}

static inline int arbint_u32_is_pow2(uint32_t v) {
  return v != 0u && (v & (v - 1u)) == 0u;
}

static inline size_t arbint_mag_ctz_or_size_max(const arbint_limb_t * p,
                                                size_t n) {
  size_t i;
  size_t ctz = 0u;

  if (n == 0u)
    return SIZE_MAX;

  for (i = 0u; i < n && p[i] == 0u; ++i)
    ctz += ARBINT_LIMB_BITS;

  if (i < n)
    ctz += arbint_ctz_limb(p[i]);

  return ctz;
}

static inline size_t arbint_mag_pow2_bit_or_size_max(const arbint_limb_t * p,
                                                     size_t n) {
  size_t i;
  size_t nonzero_idx = SIZE_MAX;

  for (i = 0u; i < n; ++i) {
    if (p[i] != 0u) {
      if (nonzero_idx != SIZE_MAX)
        return SIZE_MAX;
      nonzero_idx = i;
    }
  }

  if (nonzero_idx == SIZE_MAX)
    return SIZE_MAX;

  {
    arbint_limb_t limb = p[nonzero_idx];
    if ((limb & (limb - 1u)) != 0u)
      return SIZE_MAX;
    return nonzero_idx * ARBINT_LIMB_BITS + arbint_ctz_limb(limb);
  }
}

#endif /* ARBINT_INTERNAL_UTIL_H */
