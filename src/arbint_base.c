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

#include "arbint_base.h"

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

void * arbint_alloc(void * ud, void * ptr, size_t new_size) {
  (void) ud;
  if (new_size == 0) {
    free(ptr);
    return NULL;
  }
  return realloc(ptr, new_size);
}

arbint_err_t arbint_ctx_init_default(arbint_ctx_t * ctx) {
  arbint_alloc_t a;
  if (ctx == NULL)
    return ARBINT_EINVAL;

  a.ud = NULL;
  a.realloc = arbint_alloc;
  return arbint_ctx_init(ctx, &a, 0u);
}

arbint_err_t arbint_ctx_init(arbint_ctx_t * ctx, const arbint_alloc_t * a,
                             uint32_t flags) {
  if (ctx == NULL || a == NULL || a->realloc == NULL)
    return ARBINT_EINVAL;

  ctx->a = *a;
  ctx->flags = flags;
  ctx->rng_flags = 0u;
  return ARBINT_OK;
}

void arbint_ctx_clear(arbint_ctx_t * ctx) { (void) ctx; }

arbint_err_t arbint_init(arbint_t x, arbint_ctx_t * ctx) {
  if (x == NULL)
    return ARBINT_EINVAL;

  x[0]._cap = 0u;
  x[0]._sz = 0;
  x[0]._ptr = NULL;
  x[0]._ctx = ctx;
  return ARBINT_OK;
}

arbint_err_t arbint_init_all(arbint_ctx_t * ctx, arbint_t a, ...) {
  va_list ap, ap_saved;
  arbint_t *p, *failed;
  arbint_err_t rc;

  rc = arbint_init(a, ctx);
  if (rc != ARBINT_OK)
    return rc;

  va_start(ap, a);
  va_copy(ap_saved, ap);
  while ((p = va_arg(ap, arbint_t *)) != NULL) {
    rc = arbint_init(*p, ctx);
    if (rc != ARBINT_OK) {
      failed = p;
      va_end(ap);
      /* Clear everything initialized so far. */
      arbint_clear(a);
      while ((p = va_arg(ap_saved, arbint_t *)) != failed)
        arbint_clear(*p);
      va_end(ap_saved);
      return rc;
    }
  }
  va_end(ap);
  va_end(ap_saved);
  return ARBINT_OK;
}

void arbint_clear(arbint_t x) {
  if (x == NULL)
    return;

  if (x[0]._ptr != NULL && x[0]._ctx != NULL && x[0]._ctx->a.realloc != NULL)
    (void) x[0]._ctx->a.realloc(x[0]._ctx->a.ud, x[0]._ptr, 0u);

  x[0]._cap = 0u;
  x[0]._sz = 0;
  x[0]._ptr = NULL;
  x[0]._ctx = NULL;
}

void arbint_clear_all(arbint_t a, ...) {
  va_list ap;
  arbint_t * p;

  arbint_clear(a);
  va_start(ap, a);
  while ((p = va_arg(ap, arbint_t *)) != NULL)
    arbint_clear(*p);
  va_end(ap);
}

arbint_ctx_t * arbint_get_ctx(const arbint_t x) {
  if (x == NULL)
    return NULL;
  return x[0]._ctx;
}

arbint_err_t arbint_set_ctx(arbint_t x, arbint_ctx_t * ctx) {
  arbint_ctx_t * old_ctx;

  if (x == NULL)
    return ARBINT_EINVAL;

  old_ctx = x[0]._ctx;
  if (old_ctx == ctx)
    return ARBINT_OK;

  if (x[0]._ptr == NULL || x[0]._cap == 0u) {
    x[0]._ctx = ctx;
    return ARBINT_OK;
  }

  if (ctx == NULL || ctx->a.realloc == NULL)
    return ARBINT_EINVAL;

  if (old_ctx == NULL || old_ctx->a.realloc == NULL)
    return ARBINT_EINVAL;

  if (x[0]._cap > SIZE_MAX / sizeof(arbint_limb_t))
    return ARBINT_EOVERFLOW;

  {
    size_t used = arbint_abs_sz(x[0]._sz);
    size_t old_nbytes = x[0]._cap * sizeof(arbint_limb_t);
    size_t used_nbytes = used * sizeof(arbint_limb_t);
    void * new_ptr;

    if (used > x[0]._cap)
      return ARBINT_EINVAL;

    new_ptr = ctx->a.realloc(ctx->a.ud, NULL, old_nbytes);
    if (new_ptr == NULL)
      return ARBINT_ENOMEM;

    if (used_nbytes != 0u)
      memcpy(new_ptr, x[0]._ptr, used_nbytes);

    (void) old_ctx->a.realloc(old_ctx->a.ud, x[0]._ptr, 0u);
    x[0]._ptr = new_ptr;
    x[0]._ctx = ctx;
  }

  return ARBINT_OK;
}

void arbint_zero(arbint_t x) {
  if (x == NULL)
    return;
  x[0]._sz = 0;
}

arbint_err_t arbint_resize(arbint_t x, size_t new_cap) {
  void * new_ptr;
  size_t used;

  if (x == NULL)
    return ARBINT_EINVAL;

  if (new_cap <= x[0]._cap)
    return ARBINT_OK;

  if (new_cap > SIZE_MAX / sizeof(arbint_limb_t))
    return ARBINT_EOVERFLOW;

  if (x[0]._ctx == NULL || x[0]._ctx->a.realloc == NULL)
    return ARBINT_EINVAL;

  used = arbint_abs_sz(x[0]._sz);
  if (used > x[0]._cap)
    return ARBINT_EINVAL;

  new_ptr = x[0]._ctx->a.realloc(x[0]._ctx->a.ud, x[0]._ptr,
                                 new_cap * sizeof(arbint_limb_t));
  if (new_ptr == NULL)
    return ARBINT_ENOMEM;

  x[0]._ptr = new_ptr;
  x[0]._cap = new_cap;
  return ARBINT_OK;
}

void arbint_swap(arbint_t a, arbint_t b) {
  _arbint_struct tmp;

  if (a == NULL || b == NULL || a == b)
    return;

  tmp = a[0];
  a[0] = b[0];
  b[0] = tmp;
}

int arbint_set_signed_sz(arbint_t x, size_t used, int sign) {
  if (x == NULL)
    return 0;
  if (used == 0u) {
    x[0]._sz = 0;
    return 1;
  }
  if (used > (size_t) PTRDIFF_MAX)
    return 0;
  x[0]._sz = (sign < 0) ? -(ptrdiff_t) used : (ptrdiff_t) used;
  return 1;
}

int arbint_cmp_mag_limbs(const arbint_limb_t * a, size_t an,
                         const arbint_limb_t * b, size_t bn) {
  size_t i;

  if (an < bn)
    return -1;
  if (an > bn)
    return 1;
  for (i = an; i != 0u; --i) {
    arbint_limb_t av = a[i - 1u];
    arbint_limb_t bv = b[i - 1u];
    if (av < bv)
      return -1;
    if (av > bv)
      return 1;
  }
  return 0;
}

size_t arbint_norm_used(const arbint_limb_t * x, size_t n) {
  while (n != 0u && x[n - 1u] == 0u)
    --n;
  return n;
}

arbint_limb_t * arbint_alloc_limbs(const arbint_alloc_t * alloc, size_t n) {
  if (alloc == NULL || alloc->realloc == NULL)
    return NULL;
  if (n == 0u)
    n = 1u;
  if (n > SIZE_MAX / sizeof(arbint_limb_t))
    return NULL;
  return (arbint_limb_t *) alloc->realloc(alloc->ud, NULL,
                                          n * sizeof(arbint_limb_t));
}

void arbint_free_limbs(const arbint_alloc_t * alloc, arbint_limb_t * p) {
  if (p == NULL || alloc == NULL || alloc->realloc == NULL)
    return;
  (void) alloc->realloc(alloc->ud, p, 0u);
}
