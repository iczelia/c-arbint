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

#if ARBINT_COMPILER_MSVC && HAVE_WINDOWS_H
  #include <windows.h>
#endif

/*  Default allocator using libc malloc/realloc/free.
    Follows realloc semantics: NULL ptr = alloc, new_size = 0 = free.  */
void * arbint_alloc(void * ud, void * ptr, size_t new_size) {
  (void) ud;
  if (new_size == 0) {
    free(ptr);
    return NULL;
  }
  return realloc(ptr, new_size);
}

/*  Initialize context with default allocator (libc malloc/realloc/free).  */
arbint_err_t arbint_ctx_init_default(arbint_ctx_t * ctx) {
  arbint_alloc_t a;
  if (ctx == NULL)
    return ARBINT_EINVAL;

  a.ud = NULL;
  a.realloc = arbint_alloc;
  return arbint_ctx_init(ctx, &a, 0u);
}

/*  Initialize context with custom allocator and flags.  */
arbint_err_t arbint_ctx_init(arbint_ctx_t * ctx, const arbint_alloc_t * a,
                             uint32_t flags) {
  if (ctx == NULL || a == NULL || a->realloc == NULL)
    return ARBINT_EINVAL;

  ctx->a = *a;
  ctx->flags = flags;
  ctx->rng_flags = 0u;
  return ARBINT_OK;
}

/*  Clear context.  Currently a no-op since contexts hold no dynamic state,
    but callers should always invoke this for forward compatibility: future
    versions may cache RNG state or allocator pools that require cleanup.  */
void arbint_ctx_clear(arbint_ctx_t * ctx) { (void) ctx; }

/*  Initialize arbint to zero with given context.
    No memory allocated until first operation requiring capacity.  */
arbint_err_t arbint_init(arbint_t x, arbint_ctx_t * ctx) {
  if (x == NULL)
    return ARBINT_EINVAL;

  x[0]._cap = 0u;
  x[0]._sz = 0;
  x[0]._ptr = NULL;
  x[0]._ctx = ctx;
  return ARBINT_OK;
}

/*  Initialize multiple arbints with given context (variadic, NULL-terminated).
    On error, clears all successfully initialized arbints and returns error
    code.  */
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

/*  Free arbint memory and reset to uninitialized state.
    Safe to call multiple times or on NULL.  */
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

/*  Clear multiple arbints (variadic, NULL-terminated).  */
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

/*  Set arbint to zero without freeing memory.  */
void arbint_zero(arbint_t x) {
  if (x == NULL)
    return;
  x[0]._sz = 0;
}

/*  Grow arbint capacity to at least new_cap limbs.
    No-op if new_cap <= current capacity. Does not shrink.  */
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

/*  Swap two arbints by exchanging their internal structures.  */
void arbint_swap(arbint_t a, arbint_t b) {
  _arbint_struct tmp;

  if (a == NULL || b == NULL || a == b)
    return;

  tmp = a[0];
  a[0] = b[0];
  b[0] = tmp;
}

/*  Set arbint size field from magnitude and sign.

    Sets x->_sz to encode both the magnitude (number of used limbs) and sign
    of the integer. The sign is encoded in the sign of _sz itself:
      _sz > 0: positive number with _sz limbs
      _sz < 0: negative number with abs(_sz) limbs
      _sz = 0: exactly zero (canonical representation)

    This function enforces the normalization invariant that _sz == 0 represents
    zero (not negative zero) and abs(_sz) equals the actual used limb count.

    Parameters:
      x    - Target arbint (must be non-NULL)
      used - Number of significant limbs (magnitude)
      sign - Sign indicator: <0 for negative, 0 for zero, >0 for positive

    Returns:
      1 on success, 0 on overflow (used > PTRDIFF_MAX)

    Preconditions:
      - x != NULL
      - used <= PTRDIFF_MAX (enforced by check below)
      - If used == 0, the final _sz will be 0 regardless of sign (enforced)

    Overflow handling: Since _sz is ptrdiff_t (signed), the maximum
    representable magnitude is PTRDIFF_MAX (typically 2^62 - 1 on 64-bit). If
    used exceeds this, returns 0 to signal overflow. On 64-bit systems with
    64-bit limbs, this limits numbers to ~2^(62 * 64) = 2^3968 bits, which is
    acceptable for practical use.

    Sign handling: For negative numbers, computes -(ptrdiff_t)used safely
    without overflow since used <= PTRDIFF_MAX guarantees the negation fits in
    ptrdiff_t. The special case used==0 forces _sz=0 regardless of sign to
    maintain the invariant that zero has no sign.  */
int arbint_set_signed_sz(arbint_t x, size_t used, int sign) {
  if (x == NULL)
    return 0;
  if (used == 0u) {
    x[0]._sz = 0; /* Canonical zero: always _sz == 0, not negative zero */
    return 1;
  }
  if (used > (size_t) PTRDIFF_MAX)
    return 0; /* Magnitude too large to represent in ptrdiff_t */
  x[0]._sz = (sign < 0) ? -(ptrdiff_t) used : (ptrdiff_t) used;
  return 1;
}

/*  Compare magnitudes of two limb arrays.
    Returns -1 if |a| < |b|, 0 if equal, +1 if |a| > |b|.
    Compares from most significant limb downward for early exit.  */
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

/*  Trim leading zero limbs and return actual used count.
    Essential for maintaining normalization invariant (no leading zeros).  */
size_t arbint_norm_used(const arbint_limb_t * x, size_t n) {
  while (n != 0u && x[n - 1u] == 0u)
    --n;
  return n;
}

/*  Allocate temporary limb array for intermediate calculations.

    Allocates a temporary buffer of n limbs using the provided allocator.
    Used throughout the library for scratch space in algorithms like Karatsuba
    multiplication, division, and other multi-limb operations.

    Zero-size allocation handling: Requesting 0 limbs is treated as requesting
    1 limb. Rationale: Some allocators (including some libc implementations)
    return NULL for zero-size allocations, which creates ambiguity (NULL could
    mean "allocation failed" or "empty allocation succeeded"). Requesting at
    least 1 byte ensures non-NULL return on success, making error detection
    straightforward via NULL check.

    Parameters:
      alloc - Allocator to use (must be non-NULL with valid realloc function)
      n     - Number of limbs to allocate (0 treated as 1)

    Returns:
      Pointer to allocated limb array on success, NULL on failure (out of
      memory or overflow in size calculation).

    Preconditions:
      - alloc != NULL && alloc->realloc != NULL
      - n <= SIZE_MAX / sizeof(arbint_limb_t) (checked internally)

    Caller responsibility: Must free returned pointer via arbint_free_limbs
    using the same allocator. Failure to free causes memory leak.

    Overflow protection: Checks n * sizeof(arbint_limb_t) <= SIZE_MAX before
    allocating to prevent integer overflow in size calculation.  */
arbint_limb_t * arbint_alloc_limbs(const arbint_alloc_t * alloc, size_t n) {
  if (alloc == NULL || alloc->realloc == NULL)
    return NULL;
  if (n == 0u)
    n = 1u; /* Avoid zero-size allocation ambiguity (see comment above) */
  if (n > SIZE_MAX / sizeof(arbint_limb_t))
    return NULL;
  return (arbint_limb_t *) alloc->realloc(alloc->ud, NULL,
                                          n * sizeof(arbint_limb_t));
}

/*  Free temporary limb array allocated by arbint_alloc_limbs.

    Releases memory allocated by arbint_alloc_limbs. Safe to call with NULL
    pointer (no-op). Must use the same allocator that was used for allocation.

    Parameters:
      alloc - Same allocator used in arbint_alloc_limbs (must be non-NULL)
      p     - Pointer to limb array, or NULL

    Precondition: If p is non-NULL, it must have been allocated via
    arbint_alloc_limbs using the same allocator. Passing a pointer from a
    different source or different allocator is undefined behavior.

    Implementation note: Uses realloc(ptr, 0) to free, per the allocator
    contract (see arbint_realloc_fn typedef in arbint.h). The cast to (void)
    discards the return value since realloc(ptr, 0) may return NULL or an
    opaque marker, neither of which is meaningful for freeing.  */
void arbint_free_limbs(const arbint_alloc_t * alloc, arbint_limb_t * p) {
  if (p == NULL || alloc == NULL || alloc->realloc == NULL)
    return;
  (void) alloc->realloc(alloc->ud, p, 0u);
}

/*  Securely zero memory so the compiler cannot optimise the store away.
    Prefers explicit_bzero (glibc 2.25+, most BSDs) or memset_s (C11 Annex K),
    then SecureZeroMemory on Windows, falling back to a
    volatile-function-pointer indirection that defeats dead-store elimination
    on all known compilers.  */
void arbint_secure_zero(void * ptr, size_t len) {
#if HAVE_EXPLICIT_BZERO
  explicit_bzero(ptr, len);
#elif HAVE_MEMSET_S
  (void) memset_s(ptr, len, 0, len);
#elif ARBINT_COMPILER_MSVC && HAVE_WINDOWS_H
  SecureZeroMemory(ptr, len);
#else
  static void * (*const volatile memset_v)(void *, int, size_t) = &memset;
  (void) memset_v(ptr, 0, len);
#endif
}
