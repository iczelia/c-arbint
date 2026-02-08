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

#ifndef ARBINT_BASE_H
#define ARBINT_BASE_H

#include "arbint.h"
#include "config.h"

#include <stddef.h>
#include <stdint.h>

#if ARBINT_LIMB_BITS == 64
typedef uint64_t arbint_limb_t;
#elif ARBINT_LIMB_BITS == 32
typedef uint32_t arbint_limb_t;
#else
  #error "Unsupported ARBINT_LIMB_BITS"
#endif

#define ARBINT_HALF_BITS (ARBINT_LIMB_BITS / 2u)
#define ARBINT_HALF_MASK ((((arbint_limb_t) 1u) << ARBINT_HALF_BITS) - 1u)

#if defined(__i386__) || defined(__x86_64__) || defined(_M_IX86) ||           \
    defined(_M_X64)
  #define ARBINT_TARGET_X86_FAMILY 1
#else
  #define ARBINT_TARGET_X86_FAMILY 0
#endif

#if defined(__GNUC__) || defined(__clang__)
  #define ARBINT_COMPILER_GNU_CLANG 1
#else
  #define ARBINT_COMPILER_GNU_CLANG 0
#endif

#if defined(_MSC_VER)
  #define ARBINT_COMPILER_MSVC 1
#else
  #define ARBINT_COMPILER_MSVC 0
#endif

#define ARBINT_LIMBS(x) ((arbint_limb_t *) ((x)[0]._ptr))
#define ARBINT_CLIMBS(x) ((const arbint_limb_t *) ((x)[0]._ptr))

/*  libc-backed default allocator with realloc semantics required by
    arbint_alloc_t.  */
void * arbint_alloc(void * ud, void * ptr, size_t new_size);

/* Internal limb helpers shared by add/sub and mul implementations. */
int arbint_set_signed_sz(arbint_t x, size_t used, int sign);
int arbint_cmp_mag_limbs(const arbint_limb_t * a, size_t an,
                         const arbint_limb_t * b, size_t bn);
size_t arbint_norm_used(const arbint_limb_t * x, size_t n);

/*  Internal limb allocation/deallocation functions, primarily for
    temporaries.  */
arbint_limb_t * arbint_alloc_limbs(const arbint_alloc_t * alloc, size_t n);
void arbint_free_limbs(const arbint_alloc_t * alloc, arbint_limb_t * p);

static inline size_t arbint_abs_sz(ptrdiff_t sz) {
  if (sz >= 0)
    return (size_t) sz;
  return (size_t) (-(sz + 1)) + 1;
}

/*  Internal function to count the number of leading zeros in a limb.  */
static inline unsigned arbint_clz_limb(arbint_limb_t x) {
#if ARBINT_COMPILER_GNU_CLANG
  #if ARBINT_LIMB_BITS == 64
  return (unsigned) __builtin_clzll((unsigned long long) x);
  #else
  return (unsigned) __builtin_clz((unsigned) x);
  #endif
#elif ARBINT_COMPILER_MSVC
  unsigned long idx;
  #if ARBINT_LIMB_BITS == 64
  _BitScanReverse64(&idx, x);
  #else
  _BitScanReverse(&idx, x);
  #endif
  return (unsigned) (ARBINT_LIMB_BITS - 1u - idx);
#else
  {
    unsigned n = 0;
    arbint_limb_t mask = (arbint_limb_t) 1u << (ARBINT_LIMB_BITS - 1u);
    while (mask != 0u && (x & mask) == 0u) {
      ++n;
      mask >>= 1u;
    }
    return n;
  }
#endif
}

#endif
