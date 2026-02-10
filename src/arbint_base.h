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
#include <string.h>

#if ARBINT_LIMB_BITS == 64
typedef uint64_t arbint_limb_t;
#elif ARBINT_LIMB_BITS == 32
typedef uint32_t arbint_limb_t;
#else
  #error "Unsupported ARBINT_LIMB_BITS"
#endif /* ARBINT_LIMB_BITS */

#define ARBINT_HALF_BITS (ARBINT_LIMB_BITS / 2u)
#define ARBINT_HALF_MASK ((((arbint_limb_t) 1u) << ARBINT_HALF_BITS) - 1u)

#if defined(__i386__) || defined(__x86_64__) || defined(_M_IX86) ||           \
    defined(_M_X64)
  #define ARBINT_TARGET_X86_FAMILY 1
#else
  #define ARBINT_TARGET_X86_FAMILY 0
#endif /* x86 family detection */

#if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)
  #define ARBINT_TARGET_AARCH64 1
#else
  #define ARBINT_TARGET_AARCH64 0
#endif /* aarch64 detection */

#if defined(__GNUC__) || defined(__clang__)
  #define ARBINT_COMPILER_GNU_CLANG 1
#else
  #define ARBINT_COMPILER_GNU_CLANG 0
#endif /* defined(__GNUC__) || defined(__clang__) */

#if defined(_MSC_VER)
  #define ARBINT_COMPILER_MSVC 1
#else
  #define ARBINT_COMPILER_MSVC 0
#endif /* defined(_MSC_VER) */

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

/*  Internal function to count the number of leading zeros in a limb.
    Precondition: x != 0 (behavior is undefined for zero input).
    On GCC/Clang, __builtin_clz(0) is undefined. On MSVC, _BitScanReverse(0)
    writes an unspecified value to the output parameter. The fallback loop
    returns LIMB_BITS for zero, which is technically correct but callers
    should not rely on this. Callers MUST check for zero before calling.  */
static inline unsigned arbint_clz_limb(arbint_limb_t x) {
#if ARBINT_COMPILER_GNU_CLANG
  #if ARBINT_LIMB_BITS == 64
  return (unsigned) __builtin_clzll((unsigned long long) x);
  #else
  return (unsigned) __builtin_clz((unsigned) x);
  #endif /* ARBINT_LIMB_BITS */
#elif ARBINT_COMPILER_MSVC
  unsigned long idx;
  #if ARBINT_LIMB_BITS == 64
  _BitScanReverse64(&idx, x);
  #else
  _BitScanReverse(&idx, x);
  #endif /* ARBINT_LIMB_BITS */
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
#endif /* ARBINT_COMPILER_GNU_CLANG */
}

/*  Count trailing zeros in a single limb.
    Precondition: x != 0 (behavior is undefined for zero input).
    On GCC/Clang, __builtin_ctz(0) is undefined. On MSVC, _BitScanForward(0)
    writes an unspecified value to the output parameter. The fallback loop
    would run indefinitely. Callers MUST check for zero before calling.  */
static inline unsigned arbint_ctz_limb(arbint_limb_t x) {
#if ARBINT_COMPILER_GNU_CLANG
  #if ARBINT_LIMB_BITS == 64
  return (unsigned) __builtin_ctzll((unsigned long long) x);
  #else
  return (unsigned) __builtin_ctz((unsigned) x);
  #endif /* ARBINT_LIMB_BITS */
#elif ARBINT_COMPILER_MSVC
  {
    unsigned long idx;
  #if ARBINT_LIMB_BITS == 64
    _BitScanForward64(&idx, x);
  #else
    _BitScanForward(&idx, x);
  #endif /* ARBINT_LIMB_BITS */
    return (unsigned) idx;
  }
#else
  {
    unsigned n = 0u;
    while ((x & 1u) == 0u) {
      ++n;
      x >>= 1u;
    }
    return n;
  }
#endif /* ARBINT_COMPILER_GNU_CLANG */
}

/*  Portable popcount of a single limb.  */
static inline unsigned arbint_popcount_limb(arbint_limb_t x) {
#if ARBINT_COMPILER_GNU_CLANG
  #if ARBINT_LIMB_BITS == 64
  return (unsigned) __builtin_popcountll((unsigned long long) x);
  #else
  return (unsigned) __builtin_popcount((unsigned) x);
  #endif /* ARBINT_LIMB_BITS */
#elif ARBINT_COMPILER_MSVC
  #if ARBINT_LIMB_BITS == 64
  return (unsigned) __popcnt64(x);
  #else
  return (unsigned) __popcnt(x);
  #endif /* ARBINT_LIMB_BITS */
#else
  {
    unsigned count = 0u;
    while (x != 0u) {
      x &= x - 1u;
      ++count;
    }
    return count;
  }
#endif /* ARBINT_COMPILER_GNU_CLANG */
}

/*  Securely zero memory so the compiler cannot optimise the store away.
    Defined in arbint_base.c so that platform headers (e.g. <windows.h>)
    stay out of this header.  */
void arbint_secure_zero(void * ptr, size_t len);

#endif /* ARBINT_BASE_H */
