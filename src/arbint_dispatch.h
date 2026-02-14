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

#ifndef ARBINT_DISPATCH_H
#define ARBINT_DISPATCH_H

/*  CPU feature dispatch macros.

    These macros generate static selector functions for runtime dispatch
    based on CPU features. They implement a three-tier dispatch pattern:

    1. HAS_XXX_ALWAYS: feature available at compile time (e.g. -mbmi2),
       always use the optimized implementation.
    2. HAS_XXX: feature code compiled separately, check CPUID at runtime.
    3. Neither: only generic implementation available.

    Usage:
      ARBINT_DISPATCH_BMI2(arbint_select_foo, foo_fn_t, foo_bmi2, foo_generic)

    Generates:
      static foo_fn_t arbint_select_foo(void) { ... }  */

#include "arbint_cpu.h"
#include "config.h"

/*  Two-tier BMI2 dispatch: BMI2 > generic.  */

#if HAS_BMI2_ALWAYS

#define ARBINT_DISPATCH_BMI2(name, ret_type, fn_bmi2, fn_generic)              \
  static ret_type name(void) {                                                 \
    return fn_bmi2;                                                            \
  }

#elif HAS_BMI2

#define ARBINT_DISPATCH_BMI2(name, ret_type, fn_bmi2, fn_generic)              \
  static ret_type name(void) {                                                 \
    return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2) ? fn_bmi2           \
                                                          : fn_generic;        \
  }

#else

#define ARBINT_DISPATCH_BMI2(name, ret_type, fn_bmi2, fn_generic)              \
  static ret_type name(void) {                                                 \
    (void) fn_bmi2;                                                            \
    return fn_generic;                                                         \
  }

#endif

/*  Two-tier AVX2 dispatch: AVX2 > generic.  */

#if HAS_AVX2_ALWAYS

#define ARBINT_DISPATCH_AVX2(name, ret_type, fn_avx2, fn_generic)              \
  static ret_type name(void) {                                                 \
    return fn_avx2;                                                            \
  }

#elif HAS_AVX2

#define ARBINT_DISPATCH_AVX2(name, ret_type, fn_avx2, fn_generic)              \
  static ret_type name(void) {                                                 \
    return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_AVX2) ? fn_avx2           \
                                                          : fn_generic;        \
  }

#else

#define ARBINT_DISPATCH_AVX2(name, ret_type, fn_avx2, fn_generic)              \
  static ret_type name(void) {                                                 \
    (void) fn_avx2;                                                            \
    return fn_generic;                                                         \
  }

#endif

/*  Three-tier AVX2 > BMI2 > generic dispatch.  */

#if HAS_AVX2_ALWAYS

#define ARBINT_DISPATCH_AVX2_BMI2(name, ret_type, fn_avx2, fn_bmi2, fn_generic) \
  static ret_type name(void) {                                                 \
    (void) fn_bmi2;                                                            \
    (void) fn_generic;                                                         \
    return fn_avx2;                                                            \
  }

#elif HAS_AVX2

  #if HAS_BMI2_ALWAYS

  #define ARBINT_DISPATCH_AVX2_BMI2(name, ret_type, fn_avx2, fn_bmi2, fn_generic) \
    static ret_type name(void) {                                               \
      (void) fn_generic;                                                       \
      return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_AVX2) ? fn_avx2         \
                                                            : fn_bmi2;         \
    }

  #elif HAS_BMI2

  #define ARBINT_DISPATCH_AVX2_BMI2(name, ret_type, fn_avx2, fn_bmi2, fn_generic) \
    static ret_type name(void) {                                               \
      if (arbint_cpu_has_feature(ARBINT_CPU_FEATURE_AVX2))                     \
        return fn_avx2;                                                        \
      return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2) ? fn_bmi2         \
                                                            : fn_generic;      \
    }

  #else

  #define ARBINT_DISPATCH_AVX2_BMI2(name, ret_type, fn_avx2, fn_bmi2, fn_generic) \
    static ret_type name(void) {                                               \
      (void) fn_bmi2;                                                          \
      return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_AVX2) ? fn_avx2         \
                                                            : fn_generic;      \
    }

  #endif

#else

  /*  No AVX2, fall back to BMI2 dispatch.  */
  #define ARBINT_DISPATCH_AVX2_BMI2(name, ret_type, fn_avx2, fn_bmi2, fn_generic) \
    ARBINT_DISPATCH_BMI2(name, ret_type, fn_bmi2, fn_generic)

#endif

/*  Two-tier POPCNT dispatch: POPCNT > generic.  */

#if HAS_POPCNT_ALWAYS

#define ARBINT_DISPATCH_POPCNT(name, ret_type, fn_popcnt, fn_generic)          \
  static ret_type name(void) {                                                 \
    return fn_popcnt;                                                          \
  }

#elif HAS_POPCNT

#define ARBINT_DISPATCH_POPCNT(name, ret_type, fn_popcnt, fn_generic)          \
  static ret_type name(void) {                                                 \
    return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_POPCNT) ? fn_popcnt       \
                                                            : fn_generic;      \
  }

#else

#define ARBINT_DISPATCH_POPCNT(name, ret_type, fn_popcnt, fn_generic)          \
  static ret_type name(void) {                                                 \
    (void) fn_popcnt;                                                          \
    return fn_generic;                                                         \
  }

#endif

/*  Two-tier SHA-NI dispatch: SHA-NI > generic.  */

#if HAS_SHA_NI_ALWAYS

#define ARBINT_DISPATCH_SHA_NI(name, ret_type, fn_shani, fn_generic)           \
  static ret_type name(void) {                                                 \
    return fn_shani;                                                           \
  }

#elif HAS_SHA_NI

#define ARBINT_DISPATCH_SHA_NI(name, ret_type, fn_shani, fn_generic)           \
  static ret_type name(void) {                                                 \
    return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_SHA) ? fn_shani           \
                                                         : fn_generic;         \
  }

#else

#define ARBINT_DISPATCH_SHA_NI(name, ret_type, fn_shani, fn_generic)           \
  static ret_type name(void) {                                                 \
    (void) fn_shani;                                                           \
    return fn_generic;                                                         \
  }

#endif

/*  Two-tier SSE4.2 CRC32 dispatch: SSE4.2 > generic.  */

#if HAS_SSE42_CRC32_ALWAYS

#define ARBINT_DISPATCH_SSE42_CRC32(name, ret_type, fn_sse42, fn_generic)      \
  static ret_type name(void) {                                                 \
    return fn_sse42;                                                           \
  }

#elif HAS_SSE42_CRC32

#define ARBINT_DISPATCH_SSE42_CRC32(name, ret_type, fn_sse42, fn_generic)      \
  static ret_type name(void) {                                                 \
    return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_SSE42) ? fn_sse42         \
                                                           : fn_generic;       \
  }

#else

#define ARBINT_DISPATCH_SSE42_CRC32(name, ret_type, fn_sse42, fn_generic)      \
  static ret_type name(void) {                                                 \
    (void) fn_sse42;                                                           \
    return fn_generic;                                                         \
  }

#endif

#endif /*  ARBINT_DISPATCH_H  */
