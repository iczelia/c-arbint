AC_PREREQ([2.69])
AC_INIT([c-arbint], CA_GIT_DESCRIBE_VERSION, [k@iczelia.net], [c-arbint], [https://github.com/iczelia/c-arbint])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
m4_include([m4/arbint_x86_feature_check.m4])
m4_include([m4/arbint_arm_feature_check.m4])
AM_INIT_AUTOMAKE([foreign subdir-objects dist-xz no-dist-gzip tar-ustar])
AM_SILENT_RULES([yes])
AC_CANONICAL_HOST

AS_CASE([$host_cpu],
  [x86_64|amd64], [arbint_host_x86=yes; arbint_host_x86_64=yes; arbint_host_has_cpuid=yes; arbint_host_aarch64=no],
  [i[5-9]86|pentium*], [arbint_host_x86=yes; arbint_host_x86_64=no; arbint_host_has_cpuid=yes; arbint_host_aarch64=no],
  [i?86], [arbint_host_x86=yes; arbint_host_x86_64=no; arbint_host_has_cpuid=no; arbint_host_aarch64=no],
  [aarch64*|arm64*], [arbint_host_x86=no; arbint_host_x86_64=no; arbint_host_has_cpuid=no; arbint_host_aarch64=yes],
  [arbint_host_x86=no; arbint_host_x86_64=no; arbint_host_has_cpuid=no; arbint_host_aarch64=no])

AC_USE_SYSTEM_EXTENSIONS
AC_PROG_CC
AM_PROG_AR

dnl Save the user's explicit --enable-shared / --disable-shared choice
dnl before LT_INIT, which may override it based on platform detection.
arbint_user_enable_shared="$enable_shared"

LT_INIT

dnl Emscripten: libtool's ld_shlibs detection may fail because the default
dnl case checks for ELF targets and emcc does not always advertise ELF in
dnl its --help output.  When ld_shlibs=no, libtool forces can_build_shared
dnl and enable_shared to "no", disabling shared libraries entirely.
dnl Override here unless the user explicitly passed --disable-shared.
AS_CASE([$host_os],
  [emscripten*], [
    AS_IF([test "x$enable_shared" = "xno" && test "x$arbint_user_enable_shared" != "xno"], [
      enable_shared=yes
      can_build_shared=yes
      AC_MSG_NOTICE([emscripten detected: forcing enable_shared=yes (libtool platform detection failed)])
    ])
  ])

AC_C_INLINE
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_CHECK_SIZEOF([void *])

AS_CASE([$ac_cv_sizeof_void_p],
  [8], [arbint_limb_bits=64],
  [4], [arbint_limb_bits=32],
  [AC_MSG_ERROR([unsupported pointer size $ac_cv_sizeof_void_p: expected 4 or 8])])
AC_DEFINE_UNQUOTED([ARBINT_LIMB_BITS], [$arbint_limb_bits],
  [Number of bits in one arbint limb (autodetected by configure).])

MAN_DATE=`date "+%B %Y"`
AC_SUBST([MAN_DATE])
AC_CHECK_PROGS([NROFF], [nroff], [])
AM_CONDITIONAL([HAVE_NROFF], [test -n "$NROFF"])

AC_ARG_ENABLE([warnings],
  [AS_HELP_STRING([--enable-warnings], [Enable extra compiler warnings])],
  [],
  [enable_warnings=yes])

AC_CACHE_CHECK([whether $CC is Microsoft cl (not clang-cl)], [arbint_cv_cc_msvc_cl], [
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#if !defined(_MSC_VER) || defined(__clang__)
#error not msvc cl
#endif
  ]], [[]])],
    [arbint_cv_cc_msvc_cl=yes],
    [arbint_cv_cc_msvc_cl=no])
])
arbint_cc_msvc_cl="$arbint_cv_cc_msvc_cl"

AS_IF([test "x$arbint_cc_msvc_cl" = "xyes"], [
  AC_CACHE_CHECK([whether $CC supports /FS], [arbint_cv_cc_fs], [
    arbint_saved_cflags="$CFLAGS"
    CFLAGS="$CFLAGS /FS"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[]])],
      [arbint_cv_cc_fs=yes],
      [arbint_cv_cc_fs=no])
    CFLAGS="$arbint_saved_cflags"
  ])
  AS_IF([test "x$arbint_cv_cc_fs" = "xyes"], [
    CFLAGS="$CFLAGS /FS"
  ])
])

AS_IF([test "x$enable_warnings" = "xyes"], [
  AS_IF([test "x$arbint_cc_msvc_cl" = "xyes"], [
    CFLAGS="$CFLAGS /W4"
  ], [
    arbint_warning_cflags=""
    for arbint_warning_flag in -Wall -Wextra -pedantic -Werror=implicit-function-declaration; do
      AC_MSG_CHECKING([whether $CC supports $arbint_warning_flag])
      arbint_saved_cflags="$CFLAGS"
      CFLAGS="$CFLAGS $arbint_warning_flag"
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[]])],
        [arbint_flag_supported=yes],
        [arbint_flag_supported=no])
      CFLAGS="$arbint_saved_cflags"
      AC_MSG_RESULT([$arbint_flag_supported])
      AS_IF([test "x$arbint_flag_supported" = "xyes"], [
        arbint_warning_cflags="$arbint_warning_cflags $arbint_warning_flag"
      ])
    done
    CFLAGS="$CFLAGS$arbint_warning_cflags"
  ])
])

AC_ARG_WITH([base-optimization],
  [AS_HELP_STRING([--with-base-optimization=@<:@default|tuned@:>@],
                  [Select base optimization profile:
                   default => -O3 -g3 (default),
                   tuned => -O3 -march=... -mtune=...])],
  [arbint_base_optimization="$withval"],
  [arbint_base_optimization=default])

AC_ARG_WITH([march],
  [AS_HELP_STRING([--with-march=ARCH],
                  [Architecture/tune target to use when base optimization is tuned
                   (default: native). GCC/Clang: passed to -march/-mtune.
                   MSVC: passed to /arch: (e.g. AVX2, AVX512)])],
  [arbint_cfg_march="$withval"],
  [arbint_cfg_march=""])

AS_CASE([$arbint_base_optimization],
  [default|tuned], [],
  [AC_MSG_ERROR([invalid --with-base-optimization value: $arbint_base_optimization (expected default or tuned)])])

arbint_effective_march=""
arbint_effective_mtune=""
AS_IF([test "x$arbint_base_optimization" = "xdefault"], [
  AS_IF([test "x$arbint_cc_msvc_cl" = "xyes"], [
    CFLAGS="$CFLAGS /O2 /Zi"
  ], [
    CFLAGS="$CFLAGS -O3 -g3"
  ])
], [
  AS_IF([test "x$arbint_cc_msvc_cl" = "xyes"], [
    CFLAGS="$CFLAGS /O2 /Zi"
    AS_IF([test "x$arbint_cfg_march" != "x"], [
      CFLAGS="$CFLAGS /arch:$arbint_cfg_march"
    ])
  ], [
    arbint_effective_march="$arbint_cfg_march"
    AS_IF([test "x$arbint_effective_march" = "x"], [
      arbint_effective_march="native"
    ])
    CFLAGS="$CFLAGS -O3 -march=$arbint_effective_march -mtune=$arbint_effective_march"
  ])
])

AS_IF([test "x$arbint_cc_msvc_cl" != "xyes"], [
  AC_CACHE_CHECK([whether $CC supports -fno-semantic-interposition],
    [arbint_cv_cc_fno_semantic_interposition], [
    arbint_saved_cflags="$CFLAGS"
    CFLAGS="$CFLAGS -fno-semantic-interposition"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[]])],
      [arbint_cv_cc_fno_semantic_interposition=yes],
      [arbint_cv_cc_fno_semantic_interposition=no])
    CFLAGS="$arbint_saved_cflags"
  ])
  AS_IF([test "x$arbint_cv_cc_fno_semantic_interposition" = "xyes"], [
    CFLAGS="$CFLAGS -fno-semantic-interposition"
  ])
])

AC_ARG_ENABLE([lto],
  [AS_HELP_STRING([--enable-lto@<:@=auto|yes|no@:>@],
                  [Enable link-time optimization with -flto
                   (default: auto; auto enables LTO when supported)])],
  [arbint_enable_lto="$enableval"],
  [arbint_enable_lto=auto])

AS_CASE([$arbint_enable_lto],
  [yes|no|auto], [],
  [AC_MSG_ERROR([invalid --enable-lto value: $arbint_enable_lto (expected auto, yes, or no)])])

ARBINT_LTO_CFLAGS=""
ARBINT_LTO_LDFLAGS=""
arbint_enable_lto_effective=no
AS_IF([test "x$arbint_enable_lto" != "xno"], [
  AS_IF([test "x$arbint_cc_msvc_cl" = "xyes"], [
    AC_CACHE_CHECK([whether $CC supports LTO with /GL and /LTCG], [arbint_cv_cc_lto], [
      arbint_saved_cflags="$CFLAGS"
      arbint_saved_ldflags="$LDFLAGS"
      CFLAGS="$CFLAGS /GL"
      LDFLAGS="$LDFLAGS /LTCG"
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[return 0;]])],
        [arbint_cv_cc_lto=yes],
        [arbint_cv_cc_lto=no])
      CFLAGS="$arbint_saved_cflags"
      LDFLAGS="$arbint_saved_ldflags"
    ])
  ], [
    AC_CACHE_CHECK([whether $CC supports LTO with -flto], [arbint_cv_cc_lto], [
      arbint_saved_cflags="$CFLAGS"
      arbint_saved_ldflags="$LDFLAGS"
      CFLAGS="$CFLAGS -flto"
      LDFLAGS="$LDFLAGS -flto"
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[return 0;]])],
        [arbint_cv_cc_lto=yes],
        [arbint_cv_cc_lto=no])
      CFLAGS="$arbint_saved_cflags"
      LDFLAGS="$arbint_saved_ldflags"
    ])
  ])

  AS_IF([test "x$arbint_cv_cc_lto" = "xyes"], [
    arbint_enable_lto_effective=yes
    AS_IF([test "x$arbint_cc_msvc_cl" = "xyes"], [
      ARBINT_LTO_CFLAGS="/GL"
      ARBINT_LTO_LDFLAGS="/LTCG"
    ], [
      ARBINT_LTO_CFLAGS="-flto"
      ARBINT_LTO_LDFLAGS="-flto"
    ])
  ], [
    AS_IF([test "x$arbint_enable_lto" = "xyes"], [
      AS_IF([test "x$arbint_cc_msvc_cl" = "xyes"], [
        AC_MSG_ERROR([--enable-lto was requested, but $CC/$LD does not support /GL and /LTCG])
      ], [
        AC_MSG_ERROR([--enable-lto was requested, but $CC/$LD does not support -flto])
      ])
    ])
  ])
])
AC_SUBST([ARBINT_LTO_CFLAGS])
AC_SUBST([ARBINT_LTO_LDFLAGS])
AC_MSG_NOTICE([LTO requested=$arbint_enable_lto, enabled=$arbint_enable_lto_effective])

AC_ARG_ENABLE([assertions],
  [AS_HELP_STRING([--enable-assertions],
                  [Keep assert() checks enabled (do not define NDEBUG).
                   Default is disabled: NDEBUG is defined.])],
  [],
  [enable_assertions=no])

AS_IF([test "x$enable_assertions" != "xyes"], [
  AC_DEFINE([NDEBUG], [1], [Define to 1 to disable assert() checks.])
])

AC_CHECK_HEADERS([stdint.h stdlib.h string.h])

AS_IF([test "x$arbint_host_has_cpuid" = "xyes"], [
  arbint_has_cpuid=1
], [
  arbint_has_cpuid=0
])
AC_DEFINE_UNQUOTED([HAS_CPUID], [$arbint_has_cpuid],
  [Define to 1 if x86 CPUID is assumed available on the configured target.])

dnl BMI2 is only enabled on 64-bit platforms because our BMI2-optimized code
dnl uses _mulx_u64 which is not available on 32-bit targets.
AS_IF([test "$ac_cv_sizeof_void_p" = "8"], [
  ARBINT_CHECK_X86_FEATURE_WITH_FLAG(
    [BMI2 intrinsics], [arbint_cv_have_bmi2_plain], [arbint_cv_have_bmi2_mbmi2_v2],
    [-mbmi2],
    [[
#include <immintrin.h>
    ]],
    [[
unsigned long long hi = 0;
unsigned long long lo = _mulx_u64(3ULL, 5ULL, &hi);
return (int)(lo + hi);
    ]],
    [HAS_BMI2], [HAS_BMI2_ALWAYS],
    [Define to 1 if BMI2 intrinsics are available (with or without -mbmi2).],
    [Define to 1 if BMI2 intrinsics are available globally.]
  )

  arbint_build_bmi2=no
  ARBINT_BMI2_CFLAGS=""
  AS_IF([test "x$arbint_cv_have_bmi2_plain" = "xyes"], [
    arbint_build_bmi2=yes
  ], [test "x$arbint_cv_have_bmi2_mbmi2_v2" = "xyes"], [
    arbint_build_bmi2=yes
    ARBINT_BMI2_CFLAGS="-mbmi2"
  ])
], [
  dnl 32-bit: disable BMI2 (our code uses _mulx_u64)
  arbint_build_bmi2=no
  ARBINT_BMI2_CFLAGS=""
])
AM_CONDITIONAL([ARBINT_BUILD_BMI2], [test "x$arbint_build_bmi2" = "xyes"])
AC_SUBST([ARBINT_BMI2_CFLAGS])

ARBINT_CHECK_X86_FEATURE_WITH_FLAG(
  [POPCNT intrinsics], [arbint_cv_have_popcnt_plain], [arbint_cv_have_popcnt_mpopcnt],
  [-mpopcnt],
  [[
#include <nmmintrin.h>
  ]],
  [[
long long x = _mm_popcnt_u64(0ULL);
(void)x;
  ]],
  [HAS_POPCNT], [HAS_POPCNT_ALWAYS],
  [Define to 1 if POPCNT intrinsics are available (with or without -mpopcnt).],
  [Define to 1 if POPCNT intrinsics are available globally.]
)

arbint_build_popcnt=no
ARBINT_POPCNT_CFLAGS=""
AS_IF([test "x$arbint_cv_have_popcnt_plain" = "xyes"], [
  arbint_build_popcnt=yes
], [test "x$arbint_cv_have_popcnt_mpopcnt" = "xyes"], [
  arbint_build_popcnt=yes
  ARBINT_POPCNT_CFLAGS="-mpopcnt"
])
AM_CONDITIONAL([ARBINT_BUILD_POPCNT], [test "x$arbint_build_popcnt" = "xyes"])
AC_SUBST([ARBINT_POPCNT_CFLAGS])

dnl AVX2 is only enabled on 64-bit platforms (32-bit has limited AVX2 support and may cause issues)
AS_IF([test "$ac_cv_sizeof_void_p" = "8"], [
  ARBINT_CHECK_X86_FEATURE_WITH_FLAG(
    [AVX2 intrinsics], [arbint_cv_have_avx2_plain], [arbint_cv_have_avx2_mavx2],
    [-mavx2],
    [[
#include <immintrin.h>
    ]],
    [[
__m256i a = _mm256_setzero_si256();
__m256i b = _mm256_add_epi64(a, a);
(void)b;
    ]],
    [HAS_AVX2], [HAS_AVX2_ALWAYS],
    [Define to 1 if AVX2 intrinsics are available (with or without -mavx2).],
    [Define to 1 if AVX2 intrinsics are available globally.]
  )

  arbint_build_avx2=no
  ARBINT_AVX2_CFLAGS=""
  AS_IF([test "x$arbint_cv_have_avx2_plain" = "xyes"], [
    arbint_build_avx2=yes
  ], [test "x$arbint_cv_have_avx2_mavx2" = "xyes"], [
    arbint_build_avx2=yes
    ARBINT_AVX2_CFLAGS="-mavx2"
  ])
], [
  dnl 32-bit: disable AVX2
  arbint_build_avx2=no
  ARBINT_AVX2_CFLAGS=""
])
AM_CONDITIONAL([ARBINT_BUILD_AVX2], [test "x$arbint_build_avx2" = "xyes"])
AC_SUBST([ARBINT_AVX2_CFLAGS])

dnl Combined AVX2+BMI2 flags for files that need both (e.g., arbint_ntt_avx2.c
dnl uses AVX2 intrinsics for future vectorization and BMI2 _mulx_u64 for wide
dnl multiplication). Only set when both features are available.
ARBINT_AVX2_BMI2_CFLAGS=""
AS_IF([test "x$arbint_build_avx2" = "xyes" && test "x$arbint_build_bmi2" = "xyes"], [
  ARBINT_AVX2_BMI2_CFLAGS="$ARBINT_AVX2_CFLAGS $ARBINT_BMI2_CFLAGS"
])
AC_SUBST([ARBINT_AVX2_BMI2_CFLAGS])

ARBINT_CHECK_X86_FEATURE_WITH_FLAG(
  [SHA-NI intrinsics], [arbint_cv_have_sha_ni_plain], [arbint_cv_have_sha_ni_msha],
  [-msha -msse4.1],
  [[
#include <immintrin.h>
  ]],
  [[
__m128i a = _mm_setzero_si128();
__m128i b = _mm_sha256rnds2_epu32(a, a, a);
(void)b;
  ]],
  [HAS_SHA_NI], [HAS_SHA_NI_ALWAYS],
  [Define to 1 if SHA-NI intrinsics are available (with or without -msha).],
  [Define to 1 if SHA-NI intrinsics are available globally.]
)

arbint_build_sha_ni=no
ARBINT_SHA_NI_CFLAGS=""
AS_IF([test "x$arbint_cv_have_sha_ni_plain" = "xyes"], [
  arbint_build_sha_ni=yes
], [test "x$arbint_cv_have_sha_ni_msha" = "xyes"], [
  arbint_build_sha_ni=yes
  ARBINT_SHA_NI_CFLAGS="-msha -msse4.1"
])
AM_CONDITIONAL([ARBINT_BUILD_SHA_NI], [test "x$arbint_build_sha_ni" = "xyes"])
AC_SUBST([ARBINT_SHA_NI_CFLAGS])

dnl Use _mm_crc32_u64 in the test body: it is the widest intrinsic actually
dnl used by our implementation and is only available on x86-64, so the test
dnl naturally fails on 32-bit x86 without any extra host_cpu guard.
ARBINT_CHECK_X86_FEATURE_WITH_FLAG(
  [SSE4.2 CRC32C intrinsics (64-bit)], [arbint_cv_have_sse42_crc32_plain], [arbint_cv_have_sse42_crc32_msse42],
  [-msse4.2],
  [[
#include <nmmintrin.h>
  ]],
  [[
unsigned long long x = _mm_crc32_u64(0ULL, 0ULL);
(void)x;
  ]],
  [HAS_SSE42_CRC32], [HAS_SSE42_CRC32_ALWAYS],
  [Define to 1 if SSE4.2 CRC32C intrinsics are available (with or without -msse4.2).],
  [Define to 1 if SSE4.2 CRC32C intrinsics are available globally.]
)

arbint_build_sse42_crc32=no
ARBINT_SSE42_CRC32_CFLAGS=""
AS_IF([test "x$arbint_cv_have_sse42_crc32_plain" = "xyes"], [
  arbint_build_sse42_crc32=yes
], [test "x$arbint_cv_have_sse42_crc32_msse42" = "xyes"], [
  arbint_build_sse42_crc32=yes
  ARBINT_SSE42_CRC32_CFLAGS="-msse4.2"
])
AM_CONDITIONAL([ARBINT_BUILD_SSE42_CRC32], [test "x$arbint_build_sse42_crc32" = "xyes"])
AC_SUBST([ARBINT_SSE42_CRC32_CFLAGS])

ARBINT_CHECK_ARM_FEATURE_WITH_FLAG(
  [CRC32 intrinsics], [arbint_cv_have_arm_crc32_plain], [arbint_cv_have_arm_crc32_march],
  [-march=armv8-a+crc],
  [[
#include <arm_acle.h>
  ]],
  [[
uint32_t x = __crc32cb(0u, 0u);
(void)x;
  ]],
  [HAS_ARM_CRC32], [HAS_ARM_CRC32_ALWAYS],
  [Define to 1 if ARM CRC32 intrinsics are available (with or without -march flag).],
  [Define to 1 if ARM CRC32 intrinsics are available globally.]
)

ARBINT_CHECK_ARM_FEATURE_WITH_FLAG(
  [PMULL intrinsics], [arbint_cv_have_arm_pmull_plain], [arbint_cv_have_arm_pmull_march],
  [-march=armv8-a+crc+crypto],
  [[
#include <arm_neon.h>
  ]],
  [[
poly128_t r = vmull_p64(1, 1);
(void)r;
  ]],
  [HAS_ARM_PMULL], [HAS_ARM_PMULL_ALWAYS],
  [Define to 1 if ARM PMULL intrinsics are available (with or without -march flag).],
  [Define to 1 if ARM PMULL intrinsics are available globally.]
)

ARBINT_CHECK_ARM_FEATURE_WITH_FLAG(
  [SHA-2 crypto intrinsics], [arbint_cv_have_arm_sha2_plain], [arbint_cv_have_arm_sha2_march],
  [-march=armv8-a+crypto],
  [[
#include <arm_neon.h>
  ]],
  [[
uint32x4_t a = vdupq_n_u32(0);
uint32x4_t b = vsha256hq_u32(a, a, a);
(void)b;
  ]],
  [HAS_ARM_SHA2], [HAS_ARM_SHA2_ALWAYS],
  [Define to 1 if ARM SHA-2 crypto intrinsics are available (with or without -march flag).],
  [Define to 1 if ARM SHA-2 crypto intrinsics are available globally.]
)

arbint_build_arm_sha2=no
ARBINT_ARM_SHA2_CFLAGS=""
AS_IF([test "x$arbint_cv_have_arm_sha2_plain" = "xyes"], [
  arbint_build_arm_sha2=yes
], [test "x$arbint_cv_have_arm_sha2_march" = "xyes"], [
  arbint_build_arm_sha2=yes
  ARBINT_ARM_SHA2_CFLAGS="-march=armv8-a+crypto"
])
AM_CONDITIONAL([ARBINT_BUILD_ARM_SHA2], [test "x$arbint_build_arm_sha2" = "xyes"])
AC_SUBST([ARBINT_ARM_SHA2_CFLAGS])

dnl ARM CRC32C requires both CRC32 and PMULL for the accelerated implementation.
arbint_build_arm_crc32c=no
ARBINT_ARM_CRC32C_CFLAGS=""
arbint_have_arm_crc32=no
arbint_have_arm_pmull=no
arbint_arm_crc32_needs_flag=no
arbint_arm_pmull_needs_flag=no

AS_IF([test "x$arbint_cv_have_arm_crc32_plain" = "xyes"], [
  arbint_have_arm_crc32=yes
], [test "x$arbint_cv_have_arm_crc32_march" = "xyes"], [
  arbint_have_arm_crc32=yes
  arbint_arm_crc32_needs_flag=yes
])

AS_IF([test "x$arbint_cv_have_arm_pmull_plain" = "xyes"], [
  arbint_have_arm_pmull=yes
], [test "x$arbint_cv_have_arm_pmull_march" = "xyes"], [
  arbint_have_arm_pmull=yes
  arbint_arm_pmull_needs_flag=yes
])

AS_IF([test "x$arbint_have_arm_crc32" = "xyes" && test "x$arbint_have_arm_pmull" = "xyes"], [
  arbint_build_arm_crc32c=yes
  dnl Use the most permissive flag needed.
  AS_IF([test "x$arbint_arm_pmull_needs_flag" = "xyes"], [
    ARBINT_ARM_CRC32C_CFLAGS="-march=armv8-a+crc+crypto"
  ], [test "x$arbint_arm_crc32_needs_flag" = "xyes"], [
    ARBINT_ARM_CRC32C_CFLAGS="-march=armv8-a+crc"
  ])
])
AM_CONDITIONAL([ARBINT_BUILD_ARM_CRC32C], [test "x$arbint_build_arm_crc32c" = "xyes"])
AC_SUBST([ARBINT_ARM_CRC32C_CFLAGS])

dnl Linux getauxval for ARM runtime feature detection
AS_IF([test "x$arbint_host_aarch64" = "xyes"], [
  AC_CHECK_HEADERS([sys/auxv.h])
])

AS_IF([test "x$arbint_cc_msvc_cl" = "xyes"], [
  arbint_cv_cc_visibility_hidden=no
], [
  AC_CACHE_CHECK([whether $CC supports -fvisibility=hidden], [arbint_cv_cc_visibility_hidden], [
    arbint_saved_cflags="$CFLAGS"
    CFLAGS="$CFLAGS -fvisibility=hidden"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[]])],
      [arbint_cv_cc_visibility_hidden=yes],
      [arbint_cv_cc_visibility_hidden=no])
    CFLAGS="$arbint_saved_cflags"
  ])
])

AS_IF([test "x$arbint_cv_cc_visibility_hidden" = "xyes"], [
  ARBINT_VISIBILITY_CFLAGS="-fvisibility=hidden"
], [
  ARBINT_VISIBILITY_CFLAGS=""
])
AC_SUBST([ARBINT_VISIBILITY_CFLAGS])

AS_CASE([$host_os],
  [mingw*|msys*|cygwin*], [arbint_host_windows=yes; arbint_host_emscripten=no],
  [emscripten*], [arbint_host_windows=no; arbint_host_emscripten=yes],
  [arbint_host_windows=no; arbint_host_emscripten=no])

dnl -Wl,-Bsymbolic tells the dynamic linker to bind references to global
dnl symbols within the shared library to the definitions inside it.  This is
dnl the linker-side complement of -fno-semantic-interposition (which only
dnl affects the compiler).  Skip on Windows (PE/COFF has no interposition)
dnl and Emscripten (Wasm linker does not support it).
ARBINT_BSYMBOLIC_LDFLAGS=""
AS_IF([test "x$arbint_host_windows" != "xyes" && test "x$arbint_host_emscripten" != "xyes" && test "x$arbint_cc_msvc_cl" != "xyes"], [
  AC_CACHE_CHECK([whether the linker supports -Wl,-Bsymbolic],
    [arbint_cv_ld_bsymbolic], [
    arbint_saved_ldflags="$LDFLAGS"
    LDFLAGS="$LDFLAGS -Wl,-Bsymbolic"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[return 0;]])],
      [arbint_cv_ld_bsymbolic=yes],
      [arbint_cv_ld_bsymbolic=no])
    LDFLAGS="$arbint_saved_ldflags"
  ])
  AS_IF([test "x$arbint_cv_ld_bsymbolic" = "xyes"], [
    ARBINT_BSYMBOLIC_LDFLAGS="-Wl,-Bsymbolic"
  ])
])
AC_SUBST([ARBINT_BSYMBOLIC_LDFLAGS])

ARBINT_LIB_API_CPPFLAGS=""
AS_IF([test "x$arbint_host_windows" = "xyes" && test "x$enable_shared" = "xyes"], [
  ARBINT_LIB_API_CPPFLAGS="$ARBINT_LIB_API_CPPFLAGS -DARBINT_BUILD_DLL"
])
AS_IF([test "x$enable_shared" != "xyes"], [
  ARBINT_LIB_API_CPPFLAGS="$ARBINT_LIB_API_CPPFLAGS -DARBINT_STATIC"
])
AC_SUBST([ARBINT_LIB_API_CPPFLAGS])

ARBINT_PC_API_CFLAGS=""
AS_IF([test "x$arbint_host_windows" = "xyes"], [
  AS_IF([test "x$enable_shared" = "xyes" && test "x$enable_static" != "xyes"], [
    ARBINT_PC_API_CFLAGS="-DARBINT_USE_DLL"
  ], [test "x$enable_shared" != "xyes"], [
    ARBINT_PC_API_CFLAGS="-DARBINT_STATIC"
  ])
])
AC_SUBST([ARBINT_PC_API_CFLAGS])

# Emscripten test flags.
# Tests are MAIN_MODULE executables that dynamically load the library
# (side module).  NODERAWFS gives tests real filesystem access under
# node.js so reference data files (pi10k.txt, e10k.txt) can be read.
ARBINT_EMSCRIPTEN_TEST_LDFLAGS=""
ARBINT_EMSCRIPTEN_TEST_CFLAGS=""
AS_IF([test "x$arbint_host_emscripten" = "xyes"], [
  ARBINT_EMSCRIPTEN_TEST_LDFLAGS="-sMAIN_MODULE=2 -sNODERAWFS"
  ARBINT_EMSCRIPTEN_TEST_CFLAGS="-sMAIN_MODULE=2"
])
AC_SUBST([ARBINT_EMSCRIPTEN_TEST_LDFLAGS])
AC_SUBST([ARBINT_EMSCRIPTEN_TEST_CFLAGS])

# Detect platform entropy sources for RNG seeding
AC_CHECK_HEADERS([unistd.h fcntl.h sys/stat.h])  # Linux /dev/urandom

AS_IF([test "x$arbint_host_windows" = "xyes"], [
  AC_CHECK_HEADERS([windows.h bcrypt.h wincrypt.h])
  # Prefer BCryptGenRandom (Vista+), fall back to CryptGenRandom (XP+),
  # then RtlGenRandom (XP+, only needs windows.h)
])

AC_CHECK_HEADERS([sys/random.h])  # getentropy + arc4random_buf on macOS/BSD
AC_CHECK_FUNCS([arc4random_buf])

# getentropy: we only use it on Emscripten (routes to crypto.getRandomValues).
# It exists on native platforms too, but arc4random_buf / /dev/urandom / WinAPI
# are better and avoid header portability hassles.
AC_CHECK_FUNCS([getentropy])

# Secure memory zeroing: prefer explicit_bzero, then memset_s
AC_CHECK_FUNCS([explicit_bzero])
AC_CHECK_FUNCS([memset_s])

AC_ARG_ENABLE([entropy],
  [AS_HELP_STRING([--disable-entropy],
                  [Do not compile platform entropy sources.
                   arbint_rng_init with a NULL seed will always fail;
                   callers must supply their own seed buffer.])],
  [],
  [enable_entropy=yes])

# Define config.h macros (only when entropy is enabled)
arbint_has_getentropy=no
arbint_has_urandom=no
arbint_has_winapi_entropy=no
arbint_has_arc4random=no
AS_IF([test "x$enable_entropy" = "xyes"], [
  # getentropy: only used on Emscripten (routes to crypto.getRandomValues).
  # On native platforms arc4random_buf, /dev/urandom, or WinAPI are better
  # and avoid header portability issues (macOS wants sys/random.h, Linux
  # wants unistd.h, etc.).
  AS_IF([test "x$arbint_host_emscripten" = "xyes" && test "x$ac_cv_func_getentropy" = "xyes"], [
    AC_DEFINE([ARBINT_HAS_GETENTROPY], [1],
              [Define if getentropy() is available (Emscripten only).])
    arbint_has_getentropy=yes
  ])
  # ARBINT_HAS_URANDOM: only on non-Windows, non-Emscripten POSIX with unistd.h.
  # Emscripten has unistd.h but /dev/urandom is slow JS FS emulation;
  # getentropy is the preferred path there.
  AS_IF([test "x$ac_cv_header_unistd_h" = "xyes" && test "x$arbint_host_windows" != "xyes" && test "x$arbint_host_emscripten" != "xyes"], [
    AC_DEFINE([ARBINT_HAS_URANDOM], [1],
              [Define if /dev/urandom is likely available (non-Windows POSIX).])
    arbint_has_urandom=yes
  ])
  AS_IF([test "x$arbint_host_windows" = "xyes" && test "x$ac_cv_header_windows_h" = "xyes"], [
    AC_DEFINE([ARBINT_HAS_WINAPI_ENTROPY], [1],
              [Define if Windows entropy APIs are available])
    arbint_has_winapi_entropy=yes
  ])
  # Emscripten's arc4random_buf passes link tests but the declaration is
  # not visible from <stdlib.h>, so skip it on emscripten (getentropy is
  # the preferred path there anyway).
  AS_IF([test "x$ac_cv_func_arc4random_buf" = "xyes" && test "x$arbint_host_emscripten" != "xyes"], [
    AC_DEFINE([ARBINT_HAS_ARC4RANDOM], [1],
              [Define if arc4random_buf is available])
    arbint_has_arc4random=yes
  ])
])

# Define conditionals for Makefile.am
AM_CONDITIONAL([ARBINT_HAS_GETENTROPY],
               [test "x$arbint_has_getentropy" = "xyes"])
AM_CONDITIONAL([ARBINT_HAS_URANDOM],
               [test "x$arbint_has_urandom" = "xyes"])
AM_CONDITIONAL([ARBINT_HAS_WINAPI_ENTROPY],
               [test "x$arbint_has_winapi_entropy" = "xyes"])
AM_CONDITIONAL([ARBINT_HAS_ARC4RANDOM],
               [test "x$arbint_has_arc4random" = "xyes"])
AM_CONDITIONAL([ARBINT_HOST_EMSCRIPTEN],
               [test "x$arbint_host_emscripten" = "xyes"])

dnl GMP detection for optional benchmark tools
dnl bench_pi compares c-arbint vs GMP using the Chudnovsky algorithm.
AC_ARG_WITH([gmp],
  [AS_HELP_STRING([--with-gmp], [Enable GMP benchmark tools (default: auto)])],
  [],
  [with_gmp=auto])

arbint_have_gmp=no
AS_IF([test "x$with_gmp" != "xno"], [
  AC_CHECK_HEADER([gmp.h], [
    AC_CHECK_LIB([gmp], [__gmpz_init], [
      arbint_have_gmp=yes
      GMP_LIBS="-lgmp"
    ])
  ])
  AS_IF([test "x$with_gmp" = "xyes" && test "x$arbint_have_gmp" != "xyes"], [
    AC_MSG_ERROR([--with-gmp requested but GMP not found])
  ])
])
AM_CONDITIONAL([ARBINT_HAVE_GMP], [test "x$arbint_have_gmp" = "xyes"])
AC_SUBST([GMP_LIBS])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
  Makefile
  src/Makefile
  include/Makefile
  man/Makefile
  man/c-arbint.3
  tests/Makefile
  demos/Makefile
  tools/Makefile
  c-arbint.pc
])

AC_OUTPUT
