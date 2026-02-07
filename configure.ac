AC_PREREQ([2.69])
AC_INIT([c-arbint], CA_GIT_DESCRIBE_VERSION, [k@iczelia.net], [c-arbint], [https://github.com/iczelia/c-arbint])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
m4_include([m4/arbint_x86_feature_check.m4])
AM_INIT_AUTOMAKE([foreign subdir-objects dist-xz no-dist-gzip tar-ustar])
AM_SILENT_RULES([yes])
AC_CANONICAL_HOST

AS_CASE([$host_cpu],
  [x86_64|amd64], [arbint_host_x86=yes; arbint_host_x86_64=yes; arbint_host_has_cpuid=yes],
  [i[5-9]86|pentium*], [arbint_host_x86=yes; arbint_host_x86_64=no; arbint_host_has_cpuid=yes],
  [i?86], [arbint_host_x86=yes; arbint_host_x86_64=no; arbint_host_has_cpuid=no],
  [arbint_host_x86=no; arbint_host_x86_64=no; arbint_host_has_cpuid=no])

AC_USE_SYSTEM_EXTENSIONS
AC_PROG_CC
AM_PROG_AR
LT_INIT

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

AS_IF([test "x$enable_warnings" = "xyes"], [
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

AC_ARG_WITH([base-optimization],
  [AS_HELP_STRING([--with-base-optimization=@<:@default|tuned@:>@],
                  [Select base optimization profile:
                   default => -O2 -g3 (default),
                   tuned => -O3 -march=... -mtune=...])],
  [arbint_base_optimization="$withval"],
  [arbint_base_optimization=default])

AC_ARG_WITH([march],
  [AS_HELP_STRING([--with-march=ARCH],
                  [Architecture/tune target to use when base optimization is tuned
                   (default: native)])],
  [arbint_cfg_march="$withval"],
  [arbint_cfg_march=""])

AS_CASE([$arbint_base_optimization],
  [default|tuned], [],
  [AC_MSG_ERROR([invalid --with-base-optimization value: $arbint_base_optimization (expected default or tuned)])])

arbint_effective_march=""
arbint_effective_mtune=""
AS_IF([test "x$arbint_base_optimization" = "xdefault"], [
  CFLAGS="$CFLAGS -O2 -g3"
], [
  arbint_effective_march="$arbint_cfg_march"
  AS_IF([test "x$arbint_effective_march" = "x"], [
    arbint_effective_march="native"
  ])
  CFLAGS="$CFLAGS -O3 -march=$arbint_effective_march -mtune=$arbint_effective_march"
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

  AS_IF([test "x$arbint_cv_cc_lto" = "xyes"], [
    arbint_enable_lto_effective=yes
    ARBINT_LTO_CFLAGS="-flto"
    ARBINT_LTO_LDFLAGS="-flto"
  ], [
    AS_IF([test "x$arbint_enable_lto" = "xyes"], [
      AC_MSG_ERROR([--enable-lto was requested, but $CC/$LD does not support -flto])
    ])
  ])
])
AC_SUBST([ARBINT_LTO_CFLAGS])
AC_SUBST([ARBINT_LTO_LDFLAGS])
AC_MSG_NOTICE([LTO requested=$arbint_enable_lto, enabled=$arbint_enable_lto_effective])

AC_CHECK_HEADERS([stdint.h stdlib.h string.h])

AS_IF([test "x$arbint_host_has_cpuid" = "xyes"], [
  arbint_has_cpuid=1
], [
  arbint_has_cpuid=0
])
AC_DEFINE_UNQUOTED([HAS_CPUID], [$arbint_has_cpuid],
  [Define to 1 if x86 CPUID is assumed available on the configured target.])

ARBINT_CHECK_X86_FEATURE_WITH_FLAG(
  [BMI2 intrinsics], [arbint_cv_have_bmi2_plain], [arbint_cv_have_bmi2_mbmi2_v2],
  [-mbmi2],
  [[
#include <immintrin.h>
  ]],
  [[
unsigned int x = _pdep_u32((unsigned int) 3u, (unsigned int) 5u);
return (int) x;
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
AM_CONDITIONAL([ARBINT_BUILD_BMI2], [test "x$arbint_build_bmi2" = "xyes"])
AC_SUBST([ARBINT_BMI2_CFLAGS])

AC_CACHE_CHECK([whether $CC supports -fvisibility=hidden], [arbint_cv_cc_visibility_hidden], [
  arbint_saved_cflags="$CFLAGS"
  CFLAGS="$CFLAGS -fvisibility=hidden"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[]])],
    [arbint_cv_cc_visibility_hidden=yes],
    [arbint_cv_cc_visibility_hidden=no])
  CFLAGS="$arbint_saved_cflags"
])

AS_IF([test "x$arbint_cv_cc_visibility_hidden" = "xyes"], [
  ARBINT_VISIBILITY_CFLAGS="-fvisibility=hidden"
], [
  ARBINT_VISIBILITY_CFLAGS=""
])
AC_SUBST([ARBINT_VISIBILITY_CFLAGS])

AS_CASE([$host_os],
  [mingw*|msys*|cygwin*], [arbint_host_windows=yes],
  [arbint_host_windows=no])

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

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
  Makefile
  src/Makefile
  include/Makefile
  man/Makefile
  man/c-arbint.3
  tests/Makefile
  c-arbint.pc
])
AC_OUTPUT
