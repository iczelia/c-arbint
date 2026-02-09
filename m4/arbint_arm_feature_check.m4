dnl ARBINT_CHECK_ARM_FEATURE_WITH_FLAG
dnl ---------------------------------
dnl Probe an ARM feature by compiling and linking a test program first without
dnl extra flags, then (if needed) with a specified compiler flag.
dnl
dnl Parameters:
dnl  1: feature description (for AC_CACHE_CHECK messages)
dnl  2: cache var for plain probe result (yes/no)
dnl  3: cache var for flagged probe result (yes/no)
dnl  4: compiler flag to retry with (e.g. -march=armv8-a+crc)
dnl  5: AC_LANG_PROGRAM "includes" fragment
dnl  6: AC_LANG_PROGRAM "body" fragment
dnl  7: AC_DEFINE symbol: feature available with/without flag
dnl  8: AC_DEFINE symbol: feature always available (without extra flag)
dnl  9: description for symbol 7
dnl 10: description for symbol 8
AC_DEFUN([ARBINT_CHECK_ARM_FEATURE_WITH_FLAG], [
  AC_CACHE_CHECK([for ARM $1 without extra flags], [$2], [
    AS_IF([test "x$arbint_host_aarch64" != "xyes"], [
      AS_VAR_SET([$2], [no])
    ], [
      AC_LINK_IFELSE([AC_LANG_PROGRAM([$5], [$6])],
        [AS_VAR_SET([$2], [yes])],
        [AS_VAR_SET([$2], [no])])
    ])
  ])

  AS_VAR_IF([$2], [yes], [], [
    AC_CACHE_CHECK([for ARM $1 with $4], [$3], [
      AS_IF([test "x$arbint_host_aarch64" != "xyes"], [
        AS_VAR_SET([$3], [no])
      ], [
        arbint_saved_cflags="$CFLAGS"
        CFLAGS="$CFLAGS $4"
        AC_LINK_IFELSE([AC_LANG_PROGRAM([$5], [$6])],
          [AS_VAR_SET([$3], [yes])],
          [AS_VAR_SET([$3], [no])])
        CFLAGS="$arbint_saved_cflags"
      ])
    ])
  ])

  AS_VAR_IF([$2], [yes], [
    AC_DEFINE([$7], [1], [$9])
    AC_DEFINE([$8], [1], [$10])
  ], [
    AS_VAR_IF([$3], [yes], [
      AC_DEFINE([$7], [1], [$9])
      AC_DEFINE([$8], [0], [$10])
    ], [
      AC_DEFINE([$7], [0], [$9])
      AC_DEFINE([$8], [0], [$10])
    ])
  ])
])
