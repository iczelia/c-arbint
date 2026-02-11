# arbint - portable arbitrary-precision computation library
#
# Copyright (C) 2026 Kamila Szewczyk (k@iczelia.net)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

AM_CPPFLAGS = -I$(top_srcdir)/include

# Sources and headers shared between all platforms.
ARBINT_CORE_SOURCES = \
  arbint_base.c arbint_assign.c arbint_cmp.c arbint_addsub.c \
  arbint_div.c arbint_tdiv.c arbint_fdiv.c arbint_cdiv.c \
  arbint_iex.c arbint_cpu.c arbint_mul.c \
  arbint_mul_generic.c arbint_div_generic.c arbint_bitops.c \
  arbint_shift.c arbint_random.c arbint_hash.c arbint_hash_generic.c \
  arbint_ntt.c arbint_ntt_generic.c

noinst_HEADERS = \
  arbint_base.h arbint_assign.h arbint_cmp.h arbint_addsub.h \
  arbint_div.h arbint_iex.h arbint_cpu.h arbint_mul.h \
  arbint_mul_core.inc arbint_bitops.h arbint_shift.h arbint_random.h \
  arbint_hash.h arbint_ntt.h

# Platform-specific sources: distributed but compiled only via custom rules
# (not EXTRA_*_SOURCES, which triggers broken automake 1.18 auto-rules).
EXTRA_DIST = \
  arbint_mul_bmi2.c arbint_div_bmi2.c arbint_addsub_avx2.c \
  arbint_bitops_popcnt.c arbint_random_getentropy.c \
  arbint_random_urandom.c arbint_random_winapi.c arbint_random_arc4.c \
  arbint_hash_shani.c arbint_hash_sha2_neon.c \
  arbint_hash_crc32c_sse42.c arbint_hash_crc32c_neon.c \
  arbint_ntt_avx2.c

# -----------------------------------------------------------------------
#  Emscripten: bypass libtool's broken link mode.
#
#  Libtool's ld_shlibs detection has no emscripten case and falls through
#  to an ELF check that fails, producing a broken link command.  We use
#  libtool only for compilation (.c -> .lo/.o) via a convenience library,
#  then link the side module directly with emcc.
# -----------------------------------------------------------------------

if ARBINT_HOST_EMSCRIPTEN

noinst_LTLIBRARIES = libarbint_objs.la
libarbint_objs_la_SOURCES = $(ARBINT_CORE_SOURCES)
libarbint_objs_la_LIBADD =
libarbint_objs_la_CPPFLAGS = $(AM_CPPFLAGS) $(ARBINT_LIB_API_CPPFLAGS)
libarbint_objs_la_CFLAGS = $(ARBINT_VISIBILITY_CFLAGS) $(ARBINT_LTO_CFLAGS)

ARBINT_LTCOMPILE_EMCC = $(AM_V_CC)$(LIBTOOL) $(AM_V_lt) --tag=CC \
	$(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile \
	$(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(libarbint_objs_la_CPPFLAGS) $(CPPFLAGS) \
	$(libarbint_objs_la_CFLAGS) $(CFLAGS)

if ARBINT_HAS_GETENTROPY
libarbint_objs_la_LIBADD += libarbint_objs_la-arbint_random_getentropy.lo
endif

libarbint_objs_la-arbint_random_getentropy.lo: arbint_random_getentropy.c
	$(ARBINT_LTCOMPILE_EMCC) \
	  -c -o $@ `test -f '$<' || echo '$(srcdir)/'`$<

# Link the convenience archive into a wasm side module directly with emcc,
# bypassing libtool's broken --mode=link for emscripten.
libarbint.wasm: libarbint_objs.la
	$(AM_V_CCLD)$(CC) -sSIDE_MODULE $(ARBINT_LTO_LDFLAGS) $(LDFLAGS) \
	  -o $@ .libs/libarbint_objs.a

all-local: libarbint.wasm
clean-local:
	-rm -f libarbint.wasm

install-exec-local: libarbint.wasm
	$(MKDIR_P) '$(DESTDIR)$(libdir)'
	$(INSTALL_DATA) libarbint.wasm '$(DESTDIR)$(libdir)/libarbint.wasm'

uninstall-local:
	-rm -f '$(DESTDIR)$(libdir)/libarbint.wasm'

else
# -----------------------------------------------------------------------
#  Native platforms: standard libtool shared library build.
# -----------------------------------------------------------------------

lib_LTLIBRARIES = libarbint.la
libarbint_la_SOURCES = $(ARBINT_CORE_SOURCES)
libarbint_la_LIBADD =
libarbint_la_CPPFLAGS = $(AM_CPPFLAGS) $(ARBINT_LIB_API_CPPFLAGS)
libarbint_la_CFLAGS = $(ARBINT_VISIBILITY_CFLAGS) $(ARBINT_LTO_CFLAGS)
libarbint_la_LDFLAGS = -version-info 0:0:0 -no-undefined $(ARBINT_BSYMBOLIC_LDFLAGS) $(ARBINT_LTO_LDFLAGS)

# Libtool compilation command prefix for platform-specific object files.
ARBINT_LTCOMPILE = $(AM_V_CC)$(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile \
	$(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(libarbint_la_CPPFLAGS) $(CPPFLAGS) $(libarbint_la_CFLAGS) $(CFLAGS)

if ARBINT_BUILD_BMI2
libarbint_la_LIBADD += libarbint_la-arbint_mul_bmi2.lo libarbint_la-arbint_div_bmi2.lo
endif

if ARBINT_BUILD_AVX2
libarbint_la_LIBADD += libarbint_la-arbint_addsub_avx2.lo libarbint_la-arbint_ntt_avx2.lo
endif

if ARBINT_BUILD_POPCNT
libarbint_la_LIBADD += libarbint_la-arbint_bitops_popcnt.lo
endif

if ARBINT_HAS_GETENTROPY
libarbint_la_LIBADD += libarbint_la-arbint_random_getentropy.lo
endif

if ARBINT_HAS_URANDOM
libarbint_la_LIBADD += libarbint_la-arbint_random_urandom.lo
endif

if ARBINT_HAS_WINAPI_ENTROPY
libarbint_la_LIBADD += libarbint_la-arbint_random_winapi.lo
endif

if ARBINT_HAS_ARC4RANDOM
libarbint_la_LIBADD += libarbint_la-arbint_random_arc4.lo
endif

if ARBINT_BUILD_SHA_NI
libarbint_la_LIBADD += libarbint_la-arbint_hash_shani.lo
endif

if ARBINT_BUILD_SSE42_CRC32
libarbint_la_LIBADD += libarbint_la-arbint_hash_crc32c_sse42.lo
endif

if ARBINT_BUILD_ARM_SHA2
libarbint_la_LIBADD += libarbint_la-arbint_hash_sha2_neon.lo
endif

if ARBINT_BUILD_ARM_CRC32C
libarbint_la_LIBADD += libarbint_la-arbint_hash_crc32c_neon.lo
endif

libarbint_la-arbint_mul_bmi2.lo: arbint_mul_bmi2.c
	$(ARBINT_LTCOMPILE) $(ARBINT_BMI2_CFLAGS) \
	  -c -o $@ `test -f '$<' || echo '$(srcdir)/'`$<

libarbint_la-arbint_div_bmi2.lo: arbint_div_bmi2.c
	$(ARBINT_LTCOMPILE) $(ARBINT_BMI2_CFLAGS) \
	  -c -o $@ `test -f '$<' || echo '$(srcdir)/'`$<

libarbint_la-arbint_addsub_avx2.lo: arbint_addsub_avx2.c
	$(ARBINT_LTCOMPILE) $(ARBINT_AVX2_CFLAGS) \
	  -c -o $@ `test -f '$<' || echo '$(srcdir)/'`$<

libarbint_la-arbint_ntt_avx2.lo: arbint_ntt_avx2.c
	$(ARBINT_LTCOMPILE) $(ARBINT_AVX2_CFLAGS) \
	  -c -o $@ `test -f '$<' || echo '$(srcdir)/'`$<

libarbint_la-arbint_bitops_popcnt.lo: arbint_bitops_popcnt.c
	$(ARBINT_LTCOMPILE) $(ARBINT_POPCNT_CFLAGS) \
	  -c -o $@ `test -f '$<' || echo '$(srcdir)/'`$<

libarbint_la-arbint_random_getentropy.lo: arbint_random_getentropy.c
	$(ARBINT_LTCOMPILE) \
	  -c -o $@ `test -f '$<' || echo '$(srcdir)/'`$<

libarbint_la-arbint_random_urandom.lo: arbint_random_urandom.c
	$(ARBINT_LTCOMPILE) \
	  -c -o $@ `test -f '$<' || echo '$(srcdir)/'`$<

libarbint_la-arbint_random_winapi.lo: arbint_random_winapi.c
	$(ARBINT_LTCOMPILE) \
	  -c -o $@ `test -f '$<' || echo '$(srcdir)/'`$<

libarbint_la-arbint_random_arc4.lo: arbint_random_arc4.c
	$(ARBINT_LTCOMPILE) \
	  -c -o $@ `test -f '$<' || echo '$(srcdir)/'`$<

libarbint_la-arbint_hash_shani.lo: arbint_hash_shani.c
	$(ARBINT_LTCOMPILE) $(ARBINT_SHA_NI_CFLAGS) \
	  -c -o $@ `test -f '$<' || echo '$(srcdir)/'`$<

libarbint_la-arbint_hash_crc32c_sse42.lo: arbint_hash_crc32c_sse42.c
	$(ARBINT_LTCOMPILE) $(ARBINT_SSE42_CRC32_CFLAGS) \
	  -c -o $@ `test -f '$<' || echo '$(srcdir)/'`$<

libarbint_la-arbint_hash_sha2_neon.lo: arbint_hash_sha2_neon.c
	$(ARBINT_LTCOMPILE) $(ARBINT_ARM_SHA2_CFLAGS) \
	  -c -o $@ `test -f '$<' || echo '$(srcdir)/'`$<

libarbint_la-arbint_hash_crc32c_neon.lo: arbint_hash_crc32c_neon.c
	$(ARBINT_LTCOMPILE) $(ARBINT_ARM_CRC32C_CFLAGS) \
	  -c -o $@ `test -f '$<' || echo '$(srcdir)/'`$<

endif
