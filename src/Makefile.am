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

lib_LTLIBRARIES = libarbint.la
libarbint_la_SOURCES = arbint_base.c arbint_assign.c arbint_cmp.c arbint_addsub.c arbint_tdiv.c arbint_iex.c arbint_cpu.c arbint_mul.c arbint_mul_generic.c arbint_tdiv_generic.c arbint_bitops.c arbint_random.c arbint_hash.c arbint_hash_generic.c
EXTRA_libarbint_la_SOURCES = arbint_mul_bmi2.c arbint_tdiv_bmi2.c arbint_addsub_avx2.c arbint_random_urandom.c arbint_random_winapi.c arbint_random_arc4.c arbint_hash_shani.c arbint_hash_crc32c.c
libarbint_la_LIBADD =
libarbint_la_CPPFLAGS = $(AM_CPPFLAGS) $(ARBINT_LIB_API_CPPFLAGS)
libarbint_la_CFLAGS = $(ARBINT_VISIBILITY_CFLAGS) $(ARBINT_LTO_CFLAGS)
noinst_HEADERS = arbint_base.h arbint_assign.h arbint_cmp.h arbint_addsub.h arbint_tdiv.h arbint_iex.h arbint_cpu.h arbint_mul.h arbint_bitops.h arbint_random.h arbint_hash.h
libarbint_la_LDFLAGS = -version-info 0:0:0 -no-undefined $(ARBINT_LTO_LDFLAGS)

if ARBINT_BUILD_BMI2
libarbint_la_LIBADD += libarbint_la-arbint_mul_bmi2.lo
libarbint_la_LIBADD += libarbint_la-arbint_tdiv_bmi2.lo
endif

if ARBINT_BUILD_AVX2
libarbint_la_LIBADD += libarbint_la-arbint_addsub_avx2.lo
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
libarbint_la_LIBADD += libarbint_la-arbint_hash_crc32c.lo
endif

libarbint_la-arbint_mul_bmi2.lo: arbint_mul_bmi2.c
	$(AM_V_CC)$(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(libarbint_la_CPPFLAGS) $(CPPFLAGS) $(libarbint_la_CFLAGS) $(CFLAGS) $(ARBINT_BMI2_CFLAGS) -c -o $@ `test -f 'arbint_mul_bmi2.c' || echo '$(srcdir)/'`arbint_mul_bmi2.c

libarbint_la-arbint_tdiv_bmi2.lo: arbint_tdiv_bmi2.c
	$(AM_V_CC)$(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(libarbint_la_CPPFLAGS) $(CPPFLAGS) $(libarbint_la_CFLAGS) $(CFLAGS) $(ARBINT_BMI2_CFLAGS) -c -o $@ `test -f 'arbint_tdiv_bmi2.c' || echo '$(srcdir)/'`arbint_tdiv_bmi2.c

libarbint_la-arbint_addsub_avx2.lo: arbint_addsub_avx2.c
	$(AM_V_CC)$(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(libarbint_la_CPPFLAGS) $(CPPFLAGS) $(libarbint_la_CFLAGS) $(CFLAGS) $(ARBINT_AVX2_CFLAGS) -c -o $@ `test -f 'arbint_addsub_avx2.c' || echo '$(srcdir)/'`arbint_addsub_avx2.c

libarbint_la-arbint_random_urandom.lo: arbint_random_urandom.c
	$(AM_V_CC)$(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(libarbint_la_CPPFLAGS) $(CPPFLAGS) $(libarbint_la_CFLAGS) $(CFLAGS) -c -o $@ `test -f 'arbint_random_urandom.c' || echo '$(srcdir)/'`arbint_random_urandom.c

libarbint_la-arbint_random_winapi.lo: arbint_random_winapi.c
	$(AM_V_CC)$(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(libarbint_la_CPPFLAGS) $(CPPFLAGS) $(libarbint_la_CFLAGS) $(CFLAGS) -c -o $@ `test -f 'arbint_random_winapi.c' || echo '$(srcdir)/'`arbint_random_winapi.c

libarbint_la-arbint_random_arc4.lo: arbint_random_arc4.c
	$(AM_V_CC)$(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(libarbint_la_CPPFLAGS) $(CPPFLAGS) $(libarbint_la_CFLAGS) $(CFLAGS) -c -o $@ `test -f 'arbint_random_arc4.c' || echo '$(srcdir)/'`arbint_random_arc4.c

libarbint_la-arbint_hash_shani.lo: arbint_hash_shani.c
	$(AM_V_CC)$(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(libarbint_la_CPPFLAGS) $(CPPFLAGS) $(libarbint_la_CFLAGS) $(CFLAGS) $(ARBINT_SHA_NI_CFLAGS) -c -o $@ `test -f 'arbint_hash_shani.c' || echo '$(srcdir)/'`arbint_hash_shani.c

libarbint_la-arbint_hash_crc32c.lo: arbint_hash_crc32c.c
	$(AM_V_CC)$(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(libarbint_la_CPPFLAGS) $(CPPFLAGS) $(libarbint_la_CFLAGS) $(CFLAGS) $(ARBINT_SSE42_CRC32_CFLAGS) -c -o $@ `test -f 'arbint_hash_crc32c.c' || echo '$(srcdir)/'`arbint_hash_crc32c.c
