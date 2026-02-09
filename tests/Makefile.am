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

AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/src -I$(top_builddir)
AM_CFLAGS = $(ARBINT_LTO_CFLAGS)
AM_LDFLAGS = $(ARBINT_LTO_LDFLAGS)

check_PROGRAMS = test_basic test_mul test_tdiv test_large test_stochastic test_addsub_avx2 test_random test_hash compute_pi
TESTS = $(check_PROGRAMS)

TESTS_ENVIRONMENT = LD_LIBRARY_PATH="$(top_builddir)/src/.libs:$$LD_LIBRARY_PATH"

test_basic_SOURCES = test_basic.c
test_basic_LDADD = $(top_builddir)/src/libarbint.la
test_basic_DEPENDENCIES = $(top_builddir)/src/libarbint.la

test_mul_SOURCES = test_mul.c
test_mul_LDADD = $(top_builddir)/src/libarbint.la
test_mul_DEPENDENCIES = $(top_builddir)/src/libarbint.la

test_tdiv_SOURCES = test_tdiv.c
test_tdiv_LDADD = $(top_builddir)/src/libarbint.la
test_tdiv_DEPENDENCIES = $(top_builddir)/src/libarbint.la

test_large_SOURCES = test_large.c
test_large_LDADD = $(top_builddir)/src/libarbint.la
test_large_DEPENDENCIES = $(top_builddir)/src/libarbint.la

test_stochastic_SOURCES = test_stochastic.c
test_stochastic_LDADD = $(top_builddir)/src/libarbint.la
test_stochastic_DEPENDENCIES = $(top_builddir)/src/libarbint.la

test_addsub_avx2_SOURCES = test_addsub_avx2.c
test_addsub_avx2_LDADD = $(top_builddir)/src/libarbint.la
test_addsub_avx2_DEPENDENCIES = $(top_builddir)/src/libarbint.la

test_random_SOURCES = test_random.c
test_random_LDADD = $(top_builddir)/src/libarbint.la
test_random_DEPENDENCIES = $(top_builddir)/src/libarbint.la

test_hash_SOURCES = test_hash.c
test_hash_LDADD = $(top_builddir)/src/libarbint.la
test_hash_DEPENDENCIES = $(top_builddir)/src/libarbint.la

compute_pi_SOURCES = compute_pi.c
compute_pi_LDADD = $(top_builddir)/src/libarbint.la
compute_pi_DEPENDENCIES = $(top_builddir)/src/libarbint.la

noinst_HEADERS = test_framework.h

EXTRA_DIST = test_basic.c test_mul.c test_tdiv.c test_large.c test_stochastic.c test_addsub_avx2.c test_random.c test_hash.c compute_pi.c test_framework.h pi10k.txt
