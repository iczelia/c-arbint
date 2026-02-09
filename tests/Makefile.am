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

check_PROGRAMS = \
  test_basic \
  test_mul \
  test_tdiv \
  test_large \
  test_stochastic \
  test_addsub_avx2 \
  test_random \
  test_hash \
  test_shift \
  test_bitops \
  compute_pi \
  compute_e \
  compute_bitops

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = LD_LIBRARY_PATH="$(top_builddir)/src/.libs:$$LD_LIBRARY_PATH"

# Default linking for all test programs
# Automake automatically detects test_name.c as source for test_name binary
LDADD = $(top_builddir)/src/libarbint.la

noinst_HEADERS = test_framework.h
EXTRA_DIST = test_framework.h pi10k.txt e10k.txt
