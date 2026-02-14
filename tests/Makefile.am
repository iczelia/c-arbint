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
AM_CFLAGS = $(ARBINT_LTO_CFLAGS) $(ARBINT_EMSCRIPTEN_TEST_CFLAGS)
AM_LDFLAGS = $(ARBINT_LTO_LDFLAGS) $(ARBINT_EMSCRIPTEN_TEST_LDFLAGS)

check_PROGRAMS = \
  test_basic \
  test_mul \
  test_mul_exhaustive \
  test_sqr \
  test_tdiv \
  test_div_pow2 \
  test_div_newton \
  test_divisible \
  test_large \
  test_stochastic \
  test_addsub_avx2 \
  test_random \
  test_hash \
  test_cache \
  test_shift \
  test_bitops \
  test_pow \
  test_isqrt \
  test_root \
  test_toom3 \
  test_addmul \
  test_ntt \
  test_gcd \
  test_combin \
  test_is_power \
  test_xgcd \
  test_inv_mod \
  test_removefactor \
  test_numtheory \
  test_str \
  compute_pi \
  compute_e \
  compute_bitops

TESTS = $(check_PROGRAMS)
TESTS_ENVIRONMENT = LD_LIBRARY_PATH="$(top_builddir)/src/.libs:$$LD_LIBRARY_PATH"

# Default linking for all test programs.
# Automake automatically detects test_name.c as source for test_name binary.
# On emscripten the library is a wasm side module (no .la file); tests link
# against it directly with emcc -sMAIN_MODULE=2.
if ARBINT_HOST_EMSCRIPTEN
LDADD = $(top_builddir)/src/libarbint.wasm
else
LDADD = $(top_builddir)/src/libarbint.la
endif

# Standalone subdir builds (make -C tests ...) must refresh src first.
# Use a FORCE-driven rule on the linked library path so src is checked every
# time, but test binaries are relinked only if the library actually changed.
.PHONY: FORCE
if ARBINT_HOST_EMSCRIPTEN
$(top_builddir)/src/libarbint.wasm: FORCE
	$(MAKE) $(AM_MAKEFLAGS) -C $(top_builddir)/src libarbint.wasm
else
$(top_builddir)/src/libarbint.la: FORCE
	$(MAKE) $(AM_MAKEFLAGS) -C $(top_builddir)/src libarbint.la
endif
FORCE:

# Emscripten: libtool wrappers cannot exec wasm binaries directly.
# Copy the side module into .libs/ so the emscripten runtime can find it,
# then run tests via node on the .js files in .libs/.
if ARBINT_HOST_EMSCRIPTEN
LOG_COMPILER = $(SHELL)
AM_LOG_FLAGS = $(srcdir)/run_wasm_test.sh
endif

noinst_HEADERS = test_framework.h
EXTRA_DIST = test_framework.h pi10k.txt e10k.txt run_wasm_test.sh
