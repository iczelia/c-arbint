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

AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_builddir)
AM_CFLAGS = $(ARBINT_LTO_CFLAGS) $(ARBINT_EMSCRIPTEN_TEST_CFLAGS)
AM_LDFLAGS = $(ARBINT_LTO_LDFLAGS) $(ARBINT_EMSCRIPTEN_TEST_LDFLAGS)

# Demo programs are built but not run directly as tests
# (dc requires --test flag, handled via wrapper)
check_PROGRAMS = dc

# Use a wrapper script to pass --test to dc
TESTS = run_dc_test.sh
EXTRA_DIST = run_dc_test.sh

if ARBINT_HOST_EMSCRIPTEN
# For emscripten, the wrapper uses run_wasm_test.sh internally
TESTS_ENVIRONMENT = \
    srcdir="$(srcdir)" \
    top_srcdir="$(top_srcdir)" \
    DC_BINARY="$(builddir)/dc" \
    WASM_RUNNER="$(top_srcdir)/tests/run_wasm_test.sh"
else
TESTS_ENVIRONMENT = \
    LD_LIBRARY_PATH="$(top_builddir)/src/.libs:$$LD_LIBRARY_PATH" \
    DC_BINARY="$(builddir)/dc$(EXEEXT)"
endif

# Default linking for demo programs
if ARBINT_HOST_EMSCRIPTEN
LDADD = $(top_builddir)/src/libarbint.wasm
else
LDADD = $(top_builddir)/src/libarbint.la
endif

# Standalone subdir builds must refresh src first
.PHONY: FORCE
if ARBINT_HOST_EMSCRIPTEN
$(top_builddir)/src/libarbint.wasm: FORCE
	$(MAKE) $(AM_MAKEFLAGS) -C $(top_builddir)/src libarbint.wasm
else
$(top_builddir)/src/libarbint.la: FORCE
	$(MAKE) $(AM_MAKEFLAGS) -C $(top_builddir)/src libarbint.la
endif
FORCE:

# dc demo: a reverse-polish calculator using c-arbint
dc_SOURCES = dc.c
