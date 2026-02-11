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

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = include src man tests

# tools/ is not built by default; use 'make tools' to build
DIST_SUBDIRS = $(SUBDIRS) tools

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = c-arbint.pc

dist_doc_DATA = README.md

EXTRA_DIST = bootstrap

.PHONY: clang-format tools

clang-format:
	clang-format -i $$(find include src tests tools -type f \( -name '*.h' -o -name '*.c' \))

# Build tools (tune_thresholds, etc.)
tools:
	$(MAKE) -C tools
