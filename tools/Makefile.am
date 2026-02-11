## Process this file with automake to produce Makefile.in

## These are development/tuning tools, not installed by default.
## Build with: make -C tools tune_thresholds

noinst_PROGRAMS = tune_thresholds

tune_thresholds_SOURCES = tune_thresholds.c
tune_thresholds_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/src
tune_thresholds_LDADD = $(top_builddir)/src/libarbint.la

# Standalone subdir builds (make -C tools ...) must refresh src first.
.PHONY: FORCE
$(top_builddir)/src/libarbint.la: FORCE
	$(MAKE) $(AM_MAKEFLAGS) -C $(top_builddir)/src libarbint.la
FORCE:

## Convenience target
tune: tune_thresholds
	./tune_thresholds
