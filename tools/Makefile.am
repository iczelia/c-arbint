# Process this file with automake to produce Makefile.in

# These are development/tuning tools, not installed by default.
# Build with: make -C tools tune_thresholds
# Build with: make -C tools test_segmented_ntt

noinst_PROGRAMS = tune_thresholds test_segmented_ntt

# These tools link directly against the library object files (.o/.lo)
# rather than the shared library. This allows them to access internal
# functions (arbint_mul_mag_ntt, arbint_sqr_mag_ntt, etc.) that have
# external linkage but are not exported from libarbint.so due to
# -fvisibility=hidden.

# Extract object files from libarbint.la convenience vars.
# We link against the static .o files inside .libs/ directory.
ARBINT_OBJS = $(top_builddir)/src/.libs/*.o

tune_thresholds_SOURCES = tune_thresholds.c
tune_thresholds_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/src
tune_thresholds_LDADD = $(ARBINT_OBJS) $(LIBS)

test_segmented_ntt_SOURCES = test_segmented_ntt.c
test_segmented_ntt_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/src
test_segmented_ntt_LDADD = $(ARBINT_OBJS) $(LIBS)

# Standalone subdir builds (make -C tools ...) must refresh src first.
.PHONY: FORCE
$(top_builddir)/src/libarbint.la: FORCE
	$(MAKE) $(AM_MAKEFLAGS) -C $(top_builddir)/src libarbint.la
FORCE:

## Convenience targets
tune: tune_thresholds
	./tune_thresholds

test-segmented: test_segmented_ntt
	./test_segmented_ntt
