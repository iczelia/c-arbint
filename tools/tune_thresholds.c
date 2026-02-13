/*  arbint - portable arbitrary-precision computation library

    Copyright (C) 2026 Kamila Szewczyk (k@iczelia.net)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program. If not, see <https://www.gnu.org/licenses/>.  */

/*  Threshold tuning tool for arbint multiplication algorithms.

    This tool measures the performance of multiplication and squaring
    operations at various operand sizes to help determine optimal threshold
    values for algorithm selection (schoolbook -> Karatsuba -> Toom-3 -> NTT).

    Usage: tune_thresholds [--mul] [--sqr] [--ntt] [--ntt-sqr] [--csv]

    This tool links directly against the library object files (not the
    shared library) so it can call internal functions like arbint_mul_mag_ntt
    and arbint_sqr_mag_ntt directly. This allows accurate comparison of
    algorithms without depending on compile-time thresholds.

    The --ntt and --ntt-sqr options directly measure Toom-3 vs NTT to find
    the true crossover point where NTT becomes faster.  */

#include <arbint.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "arbint_mul.h"
#include "arbint_ntt.h"
#include "arbint_base.h"

/*  Stringify helper for printing current threshold values.  */
#define ARBINT_STR_(x) #x
#define ARBINT_STR(x) ARBINT_STR_(x)

/*  Number of iterations for timing.  */
#define WARMUP_ITERS 3
#define MIN_ITERS 5
#define TARGET_TIME_NS 100000000 /* 100ms target per size */

/*  Size ranges to test (in limbs).  */
#define MIN_SIZE 4
#define MAX_SIZE 300
#define SIZE_STEP_SMALL 2  /* step for sizes < 50 */
#define SIZE_STEP_MEDIUM 4 /* step for sizes 50-150 */
#define SIZE_STEP_LARGE 8  /* step for sizes > 150 */

/*  NTT-specific size ranges (larger operands).  */
#define NTT_MIN_SIZE 128
#define NTT_MAX_SIZE 8192
#define NTT_SIZE_STEP_SMALL 64   /* step for sizes < 512 */
#define NTT_SIZE_STEP_MEDIUM 128 /* step for sizes 512-2048 */
#define NTT_SIZE_STEP_LARGE 256  /* step for sizes > 2048 */

/*  Output format.  */
static int g_csv_mode = 0;

/*  Global RNG for generating random operands.  */
static arbint_rng_t g_rng;

/*  Simple LCG for filling limb arrays with pseudorandom data.
    Used by the direct NTT measurement functions.  */
static uint32_t g_lcg_state = 0x12345678u;

static uint32_t lcg_rand(void) {
  g_lcg_state = g_lcg_state * 1103515245u + 12345u;
  return g_lcg_state;
}

/*  Get current time in nanoseconds (platform-specific).  */
static uint64_t get_time_ns(void) {
#if defined(_WIN32)
  LARGE_INTEGER freq;
  LARGE_INTEGER count;
  QueryPerformanceFrequency(&freq);
  QueryPerformanceCounter(&count);
  return (uint64_t) ((count.QuadPart * 1000000000ull) / freq.QuadPart);
#elif defined(__APPLE__)
  return clock_gettime_nsec_np(CLOCK_MONOTONIC);
#else
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (uint64_t) ts.tv_sec * 1000000000ull + (uint64_t) ts.tv_nsec;
#endif
}

/*  Build an arbint with exactly n limbs of random data.  */
static void build_random_value(arbint_t a, size_t target_limbs) {
  /*  Generate a random value with target_limbs * LIMB_BITS bits.
      This gives us exactly target_limbs limbs of uniformly random data.  */
  arbint_urandomb(a, &g_rng, target_limbs * sizeof(arbint_limb_t) * 8);
}

/*  Measure time for multiplication at given size.
    Returns average nanoseconds per operation.  */
static double measure_mul(arbint_t a, arbint_t b, arbint_t r, size_t n) {
  uint64_t start;
  uint64_t end;
  uint64_t total_ns;
  int iters;
  int i;

  build_random_value(a, n);
  build_random_value(b, n);

  /*  Warmup.  */
  for (i = 0; i < WARMUP_ITERS; ++i)
    arbint_mul(r, a, b);

  /*  Determine iteration count for accurate timing.  */
  start = get_time_ns();
  for (i = 0; i < MIN_ITERS; ++i)
    arbint_mul(r, a, b);
  end = get_time_ns();

  total_ns = end - start;
  if (total_ns < TARGET_TIME_NS && total_ns > 0) {
    iters = (int) ((TARGET_TIME_NS * MIN_ITERS) / total_ns);
    if (iters < MIN_ITERS)
      iters = MIN_ITERS;
    if (iters > 10000)
      iters = 10000;
  } else {
    iters = MIN_ITERS;
  }

  /*  Actual measurement.  */
  start = get_time_ns();
  for (i = 0; i < iters; ++i)
    arbint_mul(r, a, b);
  end = get_time_ns();

  return (double) (end - start) / (double) iters;
}

/*  Measure time for squaring at given size.
    Returns average nanoseconds per operation.  */
static double measure_sqr(arbint_t a, arbint_t r, size_t n) {
  uint64_t start;
  uint64_t end;
  uint64_t total_ns;
  int iters;
  int i;

  build_random_value(a, n);

  /*  Warmup.  */
  for (i = 0; i < WARMUP_ITERS; ++i)
    arbint_sqr(r, a);

  /*  Determine iteration count.  */
  start = get_time_ns();
  for (i = 0; i < MIN_ITERS; ++i)
    arbint_sqr(r, a);
  end = get_time_ns();

  total_ns = end - start;
  if (total_ns < TARGET_TIME_NS && total_ns > 0) {
    iters = (int) ((TARGET_TIME_NS * MIN_ITERS) / total_ns);
    if (iters < MIN_ITERS)
      iters = MIN_ITERS;
    if (iters > 10000)
      iters = 10000;
  } else {
    iters = MIN_ITERS;
  }

  /*  Actual measurement.  */
  start = get_time_ns();
  for (i = 0; i < iters; ++i)
    arbint_sqr(r, a);
  end = get_time_ns();

  return (double) (end - start) / (double) iters;
}

/*  Compute normalized time (ns per limb^2 for schoolbook comparison).  */
static double normalize_time(double ns, size_t n) {
  return ns / ((double) n * (double) n);
}

/*  Find the crossover point where algorithm changes.
    Returns the size where the derivative of normalized time changes sign
    (indicating a transition from one algorithm to another).  */
static size_t find_crossover(double * times, size_t * sizes, size_t count,
                             size_t start_idx) {
  size_t i;
  double prev_slope = 0.0;

  for (i = start_idx + 1; i < count - 1; ++i) {
    double t_prev = normalize_time(times[i - 1], sizes[i - 1]);
    double t_curr = normalize_time(times[i], sizes[i]);
    double t_next = normalize_time(times[i + 1], sizes[i + 1]);
    double slope1 = t_curr - t_prev;
    double slope2 = t_next - t_curr;

    /*  Look for significant change in slope (algorithm transition).  */
    if (prev_slope != 0.0 && slope2 < 0.0 && prev_slope > 0.0) {
      /*  Transition from increasing to decreasing normalized time
          suggests a more efficient algorithm has kicked in.  */
      return sizes[i];
    }

    prev_slope = slope1;
  }

  return 0; /* No crossover found */
}

static void print_header(const char * op) {
  if (g_csv_mode) {
    printf("limbs,ns_per_op,ns_per_limb2\n");
  } else {
    printf("\n=== %s Timing Results ===\n", op);
    printf("%8s  %12s  %12s\n", "limbs", "ns/op", "ns/limb^2");
    printf("%8s  %12s  %12s\n", "-----", "-----", "---------");
  }
}

static void print_row(size_t n, double ns) {
  if (g_csv_mode) {
    printf("%zu,%.1f,%.4f\n", n, ns, normalize_time(ns, n));
  } else {
    printf("%8zu  %12.1f  %12.4f\n", n, ns, normalize_time(ns, n));
  }
}

static void tune_multiplication(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t r;
  double times[500];
  size_t sizes[500];
  size_t count = 0;
  size_t n;
  size_t step;
  size_t karatsuba_crossover = 0;
  size_t toom3_crossover = 0;

  if (arbint_ctx_init_default(&ctx) != ARBINT_OK) {
    fprintf(stderr, "Failed to init context\n");
    exit(1);
  }

  arbint_init(a, &ctx);
  arbint_init(b, &ctx);
  arbint_init(r, &ctx);

  print_header("Multiplication");

  for (n = MIN_SIZE; n <= MAX_SIZE;) {
    double ns = measure_mul(a, b, r, n);
    print_row(n, ns);

    if (count < 500) {
      times[count] = ns;
      sizes[count] = n;
      ++count;
    }

    /*  Variable step size.  */
    if (n < 50)
      step = SIZE_STEP_SMALL;
    else if (n < 150)
      step = SIZE_STEP_MEDIUM;
    else
      step = SIZE_STEP_LARGE;
    n += step;
  }

  /*  Analyze for crossover points.  */
  if (!g_csv_mode) {
    karatsuba_crossover = find_crossover(times, sizes, count, 0);
    if (karatsuba_crossover > 0)
      toom3_crossover = find_crossover(times, sizes, count,
                                       karatsuba_crossover / SIZE_STEP_SMALL);

    printf("\n--- Multiplication Analysis ---\n");
    if (karatsuba_crossover > 0)
      printf("Estimated Karatsuba crossover: ~%zu limbs\n",
             karatsuba_crossover);
    else
      printf("Karatsuba crossover: not detected (check smaller sizes)\n");

    if (toom3_crossover > 0)
      printf("Estimated Toom-3 crossover: ~%zu limbs\n", toom3_crossover);
    else
      printf("Toom-3 crossover: not detected (may need larger sizes)\n");

    printf("\nRecommendations:\n");
    printf("  ARBINT_KARATSUBA_THRESHOLD: %zu (current: " ARBINT_STR(
               ARBINT_KARATSUBA_THRESHOLD) ")\n",
           karatsuba_crossover > 0 ? karatsuba_crossover
                                   : ARBINT_KARATSUBA_THRESHOLD);
    printf("  ARBINT_TOOM3_THRESHOLD: %zu (current: " ARBINT_STR(
               ARBINT_TOOM3_THRESHOLD) ")\n",
           toom3_crossover > 0 ? toom3_crossover : ARBINT_TOOM3_THRESHOLD);
  }

  arbint_clear(r);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void tune_squaring(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t r;
  double times[500];
  size_t sizes[500];
  size_t count = 0;
  size_t n;
  size_t step;
  size_t karatsuba_crossover = 0;
  size_t toom3_crossover = 0;

  if (arbint_ctx_init_default(&ctx) != ARBINT_OK) {
    fprintf(stderr, "Failed to init context\n");
    exit(1);
  }

  arbint_init(a, &ctx);
  arbint_init(r, &ctx);

  print_header("Squaring");

  for (n = MIN_SIZE; n <= MAX_SIZE;) {
    double ns = measure_sqr(a, r, n);
    print_row(n, ns);

    if (count < 500) {
      times[count] = ns;
      sizes[count] = n;
      ++count;
    }

    if (n < 50)
      step = SIZE_STEP_SMALL;
    else if (n < 150)
      step = SIZE_STEP_MEDIUM;
    else
      step = SIZE_STEP_LARGE;
    n += step;
  }

  if (!g_csv_mode) {
    karatsuba_crossover = find_crossover(times, sizes, count, 0);
    if (karatsuba_crossover > 0)
      toom3_crossover = find_crossover(times, sizes, count,
                                       karatsuba_crossover / SIZE_STEP_SMALL);

    printf("\n--- Squaring Analysis ---\n");
    if (karatsuba_crossover > 0)
      printf("Estimated Karatsuba crossover: ~%zu limbs\n",
             karatsuba_crossover);
    else
      printf("Karatsuba crossover: not detected (check smaller sizes)\n");

    if (toom3_crossover > 0)
      printf("Estimated Toom-3 crossover: ~%zu limbs\n", toom3_crossover);
    else
      printf("Toom-3 crossover: not detected (may need larger sizes)\n");

    printf("\nRecommendations:\n");
    printf("  ARBINT_SQR_KARATSUBA_THRESHOLD: %zu (current: " ARBINT_STR(
               ARBINT_SQR_KARATSUBA_THRESHOLD) ")\n",
           karatsuba_crossover > 0 ? karatsuba_crossover
                                   : ARBINT_SQR_KARATSUBA_THRESHOLD);
    printf("  ARBINT_SQR_TOOM3_THRESHOLD: %zu (current: " ARBINT_STR(
               ARBINT_SQR_TOOM3_THRESHOLD) ")\n",
           toom3_crossover > 0 ? toom3_crossover : ARBINT_SQR_TOOM3_THRESHOLD);
  }

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Find crossover point in direct comparison data.
    times_a and times_b contain timings for two algorithms.
    Returns the size where algorithm B becomes faster than A.  */
static size_t find_direct_crossover(double * times_a, double * times_b,
                                    size_t * sizes, size_t count) {
  size_t i;
  int b_faster_count = 0;

  for (i = 0; i < count; ++i) {
    if (times_b[i] < times_a[i] * 0.98) {
      /*  B is at least 2% faster.  */
      ++b_faster_count;
      if (b_faster_count >= 2) {
        /*  Found sustained improvement - return first crossover.  */
        return sizes[i - b_faster_count + 1];
      }
    } else {
      b_faster_count = 0;
    }
  }

  return 0;
}

/*  Measure time for NTT multiplication at given size using internal API.
    Returns average nanoseconds per operation.  */
static double measure_ntt_mul_direct(arbint_limb_t * dst, arbint_limb_t * a,
                                     arbint_limb_t * b, size_t n,
                                     const arbint_alloc_t * alloc) {
  uint64_t start;
  uint64_t end;
  uint64_t total_ns;
  int iters;
  int i;
  size_t out_used;

  /*  Warmup.  */
  for (i = 0; i < WARMUP_ITERS; ++i)
    arbint_mul_mag_ntt(dst, &out_used, a, n, b, n, alloc);

  /*  Determine iteration count.  */
  start = get_time_ns();
  for (i = 0; i < MIN_ITERS; ++i)
    arbint_mul_mag_ntt(dst, &out_used, a, n, b, n, alloc);
  end = get_time_ns();

  total_ns = end - start;
  if (total_ns < TARGET_TIME_NS && total_ns > 0) {
    iters = (int) ((TARGET_TIME_NS * MIN_ITERS) / total_ns);
    if (iters < MIN_ITERS)
      iters = MIN_ITERS;
    if (iters > 10000)
      iters = 10000;
  } else {
    iters = MIN_ITERS;
  }

  /*  Actual measurement.  */
  start = get_time_ns();
  for (i = 0; i < iters; ++i)
    arbint_mul_mag_ntt(dst, &out_used, a, n, b, n, alloc);
  end = get_time_ns();

  return (double) (end - start) / (double) iters;
}

/*  Measure time for NTT squaring at given size using internal API.
    Returns average nanoseconds per operation.  */
static double measure_ntt_sqr_direct(arbint_limb_t * dst, arbint_limb_t * a,
                                     size_t n, const arbint_alloc_t * alloc) {
  uint64_t start;
  uint64_t end;
  uint64_t total_ns;
  int iters;
  int i;
  size_t out_used;

  /*  Warmup.  */
  for (i = 0; i < WARMUP_ITERS; ++i)
    arbint_sqr_mag_ntt(dst, &out_used, a, n, alloc);

  /*  Determine iteration count.  */
  start = get_time_ns();
  for (i = 0; i < MIN_ITERS; ++i)
    arbint_sqr_mag_ntt(dst, &out_used, a, n, alloc);
  end = get_time_ns();

  total_ns = end - start;
  if (total_ns < TARGET_TIME_NS && total_ns > 0) {
    iters = (int) ((TARGET_TIME_NS * MIN_ITERS) / total_ns);
    if (iters < MIN_ITERS)
      iters = MIN_ITERS;
    if (iters > 10000)
      iters = 10000;
  } else {
    iters = MIN_ITERS;
  }

  /*  Actual measurement.  */
  start = get_time_ns();
  for (i = 0; i < iters; ++i)
    arbint_sqr_mag_ntt(dst, &out_used, a, n, alloc);
  end = get_time_ns();

  return (double) (end - start) / (double) iters;
}

/*  Measure time for Toom-3 multiplication at given size using internal API.
    Uses arbint_mul_mag_generic which will use Toom-3 for large enough n.
    Returns average nanoseconds per operation.  */
static double measure_toom3_mul_direct(arbint_limb_t * dst, arbint_limb_t * a,
                                       arbint_limb_t * b, size_t n,
                                       const arbint_alloc_t * alloc) {
  uint64_t start;
  uint64_t end;
  uint64_t total_ns;
  int iters;
  int i;
  size_t out_used;

  /*  Warmup.  */
  for (i = 0; i < WARMUP_ITERS; ++i)
    arbint_mul_mag_generic(dst, &out_used, a, n, b, n, alloc);

  /*  Determine iteration count.  */
  start = get_time_ns();
  for (i = 0; i < MIN_ITERS; ++i)
    arbint_mul_mag_generic(dst, &out_used, a, n, b, n, alloc);
  end = get_time_ns();

  total_ns = end - start;
  if (total_ns < TARGET_TIME_NS && total_ns > 0) {
    iters = (int) ((TARGET_TIME_NS * MIN_ITERS) / total_ns);
    if (iters < MIN_ITERS)
      iters = MIN_ITERS;
    if (iters > 10000)
      iters = 10000;
  } else {
    iters = MIN_ITERS;
  }

  /*  Actual measurement.  */
  start = get_time_ns();
  for (i = 0; i < iters; ++i)
    arbint_mul_mag_generic(dst, &out_used, a, n, b, n, alloc);
  end = get_time_ns();

  return (double) (end - start) / (double) iters;
}

/*  Measure time for Toom-3 squaring at given size using internal API.
    Uses arbint_sqr via public API with threshold-based dispatch.
    Since we can't directly call Toom-3 squaring, we rely on the
    existing threshold being high enough that arbint_sqr uses Toom-3.
    Returns average nanoseconds per operation.  */
static double measure_toom3_sqr_via_public(arbint_t r, arbint_t a, size_t n) {
  /*  Just use the public API - for NTT sqr tuning we'll compare against
      direct NTT sqr calls. The public sqr will use Toom-3 if n < threshold.  */
  return measure_sqr(a, r, n);
}

static void tune_ntt(void) {
  arbint_ctx_t ctx;
  arbint_limb_t * a_limbs = NULL;
  arbint_limb_t * b_limbs = NULL;
  arbint_limb_t * dst = NULL;
  double times_toom3[500];
  double times_ntt[500];
  size_t sizes[500];
  size_t count = 0;
  size_t n;
  size_t step;
  size_t ntt_crossover = 0;
  size_t dst_cap;

  if (arbint_ctx_init_default(&ctx) != ARBINT_OK) {
    fprintf(stderr, "Failed to init context\n");
    exit(1);
  }

  if (g_csv_mode) {
    printf("limbs,toom3_ns,ntt_ns,speedup\n");
  } else {
    printf("\n=== NTT Multiplication Threshold Tuning ===\n");
    printf("%8s  %12s  %12s  %8s\n", "limbs", "Toom-3 (ns)", "NTT (ns)",
           "speedup");
    printf("%8s  %12s  %12s  %8s\n", "-----", "-----------", "--------",
           "-------");
  }

  for (n = NTT_MIN_SIZE; n <= NTT_MAX_SIZE;) {
    double ns_toom3;
    double ns_ntt;
    double speedup;

    /*  Allocate buffers for this size.  */
    dst_cap = 2u * n + 1u;
    a_limbs = arbint_alloc_limbs(&ctx.a, n);
    b_limbs = arbint_alloc_limbs(&ctx.a, n);
    dst = arbint_alloc_limbs(&ctx.a, dst_cap);

    if (!a_limbs || !b_limbs || !dst) {
      fprintf(stderr, "Failed to allocate at size %zu\n", n);
      goto cleanup;
    }

    /*  Fill with random data.  */
    for (size_t i = 0; i < n; ++i) {
      a_limbs[i] = (arbint_limb_t) lcg_rand() |
                   ((arbint_limb_t) lcg_rand() << 32);
      b_limbs[i] = (arbint_limb_t) lcg_rand() |
                   ((arbint_limb_t) lcg_rand() << 32);
    }
    /*  Ensure top limb is non-zero.  */
    a_limbs[n - 1] |= ((arbint_limb_t) 1 << 63);
    b_limbs[n - 1] |= ((arbint_limb_t) 1 << 63);

    /*  Measure both algorithms directly.  */
    ns_toom3 = measure_toom3_mul_direct(dst, a_limbs, b_limbs, n, &ctx.a);
    ns_ntt = measure_ntt_mul_direct(dst, a_limbs, b_limbs, n, &ctx.a);
    speedup = ns_toom3 / ns_ntt;

    if (g_csv_mode) {
      printf("%zu,%.1f,%.1f,%.2f\n", n, ns_toom3, ns_ntt, speedup);
    } else {
      printf("%8zu  %12.1f  %12.1f  %8.2fx\n", n, ns_toom3, ns_ntt, speedup);
    }

    if (count < 500) {
      times_toom3[count] = ns_toom3;
      times_ntt[count] = ns_ntt;
      sizes[count] = n;
      ++count;
    }

    arbint_free_limbs(&ctx.a, a_limbs);
    arbint_free_limbs(&ctx.a, b_limbs);
    arbint_free_limbs(&ctx.a, dst);
    a_limbs = b_limbs = dst = NULL;

    /*  Variable step size for NTT range.  */
    if (n < 512)
      step = NTT_SIZE_STEP_SMALL;
    else if (n < 2048)
      step = NTT_SIZE_STEP_MEDIUM;
    else
      step = NTT_SIZE_STEP_LARGE;
    n += step;
  }

  if (!g_csv_mode) {
    ntt_crossover = find_direct_crossover(times_toom3, times_ntt, sizes, count);

    printf("\n--- NTT Multiplication Analysis ---\n");
    printf("Testing range: %d to %d limbs\n", NTT_MIN_SIZE, NTT_MAX_SIZE);
    printf("Current NTT threshold: " ARBINT_STR(ARBINT_NTT_THRESHOLD) "\n");

    if (ntt_crossover > 0)
      printf("Detected NTT crossover: ~%zu limbs\n", ntt_crossover);
    else
      printf("NTT crossover: not detected in test range\n");

    printf("\nRecommendations:\n");
    printf("  ARBINT_NTT_THRESHOLD: %zu (current: " ARBINT_STR(
               ARBINT_NTT_THRESHOLD) ")\n",
           ntt_crossover > 0 ? ntt_crossover : ARBINT_NTT_THRESHOLD);
  }

cleanup:
  arbint_free_limbs(&ctx.a, a_limbs);
  arbint_free_limbs(&ctx.a, b_limbs);
  arbint_free_limbs(&ctx.a, dst);
  arbint_ctx_clear(&ctx);
}

/*  Tune NTT squaring threshold by comparing Toom-3 squaring vs NTT squaring.  */
static void tune_ntt_sqr(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t r;
  arbint_limb_t * a_limbs = NULL;
  arbint_limb_t * dst = NULL;
  double times_toom3[500];
  double times_ntt[500];
  size_t sizes[500];
  size_t count = 0;
  size_t n;
  size_t step;
  size_t ntt_crossover = 0;
  size_t dst_cap;

  if (arbint_ctx_init_default(&ctx) != ARBINT_OK) {
    fprintf(stderr, "Failed to init context\n");
    exit(1);
  }

  arbint_init(a, &ctx);
  arbint_init(r, &ctx);

  if (g_csv_mode) {
    printf("limbs,toom3_ns,ntt_ns,speedup\n");
  } else {
    printf("\n=== NTT Squaring Threshold Tuning ===\n");
    printf("%8s  %12s  %12s  %8s\n", "limbs", "Toom-3 (ns)", "NTT (ns)",
           "speedup");
    printf("%8s  %12s  %12s  %8s\n", "-----", "-----------", "--------",
           "-------");
  }

  for (n = NTT_MIN_SIZE; n <= NTT_MAX_SIZE;) {
    double ns_toom3;
    double ns_ntt;
    double speedup;

    /*  Allocate buffers for this size.  */
    dst_cap = 2u * n + 1u;
    a_limbs = arbint_alloc_limbs(&ctx.a, n);
    dst = arbint_alloc_limbs(&ctx.a, dst_cap);

    if (!a_limbs || !dst) {
      fprintf(stderr, "Failed to allocate at size %zu\n", n);
      goto cleanup;
    }

    /*  Fill with random data.  */
    for (size_t i = 0; i < n; ++i) {
      a_limbs[i] = (arbint_limb_t) lcg_rand() |
                   ((arbint_limb_t) lcg_rand() << 32);
    }
    /*  Ensure top limb is non-zero.  */
    a_limbs[n - 1] |= ((arbint_limb_t) 1 << 63);

    /*  Measure Toom-3 squaring via public API.
        (arbint_sqr will use Toom-3 since we're below the NTT threshold
        or we're directly measuring via mul_mag_generic for squaring.)  */
    build_random_value(a, n);
    ns_toom3 = measure_toom3_sqr_via_public(r, a, n);

    /*  Measure NTT squaring directly.  */
    ns_ntt = measure_ntt_sqr_direct(dst, a_limbs, n, &ctx.a);
    speedup = ns_toom3 / ns_ntt;

    if (g_csv_mode) {
      printf("%zu,%.1f,%.1f,%.2f\n", n, ns_toom3, ns_ntt, speedup);
    } else {
      printf("%8zu  %12.1f  %12.1f  %8.2fx\n", n, ns_toom3, ns_ntt, speedup);
    }

    if (count < 500) {
      times_toom3[count] = ns_toom3;
      times_ntt[count] = ns_ntt;
      sizes[count] = n;
      ++count;
    }

    arbint_free_limbs(&ctx.a, a_limbs);
    arbint_free_limbs(&ctx.a, dst);
    a_limbs = dst = NULL;

    /*  Variable step size for NTT range.  */
    if (n < 512)
      step = NTT_SIZE_STEP_SMALL;
    else if (n < 2048)
      step = NTT_SIZE_STEP_MEDIUM;
    else
      step = NTT_SIZE_STEP_LARGE;
    n += step;
  }

  if (!g_csv_mode) {
    ntt_crossover = find_direct_crossover(times_toom3, times_ntt, sizes, count);

    printf("\n--- NTT Squaring Analysis ---\n");
    printf("Testing range: %d to %d limbs\n", NTT_MIN_SIZE, NTT_MAX_SIZE);
    printf("Current NTT threshold (shared with mul): " ARBINT_STR(
        ARBINT_NTT_THRESHOLD) "\n");

    if (ntt_crossover > 0)
      printf("Detected NTT squaring crossover: ~%zu limbs\n", ntt_crossover);
    else
      printf("NTT squaring crossover: not detected in test range\n");

    printf("\nNote: NTT squaring uses the same threshold as NTT multiplication.\n");
    printf("If the squaring crossover differs significantly from multiplication,\n");
    printf("consider adding a separate ARBINT_NTT_SQR_THRESHOLD.\n");

    printf("\nRecommendations:\n");
    printf("  ARBINT_NTT_THRESHOLD (for sqr): %zu (current: " ARBINT_STR(
               ARBINT_NTT_THRESHOLD) ")\n",
           ntt_crossover > 0 ? ntt_crossover : ARBINT_NTT_THRESHOLD);
  }

cleanup:
  arbint_free_limbs(&ctx.a, a_limbs);
  arbint_free_limbs(&ctx.a, dst);
  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void print_usage(const char * prog) {
  printf("Usage: %s [OPTIONS]\n", prog);
  printf("\nOptions:\n");
  printf("  --mul      Tune multiplication thresholds (Karatsuba, Toom-3)\n");
  printf("  --sqr      Tune squaring thresholds (Karatsuba, Toom-3)\n");
  printf("  --ntt      Tune NTT multiplication threshold (Toom-3 vs NTT)\n");
  printf("  --ntt-sqr  Tune NTT squaring threshold (Toom-3 vs NTT squaring)\n");
  printf("  --csv      Output in CSV format (for plotting)\n");
  printf("  --help     Show this help message\n");
  printf("\nBy default, tunes all thresholds (mul, sqr, ntt, ntt-sqr).\n");
  printf("\nThe --ntt and --ntt-sqr options directly measure Toom-3 vs NTT\n");
  printf("by calling internal functions, providing accurate crossover detection.\n");
  printf("\nUse the recommendations to update thresholds in src/arbint_mul.h.\n");
}

int main(int argc, char ** argv) {
  int do_mul = 0;
  int do_sqr = 0;
  int do_ntt = 0;
  int do_ntt_sqr = 0;
  int i;
  arbint_err_t rc;

  /*  Initialize RNG with platform entropy.  */
  rc = arbint_rng_init(&g_rng, NULL, 0);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "Failed to initialize RNG: %d\n", rc);
    return 1;
  }

  for (i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--mul") == 0) {
      do_mul = 1;
    } else if (strcmp(argv[i], "--sqr") == 0) {
      do_sqr = 1;
    } else if (strcmp(argv[i], "--ntt") == 0) {
      do_ntt = 1;
    } else if (strcmp(argv[i], "--ntt-sqr") == 0) {
      do_ntt_sqr = 1;
    } else if (strcmp(argv[i], "--csv") == 0) {
      g_csv_mode = 1;
    } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
      print_usage(argv[0]);
      return 0;
    } else {
      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      print_usage(argv[0]);
      return 1;
    }
  }

  /*  Default: tune all.  */
  if (!do_mul && !do_sqr && !do_ntt && !do_ntt_sqr) {
    do_mul = 1;
    do_sqr = 1;
    do_ntt = 1;
    do_ntt_sqr = 1;
  }

  if (!g_csv_mode) {
    printf("arbint Threshold Tuning Tool\n");
    printf("============================\n");
    printf("\nMeasuring performance at various operand sizes...\n");
    printf("(This may take a few minutes)\n");
  }

  if (do_mul)
    tune_multiplication();

  if (do_sqr)
    tune_squaring();

  if (do_ntt)
    tune_ntt();

  if (do_ntt_sqr)
    tune_ntt_sqr();

  if (!g_csv_mode) {
    printf("\n============================\n");
    printf("Tuning complete.\n");
    printf("\nTo apply recommendations:\n");
    printf("  1. Edit src/arbint_mul.h\n");
    printf("  2. Update threshold #defines\n");
    printf("  3. Recompile: make clean && make\n");
    printf("  4. Re-run this tool to verify\n");
  }

  arbint_rng_clear(&g_rng);
  return 0;
}
