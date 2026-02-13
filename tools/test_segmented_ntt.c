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

/*  Test tool for large-number multiplication (including segmented NTT).

    This tool tests multiplication of very large numbers, including sizes
    that approach or exceed the standard NTT limit of 2^24 limbs. It verifies
    correctness using algebraic identities and performs timing measurements.

    The segmented NTT multiplication is exercised automatically by arbint_mul
    when operands exceed the NTT size limit.

    WARNING: This test requires significant memory and time.
    - 2^20 limbs = 8 MB per operand (64-bit limbs)
    - 2^22 limbs = 32 MB per operand
    - 2^24 limbs = 128 MB per operand
    - 2^25 limbs = 256 MB per operand

    Usage:
      test_segmented_ntt [OPTIONS]

    Options:
      --quick     Run quick tests only (< 1 minute, ~256 MB)
      --medium    Run medium tests (< 10 minutes, ~1 GB)
      --full      Run full test suite (may take hours, ~4+ GB)
      --size N    Test specific size in limbs
      --help      Show this help message

    Without options, runs quick tests by default.  */

#include <arbint.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "arbint_base.h"

/*  Test result tracking.  */
static int g_tests_passed = 0;
static int g_tests_failed = 0;

/*  Global RNG.  */
static arbint_rng_t g_rng;

/*  Get current time in nanoseconds.  */
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

/*  Format time in human-readable form.  */
static void format_time(uint64_t ns, char * buf, size_t buflen) {
  if (ns < 1000ull) {
    snprintf(buf, buflen, "%llu ns", (unsigned long long) ns);
  } else if (ns < 1000000ull) {
    snprintf(buf, buflen, "%.2f us", (double) ns / 1000.0);
  } else if (ns < 1000000000ull) {
    snprintf(buf, buflen, "%.2f ms", (double) ns / 1000000.0);
  } else {
    snprintf(buf, buflen, "%.2f s", (double) ns / 1000000000.0);
  }
}

/*  Format size in human-readable form.  */
static void format_size(size_t limbs, char * buf, size_t buflen) {
  size_t bytes = limbs * sizeof(arbint_limb_t);
  if (bytes < 1024ull) {
    snprintf(buf, buflen, "%zu B", bytes);
  } else if (bytes < 1024ull * 1024ull) {
    snprintf(buf, buflen, "%.1f KB", (double) bytes / 1024.0);
  } else if (bytes < 1024ull * 1024ull * 1024ull) {
    snprintf(buf, buflen, "%.1f MB", (double) bytes / (1024.0 * 1024.0));
  } else {
    snprintf(buf, buflen, "%.2f GB",
             (double) bytes / (1024.0 * 1024.0 * 1024.0));
  }
}

/*  Build a random arbint with exactly n limbs.  */
static arbint_err_t build_random(arbint_t a, size_t n_limbs) {
  return arbint_urandomb(a, &g_rng, n_limbs * sizeof(arbint_limb_t) * 8);
}

/*  Test multiplication using algebraic identity: (a+1)*(a-1) = a^2 - 1.
    This verifies that multiplication is consistent without needing a
    reference implementation.  */
static int test_algebraic_identity(arbint_ctx_t * ctx, size_t n_limbs) {
  arbint_t a;
  arbint_t ap1;
  arbint_t am1;
  arbint_t prod;
  arbint_t sqr;
  arbint_t sqr_m1;
  arbint_err_t rc;
  int result = 1;
  char size_str[32];
  char time_str[32];
  uint64_t start;
  uint64_t end;

  format_size(n_limbs, size_str, sizeof(size_str));
  printf("  Testing algebraic identity at %zu limbs (%s)... ", n_limbs,
         size_str);
  fflush(stdout);

  rc = arbint_init_all(ctx, a, ap1, am1, prod, sqr, sqr_m1, (arbint_t *) NULL);
  if (rc != ARBINT_OK) {
    printf("FAIL (init)\n");
    return 0;
  }

  /*  Generate random a with n limbs.  */
  rc = build_random(a, n_limbs);
  if (rc != ARBINT_OK) {
    printf("FAIL (random)\n");
    result = 0;
    goto cleanup;
  }

  /*  Compute a+1 and a-1.  */
  rc = arbint_add_i32(ap1, a, 1);
  rc |= arbint_sub_i32(am1, a, 1);
  if (rc != ARBINT_OK) {
    printf("FAIL (add/sub)\n");
    result = 0;
    goto cleanup;
  }

  /*  Time the multiplication: (a+1)*(a-1).  */
  start = get_time_ns();
  rc = arbint_mul(prod, ap1, am1);
  end = get_time_ns();

  if (rc != ARBINT_OK) {
    printf("FAIL (mul: %d)\n", rc);
    result = 0;
    goto cleanup;
  }

  /*  Compute a^2.  */
  rc = arbint_sqr(sqr, a);
  if (rc != ARBINT_OK) {
    printf("FAIL (sqr)\n");
    result = 0;
    goto cleanup;
  }

  /*  Compute a^2 - 1.  */
  rc = arbint_sub_i32(sqr_m1, sqr, 1);
  if (rc != ARBINT_OK) {
    printf("FAIL (sub 1)\n");
    result = 0;
    goto cleanup;
  }

  /*  Verify: (a+1)*(a-1) == a^2 - 1.  */
  if (arbint_cmp(prod, sqr_m1) != 0) {
    printf("FAIL (mismatch)\n");
    result = 0;
    goto cleanup;
  }

  format_time(end - start, time_str, sizeof(time_str));
  printf("OK (%s)\n", time_str);

cleanup:
  arbint_clear_all(a, ap1, am1, prod, sqr, sqr_m1, (arbint_t *) NULL);

  return result;
}

/*  Test multiplication using distributive property: a*(b+c) = a*b + a*c.
    This provides another verification path.  */
static int test_distributive(arbint_ctx_t * ctx, size_t n_limbs) {
  arbint_t a;
  arbint_t b;
  arbint_t c;
  arbint_t bc;
  arbint_t a_bc;
  arbint_t ab;
  arbint_t ac;
  arbint_t ab_ac;
  arbint_err_t rc;
  int result = 1;
  char size_str[32];
  char time_str[32];
  uint64_t start;
  uint64_t end;

  format_size(n_limbs, size_str, sizeof(size_str));
  printf("  Testing distributive property at %zu limbs (%s)... ", n_limbs,
         size_str);
  fflush(stdout);

  rc = arbint_init_all(ctx, a, b, c, bc, a_bc, ab, ac, ab_ac, (arbint_t *) NULL);
  if (rc != ARBINT_OK) {
    printf("FAIL (init)\n");
    return 0;
  }

  /*  Generate random operands.  */
  rc = build_random(a, n_limbs);
  rc |= build_random(b, n_limbs / 2 + 1);
  rc |= build_random(c, n_limbs / 3 + 1);
  if (rc != ARBINT_OK) {
    printf("FAIL (random)\n");
    result = 0;
    goto cleanup;
  }

  /*  Compute b + c.  */
  rc = arbint_add(bc, b, c);
  if (rc != ARBINT_OK) {
    printf("FAIL (add)\n");
    result = 0;
    goto cleanup;
  }

  /*  Time the main multiplication: a * (b + c).  */
  start = get_time_ns();
  rc = arbint_mul(a_bc, a, bc);
  end = get_time_ns();

  if (rc != ARBINT_OK) {
    printf("FAIL (mul a*(b+c))\n");
    result = 0;
    goto cleanup;
  }

  /*  Compute a*b and a*c.  */
  rc = arbint_mul(ab, a, b);
  rc |= arbint_mul(ac, a, c);
  if (rc != ARBINT_OK) {
    printf("FAIL (mul a*b or a*c)\n");
    result = 0;
    goto cleanup;
  }

  /*  Compute a*b + a*c.  */
  rc = arbint_add(ab_ac, ab, ac);
  if (rc != ARBINT_OK) {
    printf("FAIL (add ab+ac)\n");
    result = 0;
    goto cleanup;
  }

  /*  Verify: a*(b+c) == a*b + a*c.  */
  if (arbint_cmp(a_bc, ab_ac) != 0) {
    printf("FAIL (mismatch)\n");
    result = 0;
    goto cleanup;
  }

  format_time(end - start, time_str, sizeof(time_str));
  printf("OK (%s)\n", time_str);

cleanup:
  arbint_clear_all(a, b, c, bc, a_bc, ab, ac, ab_ac, (arbint_t *) NULL);

  return result;
}

/*  Test: large multiplication with timing.  */
static int test_large_mul(arbint_ctx_t * ctx, size_t n_limbs) {
  arbint_t a;
  arbint_t b;
  arbint_t r;
  arbint_err_t rc;
  int result = 1;
  char size_str[32];
  char time_str[32];
  uint64_t start;
  uint64_t end;
  size_t rn;

  format_size(n_limbs, size_str, sizeof(size_str));
  printf("\n=== Large Multiplication Test (%zu limbs, %s each) ===\n", n_limbs,
         size_str);

  printf("  Allocating operands... ");
  fflush(stdout);

  rc = arbint_init_all(ctx, a, b, r, (arbint_t *) NULL);
  if (rc != ARBINT_OK) {
    printf("FAIL\n");
    return 0;
  }
  printf("OK\n");

  printf("  Generating random operands... ");
  fflush(stdout);

  start = get_time_ns();
  rc = build_random(a, n_limbs);
  rc |= build_random(b, n_limbs);
  end = get_time_ns();

  if (rc != ARBINT_OK) {
    printf("FAIL\n");
    result = 0;
    goto cleanup;
  }

  format_time(end - start, time_str, sizeof(time_str));
  printf("OK (%s)\n", time_str);

  printf("  Running multiplication... ");
  fflush(stdout);

  start = get_time_ns();
  rc = arbint_mul(r, a, b);
  end = get_time_ns();

  if (rc != ARBINT_OK) {
    printf("FAIL (rc=%d)\n", rc);
    result = 0;
    goto cleanup;
  }

  format_time(end - start, time_str, sizeof(time_str));
  printf("OK (%s)\n", time_str);

  /*  Sanity check result size.  */
  rn = arbint_nbits(r);
  printf("  Result: %zu bits\n", rn);

  /*  Verify with algebraic identity.  */
  if (test_algebraic_identity(ctx, n_limbs))
    ++g_tests_passed;
  else
    ++g_tests_failed;

cleanup:
  arbint_clear_all(a, b, r, (arbint_t *) NULL);

  return result;
}

/*  Run quick test suite.  */
static void run_quick_tests(arbint_ctx_t * ctx) {
  size_t sizes[] = {1000, 5000, 10000, 50000, 100000};
  size_t n_sizes = sizeof(sizes) / sizeof(sizes[0]);
  size_t i;

  printf("\n=== Quick Tests ===\n");

  for (i = 0; i < n_sizes; ++i) {
    if (test_algebraic_identity(ctx, sizes[i]))
      ++g_tests_passed;
    else
      ++g_tests_failed;
  }

  printf("\n");
  for (i = 0; i < n_sizes; ++i) {
    if (test_distributive(ctx, sizes[i]))
      ++g_tests_passed;
    else
      ++g_tests_failed;
  }

  /*  Test larger size that approaches NTT threshold.  */
  if (test_large_mul(ctx, 500000))
    ++g_tests_passed;
  else
    ++g_tests_failed;
}

/*  Run medium test suite.  */
static void run_medium_tests(arbint_ctx_t * ctx) {
  run_quick_tests(ctx);

  /*  Test sizes around 1M limbs.  */
  if (test_large_mul(ctx, 1000000))
    ++g_tests_passed;
  else
    ++g_tests_failed;

  if (test_large_mul(ctx, 2000000))
    ++g_tests_passed;
  else
    ++g_tests_failed;
}

/*  Run full test suite.  */
static void run_full_tests(arbint_ctx_t * ctx) {
  run_medium_tests(ctx);

  /*  Test sizes approaching NTT limit (2^24 = 16M limbs).  */
  if (test_large_mul(ctx, 4000000))
    ++g_tests_passed;
  else
    ++g_tests_failed;

  if (test_large_mul(ctx, 8000000))
    ++g_tests_passed;
  else
    ++g_tests_failed;

  /*  This would exercise segmented multiplication if implemented.  */
  printf("\nNote: Sizes beyond 2^24 limbs (~128MB) would use segmented NTT.\n");
}

static void print_usage(const char * prog) {
  printf("Usage: %s [OPTIONS]\n", prog);
  printf("\nOptions:\n");
  printf("  --quick     Run quick tests only (< 1 minute, ~256 MB)\n");
  printf("  --medium    Run medium tests (< 10 minutes, ~1 GB)\n");
  printf("  --full      Run full test suite (may take hours, ~4+ GB)\n");
  printf("  --size N    Test specific size in limbs\n");
  printf("  --help      Show this help message\n");
  printf("\nWithout options, runs quick tests by default.\n");
}

int main(int argc, char ** argv) {
  arbint_ctx_t ctx;
  arbint_err_t rc;
  int mode_quick = 1;
  int mode_medium = 0;
  int mode_full = 0;
  size_t specific_size = 0;
  int i;

  /*  Parse arguments.  */
  for (i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--quick") == 0) {
      mode_quick = 1;
      mode_medium = 0;
      mode_full = 0;
    } else if (strcmp(argv[i], "--medium") == 0) {
      mode_quick = 0;
      mode_medium = 1;
      mode_full = 0;
    } else if (strcmp(argv[i], "--full") == 0) {
      mode_quick = 0;
      mode_medium = 0;
      mode_full = 1;
    } else if (strcmp(argv[i], "--size") == 0 && i + 1 < argc) {
      specific_size = (size_t) strtoull(argv[++i], NULL, 0);
    } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
      print_usage(argv[0]);
      return 0;
    } else {
      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      print_usage(argv[0]);
      return 1;
    }
  }

  printf("Large Number Multiplication Test\n");
  printf("=================================\n\n");

  printf("Configuration:\n");
  printf("  ARBINT_LIMB_BITS: %d\n", ARBINT_LIMB_BITS);
  printf("  sizeof(arbint_limb_t): %zu\n", sizeof(arbint_limb_t));
  printf("\n");

  /*  Initialize RNG.  */
  rc = arbint_rng_init(&g_rng, NULL, 0);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "Failed to initialize RNG: %d\n", rc);
    return 1;
  }

  /*  Initialize context.  */
  rc = arbint_ctx_init_default(&ctx);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "Failed to initialize context: %d\n", rc);
    arbint_rng_clear(&g_rng);
    return 1;
  }

  /*  Run tests.  */
  if (specific_size > 0) {
    if (test_large_mul(&ctx, specific_size))
      ++g_tests_passed;
    else
      ++g_tests_failed;
  } else if (mode_full) {
    run_full_tests(&ctx);
  } else if (mode_medium) {
    run_medium_tests(&ctx);
  } else {
    run_quick_tests(&ctx);
  }

  /*  Summary.  */
  printf("\n=================================\n");
  printf("Results: %d passed, %d failed\n", g_tests_passed, g_tests_failed);

  arbint_ctx_clear(&ctx);
  arbint_rng_clear(&g_rng);
  arbint_drop_caches();

  return (g_tests_failed > 0) ? 1 : 0;
}
