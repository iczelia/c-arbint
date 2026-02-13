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

/*  Chudnovsky pi benchmark suite: c-arbint vs GMP.

    Computes pi to a configurable number of digits using identical
    Chudnovsky implementations for both libraries, timing each phase.

    Usage: bench_pi [OPTIONS]

    Options:
      -d, --digits N       Number of digits to compute (default: 1000000)
      -r, --reference FILE Use FILE as reference for verification
      -a, --arbint-only    Run only arbint benchmark
      -g, --gmp-only       Run only GMP benchmark
      -q, --quiet          Suppress progress output
      -h, --help           Show this help message  */

#include "bench_pi_common.h"

#include <getopt.h>

/*  External benchmark functions.  */
extern int bench_pi_arbint(uint32_t digits, bench_pi_result_t * result);
extern int bench_pi_gmp(uint32_t digits, bench_pi_result_t * result);

static void print_usage(const char * prog) {
  printf("Usage: %s [OPTIONS]\n\n", prog);
  printf("Chudnovsky pi benchmark: c-arbint vs GMP\n\n");
  printf("Options:\n");
  printf("  -d, --digits N       Number of digits to compute (default: %d)\n",
         DEFAULT_PI_DIGITS);
  printf("  -r, --reference FILE Use FILE as reference for verification\n");
  printf("  -a, --arbint-only    Run only arbint benchmark\n");
  printf("  -g, --gmp-only       Run only GMP benchmark\n");
  printf("  -q, --quiet          Suppress progress output\n");
  printf("  -h, --help           Show this help message\n\n");
  printf("Examples:\n");
  printf("  %s -d 100000          Compute 100k digits with both libraries\n",
         prog);
  printf("  %s -d 1000000 -a      Compute 1M digits with arbint only\n", prog);
  printf("  %s -r pi10k.txt -d 10000  Verify against reference file\n", prog);
}

int main(int argc, char ** argv) {
  uint32_t digits = DEFAULT_PI_DIGITS;
  const char * reference_path = NULL;
  int run_arbint = 1;
  int run_gmp = 1;
  int quiet = 0;
  int opt;

  static struct option long_options[] = {{"digits", required_argument, 0, 'd'},
                                         {"reference", required_argument, 0, 'r'},
                                         {"arbint-only", no_argument, 0, 'a'},
                                         {"gmp-only", no_argument, 0, 'g'},
                                         {"quiet", no_argument, 0, 'q'},
                                         {"help", no_argument, 0, 'h'},
                                         {0, 0, 0, 0}};

  while ((opt = getopt_long(argc, argv, "d:r:agqh", long_options, NULL)) != -1) {
    switch (opt) {
    case 'd':
      digits = (uint32_t) atoi(optarg);
      if (digits < 10) {
        fprintf(stderr, "Error: digits must be at least 10\n");
        return 1;
      }
      break;
    case 'r':
      reference_path = optarg;
      break;
    case 'a':
      run_arbint = 1;
      run_gmp = 0;
      break;
    case 'g':
      run_arbint = 0;
      run_gmp = 1;
      break;
    case 'q':
      quiet = 1;
      break;
    case 'h':
      print_usage(argv[0]);
      return 0;
    default:
      print_usage(argv[0]);
      return 1;
    }
  }

  printf("================================================================\n");
  printf("  Chudnovsky Pi Benchmark: c-arbint vs GMP\n");
  printf("================================================================\n");
  printf("\nComputing %u digits of pi...\n", digits);

  bench_pi_result_t result_arbint = {0};
  bench_pi_result_t result_gmp = {0};
  char * reference = NULL;

  /*  Load reference if specified.  */
  if (reference_path) {
    reference = bench_load_reference(reference_path);
    if (!reference) {
      fprintf(stderr, "Warning: Could not load reference file: %s\n",
              reference_path);
    }
  }

  /*  Run arbint benchmark.  */
  if (run_arbint) {
    if (!quiet)
      printf("\n[arbint] Starting computation...\n");

    if (!bench_pi_arbint(digits, &result_arbint)) {
      fprintf(stderr, "Error: arbint computation failed\n");
      free(reference);
      return 1;
    }

    bench_print_timing("arbint", &result_arbint.timing);

    if (!quiet)
      printf("  Computed %zu digits\n", result_arbint.digit_count);
  }

  /*  Run GMP benchmark.  */
  if (run_gmp) {
    if (!quiet)
      printf("\n[GMP] Starting computation...\n");

    if (!bench_pi_gmp(digits, &result_gmp)) {
      fprintf(stderr, "Error: GMP computation failed\n");
      free(result_arbint.digits);
      free(reference);
      return 1;
    }

    bench_print_timing("GMP", &result_gmp.timing);

    if (!quiet)
      printf("  Computed %zu digits\n", result_gmp.digit_count);
  }

  /*  Compare results.  */
  printf("\n================================================================\n");
  printf("  Results Summary\n");
  printf("================================================================\n\n");

  if (run_arbint && run_gmp) {
    /*  Compare arbint vs GMP.  */
    size_t compare_len = digits + 1;
    if (result_arbint.digit_count < compare_len)
      compare_len = result_arbint.digit_count;
    if (result_gmp.digit_count < compare_len)
      compare_len = result_gmp.digit_count;

    printf("Comparing first %zu digits...\n", compare_len);

    if (bench_compare_digits(result_arbint.digits, result_gmp.digits,
                              compare_len, "arbint", "GMP")) {
      printf("OK: arbint and GMP results match!\n");
    } else {
      printf("FAIL: Results differ!\n");
    }

    /*  Print speedup.  */
    printf("\nPerformance comparison:\n");
    double arbint_s = (double) result_arbint.timing.total_ns / 1e9;
    double gmp_s = (double) result_gmp.timing.total_ns / 1e9;
    double ratio = arbint_s / gmp_s;

    printf("  arbint: %.3f s\n", arbint_s);
    printf("  GMP:    %.3f s\n", gmp_s);

    if (ratio > 1.0) {
      printf("  GMP is %.2fx faster than arbint\n", ratio);
    } else {
      printf("  arbint is %.2fx faster than GMP\n", 1.0 / ratio);
    }

    /*  Breakdown comparison.  */
    printf("\nPhase breakdown:\n");
    printf("  %-20s  %12s  %12s  %8s\n", "Phase", "arbint (s)", "GMP (s)",
           "Ratio");
    printf("  %-20s  %12s  %12s  %8s\n", "--------------------", "----------",
           "----------", "------");

    double a_series = (double) result_arbint.timing.series_ns / 1e9;
    double g_series = (double) result_gmp.timing.series_ns / 1e9;
    printf("  %-20s  %12.3f  %12.3f  %8.2fx\n", "Series computation", a_series,
           g_series, a_series / g_series);

    double a_sqrt = (double) result_arbint.timing.sqrt_ns / 1e9;
    double g_sqrt = (double) result_gmp.timing.sqrt_ns / 1e9;
    printf("  %-20s  %12.3f  %12.3f  %8.2fx\n", "Square root", a_sqrt, g_sqrt,
           a_sqrt / g_sqrt);

    double a_final = (double) result_arbint.timing.final_ns / 1e9;
    double g_final = (double) result_gmp.timing.final_ns / 1e9;
    printf("  %-20s  %12.3f  %12.3f  %8.2fx\n", "Final division", a_final,
           g_final, a_final / g_final);
  }

  /*  Verify against reference if available.  */
  if (reference) {
    size_t ref_len = strlen(reference);
    printf("\nReference verification:\n");

    if (run_arbint) {
      size_t check_len = (digits + 1 < ref_len) ? digits + 1 : ref_len;
      if (result_arbint.digit_count < check_len)
        check_len = result_arbint.digit_count;

      printf("  Checking arbint against reference (%zu digits)... ", check_len);
      if (bench_compare_digits(result_arbint.digits, reference, check_len,
                                "arbint", "reference")) {
        printf("OK\n");
      } else {
        printf("FAIL\n");
      }
    }

    if (run_gmp) {
      size_t check_len = (digits + 1 < ref_len) ? digits + 1 : ref_len;
      if (result_gmp.digit_count < check_len)
        check_len = result_gmp.digit_count;

      printf("  Checking GMP against reference (%zu digits)... ", check_len);
      if (bench_compare_digits(result_gmp.digits, reference, check_len, "GMP",
                                "reference")) {
        printf("OK\n");
      } else {
        printf("FAIL\n");
      }
    }
  }

  /*  Show first and last few digits.  */
  if (run_arbint && result_arbint.digit_count > 0) {
    printf("\nFirst 50 digits (arbint): ");
    size_t show = result_arbint.digit_count < 50 ? result_arbint.digit_count : 50;
    for (size_t i = 0; i < show; ++i)
      putchar(result_arbint.digits[i]);
    printf("...\n");

    if (result_arbint.digit_count > 50) {
      printf("Last 50 digits (arbint):  ...");
      size_t start = result_arbint.digit_count - 50;
      for (size_t i = start; i < result_arbint.digit_count; ++i)
        putchar(result_arbint.digits[i]);
      printf("\n");
    }
  }

  /*  Cleanup.  */
  free(result_arbint.digits);
  free(result_gmp.digits);
  free(reference);

  printf("\n================================================================\n");
  printf("  Benchmark complete.\n");
  printf("================================================================\n");

  return 0;
}
