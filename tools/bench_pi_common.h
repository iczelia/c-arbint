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

/*  Shared definitions for the Chudnovsky pi benchmark suite.

    The Chudnovsky algorithm:
      1/pi = 12 * sum(k=0 to oo) [(-1)^k (6k)! (13591409 + 545140134k)] /
                              [(3k)! (k!)^3 (640320)^(3k + 3/2)]

    Each term adds approximately 14.18 decimal digits of precision.

    Implementation uses the recurrence:
      a_{k+1} = a_k * -(6k-5)(2k-1)(6k-1) / (k^3 * C3_OVER_24)
    where C = 640320 and C3_OVER_24 = 640320^3 / 24 = 10939058860032000.  */

#ifndef BENCH_PI_COMMON_H
#define BENCH_PI_COMMON_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*  Default digits of pi to compute.  */
#define DEFAULT_PI_DIGITS 1000000

/*  Chudnovsky constants.  */
#define CHUD_A  13591409
#define CHUD_B  545140134
#define CHUD_C  426880
#define CHUD_D  10005
/*  C3_OVER_24 = 640320^3 / 24 = 10939058860032000
    Split for 32-bit construction: 10939058 * 10^9 + 860032000  */
#define CHUD_C3_24_HI  10939058u
#define CHUD_C3_24_LO  860032000u

/*  Platform-specific high-resolution timing.  */
#if defined(_WIN32)
  #define WIN32_LEAN_AND_MEAN
  #include <windows.h>
  static inline uint64_t bench_get_time_ns(void) {
    LARGE_INTEGER freq, count;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&count);
    return (uint64_t) ((count.QuadPart * 1000000000ULL) / freq.QuadPart);
  }
#elif defined(__APPLE__)
  #include <time.h>
  static inline uint64_t bench_get_time_ns(void) {
    return clock_gettime_nsec_np(CLOCK_MONOTONIC);
  }
#else
  #include <time.h>
  static inline uint64_t bench_get_time_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t) ts.tv_sec * 1000000000ULL + (uint64_t) ts.tv_nsec;
  }
#endif

typedef struct {
  uint64_t series_ns;
  uint64_t sqrt_ns;
  uint64_t final_ns;
  uint64_t total_ns;
  uint32_t terms;
} bench_pi_timing_t;

typedef struct {
  char * digits;           /*  Computed pi digits (malloc'd, caller frees).  */
  size_t digit_count;      /*  Number of digits.  */
  bench_pi_timing_t timing;
} bench_pi_result_t;

static inline void bench_print_timing(const char * name,
                                       const bench_pi_timing_t * t) {
  double series_s = (double) t->series_ns / 1e9;
  double sqrt_s = (double) t->sqrt_ns / 1e9;
  double final_s = (double) t->final_ns / 1e9;
  double total_s = (double) t->total_ns / 1e9;

  printf("\n=== %s Timing ===\n", name);
  printf("  Series computation: %8.3f s  (%u terms)\n", series_s, t->terms);
  printf("  Square root:        %8.3f s\n", sqrt_s);
  printf("  Final division:     %8.3f s\n", final_s);
  printf("  Total:              %8.3f s\n", total_s);
}

/*  Compare two digit strings and report differences.  */
static inline int bench_compare_digits(const char * a, const char * b,
                                        size_t len, const char * name_a,
                                        const char * name_b) {
  size_t i;

  for (i = 0; i < len; ++i) {
    if (a[i] != b[i]) {
      printf("MISMATCH at digit %zu: %s='%c', %s='%c'\n", i, name_a, a[i],
             name_b, b[i]);
      /*  Show context.  */
      size_t start = (i > 10) ? i - 10 : 0;
      size_t end = (i + 10 < len) ? i + 10 : len;
      printf("  %s: ", name_a);
      for (size_t j = start; j < end; ++j)
        putchar(a[j]);
      printf("\n  %s: ", name_b);
      for (size_t j = start; j < end; ++j)
        putchar(b[j]);
      printf("\n");
      return 0;
    }
  }
  return 1;
}

/*  Load reference pi digits from file.
    Returns malloc'd string with digits only (no '.'), or NULL.  */
static inline char * bench_load_reference(const char * path) {
  FILE * f = fopen(path, "r");
  if (!f)
    return NULL;

  fseek(f, 0, SEEK_END);
  long fsize = ftell(f);
  fseek(f, 0, SEEK_SET);

  if (fsize < 3) {
    fclose(f);
    return NULL;
  }

  char * buf = (char *) malloc((size_t) fsize + 1);
  if (!buf) {
    fclose(f);
    return NULL;
  }

  size_t nread = fread(buf, 1, (size_t) fsize, f);
  fclose(f);
  buf[nread] = '\0';

  /*  Strip trailing whitespace.  */
  while (nread > 0 && (buf[nread - 1] == '\n' || buf[nread - 1] == '\r' ||
                       buf[nread - 1] == ' ')) {
    buf[--nread] = '\0';
  }

  /*  Expect "3." prefix.  */
  if (buf[0] != '3' || buf[1] != '.') {
    free(buf);
    return NULL;
  }

  /*  Build digit string: "3" + fractional digits.  */
  size_t digit_len = nread - 1;
  char * digits = (char *) malloc(digit_len + 1);
  if (!digits) {
    free(buf);
    return NULL;
  }

  digits[0] = '3';
  memcpy(digits + 1, buf + 2, nread - 2);
  digits[digit_len] = '\0';

  free(buf);
  return digits;
}

#endif /*  BENCH_PI_COMMON_H  */
