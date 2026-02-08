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

/*  Calculate pi to 10000 digits using the Chudnovsky algorithm and
    compare against known-good reference digits in pi10k.txt.

    This variant purposefully relies only on the most basic functionality
    of the library to test the core features and memory management in
    particular. It might be useful to bisect issues that appear
    in the kernels of basic operations.

    The Chudnovsky algorithm:
      1/pi = 12 * sum(k=0 to oo) [(-1)^k (6k)! (13591409 + 545140134k)] /
                              [(3k)! (k!)^3 (640320)^(3k + 3/2)]

    Each term adds approximately 14.18 decimal digits of precision.  */

#include "test_framework.h"

#include <arbint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI_DIGITS 10000

/*  Known-good helpers independent of arbint.  */
static char * arbint_to_decimal_string(const arbint_t x, arbint_ctx_t * ctx) {
  arbint_t tmp, q, r;
  size_t capacity = 128;
  size_t len = 0;
  char * str = NULL;
  char * digits = NULL;
  size_t i;
  int is_negative = arbint_is_neg(x);

  if (arbint_is_zero(x)) {
    str = malloc(2);
    if (!str)
      return NULL;
    strcpy(str, "0");
    return str;
  }

  /* Allocate temporary integers */
  if (arbint_init_all(ctx, tmp, q, r, (arbint_t *) NULL) != ARBINT_OK)
    return NULL;

  /* Work with absolute value */
  if (arbint_abs(tmp, x) != ARBINT_OK)
    goto cleanup;

  /* Allocate digit buffer */
  digits = malloc(capacity);
  if (!digits)
    goto cleanup;

  /* Extract digits by repeated division by 10 */
  while (!arbint_is_zero(tmp)) {
    if (arbint_tdiv_qr_u32(q, r, tmp, 10) != ARBINT_OK)
      goto cleanup;

    uint32_t digit;
    if (arbint_get_u32(r, &digit) != ARBINT_OK)
      goto cleanup;

    /* Grow buffer if needed */
    if (len >= capacity - 1) {
      capacity *= 2;
      char * new_digits = realloc(digits, capacity);
      if (!new_digits)
        goto cleanup;
      digits = new_digits;
    }

    digits[len++] = '0' + (char) digit;

    /* tmp = q for next iteration */
    if (arbint_set(tmp, q) != ARBINT_OK)
      goto cleanup;
  }

  /* Allocate final string (digits are in reverse order) */
  str = malloc(len + (is_negative ? 2 : 1));
  if (!str)
    goto cleanup;

  size_t pos = 0;
  if (is_negative)
    str[pos++] = '-';

  /* Reverse digits into final string */
  for (i = 0; i < len; i++)
    str[pos++] = digits[len - 1 - i];
  str[pos] = '\0';

cleanup:
  free(digits);
  arbint_clear_all(r, q, tmp, (arbint_t *) NULL);
  return str;
}

/*  Make TinyCC happy: define a simple sqrt via Newton-Raphson on doubles.  */
#ifdef __TINYC__
static double sqrt(double x) {
  if (x < 0)
    return -1;
  if (x == 0)
    return 0;
  double guess = x;
  double epsilon = 1e-9;
  while ((guess * guess - x) > epsilon || (x - guess * guess) > epsilon) {
    guess = 0.5 * (guess + x / guess);
  }
  return guess;
}
#endif

/*  Fixed-point square root: result = floor(sqrt(n_val) * one).

    Uses floating-point arithmetic on n_val (a small integer) to produce
    an initial guess accurate to ~15 digits, then refines with Newton's
    method. Since Newton doubles the number of correct digits each
    iteration, only about log2(total_digits / 15) ≈ 10 iterations are
    needed for 10000+ digit precision.

    The iteration is: x_{k+1} = (x_k + n_val * one^2 / x_k) / 2
    which converges to sqrt(n_val) * one.  */
static arbint_err_t fixed_sqrt(arbint_t result, uint32_t n_val,
                               const arbint_t one, arbint_ctx_t * ctx) {
  arbint_t x, x_old, n_one, s;
  arbint_err_t rc;

  rc = arbint_init_all(ctx, x, x_old, n_one, s, (arbint_t *) NULL);
  if (rc != ARBINT_OK)
    return rc;

  /* n_one = n_val * one * one  (the value whose isqrt we want) */
  rc = arbint_mul_u32(n_one, one, n_val);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_mul(n_one, n_one, one);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Initial guess via double: x0 = round(sqrt(n_val) * one).

      We can't multiply a double by `one` directly (it's huge), so we
      build x0 = int(sqrt(n_val) * 10^16) * one / 10^16, which gives
      ~15 correct digits.  */
  {
    double sqrt_approx = sqrt((double) n_val);
    /*  fp_prec = 10^16 -- fits in a uint64, and 16 significant digits
        is about the limit of double precision.  */
    uint64_t fp_prec_hi = 10000000u;            /* 10^7  */
    uint64_t fp_prec_lo = 1000000000u;          /* 10^9  */
    uint64_t fp_prec = fp_prec_hi * fp_prec_lo; /* 10^16 */
    uint64_t scaled = (uint64_t) (sqrt_approx * (double) fp_prec);

    /* x = scaled * one */
    /* Build scaled as arbint from two 32-bit halves */
    uint32_t hi = (uint32_t) (scaled / 1000000000ULL);
    uint32_t lo = (uint32_t) (scaled % 1000000000ULL);
    arbint_set_u32(x, hi);
    arbint_mul_u32(x, x, 1000000000u);
    arbint_set_u32(s, lo);
    arbint_add(x, x, s);

    arbint_mul(x, x, one);

    /*  x = x / fp_prec. Split into two u32 divisions since fp_prec
        = 10^7 * 10^9 and the result is only an initial guess.  */
    rc = arbint_tdiv_q_u32(x, x, (uint32_t) fp_prec_hi);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_tdiv_q_u32(x, x, (uint32_t) fp_prec_lo);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

  /*  Newton iterations: x = (x + n_one / x) / 2
      Converge until x stops changing.  */
  {
    int iter = 0;
    while (1) {
      rc = arbint_set(x_old, x);
      if (rc != ARBINT_OK)
        goto cleanup;

      rc = arbint_tdiv_q(s, n_one, x);
      if (rc != ARBINT_OK)
        goto cleanup;
      rc = arbint_add(x, x, s);
      if (rc != ARBINT_OK)
        goto cleanup;
      rc = arbint_tdiv_q_u32(x, x, 2);
      if (rc != ARBINT_OK)
        goto cleanup;

      iter++;
      if (arbint_eq(x, x_old))
        break;
    }
    fprintf(stderr, "  sqrt converged in %d iterations\n", iter);
  }

  rc = arbint_set(result, x);

cleanup:
  arbint_clear_all(s, n_one, x_old, x, (arbint_t *) NULL);
  return rc;
}

/*  Compute pi using Chudnovsky algorithm with incremental term computation.

    Uses the recurrence: a_{k+1} = a_k * -(6k-5)(2k-1)(6k-1) / (k^3 * C3/24)
    where C = 640320, C3/24 = 640320^3 / 24.

    Maintains running sums a_sum and b_sum in fixed-point (scaled by 10^N),
    then computes pi = 426880 * sqrt(10005 * one) * one / total
    where total = 13591409 * a_sum + 545140134 * b_sum.  */
static arbint_err_t compute_pi_chudnovsky(arbint_t result, uint32_t digits,
                                          arbint_ctx_t * ctx) {
  arbint_t a_k, a_sum, b_sum, tmp, tmp2, total, one, sqrt_val, s;
  arbint_err_t rc;

  rc = arbint_init_all(ctx, a_k, a_sum, b_sum, tmp, tmp2, total, one, sqrt_val,
                       s, (arbint_t *) NULL);
  if (rc != ARBINT_OK)
    return rc;

  /* one = 10^(digits+10) -- extra guard digits */
  uint32_t scale = digits + 10;
  arbint_set_u32(one, 1);
  for (uint32_t i = 0; i < scale; i++) {
    rc = arbint_mul_u32(one, one, 10u);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

  /* C3_OVER_24 = 640320^3 / 24 = 10939058860032000 */
  arbint_t c3_24;
  if ((rc = arbint_init(c3_24, ctx)) != ARBINT_OK)
    goto cleanup;

  /*  Build 10939058860032000 from two 32-bit halves:
      10939058860032000 = 2546948 * 2^32 + 1502527488
      Or more simply: 10939058860032000 = 10939058 * 10^9 + 860032000.  */
  arbint_set_u32(c3_24, 10939058u);
  arbint_mul_u32(c3_24, c3_24, 1000000000u);
  arbint_add_u32(c3_24, c3_24, 860032000u);

  /* a_k = one (= 10^scale), a_sum = one, b_sum = 0 */
  arbint_set(a_k, one);
  arbint_set(a_sum, one);
  arbint_set_u32(b_sum, 0);

  fprintf(stderr, "  series: ");
  uint32_t k = 1;
  while (1) {
    if (k % 50 == 0)
      fprintf(stderr, "\r  series: term %u", k);

    /* a_k *= -(6k-5) * (2k-1) * (6k-1) */
    rc = arbint_mul_i32(tmp, a_k, -(int32_t) (6 * k - 5));
    if (rc != ARBINT_OK) {
      arbint_clear(c3_24);
      goto cleanup;
    }
    rc = arbint_mul_u32(tmp, tmp, 2 * k - 1);
    if (rc != ARBINT_OK) {
      arbint_clear(c3_24);
      goto cleanup;
    }
    rc = arbint_mul_u32(tmp, tmp, 6 * k - 1);
    if (rc != ARBINT_OK) {
      arbint_clear(c3_24);
      goto cleanup;
    }

    /* a_k /= k^3 * C3_OVER_24 */
    arbint_set_u32(tmp2, k);
    arbint_mul_u32(tmp2, tmp2, k);
    arbint_mul_u32(tmp2, tmp2, k);
    arbint_mul(tmp2, tmp2, c3_24);

    rc = arbint_tdiv_q(a_k, tmp, tmp2);
    if (rc != ARBINT_OK) {
      arbint_clear(c3_24);
      goto cleanup;
    }

    /* a_sum += a_k */
    rc = arbint_add(a_sum, a_sum, a_k);
    if (rc != ARBINT_OK) {
      arbint_clear(c3_24);
      goto cleanup;
    }

    /* b_sum += k * a_k */
    arbint_mul_u32(tmp, a_k, k);
    rc = arbint_add(b_sum, b_sum, tmp);
    if (rc != ARBINT_OK) {
      arbint_clear(c3_24);
      goto cleanup;
    }

    k++;

    /*  Terminate when a_k reaches zero (all remaining terms are zero
        at this precision).  */
    if (arbint_is_zero(a_k))
      break;
  }
  fprintf(stderr, "\r  series: %u terms computed\n", k - 1);

  arbint_clear(c3_24);

  /* total = 13591409 * a_sum + 545140134 * b_sum */
  fprintf(stderr, "  computing total...\n");
  arbint_mul_u32(total, a_sum, 13591409u);
  arbint_mul_u32(tmp, b_sum, 545140134u);
  arbint_add(total, total, tmp);

  /* sqrt_val = floor(sqrt(10005) * one) via fixed-point Newton */
  fprintf(stderr, "  computing integer sqrt...\n");
  rc = fixed_sqrt(sqrt_val, 10005, one, ctx);
  if (rc != ARBINT_OK)
    goto cleanup;

  /* pi = 426880 * sqrt_val * one / total */
  fprintf(stderr, "  final division...\n");
  arbint_mul_u32(tmp, sqrt_val, 426880u);
  arbint_mul(tmp, tmp, one);

  rc = arbint_tdiv_q(result, tmp, total);
  if (rc != ARBINT_OK)
    goto cleanup;

cleanup:
  arbint_clear_all(s, sqrt_val, one, total, tmp2, tmp, b_sum, a_sum, a_k,
                   (arbint_t *) NULL);
  return rc;
}

/* Load reference pi digits from pi10k.txt.
 * Returns a malloc'd string of just the digits (no '.' separator), or NULL. */
static char * load_reference_digits(const char * path) {
  FILE * f = fopen(path, "r");
  if (!f)
    return NULL;

  /* The file is "3." followed by 10000 digits on one line */
  fseek(f, 0, SEEK_END);
  long fsize = ftell(f);
  fseek(f, 0, SEEK_SET);

  if (fsize < 3) {
    fclose(f);
    return NULL;
  }

  char * buf = malloc((size_t) fsize + 1);
  if (!buf) {
    fclose(f);
    return NULL;
  }

  size_t nread = fread(buf, 1, (size_t) fsize, f);
  fclose(f);
  buf[nread] = '\0';

  /* Strip trailing whitespace/newline */
  while (nread > 0 && (buf[nread - 1] == '\n' || buf[nread - 1] == '\r' ||
                       buf[nread - 1] == ' ')) {
    buf[--nread] = '\0';
  }

  /* Expect "3." prefix */
  if (buf[0] != '3' || buf[1] != '.') {
    fprintf(stderr, "Reference file doesn't start with '3.'\n");
    free(buf);
    return NULL;
  }

  /* Return a copy of just the full digit string: "3" + fractional digits */
  size_t digit_len = nread - 1; /* remove the '.' */
  char * digits = malloc(digit_len + 1);
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

int main(void) {
  ARBINT_TEST_DECLARE_FAILURES();

  arbint_ctx_t ctx;
  arbint_t pi_scaled;
  arbint_err_t rc;
  char * computed_str = NULL;
  char * reference = NULL;
  char pi10k_path[4096];

  /*  Locate pi10k.txt: check srcdir env var (set by automake),
      then fall back to current directory.  */
  const char * srcdir = getenv("srcdir");
  if (srcdir)
    snprintf(pi10k_path, sizeof(pi10k_path), "%s/pi10k.txt", srcdir);
  else
    snprintf(pi10k_path, sizeof(pi10k_path), "pi10k.txt");

  /* Load reference digits */
  reference = load_reference_digits(pi10k_path);
  if (!reference) {
    /* Try fallback path relative to tests/ */
    reference = load_reference_digits("tests/pi10k.txt");
  }
  CHECK(reference != NULL);
  if (!reference) {
    fprintf(stderr, "Could not load pi10k.txt reference file\n");
    ARBINT_TEST_FINISH("compute_pi");
  }

  /* Initialize context and pi variable */
  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(pi_scaled, &ctx), ARBINT_OK);

  /* Compute pi to PI_DIGITS digits */
  fprintf(stderr, "Computing %d digits of pi...\n", PI_DIGITS);
  rc = compute_pi_chudnovsky(pi_scaled, PI_DIGITS, &ctx);
  CHECK_EQ_I(rc, ARBINT_OK);

  if (rc == ARBINT_OK) {
    fprintf(stderr, "Converting to decimal string...\n");

    /* Convert to decimal string */
    computed_str = arbint_to_decimal_string(pi_scaled, &ctx);
    CHECK(computed_str != NULL);

    if (computed_str) {
      size_t computed_len = strlen(computed_str);
      size_t reference_len = strlen(reference);

      fprintf(stderr, "Computed %zu digits, reference has %zu digits\n",
              computed_len, reference_len);

      /*  The computed result is pi * 10^(digits+10), so it should have
          digits+11 characters (the '3' plus digits+10 fractional digits).
          We compare the first digits+1 characters (the leading '3' plus
          the requested number of fractional digits).  */
      size_t compare_len = PI_DIGITS + 1; /* "3" + 10000 fractional digits */
      CHECK(computed_len >= compare_len);
      CHECK(reference_len >= compare_len);

      if (computed_len >= compare_len && reference_len >= compare_len) {
        int match = (memcmp(computed_str, reference, compare_len) == 0);
        if (!match) {
          /* Find first mismatch for diagnostic output */
          size_t i;
          for (i = 0; i < compare_len; i++) {
            if (computed_str[i] != reference[i])
              break;
          }
          fprintf(stderr,
                  "MISMATCH at digit %zu: computed '%c', expected '%c'\n", i,
                  computed_str[i], reference[i]);
          /* Show context around the mismatch */
          size_t start = (i > 10) ? i - 10 : 0;
          size_t end = (i + 10 < compare_len) ? i + 10 : compare_len;
          fprintf(stderr, "  computed:  ");
          for (size_t j = start; j < end; j++)
            fprintf(stderr, "%c", computed_str[j]);
          fprintf(stderr, "\n  expected:  ");
          for (size_t j = start; j < end; j++)
            fprintf(stderr, "%c", reference[j]);
          fprintf(stderr, "\n");
        }
        CHECK(match);
      }

      free(computed_str);
    }
  }

  free(reference);
  arbint_clear(pi_scaled);
  arbint_ctx_clear(&ctx);

  ARBINT_TEST_FINISH("compute_pi");
}
