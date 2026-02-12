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

/*  Calculate e (Euler's number) to 10000 digits using Brent's argument
    reduction method and compare against known-good reference digits in
    e10k.txt.

    This test exercises multiplication, squaring and division kernels on
    large operands (10000+ decimal digits), validating both correctness
    and precision of the library's arithmetic.

    Algorithm -- argument reduction + repeated squaring:
      1. Pick r such that x = 1/2^r is tiny (r = 32 gives x ~ 2.3e-10).
      2. Evaluate exp(x) via the Taylor series 1 + x + x^2/2! + ...
         which converges in about 35 terms at this scale.
      3. Square the result r times: exp(1) = exp(x)^(2^r).

    All arithmetic is fixed-point with one = 10^(digits + guard).  */

#include "test_framework.h"

#include <arbint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define E_DIGITS 10000

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

  if (arbint_init_all(ctx, tmp, q, r, (arbint_t *) NULL) != ARBINT_OK)
    return NULL;

  if (arbint_abs(tmp, x) != ARBINT_OK)
    goto cleanup;

  digits = malloc(capacity);
  if (!digits)
    goto cleanup;

  while (!arbint_is_zero(tmp)) {
    if (arbint_tdiv_qr_u32(q, r, tmp, 10) != ARBINT_OK)
      goto cleanup;

    uint32_t digit;
    if (arbint_get_u32(r, &digit) != ARBINT_OK)
      goto cleanup;

    if (len >= capacity - 1) {
      capacity *= 2;
      char * new_digits = realloc(digits, capacity);
      if (!new_digits)
        goto cleanup;
      digits = new_digits;
    }

    digits[len++] = '0' + (char) digit;

    if (arbint_set(tmp, q) != ARBINT_OK)
      goto cleanup;
  }

  str = malloc(len + (is_negative ? 2 : 1));
  if (!str)
    goto cleanup;

  size_t pos = 0;
  if (is_negative)
    str[pos++] = '-';

  for (i = 0; i < len; i++)
    str[pos++] = digits[len - 1 - i];
  str[pos] = '\0';

cleanup:
  free(digits);
  arbint_clear_all(r, q, tmp, (arbint_t *) NULL);
  return str;
}

/*  Compute e using Brent's argument reduction.

    1. one = 10^(digits + guard)
    2. x = one >> r              (fixed-point 1/2^r)
    3. sum = exp(x) via Taylor   (converges in ~35 terms)
    4. square r times            (exp(1) = exp(x)^(2^r))

    Each squaring: sum = sum * sum / one (fixed-point).  */
static arbint_err_t compute_e_brent(arbint_t result, uint32_t digits,
                                    arbint_ctx_t * ctx) {
  arbint_t one, x, term, sum, tmp;
  arbint_err_t rc;
  uint32_t scale = digits + 50;
  uint32_t r = 32;
  uint32_t k;
  uint32_t i;

  rc = arbint_init_all(ctx, one, x, term, sum, tmp, (arbint_t *) NULL);
  if (rc != ARBINT_OK)
    return rc;

  /*  one = 10^scale  */
  arbint_set_u32(one, 1);
  for (i = 0; i < scale; i++) {
    rc = arbint_mul_u32(one, one, 10u);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

  /*  x = one >> r  (fixed-point representation of 1/2^r)  */
  rc = arbint_shr(x, one, r);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Taylor series: sum = one + x + x^2/2! + x^3/3! + ...
      term_k = term_{k-1} * x / (one * k)  */
  rc = arbint_set(term, x);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_add(sum, one, x);
  if (rc != ARBINT_OK)
    goto cleanup;

  fprintf(stderr, "  Taylor series: ");
  for (k = 2;; k++) {
    /*  term = term * x / one  */
    rc = arbint_mul(tmp, term, x);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_tdiv_q(term, tmp, one);
    if (rc != ARBINT_OK)
      goto cleanup;

    /*  term = term / k  */
    rc = arbint_tdiv_q_u32(term, term, k);
    if (rc != ARBINT_OK)
      goto cleanup;

    if (arbint_is_zero(term))
      break;

    rc = arbint_add(sum, sum, term);
    if (rc != ARBINT_OK)
      goto cleanup;
  }
  fprintf(stderr, "%u terms\n", k - 1);

  /*  Repeated squaring: sum = sum^(2^r) in fixed-point.
      Each step: sum = sum * sum / one.  */
  fprintf(stderr, "  squaring: ");
  for (i = 0; i < r; i++) {
    if (i % 8 == 0)
      fprintf(stderr, "\r  squaring: %u / %u", i, r);

    rc = arbint_sqr(tmp, sum);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_tdiv_q(sum, tmp, one);
    if (rc != ARBINT_OK)
      goto cleanup;
  }
  fprintf(stderr, "\r  squaring: %u / %u done\n", r, r);

  rc = arbint_set(result, sum);

cleanup:
  arbint_clear_all(tmp, sum, term, x, one, (arbint_t *) NULL);
  return rc;
}

/*  Load reference e digits from e10k.txt.
    Returns a malloc'd string of just the digits (no '.' separator),
    or NULL on failure.  */
static char * load_reference_digits(const char * path) {
  FILE * f = fopen(path, "r");
  size_t nread;
  size_t digit_len;
  char * buf = NULL;
  char * digits = NULL;
  long fsize;

  if (!f)
    return NULL;

  fseek(f, 0, SEEK_END);
  fsize = ftell(f);
  fseek(f, 0, SEEK_SET);

  if (fsize < 3) {
    fclose(f);
    return NULL;
  }

  buf = malloc((size_t) fsize + 1);
  if (!buf) {
    fclose(f);
    return NULL;
  }

  nread = fread(buf, 1, (size_t) fsize, f);
  fclose(f);
  buf[nread] = '\0';

  while (nread > 0 && (buf[nread - 1] == '\n' || buf[nread - 1] == '\r' ||
                       buf[nread - 1] == ' '))
    buf[--nread] = '\0';

  /*  Expect "2." prefix  */
  if (buf[0] != '2' || buf[1] != '.') {
    fprintf(stderr, "Reference file doesn't start with '2.'\n");
    free(buf);
    return NULL;
  }

  /*  Return "2" + fractional digits (strip the '.')  */
  digit_len = nread - 1;
  digits = malloc(digit_len + 1);
  if (!digits) {
    free(buf);
    return NULL;
  }

  digits[0] = '2';
  memcpy(digits + 1, buf + 2, nread - 2);
  digits[digit_len] = '\0';

  free(buf);
  return digits;
}

int main(void) {
  ARBINT_TEST_DECLARE_FAILURES();
  ARBINT_TEST_START();

  arbint_ctx_t ctx;
  arbint_t e_scaled;
  arbint_err_t rc;
  char * computed_str = NULL;
  char * reference = NULL;
  char e10k_path[4096];

  /*  Locate e10k.txt: check srcdir env var (set by automake),
      then fall back to current directory.  */
  const char * srcdir = getenv("srcdir");
  if (srcdir)
    snprintf(e10k_path, sizeof(e10k_path), "%s/e10k.txt", srcdir);
  else
    snprintf(e10k_path, sizeof(e10k_path), "e10k.txt");

  reference = load_reference_digits(e10k_path);
  if (!reference)
    reference = load_reference_digits("tests/e10k.txt");
  CHECK(reference != NULL);
  if (!reference) {
    fprintf(stderr, "Could not load e10k.txt reference file\n");
    ARBINT_TEST_FINISH("compute_e");
  }

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(e_scaled, &ctx), ARBINT_OK);

  fprintf(stderr, "Computing %d digits of e...\n", E_DIGITS);
  rc = compute_e_brent(e_scaled, E_DIGITS, &ctx);
  CHECK_EQ_I(rc, ARBINT_OK);

  if (rc == ARBINT_OK) {
    fprintf(stderr, "Converting to decimal string...\n");

    computed_str = arbint_to_decimal_string(e_scaled, &ctx);
    CHECK(computed_str != NULL);

    if (computed_str) {
      size_t computed_len = strlen(computed_str);
      size_t reference_len = strlen(reference);

      fprintf(stderr, "Computed %zu digits, reference has %zu digits\n",
              computed_len, reference_len);

      /*  The computed result is e * 10^(digits+guard) truncated, so it
          may differ in the last reference digit if that digit was rounded.
          Drop the last digit of the reference from comparison.  */
      size_t compare_len = reference_len - 1;
      CHECK(computed_len >= compare_len);

      if (computed_len >= compare_len && reference_len >= compare_len) {
        int match = (memcmp(computed_str, reference, compare_len) == 0);
        if (!match) {
          size_t i;
          for (i = 0; i < compare_len; i++) {
            if (computed_str[i] != reference[i])
              break;
          }
          fprintf(stderr,
                  "MISMATCH at digit %zu: computed '%c', expected '%c'\n", i,
                  computed_str[i], reference[i]);
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
  arbint_clear(e_scaled);
  arbint_ctx_clear(&ctx);

  ARBINT_TEST_FINISH("compute_e");
}
