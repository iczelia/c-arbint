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

/*  Chudnovsky pi computation using GMP.

    This is the GMP implementation of the Chudnovsky algorithm for
    benchmarking against c-arbint. The algorithm is identical to ensure
    fair comparison.  */

#include "bench_pi_common.h"

#include <gmp.h>

static void fixed_sqrt_gmp(mpz_t result, uint32_t n_val, const mpz_t one) {
  mpz_t n_one;
  mpz_init(n_one);
  mpz_mul_ui(n_one, one, n_val);
  mpz_mul(n_one, n_one, one);
  mpz_sqrt(result, n_one);
  mpz_clear(n_one);
}

/*  Compute pi using Chudnovsky algorithm with incremental term computation.

    Uses the recurrence: a_{k+1} = a_k * -(6k-5)(2k-1)(6k-1) / (k^3 * C3/24)
    where C = 640320, C3/24 = 640320^3 / 24.

    Maintains running sums a_sum and b_sum in fixed-point (scaled by 10^N),
    then computes pi = 426880 * sqrt(10005 * one) * one / total
    where total = 13591409 * a_sum + 545140134 * b_sum.  */
int bench_pi_gmp(uint32_t digits, bench_pi_result_t * result) {
  mpz_t a_k, a_sum, b_sum, tmp, tmp2, total, one, sqrt_val, c3_24, pi;
  uint64_t t_start, t_series, t_sqrt, t_final;
  uint32_t k;

  memset(result, 0, sizeof(*result));

  mpz_init(a_k);
  mpz_init(a_sum);
  mpz_init(b_sum);
  mpz_init(tmp);
  mpz_init(tmp2);
  mpz_init(total);
  mpz_init(one);
  mpz_init(sqrt_val);
  mpz_init(c3_24);
  mpz_init(pi);

  t_start = bench_get_time_ns();

  {
    uint32_t scale = digits + 10;
    mpz_set_ui(one, 1);
    for (uint32_t i = 0; i < scale; ++i)
      mpz_mul_ui(one, one, 10u);
  }

  mpz_set_ui(c3_24, CHUD_C3_24_HI);
  mpz_mul_ui(c3_24, c3_24, 1000000000u);
  mpz_add_ui(c3_24, c3_24, CHUD_C3_24_LO);

  mpz_set(a_k, one);
  mpz_set(a_sum, one);
  mpz_set_ui(b_sum, 0);

  k = 1;
  while (1) {
    mpz_mul_ui(tmp, a_k, 6 * k - 5);
    mpz_mul_ui(tmp, tmp, 2 * k - 1);
    mpz_mul_ui(tmp, tmp, 6 * k - 1);
    mpz_neg(tmp, tmp);

    mpz_set_ui(tmp2, k);
    mpz_mul_ui(tmp2, tmp2, k);
    mpz_mul_ui(tmp2, tmp2, k);
    mpz_mul(tmp2, tmp2, c3_24);

    mpz_tdiv_q(a_k, tmp, tmp2);
    mpz_add(a_sum, a_sum, a_k);
    mpz_mul_ui(tmp, a_k, k);
    mpz_add(b_sum, b_sum, tmp);

    ++k;
    if (mpz_sgn(a_k) == 0)
      break;
  }

  result->timing.terms = k - 1;
  t_series = bench_get_time_ns();
  result->timing.series_ns = t_series - t_start;

  mpz_mul_ui(total, a_sum, CHUD_A);
  mpz_mul_ui(tmp, b_sum, CHUD_B);
  mpz_add(total, total, tmp);

  fixed_sqrt_gmp(sqrt_val, CHUD_D, one);

  t_sqrt = bench_get_time_ns();
  result->timing.sqrt_ns = t_sqrt - t_series;

  mpz_mul_ui(tmp, sqrt_val, CHUD_C);
  mpz_mul(tmp, tmp, one);
  mpz_tdiv_q(pi, tmp, total);

  t_final = bench_get_time_ns();
  result->timing.final_ns = t_final - t_sqrt;
  result->timing.total_ns = t_final - t_start;

  {
    char * str = mpz_get_str(NULL, 10, pi);
    if (!str)
      goto fail;

    result->digits = str;
    result->digit_count = strlen(str);
  }

  mpz_clear(pi);
  mpz_clear(c3_24);
  mpz_clear(sqrt_val);
  mpz_clear(one);
  mpz_clear(total);
  mpz_clear(tmp2);
  mpz_clear(tmp);
  mpz_clear(b_sum);
  mpz_clear(a_sum);
  mpz_clear(a_k);
  return 1;

fail:
  mpz_clear(pi);
  mpz_clear(c3_24);
  mpz_clear(sqrt_val);
  mpz_clear(one);
  mpz_clear(total);
  mpz_clear(tmp2);
  mpz_clear(tmp);
  mpz_clear(b_sum);
  mpz_clear(a_sum);
  mpz_clear(a_k);
  return 0;
}
