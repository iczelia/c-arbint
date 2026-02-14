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

/*  Chudnovsky pi computation using c-arbint.

    This is the arbint implementation of the Chudnovsky algorithm for
    benchmarking against GMP. The algorithm is identical to ensure fair
    comparison.  */

#include "bench_pi_common.h"

#include <arbint.h>
#include "arbint_div_bz.h"
#include "arbint_div_newton.h"

typedef struct {
  uint64_t calls_total;
  uint64_t calls_single_limb;
  uint64_t calls_two_limb;
  uint64_t calls_knuth_small;
  uint64_t calls_bz;
  uint64_t calls_newton;
  size_t min_nn;
  size_t max_nn;
  size_t min_dn;
  size_t max_dn;
} bench_div_stats_t;

static void bench_div_stats_record(bench_div_stats_t * stats, const arbint_t n,
                                    const arbint_t d) {
  size_t nn = arbint_abs_sz(n[0]._sz);
  size_t dn = arbint_abs_sz(d[0]._sz);

  if (dn == 0u)
    return;

  if (stats->calls_total == 0u) {
    stats->min_nn = nn;
    stats->max_nn = nn;
    stats->min_dn = dn;
    stats->max_dn = dn;
  } else {
    if (nn < stats->min_nn)
      stats->min_nn = nn;
    if (nn > stats->max_nn)
      stats->max_nn = nn;
    if (dn < stats->min_dn)
      stats->min_dn = dn;
    if (dn > stats->max_dn)
      stats->max_dn = dn;
  }

  stats->calls_total++;
  if (dn == 1u) {
    stats->calls_single_limb++;
  } else if (dn == 2u) {
    stats->calls_two_limb++;
  } else if (dn >= ARBINT_NEWTON_DIV_THRESHOLD) {
    stats->calls_newton++;
  } else if (dn >= ARBINT_BZ_THRESHOLD) {
    stats->calls_bz++;
  } else {
    stats->calls_knuth_small++;
  }
}

/*  Compute pi using Chudnovsky algorithm with incremental term computation.

    Uses the recurrence: a_{k+1} = a_k * -(6k-5)(2k-1)(6k-1) / (k^3 * C3/24)
    where C = 640320, C3/24 = 640320^3 / 24.

    Maintains running sums a_sum and b_sum in fixed-point (scaled by 10^N),
    then computes pi = 426880 * sqrt(10005 * one) * one / total
    where total = 13591409 * a_sum + 545140134 * b_sum.  */
int bench_pi_arbint(uint32_t digits, bench_pi_result_t * result) {
  arbint_ctx_t ctx;
  arbint_t a_k, a_sum, b_sum, tmp, tmp2, total, one, sqrt_val, c3_24, pi;
  arbint_err_t rc;
  uint64_t t_start, t_series, t_sqrt, t_final;
  uint64_t t_step0, t_step1, t_step2, t_step3, t_step4;
  bench_div_stats_t div_stats = {0};
  int div_stats_enabled = (getenv("ARBINT_BENCH_DIV_STATS") != NULL);
  int phase_detail_enabled = (getenv("ARBINT_BENCH_PHASE_DETAIL") != NULL);
  uint32_t k;

  memset(result, 0, sizeof(*result));

  rc = arbint_ctx_init_default(&ctx);
  if (rc != ARBINT_OK)
    return 0;

  rc = arbint_init_all(&ctx, a_k, a_sum, b_sum, tmp, tmp2, total, one, sqrt_val,
                       c3_24, pi, (arbint_t *) NULL);
  if (rc != ARBINT_OK) {
    arbint_ctx_clear(&ctx);
    return 0;
  }

  t_start = bench_get_time_ns();

  {
    uint32_t scale = digits + 10;
    arbint_set_u32(one, 1);
    for (uint32_t i = 0; i < scale; ++i) {
      rc = arbint_mul_u32(one, one, 10u);
      if (rc != ARBINT_OK)
        goto fail;
    }
  }

  arbint_set_u32(c3_24, CHUD_C3_24_HI);
  arbint_mul_u32(c3_24, c3_24, 1000000000u);
  arbint_add_u32(c3_24, c3_24, CHUD_C3_24_LO);

  arbint_set(a_k, one);
  arbint_set(a_sum, one);
  arbint_set_u32(b_sum, 0);

  k = 1;
  while (1) {
    rc = arbint_mul_i32(tmp, a_k, -(int32_t) (6 * k - 5));
    if (rc != ARBINT_OK)
      goto fail;
    rc = arbint_mul_u32(tmp, tmp, 2 * k - 1);
    if (rc != ARBINT_OK)
      goto fail;
    rc = arbint_mul_u32(tmp, tmp, 6 * k - 1);
    if (rc != ARBINT_OK)
      goto fail;

    arbint_set_u32(tmp2, k);
    arbint_mul_u32(tmp2, tmp2, k);
    arbint_mul_u32(tmp2, tmp2, k);
    arbint_mul(tmp2, tmp2, c3_24);

    if (div_stats_enabled)
      bench_div_stats_record(&div_stats, tmp, tmp2);

    rc = arbint_tdiv_q(a_k, tmp, tmp2);
    if (rc != ARBINT_OK)
      goto fail;

    rc = arbint_add(a_sum, a_sum, a_k);
    if (rc != ARBINT_OK)
      goto fail;

    arbint_mul_u32(tmp, a_k, k);
    rc = arbint_add(b_sum, b_sum, tmp);
    if (rc != ARBINT_OK)
      goto fail;

    ++k;
    if (arbint_is_zero(a_k))
      break;
  }

  result->timing.terms = k - 1;
  t_series = bench_get_time_ns();
  result->timing.series_ns = t_series - t_start;

  t_step0 = t_series;
  arbint_mul_u32(total, a_sum, CHUD_A);
  arbint_mul_u32(tmp, b_sum, CHUD_B);
  arbint_add(total, total, tmp);
  t_step1 = bench_get_time_ns();

  rc = arbint_mul_u32(tmp, one, CHUD_D);
  if (rc != ARBINT_OK)
    goto fail;
  rc = arbint_mul(tmp, tmp, one);
  if (rc != ARBINT_OK)
    goto fail;
  t_step2 = bench_get_time_ns();
  rc = arbint_isqrt(sqrt_val, tmp);
  if (rc != ARBINT_OK)
    goto fail;
  t_step3 = bench_get_time_ns();

  t_sqrt = bench_get_time_ns();
  result->timing.sqrt_ns = t_sqrt - t_series;

  arbint_mul_u32(tmp, sqrt_val, CHUD_C);
  arbint_mul(tmp, tmp, one);
  t_step4 = bench_get_time_ns();

  if (div_stats_enabled)
    bench_div_stats_record(&div_stats, tmp, total);

  rc = arbint_tdiv_q(pi, tmp, total);
  if (rc != ARBINT_OK)
    goto fail;

  t_final = bench_get_time_ns();
  result->timing.final_ns = t_final - t_sqrt;
  result->timing.total_ns = t_final - t_start;

  if (div_stats_enabled && div_stats.calls_total > 0u) {
    printf("\n[arbint] Division dispatch stats (estimated by divisor limbs):\n");
    printf("  Calls:              %llu\n",
           (unsigned long long) div_stats.calls_total);
    printf("  Dividend limbs:     min=%zu max=%zu\n", div_stats.min_nn,
           div_stats.max_nn);
    printf("  Divisor limbs:      min=%zu max=%zu\n", div_stats.min_dn,
           div_stats.max_dn);
    printf("  Single-limb path:   %llu\n",
           (unsigned long long) div_stats.calls_single_limb);
    printf("  Two-limb path:      %llu\n",
           (unsigned long long) div_stats.calls_two_limb);
    printf("  Knuth (<%u limbs):   %llu\n", (unsigned) ARBINT_BZ_THRESHOLD,
           (unsigned long long) div_stats.calls_knuth_small);
    printf("  BZ (%u..%u limbs):   %llu\n", (unsigned) ARBINT_BZ_THRESHOLD,
           (unsigned) (ARBINT_NEWTON_DIV_THRESHOLD - 1u),
           (unsigned long long) div_stats.calls_bz);
    printf("  Newton (>=%u limbs): %llu\n",
           (unsigned) ARBINT_NEWTON_DIV_THRESHOLD,
           (unsigned long long) div_stats.calls_newton);
  }

  if (phase_detail_enabled) {
    printf("\n[arbint] Phase detail:\n");
    printf("  Build total (A*a_sum + B*b_sum): %.3f s\n",
           (double) (t_step1 - t_step0) / 1e9);
    printf("  Build sqrt input (D*one*one):    %.3f s\n",
           (double) (t_step2 - t_step1) / 1e9);
    printf("  isqrt core only:                  %.3f s\n",
           (double) (t_step3 - t_step2) / 1e9);
    printf("  Build final numerator (C*s*one):  %.3f s\n",
           (double) (t_step4 - t_sqrt) / 1e9);
    printf("  Final divide only:                %.3f s\n",
           (double) (t_final - t_step4) / 1e9);
  }

  {
    char * str = NULL;
    rc = arbint_get_str(pi, &str, 10);
    if (rc != ARBINT_OK || !str)
      goto fail;

    result->digits = str;
    result->digit_count = strlen(str);
  }

  arbint_clear_all(pi, c3_24, sqrt_val, one, total, tmp2, tmp, b_sum, a_sum,
                   a_k, (arbint_t *) NULL);
  arbint_ctx_clear(&ctx);
  return 1;

fail:
  arbint_clear_all(pi, c3_24, sqrt_val, one, total, tmp2, tmp, b_sum, a_sum,
                   a_k, (arbint_t *) NULL);
  arbint_ctx_clear(&ctx);
  return 0;
}
