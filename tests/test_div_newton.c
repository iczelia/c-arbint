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

/*  Test suite for Newton-Raphson division.

    This test suite verifies:
    1. Basic correctness for various divisor sizes
    2. Threshold boundary cases (Knuth vs Newton transition)
    3. Adversarial cases (normalization boundaries, exact multiples)
    4. Large nn/dn ratios
    5. Regression against Knuth Algorithm D  */

#include "test_framework.h"

#include <arbint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*  Threshold from arbint_div_newton.h - duplicated here for test visibility.  */
#ifndef ARBINT_NEWTON_DIV_THRESHOLD
#define ARBINT_NEWTON_DIV_THRESHOLD 64u
#endif

/*  Simple LCG for deterministic random testing.  */
static uint64_t test_rng_state = 0x123456789ABCDEF0ULL;

static uint64_t test_rand64(void) {
  test_rng_state = test_rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
  return test_rng_state;
}

static void test_rand_seed(uint64_t seed) {
  test_rng_state = seed;
}

/*  Fill limb array with random data.  */
static void fill_random_limbs(arbint_limb_t * p, size_t n) {
  size_t i;
  for (i = 0; i < n; ++i) {
#if ARBINT_LIMB_BITS == 64
    p[i] = test_rand64();
#else
    p[i] = (arbint_limb_t) test_rand64();
#endif
  }
  /*  Ensure top limb is nonzero.  */
  if (n > 0 && p[n - 1] == 0)
    p[n - 1] = 1;
}

/*  Set arbint from raw limb array using public API.  */
static arbint_err_t set_from_limbs(arbint_t x, const arbint_limb_t * limbs,
                                    size_t n) {
  /*  Use arbint_import to set from bytes (little-endian).
      The import function normalizes leading zeros automatically.  */
  return arbint_import(x, limbs, n * sizeof(arbint_limb_t));
}

/*  Verify division: check that n = q*d + r and 0 <= r < d.  */
static int verify_division(const arbint_t n, const arbint_t d,
                            const arbint_t q, const arbint_t r,
                            int * failures) {
  arbint_t check, qd;
  arbint_err_t rc;
  int ok = 1;

  rc = arbint_init(check, n[0]._ctx);
  if (rc != ARBINT_OK) {
    ++(*failures);
    return 0;
  }
  rc = arbint_init(qd, n[0]._ctx);
  if (rc != ARBINT_OK) {
    arbint_clear(check);
    ++(*failures);
    return 0;
  }

  /*  check = q * d + r  */
  rc = arbint_mul(qd, q, d);
  if (rc != ARBINT_OK) {
    ok = 0;
    goto cleanup;
  }
  rc = arbint_add(check, qd, r);
  if (rc != ARBINT_OK) {
    ok = 0;
    goto cleanup;
  }

  /*  Verify n == q*d + r.  */
  if (arbint_cmp(check, n) != 0) {
    fprintf(stderr, "FAIL: n != q*d + r\n");
    ok = 0;
    ++(*failures);
  }

  /*  Verify 0 <= r < |d|.  */
  if (arbint_signum(r) < 0) {
    fprintf(stderr, "FAIL: r < 0\n");
    ok = 0;
    ++(*failures);
  }
  if (arbint_cmpabs(r, d) >= 0) {
    fprintf(stderr, "FAIL: r >= |d|\n");
    ok = 0;
    ++(*failures);
  }

cleanup:
  arbint_clear(qd);
  arbint_clear(check);
  return ok;
}

int main(void) {
  ARBINT_TEST_DECLARE_FAILURES();

  arbint_ctx_t ctx;
  arbint_t n, d, q, r;
  arbint_limb_t * limbs_n = NULL;
  arbint_limb_t * limbs_d = NULL;
  size_t trial;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  test_rand_seed(42);

  printf("Test 1: Basic small divisions (sanity check)\n");
  {
    /*  7 / 3 = 2 remainder 1  */
    CHECK_EQ_I(arbint_set_i32(n, 7), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(d, 3), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
    check_i32_value(q, 2);
    check_i32_value(r, 1);

    /*  100 / 7 = 14 remainder 2  */
    CHECK_EQ_I(arbint_set_i32(n, 100), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(d, 7), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
    check_i32_value(q, 14);
    check_i32_value(r, 2);
  }

  printf("Test 2: Threshold boundary tests (Knuth vs Newton transition)\n");
  {
    /*  Test divisors around the threshold.  */
    size_t dn_vals[] = {
        ARBINT_NEWTON_DIV_THRESHOLD - 2,
        ARBINT_NEWTON_DIV_THRESHOLD - 1,
        ARBINT_NEWTON_DIV_THRESHOLD,
        ARBINT_NEWTON_DIV_THRESHOLD + 1,
        ARBINT_NEWTON_DIV_THRESHOLD + 2
    };
    size_t num_dn = sizeof(dn_vals) / sizeof(dn_vals[0]);
    size_t i;

    for (i = 0; i < num_dn; ++i) {
      size_t dn = dn_vals[i];
      size_t nn = dn + 10;

      if (dn < 2)
        continue;

      limbs_n = (arbint_limb_t *) malloc(nn * sizeof(arbint_limb_t));
      limbs_d = (arbint_limb_t *) malloc(dn * sizeof(arbint_limb_t));
      CHECK(limbs_n != NULL && limbs_d != NULL);

      fill_random_limbs(limbs_n, nn);
      fill_random_limbs(limbs_d, dn);

      CHECK_EQ_I(set_from_limbs(n, limbs_n, nn), ARBINT_OK);
      CHECK_EQ_I(set_from_limbs(d, limbs_d, dn), ARBINT_OK);
      CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);

      verify_division(n, d, q, r, &g_failures);

      free(limbs_n);
      free(limbs_d);
      limbs_n = NULL;
      limbs_d = NULL;
    }
  }

  printf("Test 3: Adversarial - divisor top limb near boundaries\n");
  {
    /*  Test with top limb values that stress normalization.  */
    arbint_limb_t top_vals[] = {
#if ARBINT_LIMB_BITS == 64
        0x8000000000000000ULL, /*  Exact MSB  */
        0x8000000000000001ULL, /*  Just above MSB  */
        0xFFFFFFFFFFFFFFFFULL, /*  All ones  */
        0xFFFFFFFFFFFFFFFEULL, /*  All ones minus 1  */
        0x8123456789ABCDEFULL  /*  Random normalized  */
#else
        0x80000000u,
        0x80000001u,
        0xFFFFFFFFu,
        0xFFFFFFFEu,
        0x89ABCDEFu
#endif
    };
    size_t num_tops = sizeof(top_vals) / sizeof(top_vals[0]);
    size_t i;

    for (i = 0; i < num_tops; ++i) {
      size_t dn = 4;
      size_t nn = 8;

      limbs_n = (arbint_limb_t *) malloc(nn * sizeof(arbint_limb_t));
      limbs_d = (arbint_limb_t *) malloc(dn * sizeof(arbint_limb_t));
      CHECK(limbs_n != NULL && limbs_d != NULL);

      fill_random_limbs(limbs_n, nn);
      fill_random_limbs(limbs_d, dn);
      limbs_d[dn - 1] = top_vals[i];

      CHECK_EQ_I(set_from_limbs(n, limbs_n, nn), ARBINT_OK);
      CHECK_EQ_I(set_from_limbs(d, limbs_d, dn), ARBINT_OK);
      CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);

      verify_division(n, d, q, r, &g_failures);

      free(limbs_n);
      free(limbs_d);
      limbs_n = NULL;
      limbs_d = NULL;
    }
  }

  printf("Test 4: Adversarial - exact multiples and near-multiples\n");
  {
    size_t dn = 4;

    limbs_d = (arbint_limb_t *) malloc(dn * sizeof(arbint_limb_t));
    CHECK(limbs_d != NULL);
    fill_random_limbs(limbs_d, dn);
    CHECK_EQ_I(set_from_limbs(d, limbs_d, dn), ARBINT_OK);

    /*  n = q * d (exact multiple, remainder = 0)  */
    CHECK_EQ_I(arbint_set_i32(q, 12345), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(n, q, d), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
    CHECK(arbint_is_zero(r));
    verify_division(n, d, q, r, &g_failures);

    /*  n = q * d + 1 (remainder = 1)  */
    CHECK_EQ_I(arbint_set_i32(q, 12345), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(n, q, d), ARBINT_OK);
    CHECK_EQ_I(arbint_add_i32(n, n, 1), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
    check_i32_value(r, 1);
    verify_division(n, d, q, r, &g_failures);

    /*  n = q * d + (d - 1) (remainder = d - 1)  */
    CHECK_EQ_I(arbint_set_i32(q, 12345), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(n, q, d), ARBINT_OK);
    CHECK_EQ_I(arbint_add(n, n, d), ARBINT_OK);
    CHECK_EQ_I(arbint_sub_i32(n, n, 1), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
    {
      arbint_t d_minus_1;
      CHECK_EQ_I(arbint_init(d_minus_1, &ctx), ARBINT_OK);
      CHECK_EQ_I(arbint_sub_i32(d_minus_1, d, 1), ARBINT_OK);
      CHECK(arbint_cmp(r, d_minus_1) == 0);
      arbint_clear(d_minus_1);
    }
    verify_division(n, d, q, r, &g_failures);

    free(limbs_d);
    limbs_d = NULL;
  }

  printf("Test 5: Large nn/dn ratios\n");
  {
    size_t dn = 4;
    size_t nn_vals[] = {dn, 2 * dn, 5 * dn, 10 * dn, 20 * dn};
    size_t num_nn = sizeof(nn_vals) / sizeof(nn_vals[0]);
    size_t i;

    limbs_d = (arbint_limb_t *) malloc(dn * sizeof(arbint_limb_t));
    CHECK(limbs_d != NULL);
    fill_random_limbs(limbs_d, dn);
    CHECK_EQ_I(set_from_limbs(d, limbs_d, dn), ARBINT_OK);

    for (i = 0; i < num_nn; ++i) {
      size_t nn = nn_vals[i];

      limbs_n = (arbint_limb_t *) malloc(nn * sizeof(arbint_limb_t));
      CHECK(limbs_n != NULL);
      fill_random_limbs(limbs_n, nn);

      CHECK_EQ_I(set_from_limbs(n, limbs_n, nn), ARBINT_OK);
      CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);

      verify_division(n, d, q, r, &g_failures);

      free(limbs_n);
      limbs_n = NULL;
    }

    free(limbs_d);
    limbs_d = NULL;
  }

  printf("Test 6: Stochastic testing\n");
  {
    for (trial = 0; trial < 200; ++trial) {
      size_t dn = 2 + (trial % 20);
      size_t nn = dn + (trial % (dn + 1));

      limbs_n = (arbint_limb_t *) malloc(nn * sizeof(arbint_limb_t));
      limbs_d = (arbint_limb_t *) malloc(dn * sizeof(arbint_limb_t));
      CHECK(limbs_n != NULL && limbs_d != NULL);

      fill_random_limbs(limbs_n, nn);
      fill_random_limbs(limbs_d, dn);

      CHECK_EQ_I(set_from_limbs(n, limbs_n, nn), ARBINT_OK);
      CHECK_EQ_I(set_from_limbs(d, limbs_d, dn), ARBINT_OK);
      CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);

      verify_division(n, d, q, r, &g_failures);

      free(limbs_n);
      free(limbs_d);
      limbs_n = NULL;
      limbs_d = NULL;
    }
  }

  printf("Test 7: Newton-specific sizes (large divisors)\n");
  {
    /*  Test at sizes where Newton should be used.  */
    size_t dn_vals[] = {64, 100, 128, 200, 256};
    size_t num_dn = sizeof(dn_vals) / sizeof(dn_vals[0]);
    size_t i;

    for (i = 0; i < num_dn; ++i) {
      size_t dn = dn_vals[i];
      size_t nn = dn + dn / 2;

      limbs_n = (arbint_limb_t *) malloc(nn * sizeof(arbint_limb_t));
      limbs_d = (arbint_limb_t *) malloc(dn * sizeof(arbint_limb_t));
      CHECK(limbs_n != NULL && limbs_d != NULL);

      fill_random_limbs(limbs_n, nn);
      fill_random_limbs(limbs_d, dn);

      CHECK_EQ_I(set_from_limbs(n, limbs_n, nn), ARBINT_OK);
      CHECK_EQ_I(set_from_limbs(d, limbs_d, dn), ARBINT_OK);
      CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);

      verify_division(n, d, q, r, &g_failures);

      free(limbs_n);
      free(limbs_d);
      limbs_n = NULL;
      limbs_d = NULL;
    }
  }

  printf("Test 8: All-ones dividend and divisor\n");
  {
    size_t dn = 8;
    size_t nn = 16;
    size_t i;

    limbs_n = (arbint_limb_t *) malloc(nn * sizeof(arbint_limb_t));
    limbs_d = (arbint_limb_t *) malloc(dn * sizeof(arbint_limb_t));
    CHECK(limbs_n != NULL && limbs_d != NULL);

    for (i = 0; i < nn; ++i)
      limbs_n[i] = ~(arbint_limb_t) 0;
    for (i = 0; i < dn; ++i)
      limbs_d[i] = ~(arbint_limb_t) 0;

    CHECK_EQ_I(set_from_limbs(n, limbs_n, nn), ARBINT_OK);
    CHECK_EQ_I(set_from_limbs(d, limbs_d, dn), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);

    verify_division(n, d, q, r, &g_failures);

    free(limbs_n);
    free(limbs_d);
    limbs_n = NULL;
    limbs_d = NULL;
  }

  arbint_clear(r);
  arbint_clear(q);
  arbint_clear(d);
  arbint_clear(n);
  arbint_ctx_clear(&ctx);

  ARBINT_TEST_FINISH("test_div_newton");
}
