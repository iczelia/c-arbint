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

#include "test_framework.h"

ARBINT_TEST_DECLARE_FAILURES();

/*  Helper: verify that a == p^k * rest.  */
static void verify_factorization(const arbint_t a, uint32_t p, uint32_t k,
                                 const arbint_t rest, arbint_t tmp,
                                 const char * desc) {
  arbint_err_t rc;

  /*  tmp = p^k.  */
  rc = arbint_set_u32(tmp, p);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "FAIL %s: set p failed with %d\n", desc, (int) rc);
    ++g_failures;
    return;
  }

  rc = arbint_pow_u32(tmp, tmp, k);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "FAIL %s: p^k failed with %d\n", desc, (int) rc);
    ++g_failures;
    return;
  }

  /*  tmp = p^k * rest.  */
  rc = arbint_mul(tmp, tmp, rest);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "FAIL %s: p^k * rest failed with %d\n", desc, (int) rc);
    ++g_failures;
    return;
  }

  /*  Verify tmp == a.  */
  if (!arbint_eq(tmp, a)) {
    fprintf(stderr, "FAIL %s: a != p^k * rest\n", desc);
    ++g_failures;
  }
}

/*  Helper: verify that rest is not divisible by p.  */
static void verify_coprime(const arbint_t rest, uint32_t p,
                           const char * desc) {
  int divisible;
  arbint_err_t rc;

  if (arbint_is_zero(rest))
    return;

  rc = arbint_divisible_u32(rest, p, &divisible);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "FAIL %s: divisible check failed with %d\n", desc,
            (int) rc);
    ++g_failures;
    return;
  }

  if (divisible) {
    fprintf(stderr, "FAIL %s: rest is still divisible by p\n", desc);
    ++g_failures;
  }
}

int main(void) {
  ARBINT_TEST_START();
  arbint_ctx_t ctx;
  arbint_t a, rest, tmp;
  uint32_t k;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init_all(&ctx, a, rest, tmp, (arbint_t *) NULL), ARBINT_OK);

  /*  ================================================================
      Basic functionality tests
      ================================================================  */

  /*  12 = 2^2 * 3.  */
  CHECK_EQ_I(arbint_set_i32(a, 12), ARBINT_OK);
  CHECK_EQ_I(arbint_removefactor_u32(rest, a, 2u, &k), ARBINT_OK);
  CHECK_EQ_I(k, 2u);
  check_i32_value(rest, 3);
  verify_factorization(a, 2u, k, rest, tmp, "12 = 2^k * rest");
  verify_coprime(rest, 2u, "12 = 2^k * rest");

  /*  12 = 3^1 * 4.  */
  CHECK_EQ_I(arbint_removefactor_u32(rest, a, 3u, &k), ARBINT_OK);
  CHECK_EQ_I(k, 1u);
  check_i32_value(rest, 4);
  verify_factorization(a, 3u, k, rest, tmp, "12 = 3^k * rest");
  verify_coprime(rest, 3u, "12 = 3^k * rest");

  /*  12 = 5^0 * 12 (5 does not divide 12).  */
  CHECK_EQ_I(arbint_removefactor_u32(rest, a, 5u, &k), ARBINT_OK);
  CHECK_EQ_I(k, 0u);
  check_i32_value(rest, 12);

  /*  1000 = 2^3 * 125.  */
  CHECK_EQ_I(arbint_set_i32(a, 1000), ARBINT_OK);
  CHECK_EQ_I(arbint_removefactor_u32(rest, a, 2u, &k), ARBINT_OK);
  CHECK_EQ_I(k, 3u);
  check_i32_value(rest, 125);
  verify_factorization(a, 2u, k, rest, tmp, "1000 = 2^k * rest");

  /*  1000 = 5^3 * 8.  */
  CHECK_EQ_I(arbint_removefactor_u32(rest, a, 5u, &k), ARBINT_OK);
  CHECK_EQ_I(k, 3u);
  check_i32_value(rest, 8);
  verify_factorization(a, 5u, k, rest, tmp, "1000 = 5^k * rest");

  /*  ================================================================
      Edge cases
      ================================================================  */

  /*  a == 0: no factors to remove.  */
  arbint_zero(a);
  CHECK_EQ_I(arbint_removefactor_u32(rest, a, 2u, &k), ARBINT_OK);
  CHECK_EQ_I(k, 0u);
  CHECK(arbint_is_zero(rest));

  /*  a == 1: no factors of any p > 1.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_removefactor_u32(rest, a, 2u, &k), ARBINT_OK);
  CHECK_EQ_I(k, 0u);
  check_i32_value(rest, 1);

  CHECK_EQ_I(arbint_removefactor_u32(rest, a, 7u, &k), ARBINT_OK);
  CHECK_EQ_I(k, 0u);
  check_i32_value(rest, 1);

  /*  p == 1: divides everything, but we treat as removing 0 factors.  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_removefactor_u32(rest, a, 1u, &k), ARBINT_OK);
  CHECK_EQ_I(k, 0u);
  check_i32_value(rest, 42);

  /*  a is a power of p: rest should be 1.  */
  CHECK_EQ_I(arbint_set_i32(a, 64), ARBINT_OK);  /*  2^6  */
  CHECK_EQ_I(arbint_removefactor_u32(rest, a, 2u, &k), ARBINT_OK);
  CHECK_EQ_I(k, 6u);
  check_i32_value(rest, 1);
  verify_factorization(a, 2u, k, rest, tmp, "64 = 2^k * rest");

  CHECK_EQ_I(arbint_set_i32(a, 243), ARBINT_OK);  /*  3^5  */
  CHECK_EQ_I(arbint_removefactor_u32(rest, a, 3u, &k), ARBINT_OK);
  CHECK_EQ_I(k, 5u);
  check_i32_value(rest, 1);
  verify_factorization(a, 3u, k, rest, tmp, "243 = 3^k * rest");

  /*  ================================================================
      Negative inputs
      ================================================================  */

  /*  -12 = 2^2 * (-3).  */
  CHECK_EQ_I(arbint_set_i32(a, -12), ARBINT_OK);
  CHECK_EQ_I(arbint_removefactor_u32(rest, a, 2u, &k), ARBINT_OK);
  CHECK_EQ_I(k, 2u);
  check_i32_value(rest, -3);
  verify_factorization(a, 2u, k, rest, tmp, "-12 = 2^k * rest");

  /*  -1000 = 5^3 * (-8).  */
  CHECK_EQ_I(arbint_set_i32(a, -1000), ARBINT_OK);
  CHECK_EQ_I(arbint_removefactor_u32(rest, a, 5u, &k), ARBINT_OK);
  CHECK_EQ_I(k, 3u);
  check_i32_value(rest, -8);
  verify_factorization(a, 5u, k, rest, tmp, "-1000 = 5^k * rest");

  /*  ================================================================
      Error conditions
      ================================================================  */

  CHECK_EQ_I(arbint_removefactor_u32(NULL, a, 2u, &k), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_removefactor_u32(rest, NULL, 2u, &k), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_removefactor_u32(rest, a, 2u, NULL), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_removefactor_u32(rest, a, 0u, &k), ARBINT_EZERO);

  /*  ================================================================
      Aliasing tests
      ================================================================  */

  /*  rest == a: in-place update.  */
  CHECK_EQ_I(arbint_set_i32(rest, 1000), ARBINT_OK);
  CHECK_EQ_I(arbint_removefactor_u32(rest, rest, 2u, &k), ARBINT_OK);
  CHECK_EQ_I(k, 3u);
  check_i32_value(rest, 125);

  /*  ================================================================
      Large number tests
      ================================================================  */

  /*  a = 2^100.  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_removefactor_u32(rest, a, 2u, &k), ARBINT_OK);
  CHECK_EQ_I(k, 100u);
  check_i32_value(rest, 1);
  verify_factorization(a, 2u, k, rest, tmp, "2^100 = 2^k * rest");

  /*  a = 2^50 * 3^30.  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 50), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(tmp, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(tmp, tmp, 30), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(a, a, tmp), ARBINT_OK);

  CHECK_EQ_I(arbint_removefactor_u32(rest, a, 2u, &k), ARBINT_OK);
  CHECK_EQ_I(k, 50u);
  /*  rest should be 3^30.  */
  CHECK_EQ_I(arbint_set_i32(tmp, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(tmp, tmp, 30), ARBINT_OK);
  CHECK(arbint_eq(rest, tmp));
  verify_coprime(rest, 2u, "2^50 * 3^30, remove 2s");

  CHECK_EQ_I(arbint_removefactor_u32(rest, a, 3u, &k), ARBINT_OK);
  CHECK_EQ_I(k, 30u);
  /*  rest should be 2^50.  */
  CHECK_EQ_I(arbint_set_i32(tmp, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(tmp, tmp, 50), ARBINT_OK);
  CHECK(arbint_eq(rest, tmp));
  verify_coprime(rest, 3u, "2^50 * 3^30, remove 3s");

  /*  ================================================================
      Stochastic tests
      ================================================================  */
  {
    uint32_t primes[] = {2, 3, 5, 7, 11, 13};
    size_t num_primes = sizeof(primes) / sizeof(primes[0]);

    for (size_t i = 0; i < num_primes; ++i) {
      uint32_t p = primes[i];

      for (uint32_t exp = 0; exp <= 10; ++exp) {
        char desc[64];

        /*  a = p^exp * 17 (17 is coprime to all test primes except itself).  */
        CHECK_EQ_I(arbint_set_u32(a, p), ARBINT_OK);
        CHECK_EQ_I(arbint_pow_u32(a, a, exp), ARBINT_OK);
        CHECK_EQ_I(arbint_mul_u32(a, a, 17u), ARBINT_OK);

        snprintf(desc, sizeof(desc), "%u^%u * 17", p, exp);

        CHECK_EQ_I(arbint_removefactor_u32(rest, a, p, &k), ARBINT_OK);
        if (p == 17u) {
          CHECK_EQ_I(k, exp + 1u);
          check_i32_value(rest, 1);
        } else {
          CHECK_EQ_I(k, exp);
          check_i32_value(rest, 17);
        }
        verify_factorization(a, p, k, rest, tmp, desc);
        verify_coprime(rest, p, desc);
      }
    }
  }

  arbint_clear_all(a, rest, tmp, (arbint_t *) NULL);
  arbint_ctx_clear(&ctx);

  ARBINT_TEST_FINISH("test_removefactor");
}
