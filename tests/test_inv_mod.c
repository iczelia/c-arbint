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

/*  Helper: verify that rop is the modular inverse of a mod |m|.
    Checks that (rop * a) mod |m| == 1.  */
static void verify_inverse(const arbint_t rop, const arbint_t a,
                           const arbint_t m, arbint_t tmp1, arbint_t tmp2,
                           const char * desc) {
  arbint_err_t rc;

  /*  tmp1 = rop * a.  */
  rc = arbint_mul(tmp1, rop, a);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "FAIL %s: rop * a failed with %d\n", desc, (int) rc);
    ++g_failures;
    return;
  }

  /*  Use |m| for verification to get positive canonical remainder.  */
  rc = arbint_abs(tmp2, m);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "FAIL %s: abs(m) failed with %d\n", desc, (int) rc);
    ++g_failures;
    return;
  }

  /*  tmp1 = tmp1 mod |m| (floor division remainder for canonical result).  */
  rc = arbint_fdiv_r(tmp1, tmp1, tmp2);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "FAIL %s: (rop * a) mod |m| failed with %d\n", desc,
            (int) rc);
    ++g_failures;
    return;
  }

  /*  Verify tmp1 == 1.  */
  if (!arbint_is_one(tmp1)) {
    fprintf(stderr, "FAIL %s: (rop * a) mod |m| != 1\n", desc);
    ++g_failures;
  }
}

/*  Helper: verify inverse with u32 modulus.  */
static void verify_inverse_u32(const arbint_t rop, const arbint_t a,
                               uint32_t m, arbint_t tmp1, arbint_t tmp2,
                               const char * desc) {
  arbint_err_t rc;

  rc = arbint_mul(tmp1, rop, a);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "FAIL %s: rop * a failed with %d\n", desc, (int) rc);
    ++g_failures;
    return;
  }

  rc = arbint_fdiv_r_u32(tmp2, tmp1, m);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "FAIL %s: (rop * a) mod m failed with %d\n", desc,
            (int) rc);
    ++g_failures;
    return;
  }

  if (!arbint_is_one(tmp2)) {
    fprintf(stderr, "FAIL %s: (rop * a) mod m != 1\n", desc);
    ++g_failures;
  }
}

/*  Test arbint_inv_mod with i32 inputs.  */
static void test_inv_mod_i32(arbint_t rop, arbint_t a, arbint_t m, arbint_t tmp1,
                             arbint_t tmp2, int32_t av, int32_t mv,
                             const char * desc) {
  arbint_err_t rc;

  CHECK_EQ_I(arbint_set_i32(a, av), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(m, mv), ARBINT_OK);

  rc = arbint_inv_mod(rop, a, m);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "FAIL %s: arbint_inv_mod(%d, %d) failed with %d\n", desc,
            av, mv, (int) rc);
    ++g_failures;
    return;
  }

  /*  Result should be in [0, |m|).  */
  if (arbint_is_neg(rop)) {
    fprintf(stderr, "FAIL %s: result is negative\n", desc);
    ++g_failures;
    return;
  }

  verify_inverse(rop, a, m, tmp1, tmp2, desc);
}

/*  Test arbint_inv_mod_u32 with i32 input and u32 modulus.  */
static void test_inv_mod_u32_i32(arbint_t rop, arbint_t a, arbint_t tmp1,
                                 arbint_t tmp2, int32_t av, uint32_t mv,
                                 const char * desc) {
  arbint_err_t rc;

  CHECK_EQ_I(arbint_set_i32(a, av), ARBINT_OK);

  rc = arbint_inv_mod_u32(rop, a, mv);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "FAIL %s: arbint_inv_mod_u32(%d, %u) failed with %d\n", desc,
            av, mv, (int) rc);
    ++g_failures;
    return;
  }

  if (arbint_is_neg(rop)) {
    fprintf(stderr, "FAIL %s: result is negative\n", desc);
    ++g_failures;
    return;
  }

  verify_inverse_u32(rop, a, mv, tmp1, tmp2, desc);
}

int main(void) {
  arbint_ctx_t ctx;
  arbint_t rop, a, m, tmp1, tmp2;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init_all(&ctx, rop, a, m, tmp1, tmp2, (arbint_t *) NULL),
             ARBINT_OK);

  /*  ================================================================
      Basic functionality tests for arbint_inv_mod
      ================================================================  */

  /*  inv(3, 7) = 5 since 3 * 5 = 15 = 1 (mod 7).  */
  test_inv_mod_i32(rop, a, m, tmp1, tmp2, 3, 7, "inv(3, 7)");
  check_i32_value(rop, 5);

  /*  inv(2, 5) = 3 since 2 * 3 = 6 = 1 (mod 5).  */
  test_inv_mod_i32(rop, a, m, tmp1, tmp2, 2, 5, "inv(2, 5)");
  check_i32_value(rop, 3);

  /*  inv(17, 31): verify result via multiplication.  */
  test_inv_mod_i32(rop, a, m, tmp1, tmp2, 17, 31, "inv(17, 31)");

  /*  inv(7, 11) = 8 since 7 * 8 = 56 = 1 (mod 11).  */
  test_inv_mod_i32(rop, a, m, tmp1, tmp2, 7, 11, "inv(7, 11)");
  check_i32_value(rop, 8);

  /*  inv(13, 17) = 4 since 13 * 4 = 52 = 3 * 17 + 1 = 1 (mod 17).  */
  test_inv_mod_i32(rop, a, m, tmp1, tmp2, 13, 17, "inv(13, 17)");
  check_i32_value(rop, 4);

  /*  ================================================================
      Edge cases
      ================================================================  */

  /*  inv(any, 1) = 0: any integer is 0 mod 1.  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(m, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_inv_mod(rop, a, m), ARBINT_OK);
  check_i32_value(rop, 0);

  /*  inv(1, m) = 1 for any m > 1.  */
  test_inv_mod_i32(rop, a, m, tmp1, tmp2, 1, 17, "inv(1, 17)");
  check_i32_value(rop, 1);

  /*  inv(-1, m) = m - 1 since (-1) * (m-1) = -(m-1) = 1 - m = 1 (mod m).  */
  test_inv_mod_i32(rop, a, m, tmp1, tmp2, -1, 7, "inv(-1, 7)");
  check_i32_value(rop, 6);

  /*  inv(m-1, m) = m - 1 since (m-1)^2 = m^2 - 2m + 1 = 1 (mod m).  */
  test_inv_mod_i32(rop, a, m, tmp1, tmp2, 10, 11, "inv(10, 11)");
  check_i32_value(rop, 10);

  /*  ================================================================
      Error conditions
      ================================================================  */

  /*  gcd(a, m) != 1: inverse does not exist.  */
  CHECK_EQ_I(arbint_set_i32(a, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(m, 6), ARBINT_OK);
  CHECK_EQ_I(arbint_inv_mod(rop, a, m), ARBINT_EDOM);

  CHECK_EQ_I(arbint_set_i32(a, 6), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(m, 9), ARBINT_OK);
  CHECK_EQ_I(arbint_inv_mod(rop, a, m), ARBINT_EDOM);

  /*  a == 0: gcd(0, m) = m != 1 (for m > 1).  */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(m, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_inv_mod(rop, a, m), ARBINT_EDOM);

  /*  mod == 0: division by zero.  */
  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  arbint_zero(m);
  CHECK_EQ_I(arbint_inv_mod(rop, a, m), ARBINT_EZERO);

  /*  NULL inputs.  */
  CHECK_EQ_I(arbint_inv_mod(NULL, a, m), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_set_i32(m, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_inv_mod(rop, NULL, m), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_inv_mod(rop, a, NULL), ARBINT_EINVAL);

  /*  ================================================================
      Negative input tests
      ================================================================  */

  /*  inv(-3, 7): -3 = 4 (mod 7), inv(4, 7) = 2 since 4 * 2 = 8 = 1 (mod 7).  */
  test_inv_mod_i32(rop, a, m, tmp1, tmp2, -3, 7, "inv(-3, 7)");
  check_i32_value(rop, 2);

  /*  inv(-2, 5): -2 = 3 (mod 5), inv(3, 5) = 2 since 3 * 2 = 6 = 1 (mod 5).  */
  test_inv_mod_i32(rop, a, m, tmp1, tmp2, -2, 5, "inv(-2, 5)");
  check_i32_value(rop, 2);

  /*  Negative modulus: result should be in [0, |mod|).  */
  test_inv_mod_i32(rop, a, m, tmp1, tmp2, 3, -7, "inv(3, -7)");
  check_i32_value(rop, 5);

  test_inv_mod_i32(rop, a, m, tmp1, tmp2, -3, -7, "inv(-3, -7)");
  check_i32_value(rop, 2);

  /*  ================================================================
      Aliasing tests
      ================================================================  */

  /*  rop == a: in-place update.  */
  CHECK_EQ_I(arbint_set_i32(rop, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(m, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_inv_mod(rop, rop, m), ARBINT_OK);
  check_i32_value(rop, 5);

  /*  rop == m: in-place update.  */
  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(rop, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_inv_mod(rop, a, rop), ARBINT_OK);
  check_i32_value(rop, 5);

  /*  ================================================================
      Large number tests
      ================================================================  */

  /*  a = 2^64 + 1, m = large prime.  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 64), ARBINT_OK);
  CHECK_EQ_I(arbint_add_i32(a, a, 1), ARBINT_OK);

  /*  Use 65537 (Fermat prime F4) as modulus.  */
  CHECK_EQ_I(arbint_set_i32(m, 65537), ARBINT_OK);
  CHECK_EQ_I(arbint_inv_mod(rop, a, m), ARBINT_OK);
  CHECK(!arbint_is_neg(rop));
  verify_inverse(rop, a, m, tmp1, tmp2, "inv(2^64+1, 65537)");

  /*  Large modulus: 2^128 - 159 is prime.  */
  CHECK_EQ_I(arbint_set_i32(m, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(m, m, 128), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_i32(m, m, 159), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(a, 12345), ARBINT_OK);
  CHECK_EQ_I(arbint_inv_mod(rop, a, m), ARBINT_OK);
  CHECK(!arbint_is_neg(rop));
  verify_inverse(rop, a, m, tmp1, tmp2, "inv(12345, 2^128-159)");

  /*  ================================================================
      arbint_inv_mod_u32 tests
      ================================================================  */

  /*  Basic functionality.  */
  test_inv_mod_u32_i32(rop, a, tmp1, tmp2, 3, 7u, "inv_u32(3, 7)");
  check_i32_value(rop, 5);

  test_inv_mod_u32_i32(rop, a, tmp1, tmp2, 2, 5u, "inv_u32(2, 5)");
  check_i32_value(rop, 3);

  test_inv_mod_u32_i32(rop, a, tmp1, tmp2, 17, 31u, "inv_u32(17, 31)");

  /*  Edge cases.  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_inv_mod_u32(rop, a, 1u), ARBINT_OK);
  check_i32_value(rop, 0);

  test_inv_mod_u32_i32(rop, a, tmp1, tmp2, 1, 17u, "inv_u32(1, 17)");
  check_i32_value(rop, 1);

  /*  Error conditions.  */
  CHECK_EQ_I(arbint_set_i32(a, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_inv_mod_u32(rop, a, 6u), ARBINT_EDOM);

  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_inv_mod_u32(rop, a, 7u), ARBINT_EDOM);

  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_inv_mod_u32(rop, a, 0u), ARBINT_EZERO);

  CHECK_EQ_I(arbint_inv_mod_u32(NULL, a, 7u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_inv_mod_u32(rop, NULL, 7u), ARBINT_EINVAL);

  /*  Negative a with u32 modulus.  */
  test_inv_mod_u32_i32(rop, a, tmp1, tmp2, -3, 7u, "inv_u32(-3, 7)");
  check_i32_value(rop, 2);

  /*  Large a with small mod.  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 256), ARBINT_OK);
  CHECK_EQ_I(arbint_inv_mod_u32(rop, a, 17u), ARBINT_OK);
  CHECK(!arbint_is_neg(rop));
  verify_inverse_u32(rop, a, 17u, tmp1, tmp2, "inv_u32(2^256, 17)");

  /*  Aliasing: rop == a.  */
  CHECK_EQ_I(arbint_set_i32(rop, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_inv_mod_u32(rop, rop, 7u), ARBINT_OK);
  check_i32_value(rop, 5);

  /*  ================================================================
      Stochastic tests: verify random coprime pairs
      ================================================================  */
  {
    uint32_t primes[] = {3,  5,  7,  11, 13, 17, 19, 23, 29, 31,
                         37, 41, 43, 47, 53, 59, 61, 67, 71, 73};
    size_t num_primes = sizeof(primes) / sizeof(primes[0]);

    for (size_t i = 0; i < num_primes; ++i) {
      uint32_t p = primes[i];

      for (int32_t av = 1; av < (int32_t) p; ++av) {
        char desc[64];
        snprintf(desc, sizeof(desc), "stochastic inv(%d, %u)", av, p);

        CHECK_EQ_I(arbint_set_i32(a, av), ARBINT_OK);
        CHECK_EQ_I(arbint_set_u32(m, p), ARBINT_OK);
        CHECK_EQ_I(arbint_inv_mod(rop, a, m), ARBINT_OK);
        verify_inverse(rop, a, m, tmp1, tmp2, desc);

        CHECK_EQ_I(arbint_inv_mod_u32(rop, a, p), ARBINT_OK);
        verify_inverse_u32(rop, a, p, tmp1, tmp2, desc);
      }
    }
  }

  arbint_clear_all(rop, a, m, tmp1, tmp2, (arbint_t *) NULL);
  arbint_ctx_clear(&ctx);

  ARBINT_TEST_FINISH("test_inv_mod");
}
