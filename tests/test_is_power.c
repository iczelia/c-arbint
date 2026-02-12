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

/*  Helper: test that arbint_is_power(a, &out) returns expected result.  */
static void test_is_power_value(const arbint_t a, int expected_out,
                                const char * desc) {
  int out = -1;
  arbint_err_t rc = arbint_is_power(a, &out);
  if (rc != ARBINT_OK) {
    fprintf(stderr, "FAIL %s: arbint_is_power returned %d\n", desc, (int) rc);
    ++g_failures;
  } else if (out != expected_out) {
    fprintf(stderr, "FAIL %s: expected out=%d, got out=%d\n", desc,
            expected_out, out);
    ++g_failures;
  }
}

/*  Helper: test is_power for an i32 value.  */
static void test_is_power_i32(arbint_t a, int32_t val, int expected_out,
                              const char * desc) {
  arbint_err_t rc = arbint_set_i32(a, val);
  CHECK_EQ_I(rc, ARBINT_OK);
  test_is_power_value(a, expected_out, desc);
}

int main(void) {
  arbint_ctx_t ctx;
  arbint_t a, b;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);

  /*  Test NULL inputs.  */
  {
    int out;
    CHECK_EQ_I(arbint_is_power(NULL, &out), ARBINT_EINVAL);
    CHECK_EQ_I(arbint_is_power(a, NULL), ARBINT_EINVAL);
  }

  /*  Test boundary values.  */
  test_is_power_i32(a, 0, 1, "0 = 0^k");
  test_is_power_i32(a, 1, 1, "1 = 1^k");
  test_is_power_i32(a, -1, 1, "-1 = (-1)^3");

  /*  Test small perfect squares.  */
  test_is_power_i32(a, 4, 1, "4 = 2^2");
  test_is_power_i32(a, 9, 1, "9 = 3^2");
  test_is_power_i32(a, 16, 1, "16 = 2^4 = 4^2");
  test_is_power_i32(a, 25, 1, "25 = 5^2");
  test_is_power_i32(a, 36, 1, "36 = 6^2");
  test_is_power_i32(a, 49, 1, "49 = 7^2");
  test_is_power_i32(a, 64, 1, "64 = 2^6 = 4^3 = 8^2");
  test_is_power_i32(a, 81, 1, "81 = 3^4 = 9^2");
  test_is_power_i32(a, 100, 1, "100 = 10^2");
  test_is_power_i32(a, 121, 1, "121 = 11^2");
  test_is_power_i32(a, 144, 1, "144 = 12^2");

  /*  Test small perfect cubes.  */
  test_is_power_i32(a, 8, 1, "8 = 2^3");
  test_is_power_i32(a, 27, 1, "27 = 3^3");
  test_is_power_i32(a, 125, 1, "125 = 5^3");
  test_is_power_i32(a, 216, 1, "216 = 6^3");
  test_is_power_i32(a, 343, 1, "343 = 7^3");
  test_is_power_i32(a, 512, 1, "512 = 8^3 = 2^9");
  test_is_power_i32(a, 729, 1, "729 = 3^6 = 9^3 = 27^2");
  test_is_power_i32(a, 1000, 1, "1000 = 10^3");

  /*  Test negative odd powers.  */
  test_is_power_i32(a, -8, 1, "-8 = (-2)^3");
  test_is_power_i32(a, -27, 1, "-27 = (-3)^3");
  test_is_power_i32(a, -125, 1, "-125 = (-5)^3");
  test_is_power_i32(a, -216, 1, "-216 = (-6)^3");
  test_is_power_i32(a, -343, 1, "-343 = (-7)^3");
  test_is_power_i32(a, -512, 1, "-512 = (-8)^3");
  test_is_power_i32(a, -729, 1, "-729 = (-9)^3");
  test_is_power_i32(a, -1000, 1, "-1000 = (-10)^3");

  /*  Test negative non-powers (no odd root).  */
  test_is_power_i32(a, -4, 0, "-4 not a perfect power");
  test_is_power_i32(a, -9, 0, "-9 not a perfect power");
  test_is_power_i32(a, -16, 0, "-16 not a perfect power");
  test_is_power_i32(a, -25, 0, "-25 not a perfect power");
  test_is_power_i32(a, -100, 0, "-100 not a perfect power");

  /*  Test higher powers.  */
  test_is_power_i32(a, 32, 1, "32 = 2^5");
  test_is_power_i32(a, 243, 1, "243 = 3^5");
  test_is_power_i32(a, 128, 1, "128 = 2^7");
  test_is_power_i32(a, 256, 1, "256 = 2^8 = 4^4 = 16^2");

  /*  Test negative higher odd powers.  */
  test_is_power_i32(a, -32, 1, "-32 = (-2)^5");
  test_is_power_i32(a, -243, 1, "-243 = (-3)^5");
  test_is_power_i32(a, -128, 1, "-128 = (-2)^7");

  /*  Test non-powers (primes and composite non-powers).  */
  test_is_power_i32(a, 2, 0, "2 is prime, not a perfect power");
  test_is_power_i32(a, 3, 0, "3 is prime, not a perfect power");
  test_is_power_i32(a, 5, 0, "5 is prime, not a perfect power");
  test_is_power_i32(a, 6, 0, "6 = 2*3, not a perfect power");
  test_is_power_i32(a, 7, 0, "7 is prime, not a perfect power");
  test_is_power_i32(a, 10, 0, "10 = 2*5, not a perfect power");
  test_is_power_i32(a, 11, 0, "11 is prime, not a perfect power");
  test_is_power_i32(a, 12, 0, "12 = 4*3, not a perfect power");
  test_is_power_i32(a, 13, 0, "13 is prime, not a perfect power");
  test_is_power_i32(a, 14, 0, "14 = 2*7, not a perfect power");
  test_is_power_i32(a, 15, 0, "15 = 3*5, not a perfect power");
  test_is_power_i32(a, 17, 0, "17 is prime, not a perfect power");
  test_is_power_i32(a, 18, 0, "18 = 2*9, not a perfect power");
  test_is_power_i32(a, 19, 0, "19 is prime, not a perfect power");
  test_is_power_i32(a, 20, 0, "20 = 4*5, not a perfect power");

  /*  Test large prime exponents (within static table).  */
  /*  2^127: use repeated squaring to build it.  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 127), ARBINT_OK);
  test_is_power_value(a, 1, "2^127");

  /*  3^67.  */
  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 67), ARBINT_OK);
  test_is_power_value(a, 1, "3^67");

  /*  2^257.  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 257), ARBINT_OK);
  test_is_power_value(a, 1, "2^257");

  /*  5^101.  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 101), ARBINT_OK);
  test_is_power_value(a, 1, "5^101");

  /*  Test large prime exponents (requiring dynamic sieve).  */
  /*  2^523 (523 is prime, beyond static table which ends at 521).  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 523), ARBINT_OK);
  test_is_power_value(a, 1, "2^523");

  /*  3^541 (541 is prime).  */
  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 541), ARBINT_OK);
  test_is_power_value(a, 1, "3^541");

  /*  2^1009 (1009 is prime).  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 1009), ARBINT_OK);
  test_is_power_value(a, 1, "2^1009");

  /*  Drop caches and verify dynamic sieve path still works after rebuild.  */
  arbint_drop_caches();
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 541), ARBINT_OK);
  test_is_power_value(a, 1, "2^541 after cache drop");

  /*  Test large non-powers.  */
  /*  2^127 + 1.  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 127), ARBINT_OK);
  CHECK_EQ_I(arbint_add_i32(a, a, 1), ARBINT_OK);
  test_is_power_value(a, 0, "2^127 + 1");

  /*  2^127 - 1 (Mersenne number, but 127-bit, not a perfect power).  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 127), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_i32(a, a, 1), ARBINT_OK);
  test_is_power_value(a, 0, "2^127 - 1");

  /*  2^523 + 1.  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 523), ARBINT_OK);
  CHECK_EQ_I(arbint_add_i32(a, a, 1), ARBINT_OK);
  test_is_power_value(a, 0, "2^523 + 1");

  /*  Test negative large odd powers.  */
  CHECK_EQ_I(arbint_set_i32(a, -2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(b, a, 127), ARBINT_OK);
  /*  (-2)^127 = -(2^127) because 127 is odd.  */
  test_is_power_value(b, 1, "(-2)^127");

  CHECK_EQ_I(arbint_set_i32(a, -3), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(b, a, 67), ARBINT_OK);
  test_is_power_value(b, 1, "(-3)^67");

  /*  Test 4th, 5th, 6th, 7th powers.  */
  test_is_power_i32(a, 16, 1, "16 = 2^4");
  test_is_power_i32(a, 81, 1, "81 = 3^4");
  test_is_power_i32(a, 625, 1, "625 = 5^4");
  test_is_power_i32(a, 1296, 1, "1296 = 6^4");
  test_is_power_i32(a, 2401, 1, "2401 = 7^4");
  test_is_power_i32(a, 3125, 1, "3125 = 5^5");
  test_is_power_i32(a, 7776, 1, "7776 = 6^5");
  test_is_power_i32(a, 46656, 1, "46656 = 6^6");
  test_is_power_i32(a, 78125, 1, "78125 = 5^7");
  test_is_power_i32(a, 2187, 1, "2187 = 3^7");

  /*  Negative 5th and 7th powers.  */
  test_is_power_i32(a, -3125, 1, "-3125 = (-5)^5");
  test_is_power_i32(a, -7776, 1, "-7776 = (-6)^5");
  test_is_power_i32(a, -78125, 1, "-78125 = (-5)^7");
  test_is_power_i32(a, -2187, 1, "-2187 = (-3)^7");

  /*  Composite exponent detection (detected via prime factor).  */
  /*  2^6 = (2^2)^3 = (2^3)^2, detected as square or cube.  */
  test_is_power_i32(a, 64, 1, "64 = 2^6 detected via k=2 or k=3");
  /*  2^12 = (2^2)^6 = (2^3)^4 = (2^4)^3 = (2^6)^2.  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 12), ARBINT_OK);
  test_is_power_value(a, 1, "2^12 = 4096");
  /*  3^15 = (3^3)^5 = (3^5)^3.  */
  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 15), ARBINT_OK);
  test_is_power_value(a, 1, "3^15");

  /*  Near-miss values (off by 1 from perfect powers).  */
  test_is_power_i32(a, 3, 0, "3 = 4-1, not a power");
  test_is_power_i32(a, 5, 0, "5 = 4+1, not a power");
  test_is_power_i32(a, 7, 0, "7 = 8-1, not a power");
  test_is_power_i32(a, 26, 0, "26 = 27-1, not a power");
  test_is_power_i32(a, 28, 0, "28 = 27+1, not a power");
  test_is_power_i32(a, 63, 0, "63 = 64-1, not a power");
  test_is_power_i32(a, 65, 0, "65 = 64+1, not a power");
  test_is_power_i32(a, 124, 0, "124 = 125-1, not a power");
  test_is_power_i32(a, 126, 0, "126 = 125+1, not a power");
  test_is_power_i32(a, 255, 0, "255 = 256-1, not a power");
  test_is_power_i32(a, 257, 0, "257 = 256+1, not a power");

  /*  Near-miss for large values.  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 64), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_i32(a, a, 1), ARBINT_OK);
  test_is_power_value(a, 0, "2^64 - 1");

  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 64), ARBINT_OK);
  CHECK_EQ_I(arbint_add_i32(a, a, 1), ARBINT_OK);
  test_is_power_value(a, 0, "2^64 + 1");

  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_i32(a, a, 1), ARBINT_OK);
  test_is_power_value(a, 0, "3^100 - 1");

  /*  Semiprimes and products that aren't powers.  */
  test_is_power_i32(a, 21, 0, "21 = 3*7");
  test_is_power_i32(a, 22, 0, "22 = 2*11");
  test_is_power_i32(a, 33, 0, "33 = 3*11");
  test_is_power_i32(a, 35, 0, "35 = 5*7");
  test_is_power_i32(a, 39, 0, "39 = 3*13");
  test_is_power_i32(a, 51, 0, "51 = 3*17");
  test_is_power_i32(a, 55, 0, "55 = 5*11");
  test_is_power_i32(a, 57, 0, "57 = 3*19");

  /*  Large base with small exponent.  */
  CHECK_EQ_I(arbint_set_i32(a, 1000000), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
  test_is_power_value(a, 1, "1000000^2 = 10^12");

  CHECK_EQ_I(arbint_set_i32(a, 12345), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
  test_is_power_value(a, 1, "12345^2");

  CHECK_EQ_I(arbint_set_i32(a, 12345), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 3), ARBINT_OK);
  test_is_power_value(a, 1, "12345^3");

  CHECK_EQ_I(arbint_set_i32(a, 999), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 5), ARBINT_OK);
  test_is_power_value(a, 1, "999^5");

  /*  Large base, not a power after modification.  */
  CHECK_EQ_I(arbint_set_i32(a, 12345), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
  CHECK_EQ_I(arbint_add_i32(a, a, 1), ARBINT_OK);
  test_is_power_value(a, 0, "12345^2 + 1");

  /*  Test prime base powers.  */
  CHECK_EQ_I(arbint_set_i32(a, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 11), ARBINT_OK);
  test_is_power_value(a, 1, "7^11");

  CHECK_EQ_I(arbint_set_i32(a, 11), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 13), ARBINT_OK);
  test_is_power_value(a, 1, "11^13");

  CHECK_EQ_I(arbint_set_i32(a, 13), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 17), ARBINT_OK);
  test_is_power_value(a, 1, "13^17");

  /*  Test at static table boundary (521 is last prime in table).  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 521), ARBINT_OK);
  test_is_power_value(a, 1, "2^521 (at static table boundary)");

  /*  Test just past static table boundary.  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 523), ARBINT_OK);
  test_is_power_value(a, 1, "2^523 (past static table)");

  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 541), ARBINT_OK);
  test_is_power_value(a, 1, "2^541 (needs dynamic sieve)");

  /*  Very high exponent test (requires larger sieve).  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 2003), ARBINT_OK);
  test_is_power_value(a, 1, "2^2003 (prime exponent)");

  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 3001), ARBINT_OK);
  test_is_power_value(a, 1, "2^3001 (prime exponent)");

  /*  Negative large odd prime exponent.  */
  CHECK_EQ_I(arbint_set_i32(a, -2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(b, a, 523), ARBINT_OK);
  test_is_power_value(b, 1, "(-2)^523");

  CHECK_EQ_I(arbint_set_i32(a, -3), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(b, a, 541), ARBINT_OK);
  test_is_power_value(b, 1, "(-3)^541");

  /*  Test that negative even exponents are rejected.  */
  /*  -64 is NOT a perfect power: 64 = 2^6 = 4^3 = 8^2, all even bases.
      For -64 = b^k, k must be odd. -64 = (-4)^3? No, (-4)^3 = -64. Yes!  */
  test_is_power_i32(a, -64, 1, "-64 = (-4)^3");

  /*  -256 = (-4)^4? No, that's positive. -256 = b^k for odd k?
      256 = 2^8 = 4^4 = 16^2. No odd power form. -256 = (-b)^k for odd k
      requires 256 = b^k for some b, which needs k|8. k=1 trivial.
      Actually -256 is not a perfect power.  */
  test_is_power_i32(a, -256, 0, "-256 not a perfect power");

  /*  -1024 = 2^10. No odd divisors of 10 except 5? 10 = 2*5.
      Is -1024 = b^5 for some b? 1024^(1/5) = 4. (-4)^5 = -1024. Yes!  */
  test_is_power_i32(a, -1024, 1, "-1024 = (-4)^5");

  /*  -4096 = 2^12. 12 = 4*3. Is -4096 = b^3? 4096^(1/3) = 16. (-16)^3 = -4096.
   * Yes!  */
  CHECK_EQ_I(arbint_set_i32(a, -4096), ARBINT_OK);
  test_is_power_value(a, 1, "-4096 = (-16)^3");

  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);

  ARBINT_TEST_FINISH("test_is_power");
}
