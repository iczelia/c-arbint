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

#include <stdlib.h>
#include <time.h>

int main(void) {
  ARBINT_TEST_DECLARE_FAILURES();

  arbint_ctx_t ctx;
  arbint_t a, b, expected, result;
  size_t i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(result, &ctx), ARBINT_OK);

  /*  Test 1: Basic small values across threshold boundary.  */
  {
    int32_t test_values[] = {1,      2,       100, 1000, 10000,
                             100000, 1000000, -1,  -42,  -999};

    for (i = 0; i < sizeof(test_values) / sizeof(test_values[0]); ++i) {
      int32_t val = test_values[i];

      CHECK_EQ_I(arbint_set_i32(a, val), ARBINT_OK);

      /*  Double using addition.  */
      CHECK_EQ_I(arbint_add(result, a, a), ARBINT_OK);

      /*  Compute expected via multiplication.  */
      CHECK_EQ_I(arbint_mul_i32(expected, a, 2), ARBINT_OK);

      /*  Results should match.  */
      CHECK_EQ_I(arbint_cmp(result, expected), 0);
    }
  }

  /*  Test 2: Edge cases.  */
  {
    /*  Zero.  */
    CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
    CHECK_EQ_I(arbint_add(result, a, a), ARBINT_OK);
    check_i32_value(result, 0);

    /*  Powers of 2 (trigger specific carry patterns).  */
    for (i = 0; i < 20; ++i) {
      int32_t val = 1 << i;

      CHECK_EQ_I(arbint_set_i32(a, val), ARBINT_OK);
      CHECK_EQ_I(arbint_add(result, a, a), ARBINT_OK);
      CHECK_EQ_I(arbint_mul_i32(expected, a, 2), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(result, expected), 0);
    }

    /*  Near-max i32 values (test overflow to larger sizes).  */
    CHECK_EQ_I(arbint_set_i32(a, 1073741824), ARBINT_OK); /*  2^30  */
    CHECK_EQ_I(arbint_add(result, a, a), ARBINT_OK);
    CHECK_EQ_I(arbint_mul_i32(expected, a, 2), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(result, expected), 0);
  }

  /*  Test 3: In-place doubling (aliasing).  */
  {
    CHECK_EQ_I(arbint_set_i32(a, 12345), ARBINT_OK);
    CHECK_EQ_I(arbint_set_i32(b, 12345), ARBINT_OK);

    /*  Double in-place.  */
    CHECK_EQ_I(arbint_add(a, a, a), ARBINT_OK);

    /*  Double via separate variable.  */
    CHECK_EQ_I(arbint_add(b, b, b), ARBINT_OK);

    /*  Results should match.  */
    CHECK_EQ_I(arbint_cmp(a, b), 0);
    check_i32_value(a, 24690);
  }

  /*  Test 4: Large number doubling (exercises AVX2 path for threshold >= 16
   * limbs).  */
  {
    /*  Build a large number by repeated squaring.  */
    CHECK_EQ_I(arbint_set_i32(a, 123456789), ARBINT_OK);

    for (i = 0; i < 10; ++i) {
      CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
    }

    /*  Now double it.  */
    CHECK_EQ_I(arbint_add(result, a, a), ARBINT_OK);

    /*  Verify via multiplication.  */
    CHECK_EQ_I(arbint_mul_i32(expected, a, 2), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(result, expected), 0);
  }

  /*  Test 5: Random stochastic testing with public API.  */
  {
    unsigned int seed = (unsigned int) time(NULL);
    size_t test_count = 200;

    srand(seed);

    for (i = 0; i < test_count; ++i) {
      /*  Generate random i32 value.  */
      int32_t val = (int32_t) rand();
      if (rand() % 2)
        val = -val;

      CHECK_EQ_I(arbint_set_i32(a, val), ARBINT_OK);

      /*  Build larger number by squaring a few times.  */
      size_t squares = (size_t) (rand() % 5);
      size_t j;
      for (j = 0; j < squares; ++j) {
        CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
      }

      /*  Double using addition.  */
      CHECK_EQ_I(arbint_add(result, a, a), ARBINT_OK);

      /*  Verify via multiplication.  */
      CHECK_EQ_I(arbint_mul_i32(expected, a, 2), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(result, expected), 0);
    }
  }

  /*  Test 6: Negative numbers.  */
  {
    CHECK_EQ_I(arbint_set_i32(a, -42), ARBINT_OK);
    CHECK_EQ_I(arbint_add(result, a, a), ARBINT_OK);
    check_i32_value(result, -84);

    CHECK_EQ_I(arbint_set_i32(a, -1000000), ARBINT_OK);
    CHECK_EQ_I(arbint_add(result, a, a), ARBINT_OK);
    CHECK_EQ_I(arbint_mul_i32(expected, a, 2), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(result, expected), 0);
  }

  arbint_clear(result);
  arbint_clear(expected);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);

  ARBINT_TEST_FINISH("test_addsub_avx2");
}
