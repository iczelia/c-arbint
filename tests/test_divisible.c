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

#include <stdint.h>
#include <string.h>

ARBINT_TEST_DECLARE_FAILURES();

/*  arbint_divisible_u32 tests.  */

static void test_divisible_u32_basic(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  int out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);

  /* Error cases. */
  CHECK_EQ_I(arbint_divisible_u32(NULL, 5u, &out), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_divisible_u32(n, 5u, NULL), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_divisible_u32(n, 0u, &out), ARBINT_EZERO);

  /* Zero is divisible by anything nonzero. */
  CHECK_EQ_I(arbint_set_i32(n, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible_u32(n, 1u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 7u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 1000000u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /* Simple divisible cases. */
  CHECK_EQ_I(arbint_set_i32(n, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible_u32(n, 1u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 2u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 4u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 5u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 10u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 20u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 25u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 50u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 100u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /* Simple non-divisible cases. */
  CHECK_EQ_I(arbint_divisible_u32(n, 3u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);
  CHECK_EQ_I(arbint_divisible_u32(n, 7u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);
  CHECK_EQ_I(arbint_divisible_u32(n, 11u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);
  CHECK_EQ_I(arbint_divisible_u32(n, 101u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  /* Negative numbers: divisibility ignores sign. */
  CHECK_EQ_I(arbint_set_i32(n, -100), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible_u32(n, 2u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 5u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 10u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 3u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  arbint_clear(n);
  arbint_ctx_clear(&ctx);
}

static void test_divisible_u32_powers_of_two(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  int out;
  uint32_t pow2;
  unsigned i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);

  /* Test powers of 2 divisors (uses fast path). */
  CHECK_EQ_I(arbint_set_i32(n, 256), ARBINT_OK);

  for (i = 0u; i <= 8u; ++i) {
    pow2 = 1u << i;
    CHECK_EQ_I(arbint_divisible_u32(n, pow2, &out), ARBINT_OK);
    CHECK_EQ_I(out, 1);
  }

  /* 256 is not divisible by 512. */
  CHECK_EQ_I(arbint_divisible_u32(n, 512u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  /* Test odd number with powers of 2. */
  CHECK_EQ_I(arbint_set_i32(n, 255), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible_u32(n, 1u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 2u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);
  CHECK_EQ_I(arbint_divisible_u32(n, 4u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  /* Test with number that has some trailing zeros. */
  CHECK_EQ_I(arbint_set_i32(n, 96), ARBINT_OK); /* 96 = 32 * 3 = 2^5 * 3 */
  CHECK_EQ_I(arbint_divisible_u32(n, 1u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 2u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 4u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 8u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 16u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 32u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 64u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  arbint_clear(n);
  arbint_ctx_clear(&ctx);
}

static void test_divisible_u32_large_divisors(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  int out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);

  /* Large divisor that fits in u32. */
  CHECK_EQ_I(arbint_set_u32(n, 1000000000u), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible_u32(n, 1000000u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 1000000000u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 7u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  /* Prime divisor. */
  CHECK_EQ_I(arbint_set_u32(n, 2147483647u), ARBINT_OK); /* Mersenne prime */
  CHECK_EQ_I(arbint_divisible_u32(n, 2147483647u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 2u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);
  CHECK_EQ_I(arbint_divisible_u32(n, 3u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  arbint_clear(n);
  arbint_ctx_clear(&ctx);
}

static void test_divisible_u32_multilimb(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t factor;
  int out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(factor, &ctx), ARBINT_OK);

  /* Build a large number via repeated squaring: 2^64 * 12345. */
  CHECK_EQ_I(arbint_set_u32(n, 2u), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK);             /* 4 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK);             /* 16 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK);             /* 256 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK);             /* 65536 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK);             /* 2^32 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK);             /* 2^64 */
  CHECK_EQ_I(arbint_mul_u32(n, n, 12345u), ARBINT_OK); /* 2^64 * 12345 */

  /* n = 2^64 * 12345 is divisible by 3, 5, 15, 823, etc. */
  CHECK_EQ_I(arbint_divisible_u32(n, 1u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 3u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 5u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 15u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 12345u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /* Not divisible by 7. */
  CHECK_EQ_I(arbint_divisible_u32(n, 7u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  /* Power of 2 divisor on multi-limb number. */
  CHECK_EQ_I(arbint_divisible_u32(n, 2u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 4u, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);
  CHECK_EQ_I(arbint_divisible_u32(n, 1u << 20, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  arbint_clear(factor);
  arbint_clear(n);
  arbint_ctx_clear(&ctx);
}

/*  arbint_divisible tests.  */

static void test_divisible_basic(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t d;
  int out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);

  /* Error cases. */
  CHECK_EQ_I(arbint_divisible(NULL, d, &out), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_divisible(n, NULL, &out), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_divisible(n, d, NULL), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_set_i32(d, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_EZERO);

  /* Zero is divisible by any nonzero. */
  CHECK_EQ_I(arbint_set_i32(n, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(d, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /* Basic divisibility. */
  CHECK_EQ_I(arbint_set_i32(n, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(d, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(d, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  /* Negative divisor. */
  CHECK_EQ_I(arbint_set_i32(d, -10), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /* Negative dividend. */
  CHECK_EQ_I(arbint_set_i32(n, -100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(d, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /* Both negative. */
  CHECK_EQ_I(arbint_set_i32(n, -100), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(d, -10), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /* |n| < |d|: not divisible (unless n == 0). */
  CHECK_EQ_I(arbint_set_i32(n, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(d, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  arbint_clear(d);
  arbint_clear(n);
  arbint_ctx_clear(&ctx);
}

static void test_divisible_powers_of_two(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t d;
  int out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);

  /* Test power-of-2 divisor path. */
  CHECK_EQ_I(arbint_set_i32(n, 1024), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(d, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(d, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(d, 512), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(d, 1024), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(d, 2048), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  /* Large power of 2 divisor (multi-limb on 32-bit). */
  CHECK_EQ_I(arbint_set_u32(n, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(n, n, 100), ARBINT_OK); /* n = 2^100 */

  CHECK_EQ_I(arbint_set_u32(d, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(d, d, 50), ARBINT_OK); /* d = 2^50 */
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_shl(d, d, 60), ARBINT_OK); /* d = 2^110 */
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  arbint_clear(d);
  arbint_clear(n);
  arbint_ctx_clear(&ctx);
}

static void test_divisible_multilimb(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t d;
  arbint_t q;
  int out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);

  /* Build large multi-limb numbers. */
  /* n = 2^128 * 999983 (a prime) */
  CHECK_EQ_I(arbint_set_u32(n, 2u), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 4 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 16 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 256 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 65536 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 2^32 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 2^64 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 2^128 */
  CHECK_EQ_I(arbint_mul_u32(n, n, 999983u), ARBINT_OK);

  /* d = 2^64 (multi-limb on 32-bit, single on 64-bit). */
  CHECK_EQ_I(arbint_set_u32(d, 2u), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 2^64 */

  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /* d = 999983 (prime factor of n). */
  CHECK_EQ_I(arbint_set_u32(d, 999983u), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /* d = 999979 (different prime, not a factor). */
  CHECK_EQ_I(arbint_set_u32(d, 999979u), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  /* Test with multi-limb divisor: d = 2^64 * 3. */
  CHECK_EQ_I(arbint_set_u32(d, 2u), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK);         /* 2^64 */
  CHECK_EQ_I(arbint_mul_u32(d, d, 3u), ARBINT_OK); /* 2^64 * 3 */

  /* n = 2^128 * 999983 is not divisible by 2^64 * 3 (999983 is not div by 3).
   */
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  /* n = 2^128 * 3. */
  CHECK_EQ_I(arbint_set_u32(n, 2u), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 2^128 */
  CHECK_EQ_I(arbint_mul_u32(n, n, 3u), ARBINT_OK);

  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  arbint_clear(q);
  arbint_clear(d);
  arbint_clear(n);
  arbint_ctx_clear(&ctx);
}

static void test_divisible_consistency(void) {
  /* Verify divisible(n, d) == (tdiv_r(n, d) == 0). */
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t d;
  arbint_t r;
  int div_out;
  int rem_is_zero;
  uint32_t test_values[] = {1u,     2u,     3u,     5u,       7u,
                            10u,    13u,    17u,    100u,     127u,
                            128u,   255u,   256u,   1000u,    10000u,
                            65535u, 65536u, 99999u, 1000000u, 2147483647u};
  size_t nvals = sizeof(test_values) / sizeof(test_values[0]);
  size_t i;
  size_t j;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  for (i = 0u; i < nvals; ++i) {
    CHECK_EQ_I(arbint_set_u32(n, test_values[i]), ARBINT_OK);
    for (j = 0u; j < nvals; ++j) {
      if (test_values[j] == 0u)
        continue;

      CHECK_EQ_I(arbint_set_u32(d, test_values[j]), ARBINT_OK);

      CHECK_EQ_I(arbint_divisible(n, d, &div_out), ARBINT_OK);
      CHECK_EQ_I(arbint_tdiv_r(r, n, d), ARBINT_OK);
      rem_is_zero = arbint_is_zero(r);

      CHECK_EQ_I(div_out, rem_is_zero);
    }
  }

  arbint_clear(r);
  arbint_clear(d);
  arbint_clear(n);
  arbint_ctx_clear(&ctx);
}

static void test_divisible_u32_consistency(void) {
  /* Verify divisible_u32(n, d) == (tdiv_r_u32(n, d) == 0). */
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t r;
  int div_out;
  uint32_t divisors[] = {
      1u,   2u,   3u,    4u,     5u,     7u,       8u,          9u,
      10u,  16u,  17u,   32u,    64u,    100u,     127u,        128u,
      255u, 256u, 1000u, 65535u, 65536u, 1000000u, 2147483647u, 0xFFFFFFFFu};
  size_t ndivs = sizeof(divisors) / sizeof(divisors[0]);
  size_t i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /* Build a multi-limb number. */
  CHECK_EQ_I(arbint_set_u32(n, 2u), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 2^64 */
  CHECK_EQ_I(arbint_mul_u32(n, n, 123456789u), ARBINT_OK);
  CHECK_EQ_I(arbint_add_u32(n, n, 42u), ARBINT_OK);

  for (i = 0u; i < ndivs; ++i) {
    CHECK_EQ_I(arbint_divisible_u32(n, divisors[i], &div_out), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_r_u32(r, n, divisors[i]), ARBINT_OK);
    CHECK_EQ_I(div_out, arbint_is_zero(r));
  }

  arbint_clear(r);
  arbint_clear(n);
  arbint_ctx_clear(&ctx);
}

static void test_divisible_edge_cases(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t d;
  int out;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);

  /* n == d: always divisible. */
  CHECK_EQ_I(arbint_set_i32(n, 12345), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(d, 12345), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /* d == 1: always divisible. */
  CHECK_EQ_I(arbint_set_i32(d, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /* d == -1: always divisible. */
  CHECK_EQ_I(arbint_set_i32(d, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /* n == 1: only divisible by 1 or -1. */
  CHECK_EQ_I(arbint_set_i32(n, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(d, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(d, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(d, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_divisible(n, d, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  arbint_clear(d);
  arbint_clear(n);
  arbint_ctx_clear(&ctx);
}

int main(void) {
  ARBINT_TEST_START();
  test_divisible_u32_basic();
  test_divisible_u32_powers_of_two();
  test_divisible_u32_large_divisors();
  test_divisible_u32_multilimb();
  test_divisible_basic();
  test_divisible_powers_of_two();
  test_divisible_multilimb();
  test_divisible_consistency();
  test_divisible_u32_consistency();
  test_divisible_edge_cases();

  ARBINT_TEST_FINISH("test_divisible");
}
