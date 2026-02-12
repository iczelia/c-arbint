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

/*  Helper: verify that r = root_k(a), i.e. r^k <= a < (r+1)^k.  */
static void check_root_invariant(const arbint_t r, const arbint_t a,
                                 uint32_t k) {
  arbint_ctx_t * ctx = arbint_get_ctx(r);
  arbint_t rk, rp1, rp1k;

  CHECK_EQ_I(arbint_init(rk, ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(rp1, ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(rp1k, ctx), ARBINT_OK);

  /*  r^k <= a  */
  CHECK_EQ_I(arbint_pow_u32(rk, r, k), ARBINT_OK);
  CHECK(arbint_cmp(rk, a) <= 0);

  /*  (r+1)^k > a  */
  CHECK_EQ_I(arbint_add_i32(rp1, r, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(rp1k, rp1, k), ARBINT_OK);
  CHECK(arbint_cmp(rp1k, a) > 0);

  arbint_clear(rp1k);
  arbint_clear(rp1);
  arbint_clear(rk);
}

static void test_root_edge_cases(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  k=0: EDOM (undefined).  */
  CHECK_EQ_I(arbint_set_i32(a, 8), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 0u), ARBINT_EDOM);

  /*  k=1: rop = a.  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 1u), ARBINT_OK);
  check_i32_value(r, 42);

  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 1u), ARBINT_OK);
  check_i32_value(r, 0);

  /*  a=0: rop = 0 for any k.  */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 3u), ARBINT_OK);
  check_i32_value(r, 0);

  CHECK_EQ_I(arbint_root(r, a, 10u), ARBINT_OK);
  check_i32_value(r, 0);

  /*  a=1: rop = 1 for any k.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 3u), ARBINT_OK);
  check_i32_value(r, 1);

  CHECK_EQ_I(arbint_root(r, a, 100u), ARBINT_OK);
  check_i32_value(r, 1);

  /*  Negative input: EDOM.  */
  CHECK_EQ_I(arbint_set_i32(a, -8), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 3u), ARBINT_EDOM);

  CHECK_EQ_I(arbint_set_i32(a, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 5u), ARBINT_EDOM);

  /*  NULL inputs: EINVAL.  */
  CHECK_EQ_I(arbint_root(NULL, a, 3u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_root(r, NULL, 3u), ARBINT_EINVAL);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_root_perfect_powers(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Cube roots (k=3).  */
  CHECK_EQ_I(arbint_set_i32(a, 8), ARBINT_OK); /*  2^3 = 8  */
  CHECK_EQ_I(arbint_root(r, a, 3u), ARBINT_OK);
  check_i32_value(r, 2);

  CHECK_EQ_I(arbint_set_i32(a, 27), ARBINT_OK); /*  3^3 = 27  */
  CHECK_EQ_I(arbint_root(r, a, 3u), ARBINT_OK);
  check_i32_value(r, 3);

  CHECK_EQ_I(arbint_set_i32(a, 64), ARBINT_OK); /*  4^3 = 64  */
  CHECK_EQ_I(arbint_root(r, a, 3u), ARBINT_OK);
  check_i32_value(r, 4);

  CHECK_EQ_I(arbint_set_i32(a, 125), ARBINT_OK); /*  5^3 = 125  */
  CHECK_EQ_I(arbint_root(r, a, 3u), ARBINT_OK);
  check_i32_value(r, 5);

  CHECK_EQ_I(arbint_set_i32(a, 1000), ARBINT_OK); /*  10^3 = 1000  */
  CHECK_EQ_I(arbint_root(r, a, 3u), ARBINT_OK);
  check_i32_value(r, 10);

  /*  Fourth roots (k=4).  */
  CHECK_EQ_I(arbint_set_i32(a, 16), ARBINT_OK); /*  2^4 = 16  */
  CHECK_EQ_I(arbint_root(r, a, 4u), ARBINT_OK);
  check_i32_value(r, 2);

  CHECK_EQ_I(arbint_set_i32(a, 81), ARBINT_OK); /*  3^4 = 81  */
  CHECK_EQ_I(arbint_root(r, a, 4u), ARBINT_OK);
  check_i32_value(r, 3);

  CHECK_EQ_I(arbint_set_i32(a, 256), ARBINT_OK); /*  4^4 = 256  */
  CHECK_EQ_I(arbint_root(r, a, 4u), ARBINT_OK);
  check_i32_value(r, 4);

  CHECK_EQ_I(arbint_set_i32(a, 10000), ARBINT_OK); /*  10^4 = 10000  */
  CHECK_EQ_I(arbint_root(r, a, 4u), ARBINT_OK);
  check_i32_value(r, 10);

  /*  Fifth roots (k=5).  */
  CHECK_EQ_I(arbint_set_i32(a, 32), ARBINT_OK); /*  2^5 = 32  */
  CHECK_EQ_I(arbint_root(r, a, 5u), ARBINT_OK);
  check_i32_value(r, 2);

  CHECK_EQ_I(arbint_set_i32(a, 243), ARBINT_OK); /*  3^5 = 243  */
  CHECK_EQ_I(arbint_root(r, a, 5u), ARBINT_OK);
  check_i32_value(r, 3);

  /*  Tenth root (k=10).  */
  CHECK_EQ_I(arbint_set_i32(a, 1024), ARBINT_OK); /*  2^10 = 1024  */
  CHECK_EQ_I(arbint_root(r, a, 10u), ARBINT_OK);
  check_i32_value(r, 2);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_root_non_perfect(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Non-perfect cube roots: floor behavior.  */
  /*  root_3(9) = 2 (since 2^3=8 < 9 < 27=3^3)  */
  CHECK_EQ_I(arbint_set_i32(a, 9), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 3u), ARBINT_OK);
  check_i32_value(r, 2);
  check_root_invariant(r, a, 3u);

  /*  root_3(26) = 2 (since 2^3=8 < 26 < 27=3^3)  */
  CHECK_EQ_I(arbint_set_i32(a, 26), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 3u), ARBINT_OK);
  check_i32_value(r, 2);
  check_root_invariant(r, a, 3u);

  /*  root_3(28) = 3 (since 3^3=27 < 28 < 64=4^3)  */
  CHECK_EQ_I(arbint_set_i32(a, 28), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 3u), ARBINT_OK);
  check_i32_value(r, 3);
  check_root_invariant(r, a, 3u);

  /*  Non-perfect fourth roots.  */
  /*  root_4(80) = 2 (since 2^4=16 < 80 < 81=3^4)  */
  CHECK_EQ_I(arbint_set_i32(a, 80), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 4u), ARBINT_OK);
  check_i32_value(r, 2);
  check_root_invariant(r, a, 4u);

  /*  root_4(82) = 3 (since 3^4=81 < 82 < 256=4^4)  */
  CHECK_EQ_I(arbint_set_i32(a, 82), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 4u), ARBINT_OK);
  check_i32_value(r, 3);
  check_root_invariant(r, a, 4u);

  /*  Test boundary values: just below a perfect power.  */
  /*  root_3(7) = 1 (since 1^3=1 < 7 < 8=2^3)  */
  CHECK_EQ_I(arbint_set_i32(a, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 3u), ARBINT_OK);
  check_i32_value(r, 1);
  check_root_invariant(r, a, 3u);

  /*  root_5(31) = 1 (since 1^5=1 < 31 < 32=2^5)  */
  CHECK_EQ_I(arbint_set_i32(a, 31), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 5u), ARBINT_OK);
  check_i32_value(r, 1);
  check_root_invariant(r, a, 5u);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_root_k2_delegation(void) {
  arbint_ctx_t ctx;
  arbint_t a, r_root, r_isqrt;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r_root, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r_isqrt, &ctx), ARBINT_OK);

  /*  Verify root(a, 2) == isqrt(a) for various values.  */
  int32_t test_values[] = {0,  1,  2,  3,   4,   5,    9,     10,    15,
                           16, 17, 99, 100, 101, 1000, 10000, 65536, 1000000};
  size_t i;

  for (i = 0; i < sizeof(test_values) / sizeof(test_values[0]); ++i) {
    CHECK_EQ_I(arbint_set_i32(a, test_values[i]), ARBINT_OK);
    CHECK_EQ_I(arbint_root(r_root, a, 2u), ARBINT_OK);
    CHECK_EQ_I(arbint_isqrt(r_isqrt, a), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(r_root, r_isqrt), 0);
  }

  arbint_clear(r_isqrt);
  arbint_clear(r_root);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_root_aliasing(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  In-place: rop == a.  */
  CHECK_EQ_I(arbint_set_i32(a, 27), ARBINT_OK);
  CHECK_EQ_I(arbint_root(a, a, 3u), ARBINT_OK);
  check_i32_value(a, 3);

  CHECK_EQ_I(arbint_set_i32(a, 1000), ARBINT_OK);
  CHECK_EQ_I(arbint_root(a, a, 3u), ARBINT_OK);
  check_i32_value(a, 10);

  /*  In-place k=1.  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_root(a, a, 1u), ARBINT_OK);
  check_i32_value(a, 42);

  /*  In-place k=2 (delegates to isqrt).  */
  CHECK_EQ_I(arbint_set_i32(a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_root(a, a, 2u), ARBINT_OK);
  check_i32_value(a, 10);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_root_large(void) {
  arbint_ctx_t ctx;
  arbint_t base, power, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(base, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(power, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Test: root_3(2^192) = 2^64.  */
  CHECK_EQ_I(arbint_set_i32(base, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(power, base, 192u), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, power, 3u), ARBINT_OK);
  /*  Verify r = 2^64 by checking r^3 = 2^192.  */
  CHECK_EQ_I(arbint_pow_u32(base, r, 3u), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(base, power), 0);

  /*  Test: root_4(2^256) = 2^64.  */
  CHECK_EQ_I(arbint_set_i32(base, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(power, base, 256u), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, power, 4u), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(base, r, 4u), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(base, power), 0);

  /*  Test: root_10(2^640) = 2^64.  */
  CHECK_EQ_I(arbint_set_i32(base, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(power, base, 640u), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, power, 10u), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(base, r, 10u), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(base, power), 0);

  /*  Test with a large base: root_3(100^3) = 100.  */
  CHECK_EQ_I(arbint_set_i32(base, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(power, base, 3u), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, power, 3u), ARBINT_OK);
  check_i32_value(r, 100);

  /*  Test: root_5(1000^5) = 1000.  */
  CHECK_EQ_I(arbint_set_i32(base, 1000), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(power, base, 5u), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, power, 5u), ARBINT_OK);
  check_i32_value(r, 1000);

  arbint_clear(r);
  arbint_clear(power);
  arbint_clear(base);
  arbint_ctx_clear(&ctx);
}

static void test_root_pow_roundtrip(void) {
  arbint_ctx_t ctx;
  arbint_t base, power, r, expected;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(base, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(power, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  /*  For various base and k: root_k(base^k) = base.  */
  int32_t bases[] = {2, 3, 5, 7, 10, 13, 17, 100, 127, 255, 1000};
  uint32_t ks[] = {3, 4, 5, 6, 7, 8, 10, 12, 15, 20};
  size_t bi, ki;

  for (bi = 0; bi < sizeof(bases) / sizeof(bases[0]); ++bi) {
    for (ki = 0; ki < sizeof(ks) / sizeof(ks[0]); ++ki) {
      CHECK_EQ_I(arbint_set_i32(base, bases[bi]), ARBINT_OK);
      CHECK_EQ_I(arbint_pow_u32(power, base, ks[ki]), ARBINT_OK);
      CHECK_EQ_I(arbint_root(r, power, ks[ki]), ARBINT_OK);
      CHECK_EQ_I(arbint_cmp(r, base), 0);
    }
  }

  /*  root_k(base^k - 1) = base - 1.  */
  CHECK_EQ_I(arbint_set_i32(base, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(power, base, 3u), ARBINT_OK); /*  1000  */
  CHECK_EQ_I(arbint_sub_i32(power, power, 1), ARBINT_OK); /*  999  */
  CHECK_EQ_I(arbint_root(r, power, 3u), ARBINT_OK);
  check_i32_value(r, 9);
  check_root_invariant(r, power, 3u);

  /*  root_k(base^k + 1) = base (if (base+1)^k > base^k + 1).  */
  CHECK_EQ_I(arbint_set_i32(base, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(power, base, 3u), ARBINT_OK); /*  1000  */
  CHECK_EQ_I(arbint_add_i32(power, power, 1), ARBINT_OK); /*  1001  */
  CHECK_EQ_I(arbint_root(r, power, 3u), ARBINT_OK);
  check_i32_value(r, 10);
  check_root_invariant(r, power, 3u);

  arbint_clear(expected);
  arbint_clear(r);
  arbint_clear(power);
  arbint_clear(base);
  arbint_ctx_clear(&ctx);
}

static void test_root_stochastic(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;
  int i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Test random values and verify invariant.  */
  for (i = 0; i < 200; ++i) {
    /*  Generate a random value by squaring repeatedly.  */
    int32_t seed = 1 + (i * 7) % 100;
    uint32_t k = 3u + (uint32_t) (i % 8); /*  k in [3, 10]  */
    int squarings = 1 + (i % 5);          /*  1-5 squarings  */
    int j;

    CHECK_EQ_I(arbint_set_i32(a, seed), ARBINT_OK);
    for (j = 0; j < squarings; ++j) {
      CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
    }

    CHECK_EQ_I(arbint_root(r, a, k), ARBINT_OK);
    check_root_invariant(r, a, k);
  }

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

static void test_root_large_k(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Large k with small a: root_100(50) = 1.  */
  CHECK_EQ_I(arbint_set_i32(a, 50), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 100u), ARBINT_OK);
  check_i32_value(r, 1);
  check_root_invariant(r, a, 100u);

  /*  Large k with medium a: root_32(2^32) = 2.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 32u), ARBINT_OK); /*  2^32  */
  CHECK_EQ_I(arbint_root(r, a, 32u), ARBINT_OK);
  check_i32_value(r, 2);

  /*  root_64(2^64) = 2.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 64u), ARBINT_OK); /*  2^64  */
  CHECK_EQ_I(arbint_root(r, a, 64u), ARBINT_OK);
  check_i32_value(r, 2);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test consecutive values 2..100 for cube root to verify floor behavior.  */
static void test_root_consecutive(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;
  int32_t i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  for (i = 2; i <= 100; ++i) {
    CHECK_EQ_I(arbint_set_i32(a, i), ARBINT_OK);
    CHECK_EQ_I(arbint_root(r, a, 3u), ARBINT_OK);
    check_root_invariant(r, a, 3u);
  }

  /*  Also test k=4 for some values.  */
  for (i = 2; i <= 100; ++i) {
    CHECK_EQ_I(arbint_set_i32(a, i), ARBINT_OK);
    CHECK_EQ_I(arbint_root(r, a, 4u), ARBINT_OK);
    check_root_invariant(r, a, 4u);
  }

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test boundary values around perfect powers to verify Newton convergence. */
static void test_root_boundary(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Values around 2^30 (for k=3: ~1024^3).  */
  /*  1024^3 = 1073741824.  */
  CHECK_EQ_I(arbint_set_i32(a, 1073741823), ARBINT_OK); /*  1024^3 - 1  */
  CHECK_EQ_I(arbint_root(r, a, 3u), ARBINT_OK);
  check_i32_value(r, 1023);
  check_root_invariant(r, a, 3u);

  CHECK_EQ_I(arbint_set_i32(a, 1073741824), ARBINT_OK); /*  1024^3  */
  CHECK_EQ_I(arbint_root(r, a, 3u), ARBINT_OK);
  check_i32_value(r, 1024);
  check_root_invariant(r, a, 3u);

  CHECK_EQ_I(arbint_set_i32(a, 1073741825), ARBINT_OK); /*  1024^3 + 1  */
  CHECK_EQ_I(arbint_root(r, a, 3u), ARBINT_OK);
  check_i32_value(r, 1024);
  check_root_invariant(r, a, 3u);

  /*  Values around 256^4 = 4294967296 (for k=4).  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 32u), ARBINT_OK); /*  2^32 = 256^4 = 65536^2  */
  CHECK_EQ_I(arbint_root(r, a, 4u), ARBINT_OK);
  check_i32_value(r, 256);
  check_root_invariant(r, a, 4u);

  /*  2^32 - 1  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(a, a, 32u), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_i32(a, a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 4u), ARBINT_OK);
  check_i32_value(r, 255);
  check_root_invariant(r, a, 4u);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test specific k values with known results.  */
static void test_root_specific_k(void) {
  arbint_ctx_t ctx;
  arbint_t a, r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  k=6: 64 = 2^6.  */
  CHECK_EQ_I(arbint_set_i32(a, 64), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 6u), ARBINT_OK);
  check_i32_value(r, 2);

  /*  k=7: 128 = 2^7.  */
  CHECK_EQ_I(arbint_set_i32(a, 128), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 7u), ARBINT_OK);
  check_i32_value(r, 2);

  /*  k=8: 256 = 2^8, 6561 = 3^8.  */
  CHECK_EQ_I(arbint_set_i32(a, 256), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 8u), ARBINT_OK);
  check_i32_value(r, 2);

  CHECK_EQ_I(arbint_set_i32(a, 6561), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 8u), ARBINT_OK);
  check_i32_value(r, 3);

  /*  k=9: 512 = 2^9.  */
  CHECK_EQ_I(arbint_set_i32(a, 512), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 9u), ARBINT_OK);
  check_i32_value(r, 2);

  /*  k=12: 4096 = 2^12.  */
  CHECK_EQ_I(arbint_set_i32(a, 4096), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 12u), ARBINT_OK);
  check_i32_value(r, 2);

  /*  k=16: 65536 = 2^16.  */
  CHECK_EQ_I(arbint_set_i32(a, 65536), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 16u), ARBINT_OK);
  check_i32_value(r, 2);

  /*  k=20: 1048576 = 2^20.  */
  CHECK_EQ_I(arbint_set_i32(a, 1048576), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, a, 20u), ARBINT_OK);
  check_i32_value(r, 2);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/*  Test very large operands with multi-limb results.  */
static void test_root_multilimb_result(void) {
  arbint_ctx_t ctx;
  arbint_t base, power, r, check;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(base, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(power, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(check, &ctx), ARBINT_OK);

  /*  Build a large base (2^100) and compute root_3((2^100)^3) = 2^100.  */
  CHECK_EQ_I(arbint_set_i32(base, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(base, base, 100u), ARBINT_OK);    /*  base = 2^100  */
  CHECK_EQ_I(arbint_pow_u32(power, base, 3u), ARBINT_OK); /*  power = 2^300  */
  CHECK_EQ_I(arbint_root(r, power, 3u), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, base), 0);

  /*  root_5((2^50)^5) = 2^50.  */
  CHECK_EQ_I(arbint_set_i32(base, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(base, base, 50u), ARBINT_OK);     /*  base = 2^50  */
  CHECK_EQ_I(arbint_pow_u32(power, base, 5u), ARBINT_OK); /*  power = 2^250  */
  CHECK_EQ_I(arbint_root(r, power, 5u), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(r, base), 0);

  /*  Test with a non-power-of-two base: root_3((12345)^3).  */
  CHECK_EQ_I(arbint_set_i32(base, 12345), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(power, base, 3u), ARBINT_OK);
  CHECK_EQ_I(arbint_root(r, power, 3u), ARBINT_OK);
  check_i32_value(r, 12345);

  /*  Non-perfect power near a large perfect power.  */
  CHECK_EQ_I(arbint_add_i32(power, power, 1), ARBINT_OK); /*  12345^3 + 1  */
  CHECK_EQ_I(arbint_root(r, power, 3u), ARBINT_OK);
  check_i32_value(r, 12345);
  check_root_invariant(r, power, 3u);

  arbint_clear(check);
  arbint_clear(r);
  arbint_clear(power);
  arbint_clear(base);
  arbint_ctx_clear(&ctx);
}

int main(void) {
  test_root_edge_cases();
  test_root_perfect_powers();
  test_root_non_perfect();
  test_root_k2_delegation();
  test_root_aliasing();
  test_root_large();
  test_root_pow_roundtrip();
  test_root_stochastic();
  test_root_large_k();
  test_root_consecutive();
  test_root_boundary();
  test_root_specific_k();
  test_root_multilimb_result();

  ARBINT_TEST_FINISH("test_root");
}
