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

/*  Comprehensive tests for specialized squaring kernels (schoolbook and
    Karatsuba) and the power-of-two fast path in arbint_mul_u32/i32.

    Strategy:
      - Cross-check arbint_sqr(a) == arbint_mul(a, a) at many sizes
      - Algebraic identities: (a+b)^2 = a^2 + 2ab + b^2
      - Threshold boundary tests (sizes near KARATSUBA_THRESHOLD)
      - Aliasing (rop == a)
      - Edge cases: 0, 1, -1, single-limb, max-limb values
      - Power-of-two mul_u32 / mul_i32 vs shl cross-check
      - Stochastic tests at varied sizes  */

#include "test_framework.h"

#include <stdint.h>
#include <string.h>

ARBINT_TEST_DECLARE_FAILURES();

/* ===== Deterministic RNG (xorshift64) ===== */

typedef struct {
  uint64_t state;
} xrng_t;

static void xrng_seed(xrng_t * rng, uint64_t seed) {
  if (seed == 0u)
    seed = 0x9e3779b97f4a7c15ull;
  rng->state = seed;
}

static uint64_t xrng_u64(xrng_t * rng) {
  uint64_t x = rng->state;
  x ^= x >> 12;
  x ^= x << 25;
  x ^= x >> 27;
  rng->state = x;
  return x * 0x2545f4914f6cdd1dull;
}

/*  Fill n limbs of an arbint with random data, then normalize.
    sign: +1 or -1.  n must be > 0.  */
static void fill_random(arbint_t x, int sign, size_t n, xrng_t * rng) {
  arbint_limb_t limbs[256];
  size_t i;

  if (n > 256u)
    n = 256u;

  for (i = 0u; i < n; ++i) {
    uint64_t v = xrng_u64(rng);
#if ARBINT_LIMB_BITS == 64
    limbs[i] = (arbint_limb_t) v;
#else
    limbs[i] = (arbint_limb_t) (v & 0xFFFFFFFFu);
#endif
  }

  /*  Ensure top limb is nonzero so set_mag_limbs succeeds.  */
  if (limbs[n - 1u] == 0u)
    limbs[n - 1u] = 1u;

  arbint_test_set_mag_limbs_raw(x, sign, limbs, n);
}

/*  Build a random arbint with exactly n limbs via repeated squaring
    starting from a random seed value.  Returns actual limb count.  */
static size_t build_large(arbint_t x, size_t target_limbs, xrng_t * rng) {
  size_t cur;

  arbint_set_i32(x, (int32_t) (xrng_u64(rng) % 1000000u + 3u));

  while (1) {
    arbint_sqr(x, x);
    cur = arbint_abs_sz(x[0]._sz);
    if (cur >= target_limbs)
      break;
  }

  return cur;
}

/* ===== Test: NULL pointer checks ===== */

static void test_sqr_null(void) {
  arbint_ctx_t ctx;
  arbint_t a;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_sqr(NULL, a), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_sqr(a, NULL), ARBINT_EINVAL);

  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ===== Test: zero, one, minus one ===== */

static void test_sqr_trivial(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  0^2 = 0  */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
  CHECK(arbint_is_zero(r));

  /*  1^2 = 1  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
  check_i32_value(r, 1);

  /*  (-1)^2 = 1  */
  CHECK_EQ_I(arbint_set_i32(a, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
  check_i32_value(r, 1);

  /*  2^2 = 4  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
  check_i32_value(r, 4);

  /*  (-7)^2 = 49  */
  CHECK_EQ_I(arbint_set_i32(a, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
  check_i32_value(r, 49);

  /*  12345^2 = 152399025  */
  CHECK_EQ_I(arbint_set_i32(a, 12345), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
  check_i32_value(r, 152399025);

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ===== Test: sqr == mul(a,a) cross-check for single-limb values ===== */

static void test_sqr_vs_mul_small(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t sqr_r;
  arbint_t mul_r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(sqr_r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(mul_r, &ctx), ARBINT_OK);

  {
    int32_t vals[] = {0,    1,     -1,     2,         -2,           127,
                      -128, 32767, -32768, INT32_MAX, INT32_MIN + 1};
    size_t i;

    for (i = 0u; i < sizeof(vals) / sizeof(vals[0]); ++i) {
      CHECK_EQ_I(arbint_set_i32(a, vals[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_sqr(sqr_r, a), ARBINT_OK);
      CHECK_EQ_I(arbint_mul(mul_r, a, a), ARBINT_OK);
      CHECK(arbint_eq(sqr_r, mul_r));
      /*  Result must be non-negative.  */
      CHECK(arbint_is_zero(sqr_r) || arbint_cmp_i32(sqr_r, 0) > 0);
    }
  }

  arbint_clear(mul_r);
  arbint_clear(sqr_r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ===== Test: aliasing (rop == a) ===== */

static void test_sqr_aliasing(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t expected;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  /*  Small value aliasing.  */
  CHECK_EQ_I(arbint_set_i32(a, -19), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(expected, -19), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(expected, expected, expected), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
  CHECK(arbint_eq(a, expected));

  /*  Multi-limb aliasing: build a value > 1 limb.  */
  {
    arbint_limb_t limbs[4];
    limbs[0] = (arbint_limb_t) ~((arbint_limb_t) 0u);
    limbs[1] = (arbint_limb_t) 42u;
    limbs[2] = (arbint_limb_t) 0x123u;
    limbs[3] = (arbint_limb_t) 1u;
    set_mag_limbs(a, 1, limbs, 4u);
    set_mag_limbs(expected, 1, limbs, 4u);

    CHECK_EQ_I(arbint_mul(expected, expected, expected), ARBINT_OK);
    CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
    CHECK(arbint_eq(a, expected));
  }

  /*  Negative multi-limb aliasing.  */
  {
    arbint_limb_t limbs[4];
    limbs[0] = (arbint_limb_t) ~((arbint_limb_t) 0u);
    limbs[1] = (arbint_limb_t) 42u;
    limbs[2] = (arbint_limb_t) 0x123u;
    limbs[3] = (arbint_limb_t) 1u;
    set_mag_limbs(a, -1, limbs, 4u);
    set_mag_limbs(expected, -1, limbs, 4u);

    CHECK_EQ_I(arbint_mul(expected, expected, expected), ARBINT_OK);
    CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
    CHECK(arbint_eq(a, expected));
  }

  arbint_clear(expected);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ===== Test: max-limb value (all bits set) ===== */

static void test_sqr_max_limb(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t sqr_r;
  arbint_t mul_r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(sqr_r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(mul_r, &ctx), ARBINT_OK);

  /*  Single limb with all bits set.  */
  {
    arbint_limb_t limbs[1];
    limbs[0] = (arbint_limb_t) ~((arbint_limb_t) 0u);
    set_mag_limbs(a, 1, limbs, 1u);

    CHECK_EQ_I(arbint_sqr(sqr_r, a), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(mul_r, a, a), ARBINT_OK);
    CHECK(arbint_eq(sqr_r, mul_r));
  }

  /*  Two limbs with all bits set.  */
  {
    arbint_limb_t limbs[2];
    limbs[0] = (arbint_limb_t) ~((arbint_limb_t) 0u);
    limbs[1] = (arbint_limb_t) ~((arbint_limb_t) 0u);
    set_mag_limbs(a, 1, limbs, 2u);

    CHECK_EQ_I(arbint_sqr(sqr_r, a), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(mul_r, a, a), ARBINT_OK);
    CHECK(arbint_eq(sqr_r, mul_r));
  }

  arbint_clear(mul_r);
  arbint_clear(sqr_r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ===== Test: multi-limb cross-check sqr == mul at specific sizes ===== */

static void test_sqr_vs_mul_multi_limb(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t sqr_r;
  arbint_t mul_r;
  xrng_t rng;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(sqr_r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(mul_r, &ctx), ARBINT_OK);
  xrng_seed(&rng, 0xdeadbeef12345678ull);

  {
    /*  Test at sizes that exercise schoolbook, threshold boundary,
        and Karatsuba.  KARATSUBA_THRESHOLD is 32.  */
    static const size_t sizes[] = {1u,  2u,  3u,  4u,  7u,   8u,  15u,
                                   16u, 30u, 31u, 32u, 33u,  34u, 48u,
                                   63u, 64u, 65u, 96u, 100u, 128u};
    size_t si;

    for (si = 0u; si < sizeof(sizes) / sizeof(sizes[0]) && g_failures == 0;
         ++si) {
      size_t n = sizes[si];
      int sign = (xrng_u64(&rng) & 1u) ? 1 : -1;

      fill_random(a, sign, n, &rng);

      CHECK_EQ_I(arbint_sqr(sqr_r, a), ARBINT_OK);
      CHECK_EQ_I(arbint_mul(mul_r, a, a), ARBINT_OK);

      if (!arbint_eq(sqr_r, mul_r)) {
        fprintf(stderr, "FAIL sqr vs mul mismatch at size %zu (sign %d)\n", n,
                sign);
        ++g_failures;
      }

      /*  Result always positive.  */
      CHECK(arbint_is_zero(sqr_r) || arbint_cmp_i32(sqr_r, 0) > 0);
    }
  }

  arbint_clear(mul_r);
  arbint_clear(sqr_r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ===== Test: algebraic identity (a+b)^2 = a^2 + 2ab + b^2 ===== */

static void test_sqr_algebraic_identity(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t apb;
  arbint_t lhs;
  arbint_t a2;
  arbint_t b2;
  arbint_t ab;
  arbint_t rhs;
  xrng_t rng;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(apb, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(lhs, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a2, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b2, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(ab, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(rhs, &ctx), ARBINT_OK);
  xrng_seed(&rng, 0xcafebabe42ull);

  {
    static const size_t sizes[] = {1u, 4u, 16u, 31u, 32u, 33u, 50u};
    size_t si;

    for (si = 0u; si < sizeof(sizes) / sizeof(sizes[0]) && g_failures == 0;
         ++si) {
      size_t n = sizes[si];

      fill_random(a, 1, n, &rng);
      fill_random(b, 1, n, &rng);

      /*  LHS: (a + b)^2  */
      CHECK_EQ_I(arbint_add(apb, a, b), ARBINT_OK);
      CHECK_EQ_I(arbint_sqr(lhs, apb), ARBINT_OK);

      /*  RHS: a^2 + 2*a*b + b^2  */
      CHECK_EQ_I(arbint_sqr(a2, a), ARBINT_OK);
      CHECK_EQ_I(arbint_sqr(b2, b), ARBINT_OK);
      CHECK_EQ_I(arbint_mul(ab, a, b), ARBINT_OK);
      CHECK_EQ_I(arbint_add(rhs, a2, b2), ARBINT_OK);
      CHECK_EQ_I(arbint_add(rhs, rhs, ab), ARBINT_OK);
      CHECK_EQ_I(arbint_add(rhs, rhs, ab), ARBINT_OK);

      if (!arbint_eq(lhs, rhs)) {
        fprintf(stderr, "FAIL (a+b)^2 identity at size %zu\n", n);
        ++g_failures;
      }
    }
  }

  arbint_clear(rhs);
  arbint_clear(ab);
  arbint_clear(b2);
  arbint_clear(a2);
  arbint_clear(lhs);
  arbint_clear(apb);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ===== Test: difference of squares: (a-b)(a+b) = a^2 - b^2 ===== */

static void test_sqr_diff_of_squares(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t apb;
  arbint_t amb;
  arbint_t lhs;
  arbint_t a2;
  arbint_t b2;
  arbint_t rhs;
  xrng_t rng;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(apb, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(amb, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(lhs, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a2, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b2, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(rhs, &ctx), ARBINT_OK);
  xrng_seed(&rng, 0xfeedface99ull);

  {
    static const size_t sizes[] = {2u, 8u, 31u, 32u, 33u, 64u};
    size_t si;

    for (si = 0u; si < sizeof(sizes) / sizeof(sizes[0]) && g_failures == 0;
         ++si) {
      size_t n = sizes[si];

      fill_random(a, 1, n, &rng);
      fill_random(b, 1, n, &rng);

      /*  LHS: (a+b) * (a-b)  */
      CHECK_EQ_I(arbint_add(apb, a, b), ARBINT_OK);
      CHECK_EQ_I(arbint_sub(amb, a, b), ARBINT_OK);
      CHECK_EQ_I(arbint_mul(lhs, apb, amb), ARBINT_OK);

      /*  RHS: a^2 - b^2  */
      CHECK_EQ_I(arbint_sqr(a2, a), ARBINT_OK);
      CHECK_EQ_I(arbint_sqr(b2, b), ARBINT_OK);
      CHECK_EQ_I(arbint_sub(rhs, a2, b2), ARBINT_OK);

      if (!arbint_eq(lhs, rhs)) {
        fprintf(stderr, "FAIL diff-of-squares identity at size %zu\n", n);
        ++g_failures;
      }
    }
  }

  arbint_clear(rhs);
  arbint_clear(b2);
  arbint_clear(a2);
  arbint_clear(lhs);
  arbint_clear(amb);
  arbint_clear(apb);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ===== Test: large Karatsuba squaring via repeated squaring ===== */

static void test_sqr_large_karatsuba(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t sqr_r;
  arbint_t mul_r;
  xrng_t rng;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(sqr_r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(mul_r, &ctx), ARBINT_OK);
  xrng_seed(&rng, 0xbaadf00d42ull);

  /*  Build a number with ~64 limbs by repeated squaring, then verify
      one more squaring step matches mul(a,a).  */
  build_large(a, 64u, &rng);
  CHECK_EQ_I(arbint_sqr(sqr_r, a), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(mul_r, a, a), ARBINT_OK);
  CHECK(arbint_eq(sqr_r, mul_r));

  /*  Build ~200 limbs.  */
  build_large(a, 200u, &rng);
  CHECK_EQ_I(arbint_sqr(sqr_r, a), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(mul_r, a, a), ARBINT_OK);
  CHECK(arbint_eq(sqr_r, mul_r));

  arbint_clear(mul_r);
  arbint_clear(sqr_r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ===== Test: large aliasing (rop == a) with Karatsuba sizes ===== */

static void test_sqr_large_aliasing(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t expected;
  xrng_t rng;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);
  xrng_seed(&rng, 0x1234abcd5678ull);

  build_large(a, 64u, &rng);
  CHECK_EQ_I(arbint_set(b, a), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(expected, a), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
  CHECK(arbint_eq(a, expected));

  arbint_clear(expected);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ===== Test: stochastic sqr == mul at random sizes ===== */

static void test_sqr_stochastic(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t sqr_r;
  arbint_t mul_r;
  xrng_t rng;
  enum { ITERS = 200 };
  int i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(sqr_r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(mul_r, &ctx), ARBINT_OK);
  xrng_seed(&rng, 0x5555aaaa3333ull);

  for (i = 0; i < ITERS && g_failures == 0; ++i) {
    size_t n = (size_t) (xrng_u64(&rng) % 120u) + 1u;
    int sign = (xrng_u64(&rng) & 1u) ? 1 : -1;

    fill_random(a, sign, n, &rng);

    CHECK_EQ_I(arbint_sqr(sqr_r, a), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(mul_r, a, a), ARBINT_OK);

    if (!arbint_eq(sqr_r, mul_r)) {
      fprintf(stderr, "FAIL stochastic sqr vs mul at iter %d (size %zu)\n", i,
              n);
      ++g_failures;
    }
  }

  arbint_clear(mul_r);
  arbint_clear(sqr_r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ===== Test: power-of-two mul_u32 fast path ===== */

static void test_mul_u32_power_of_two(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t mul_r;
  arbint_t shl_r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(mul_r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(shl_r, &ctx), ARBINT_OK);

  /*  Test various power-of-two multipliers against shl.  */
  {
    uint32_t pows[] = {2u,         4u,         8u,         16u,
                       32u,        64u,        256u,       1024u,
                       (1u << 16), (1u << 20), (1u << 30), (1u << 31)};
    uint32_t shifts[] = {1u, 2u, 3u, 4u, 5u, 6u, 8u, 10u, 16u, 20u, 30u, 31u};
    size_t i;

    for (i = 0u; i < sizeof(pows) / sizeof(pows[0]); ++i) {
      /*  Positive.  */
      CHECK_EQ_I(arbint_set_i32(a, 12345), ARBINT_OK);
      CHECK_EQ_I(arbint_mul_u32(mul_r, a, pows[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_shl(shl_r, a, shifts[i]), ARBINT_OK);
      CHECK(arbint_eq(mul_r, shl_r));

      /*  Negative.  */
      CHECK_EQ_I(arbint_set_i32(a, -9999), ARBINT_OK);
      CHECK_EQ_I(arbint_mul_u32(mul_r, a, pows[i]), ARBINT_OK);
      CHECK_EQ_I(arbint_shl(shl_r, a, shifts[i]), ARBINT_OK);
      CHECK(arbint_eq(mul_r, shl_r));
    }
  }

  /*  Multi-limb power-of-two.  */
  {
    arbint_limb_t limbs[4];
    limbs[0] = (arbint_limb_t) ~((arbint_limb_t) 0u);
    limbs[1] = (arbint_limb_t) ~((arbint_limb_t) 0u);
    limbs[2] = (arbint_limb_t) 42u;
    limbs[3] = (arbint_limb_t) 1u;
    set_mag_limbs(a, 1, limbs, 4u);

    CHECK_EQ_I(arbint_mul_u32(mul_r, a, 16u), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(shl_r, a, 4u), ARBINT_OK);
    CHECK(arbint_eq(mul_r, shl_r));

    CHECK_EQ_I(arbint_mul_u32(mul_r, a, (1u << 31)), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(shl_r, a, 31u), ARBINT_OK);
    CHECK(arbint_eq(mul_r, shl_r));
  }

  /*  Aliasing: rop == a with power-of-two.  */
  {
    CHECK_EQ_I(arbint_set_i32(a, 777), ARBINT_OK);
    CHECK_EQ_I(arbint_shl(shl_r, a, 3u), ARBINT_OK);
    CHECK_EQ_I(arbint_mul_u32(a, a, 8u), ARBINT_OK);
    CHECK(arbint_eq(a, shl_r));
  }

  arbint_clear(shl_r);
  arbint_clear(mul_r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ===== Test: power-of-two mul_i32 fast path ===== */

static void test_mul_i32_power_of_two(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t mul_r;
  arbint_t shl_r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(mul_r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(shl_r, &ctx), ARBINT_OK);

  /*  Positive power-of-two via i32.  */
  CHECK_EQ_I(arbint_set_i32(a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_i32(mul_r, a, 8), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(shl_r, a, 3u), ARBINT_OK);
  CHECK(arbint_eq(mul_r, shl_r));

  /*  Negative power-of-two: mul_i32(a, -8) should equal -shl(a, 3).  */
  CHECK_EQ_I(arbint_set_i32(a, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_i32(mul_r, a, -8), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(shl_r, a, 3u), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(shl_r, shl_r), ARBINT_OK);
  CHECK(arbint_eq(mul_r, shl_r));

  /*  -1 * a (edge case for negative power-of-two path).  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_i32(mul_r, a, -1), ARBINT_OK);
  check_i32_value(mul_r, -42);

  /*  INT32_MIN = -2^31 (negative power-of-two).  */
  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_i32(mul_r, a, INT32_MIN), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(shl_r, a, 31u), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(shl_r, shl_r), ARBINT_OK);
  CHECK(arbint_eq(mul_r, shl_r));

  arbint_clear(shl_r);
  arbint_clear(mul_r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ===== Test: mul_u32 non-power-of-two still works ===== */

static void test_mul_u32_non_power(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t r;
  arbint_t expected;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  /*  Make sure non-power-of-two still goes through the limb multiply path.  */
  CHECK_EQ_I(arbint_set_i32(a, 1234), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_u32(r, a, 567u), ARBINT_OK);
  check_i32_value(r, 699678);

  CHECK_EQ_I(arbint_mul_u32(r, a, 3u), ARBINT_OK);
  check_i32_value(r, 3702);

  CHECK_EQ_I(arbint_mul_u32(r, a, UINT32_MAX), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(expected, UINT32_MAX), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(expected, a, expected), ARBINT_OK);
  CHECK(arbint_eq(r, expected));

  arbint_clear(expected);
  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ===== Test: known exact squares ===== */

static void test_sqr_known_values(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  (B+1)^2 = B^2 + 2B + 1 where B = 2^LIMB_BITS.  */
  {
    arbint_limb_t a_limbs[2];
    arbint_limb_t e_limbs[3];

    a_limbs[0] = 1u;
    a_limbs[1] = 1u;
    set_mag_limbs(a, 1, a_limbs, 2u);

    CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);

    e_limbs[0] = 1u;
    e_limbs[1] = 2u;
    e_limbs[2] = 1u;
    CHECK(limbs_equal(r, 1, e_limbs, 3u));
  }

  /*  (B-1)^2 = B^2 - 2B + 1.  */
  {
    arbint_limb_t max = (arbint_limb_t) ~((arbint_limb_t) 0u);
    arbint_limb_t a_limbs[1];
    arbint_limb_t e_limbs[2];

    a_limbs[0] = max;
    set_mag_limbs(a, 1, a_limbs, 1u);

    CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);

    /*  max^2 = (B-1)^2 = B^2 - 2B + 1
        In limbs: [1, B-2] (since B^2 - 2B + 1 = (B-2)*B + 1).  */
    e_limbs[0] = 1u;
    e_limbs[1] = max - 1u;
    CHECK(limbs_equal(r, 1, e_limbs, 2u));
  }

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

/* ===== Test: negative squaring always produces positive ===== */

static void test_sqr_sign(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t r;
  xrng_t rng;
  int i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  xrng_seed(&rng, 0x7777888899ull);

  for (i = 0; i < 50 && g_failures == 0; ++i) {
    size_t n = (size_t) (xrng_u64(&rng) % 40u) + 1u;
    fill_random(a, -1, n, &rng);

    CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
    CHECK(arbint_cmp_i32(r, 0) > 0);
  }

  arbint_clear(r);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);
}

int main(void) {
  ARBINT_TEST_START();
  test_sqr_null();
  test_sqr_trivial();
  test_sqr_vs_mul_small();
  test_sqr_aliasing();
  test_sqr_max_limb();
  test_sqr_vs_mul_multi_limb();
  test_sqr_algebraic_identity();
  test_sqr_diff_of_squares();
  test_sqr_large_karatsuba();
  test_sqr_large_aliasing();
  test_sqr_stochastic();
  test_mul_u32_power_of_two();
  test_mul_i32_power_of_two();
  test_mul_u32_non_power();
  test_sqr_known_values();
  test_sqr_sign();
  ARBINT_TEST_FINISH("test_sqr");
}
