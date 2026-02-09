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

/*  Exhaustive multiplication tests covering schoolbook, Karatsuba,
    threshold boundaries, aliasing, algebraic properties, and edge cases.

    Tests are identified by the audit plan IDs A1-L9 in comments.
    Only tests NOT already covered by test_mul.c / test_stochastic.c /
    test_large.c are implemented here.  */

#include "test_framework.h"

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

ARBINT_TEST_DECLARE_FAILURES();

#define MUL_EXHAUST_MAX_LIMBS 2200u

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

static size_t xrng_range(xrng_t * rng, size_t hi_exclusive) {
  if (hi_exclusive == 0u)
    return 0u;
  return (size_t) (xrng_u64(rng) % (uint64_t) hi_exclusive);
}

static arbint_limb_t xrng_limb(xrng_t * rng) {
#if ARBINT_LIMB_BITS == 64
  return (arbint_limb_t) xrng_u64(rng);
#else
  return (arbint_limb_t) (xrng_u64(rng) & 0xffffffffu);
#endif /* ARBINT_LIMB_BITS */
}

/* ===== Helper: build arbint from limb pattern ===== */

/*  Fill x with n limbs of the given pattern, sign +1 or -1.
    Ensures top limb is nonzero (sets bit 0 if needed).  */
static void fill_pattern(arbint_t x, size_t n, int sign, int pattern,
                         xrng_t * rng) {
  size_t i;
  arbint_limb_t maxv = (arbint_limb_t) ~((arbint_limb_t) 0u);

  if (n == 0u) {
    arbint_zero(x);
    return;
  }

  CHECK_EQ_I(arbint_resize(x, n), ARBINT_OK);
  switch (pattern) {
  case 0: /* all ones */
    for (i = 0u; i < n; ++i)
      ARBINT_LIMBS(x)[i] = maxv;
    break;
  case 1: /* single bit in top limb only */
    memset(ARBINT_LIMBS(x), 0, n * sizeof(arbint_limb_t));
    ARBINT_LIMBS(x)[n - 1u] = (arbint_limb_t) 1u;
    break;
  case 2: /* alternating max/0 */
    for (i = 0u; i < n; ++i)
      ARBINT_LIMBS(x)[i] = (i & 1u) ? maxv : (arbint_limb_t) 0u;
    if (ARBINT_LIMBS(x)[n - 1u] == 0u)
      ARBINT_LIMBS(x)[n - 1u] = (arbint_limb_t) 1u;
    break;
  case 3: /* random */
    for (i = 0u; i < n; ++i)
      ARBINT_LIMBS(x)[i] = xrng_limb(rng);
    if (ARBINT_LIMBS(x)[n - 1u] == 0u)
      ARBINT_LIMBS(x)[n - 1u] = (arbint_limb_t) 1u;
    break;
  case 4: /* low half max, high half 1 */
    for (i = 0u; i < n; ++i)
      ARBINT_LIMBS(x)[i] = (i < n / 2u) ? maxv : (arbint_limb_t) 1u;
    break;
  default: /* all 1 */
    for (i = 0u; i < n; ++i)
      ARBINT_LIMBS(x)[i] = (arbint_limb_t) 1u;
    break;
  }

  x[0]._sz = (sign < 0) ? -(ptrdiff_t) n : (ptrdiff_t) n;
}

/* ===== A: mul_wide_limb (tested indirectly via 1x1 mul) ===== */

static void test_wide_limb_via_mul(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r;
  arbint_limb_t maxv = (arbint_limb_t) ~((arbint_limb_t) 0u);
  arbint_limb_t half_bit = (arbint_limb_t) 1u << (ARBINT_LIMB_BITS / 2u);
  arbint_limb_t half_max = half_bit - 1u;
  arbint_limb_t am[1];
  arbint_limb_t bm[1];
  arbint_limb_t em[2];

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  A1: 0 * 0 = 0  */
  CHECK_EQ_I(arbint_set_u32(a, 0u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 0u), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK(arbint_is_zero(r));

  /*  A2: 1 * 1 = 1  */
  CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  check_u32_value(r, 1u);

  /*  A3: LIMB_MAX^2 = [1, LIMB_MAX - 1].  */
  am[0] = maxv;
  bm[0] = maxv;
  set_mag_limbs(a, 1, am, 1u);
  set_mag_limbs(b, 1, bm, 1u);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  em[0] = (arbint_limb_t) 1u;
  em[1] = maxv - (arbint_limb_t) 1u;
  CHECK(limbs_equal(r, 1, em, 2u));

  /*  A4: LIMB_MAX * 1 = LIMB_MAX  */
  am[0] = maxv;
  bm[0] = (arbint_limb_t) 1u;
  set_mag_limbs(a, 1, am, 1u);
  set_mag_limbs(b, 1, bm, 1u);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK(limbs_equal(r, 1, am, 1u));

  /*  A5: half_bit * half_bit = (1, 0) -- product crosses into hi  */
  am[0] = half_bit;
  bm[0] = half_bit;
  set_mag_limbs(a, 1, am, 1u);
  set_mag_limbs(b, 1, bm, 1u);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  em[0] = (arbint_limb_t) 0u;
  em[1] = (arbint_limb_t) 1u;
  CHECK(limbs_equal(r, 1, em, 2u));

  /*  A6: half_max * half_max -- fits in one limb since (2^(N/2)-1)^2 < 2^N.  */
  am[0] = half_max;
  bm[0] = half_max;
  set_mag_limbs(a, 1, am, 1u);
  set_mag_limbs(b, 1, bm, 1u);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  {
    arbint_limb_t expected_lo = half_max * half_max;
    em[0] = expected_lo;
    CHECK(limbs_equal(r, 1, em, 1u));
  }

  arbint_clear(r);
  arbint_clear(b);
  arbint_clear(a);
}

/* ===== B: muladd_limb (tested indirectly) ===== */

static void test_muladd_via_mul(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r, expected;
  arbint_limb_t maxv = (arbint_limb_t) ~((arbint_limb_t) 0u);
  arbint_limb_t am[2];
  arbint_limb_t bm[1];
  arbint_limb_t em[3];

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(expected, &ctx), ARBINT_OK);

  /*  B1: 0*0 + 0 + 0 tested via 0 * 0 = 0, already in A1.  */

  /*  B2: (B^2 - 1) * (B - 1) = B^3 - B^2 - B + 1 = [1, B-1, B-2].  */
  am[0] = maxv;
  am[1] = maxv;
  bm[0] = maxv;
  set_mag_limbs(a, 1, am, 2u);
  set_mag_limbs(b, 1, bm, 1u);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  em[0] = (arbint_limb_t) 1u;
  em[1] = maxv;
  em[2] = maxv - (arbint_limb_t) 1u;
  CHECK(limbs_equal(r, 1, em, 3u));

  /*  B3: (2B - 1) * 1 = [LIMB_MAX, 1].  */
  am[0] = maxv;
  am[1] = (arbint_limb_t) 1u;
  bm[0] = (arbint_limb_t) 1u;
  set_mag_limbs(a, 1, am, 2u);
  set_mag_limbs(b, 1, bm, 1u);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  em[0] = maxv;
  em[1] = (arbint_limb_t) 1u;
  CHECK(limbs_equal(r, 1, em, 2u));

  /*  B4/B5: (B^2 - 1)^2 = B^4 - 2B^2 + 1 = [1, 0, B-2, B-1].  */
  {
    arbint_limb_t am2[2];
    arbint_limb_t bm2[2];
    arbint_limb_t em4[4];

    am2[0] = maxv;
    am2[1] = maxv;
    bm2[0] = maxv;
    bm2[1] = maxv;
    set_mag_limbs(a, 1, am2, 2u);
    set_mag_limbs(b, 1, bm2, 2u);
    CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
    em4[0] = (arbint_limb_t) 1u;
    em4[1] = (arbint_limb_t) 0u;
    em4[2] = maxv - (arbint_limb_t) 1u;
    em4[3] = maxv;
    CHECK(limbs_equal(r, 1, em4, 4u));
  }

  arbint_clear(expected);
  arbint_clear(r);
  arbint_clear(b);
  arbint_clear(a);
}

/* ===== C: schoolbook-specific patterns (below Karatsuba threshold) ===== */

static void test_schoolbook_patterns(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r, t;
  arbint_limb_t maxv = (arbint_limb_t) ~((arbint_limb_t) 0u);
  size_t n;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t, &ctx), ARBINT_OK);

  /*  C1: 0 * N-limb = 0 (both directions).  */
  {
    xrng_t rng;
    xrng_seed(&rng, 0xc1c1c1c1ull);
    fill_pattern(a, 20u, 1, 3, &rng);
    arbint_zero(b);
    CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
    CHECK(arbint_is_zero(r));
    CHECK_EQ_I(arbint_mul(r, b, a), ARBINT_OK);
    CHECK(arbint_is_zero(r));
  }

  /*  C2: 1-limb * 1-limb small values.  */
  CHECK_EQ_I(arbint_set_i32(a, 12345), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 67890), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  /*  12345 * 67890 = 838102050  */
  {
    int32_t got;
    CHECK_EQ_I(arbint_get_i32(r, &got), ARBINT_OK);
    CHECK_EQ_I(got, 838102050);
  }

  /*  C4: 1-limb * N-limb for several N below threshold.  */
  {
    static const size_t sizes[] = {2u, 3u, 4u, 15u, 31u};
    size_t si;
    xrng_t rng;
    xrng_seed(&rng, 0xc4c4c4c4ull);

    for (si = 0u; si < sizeof(sizes) / sizeof(sizes[0]); ++si) {
      n = sizes[si];
      fill_pattern(a, n, 1, 3, &rng);
      fill_pattern(b, 1u, 1, 3, &rng);
      CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
      CHECK_EQ_I(arbint_mul(t, b, a), ARBINT_OK);
      CHECK(arbint_eq(r, t));
    }
  }

  /*  C5: N-limb * N-limb symmetric, below threshold.  */
  {
    static const size_t sizes[] = {2u, 4u, 8u, 16u, 31u};
    size_t si;
    xrng_t rng;
    xrng_seed(&rng, 0xc5c5c5c5ull);

    for (si = 0u; si < sizeof(sizes) / sizeof(sizes[0]); ++si) {
      n = sizes[si];
      fill_pattern(a, n, 1, 3, &rng);
      fill_pattern(b, n, 1, 3, &rng);
      CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
      CHECK_EQ_I(arbint_mul(t, b, a), ARBINT_OK);
      CHECK(arbint_eq(r, t));
    }
  }

  /*  C6: all-ones * all-ones (max carry).
      (B^n - 1)^2 = B^(2n) - 2*B^n + 1 = [1, 0,..,0, B-2, B-1,..,B-1] (2n limbs).  */
  {
    static const size_t sizes[] = {2u, 4u, 8u, 16u, 31u};
    size_t si;

    for (si = 0u; si < sizeof(sizes) / sizeof(sizes[0]); ++si) {
      n = sizes[si];
      fill_pattern(a, n, 1, 0, NULL);
      fill_pattern(b, n, 1, 0, NULL);
      CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
      {
        arbint_limb_t * ep;
        size_t rn = 2u * n;
        size_t j;

        ep = (arbint_limb_t *) calloc(rn, sizeof(arbint_limb_t));
        CHECK(ep != NULL);
        ep[0] = (arbint_limb_t) 1u;
        for (j = 1u; j < n; ++j)
          ep[j] = (arbint_limb_t) 0u;
        ep[n] = maxv - (arbint_limb_t) 1u;
        for (j = n + 1u; j < rn; ++j)
          ep[j] = maxv;
        CHECK(limbs_equal(r, 1, ep, rn));
        free(ep);
      }
    }
  }

  /*  C7: B^(n-1) * B^(n-1) = B^(2n-2).  */
  {
    static const size_t sizes[] = {2u, 5u, 16u, 31u};
    size_t si;

    for (si = 0u; si < sizeof(sizes) / sizeof(sizes[0]); ++si) {
      arbint_limb_t * ep;
      size_t rn;

      n = sizes[si];
      fill_pattern(a, n, 1, 1, NULL);
      fill_pattern(b, n, 1, 1, NULL);
      CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);

      rn = 2u * n - 1u;
      ep = (arbint_limb_t *) calloc(rn, sizeof(arbint_limb_t));
      CHECK(ep != NULL);
      ep[rn - 1u] = (arbint_limb_t) 1u;
      CHECK(limbs_equal(r, 1, ep, rn));
      free(ep);
    }
  }

  /*  C8: alternating LIMB_MAX and 0 limbs.  */
  {
    xrng_t rng;
    xrng_seed(&rng, 0xc8c8c8c8ull);

    fill_pattern(a, 16u, 1, 2, NULL);
    fill_pattern(b, 16u, 1, 2, NULL);
    CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
    /*  Just verify commutativity and nonzero.  */
    CHECK_EQ_I(arbint_mul(t, b, a), ARBINT_OK);
    CHECK(arbint_eq(r, t));
    CHECK(!arbint_is_zero(r));
  }

  /*  C10: LIMB_MAX * LIMB_MAX = [1, LIMB_MAX - 1] (n=1 case of C6).  */
  {
    arbint_limb_t am1[1];
    arbint_limb_t em2[2];
    am1[0] = maxv;
    set_mag_limbs(a, 1, am1, 1u);
    set_mag_limbs(b, 1, am1, 1u);
    CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
    em2[0] = (arbint_limb_t) 1u;
    em2[1] = maxv - (arbint_limb_t) 1u;
    CHECK(limbs_equal(r, 1, em2, 2u));
  }

  /*  C11: (B^n - 1) * 1 = B^n - 1.  */
  {
    static const size_t sizes[] = {1u, 5u, 16u, 31u};
    size_t si;

    for (si = 0u; si < sizeof(sizes) / sizeof(sizes[0]); ++si) {
      n = sizes[si];
      fill_pattern(a, n, 1, 0, NULL);
      CHECK_EQ_I(arbint_set_u32(b, 1u), ARBINT_OK);
      CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
      CHECK(arbint_eq(r, a));
    }
  }

  /*  C12: (B^n - 1) * (B - 1) = B^(n+1) - B^n - B + 1
      = [1, B-1, .., B-1, B-2] (n+1 limbs).  */
  {
    static const size_t sizes[] = {1u, 2u, 5u, 16u, 31u};
    size_t si;

    for (si = 0u; si < sizeof(sizes) / sizeof(sizes[0]); ++si) {
      arbint_limb_t * ep;
      arbint_limb_t bm1[1];
      size_t rn;
      size_t j;

      n = sizes[si];
      fill_pattern(a, n, 1, 0, NULL);
      bm1[0] = maxv;
      set_mag_limbs(b, 1, bm1, 1u);
      CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);

      rn = n + 1u;
      ep = (arbint_limb_t *) calloc(rn, sizeof(arbint_limb_t));
      CHECK(ep != NULL);
      ep[0] = (arbint_limb_t) 1u;
      for (j = 1u; j < n; ++j)
        ep[j] = maxv;
      ep[n] = maxv - (arbint_limb_t) 1u;
      CHECK(limbs_equal(r, 1, ep, rn));
      free(ep);
    }
  }

  arbint_clear(t);
  arbint_clear(r);
  arbint_clear(b);
  arbint_clear(a);
}

/* ===== D: Karatsuba threshold and split tests ===== */

static void test_karatsuba_threshold(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r, t;
  xrng_t rng;
  size_t i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t, &ctx), ARBINT_OK);
  xrng_seed(&rng, 0xd1d2d3d4ull);

  /*  D1-D3: threshold boundary (31, 32, 33 limbs).
      For each, multiply random operands and cross-check commutativity.  */
  {
    static const size_t sizes[] = {31u, 32u, 33u};
    size_t si;

    for (si = 0u; si < 3u; ++si) {
      size_t n = sizes[si];
      for (i = 0u; i < 20u && g_failures == 0; ++i) {
        fill_pattern(a, n, 1, 3, &rng);
        fill_pattern(b, n, 1, 3, &rng);
        CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
        CHECK_EQ_I(arbint_mul(t, b, a), ARBINT_OK);
        CHECK(arbint_eq(r, t));
      }
    }
  }

  /*  D4: cross-validate Karatsuba at 32 limbs via (a+b)^2 = a^2 + 2ab + b^2.  */
  {
    arbint_t a2, b2, ab2, sum, lhs;

    CHECK_EQ_I(arbint_init(a2, &ctx), ARBINT_OK);
    CHECK_EQ_I(arbint_init(b2, &ctx), ARBINT_OK);
    CHECK_EQ_I(arbint_init(ab2, &ctx), ARBINT_OK);
    CHECK_EQ_I(arbint_init(sum, &ctx), ARBINT_OK);
    CHECK_EQ_I(arbint_init(lhs, &ctx), ARBINT_OK);

    for (i = 0u; i < 10u && g_failures == 0; ++i) {
      fill_pattern(a, 32u, 1, 3, &rng);
      fill_pattern(b, 32u, 1, 3, &rng);

      /*  lhs = (a + b)^2  */
      CHECK_EQ_I(arbint_add(lhs, a, b), ARBINT_OK);
      CHECK_EQ_I(arbint_sqr(lhs, lhs), ARBINT_OK);

      /*  rhs = a^2 + 2ab + b^2  */
      CHECK_EQ_I(arbint_sqr(a2, a), ARBINT_OK);
      CHECK_EQ_I(arbint_sqr(b2, b), ARBINT_OK);
      CHECK_EQ_I(arbint_mul(ab2, a, b), ARBINT_OK);
      CHECK_EQ_I(arbint_add(ab2, ab2, ab2), ARBINT_OK);
      CHECK_EQ_I(arbint_add(sum, a2, b2), ARBINT_OK);
      CHECK_EQ_I(arbint_add(sum, sum, ab2), ARBINT_OK);

      CHECK(arbint_eq(lhs, sum));
    }

    arbint_clear(lhs);
    arbint_clear(sum);
    arbint_clear(ab2);
    arbint_clear(b2);
    arbint_clear(a2);
  }

  /*  D5: an == bn == 32, all-ones pattern.  */
  fill_pattern(a, 32u, 1, 0, NULL);
  fill_pattern(b, 32u, 1, 0, NULL);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(t, a), ARBINT_OK);
  CHECK(arbint_eq(r, t));

  /*  D6: balanced odd split: an == bn == 33.  */
  fill_pattern(a, 33u, 1, 0, NULL);
  fill_pattern(b, 33u, 1, 0, NULL);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(t, a), ARBINT_OK);
  CHECK(arbint_eq(r, t));

  /*  D7: an=63, bn=32 (ratio ~2:1, Karatsuba path).  */
  fill_pattern(a, 63u, 1, 3, &rng);
  fill_pattern(b, 32u, 1, 3, &rng);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(t, b, a), ARBINT_OK);
  CHECK(arbint_eq(r, t));

  /*  D8: an=64, bn=32 (ratio exactly 2:1).  */
  fill_pattern(a, 64u, 1, 3, &rng);
  fill_pattern(b, 32u, 1, 3, &rng);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(t, b, a), ARBINT_OK);
  CHECK(arbint_eq(r, t));

  /*  D9: an=65, bn=32 (ratio > 2:1, schoolbook fallback).  */
  fill_pattern(a, 65u, 1, 3, &rng);
  fill_pattern(b, 32u, 1, 3, &rng);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(t, b, a), ARBINT_OK);
  CHECK(arbint_eq(r, t));

  /*  D10: an=100, bn=33 (ratio ~3:1, schoolbook).  */
  fill_pattern(a, 100u, 1, 3, &rng);
  fill_pattern(b, 33u, 1, 3, &rng);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(t, b, a), ARBINT_OK);
  CHECK(arbint_eq(r, t));

  /*  D11: an=64, bn=33 (m1 = 1 in Karatsuba split).  */
  fill_pattern(a, 64u, 1, 3, &rng);
  fill_pattern(b, 33u, 1, 3, &rng);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(t, b, a), ARBINT_OK);
  CHECK(arbint_eq(r, t));

  /*  D12: an=64, bn=31 -- bn below threshold, schoolbook.  */
  fill_pattern(a, 64u, 1, 3, &rng);
  fill_pattern(b, 31u, 1, 3, &rng);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(t, b, a), ARBINT_OK);
  CHECK(arbint_eq(r, t));

  /*  D13-D16: recursive Karatsuba at various depths.  */
  {
    static const size_t sizes[] = {64u, 128u, 256u, 512u, 1024u};
    size_t si;

    for (si = 0u; si < sizeof(sizes) / sizeof(sizes[0]) && g_failures == 0;
         ++si) {
      size_t n = sizes[si];
      fill_pattern(a, n, 1, 3, &rng);
      fill_pattern(b, n, 1, 3, &rng);
      CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
      CHECK_EQ_I(arbint_mul(t, b, a), ARBINT_OK);
      CHECK(arbint_eq(r, t));
    }
  }

  /*  D17-D18: z0/z2 non-overlap via (B^34 - 1)^2 limb pattern.  */
  {
    arbint_limb_t maxl = (arbint_limb_t) ~((arbint_limb_t) 0u);
    arbint_limb_t * ep;
    size_t rn;
    size_t j;

    fill_pattern(a, 34u, 1, 0, NULL);
    CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
    rn = 68u;
    ep = (arbint_limb_t *) calloc(rn, sizeof(arbint_limb_t));
    CHECK(ep != NULL);
    ep[0] = (arbint_limb_t) 1u;
    for (j = 1u; j < 34u; ++j)
      ep[j] = (arbint_limb_t) 0u;
    ep[34u] = maxl - (arbint_limb_t) 1u;
    for (j = 35u; j < rn; ++j)
      ep[j] = maxl;
    CHECK(limbs_equal(r, 1, ep, rn));
    free(ep);
  }

  /*  D19: all-LIMB_MAX * all-LIMB_MAX at n=32 (max carry in add_shifted).  */
  fill_pattern(a, 32u, 1, 0, NULL);
  fill_pattern(b, 32u, 1, 0, NULL);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  {
    arbint_limb_t maxl = (arbint_limb_t) ~((arbint_limb_t) 0u);
    arbint_limb_t * ep;
    size_t rn = 64u;
    size_t j;

    ep = (arbint_limb_t *) calloc(rn, sizeof(arbint_limb_t));
    CHECK(ep != NULL);
    ep[0] = (arbint_limb_t) 1u;
    for (j = 1u; j < 32u; ++j)
      ep[j] = (arbint_limb_t) 0u;
    ep[32u] = maxl - (arbint_limb_t) 1u;
    for (j = 33u; j < rn; ++j)
      ep[j] = maxl;
    CHECK(limbs_equal(r, 1, ep, rn));
    free(ep);
  }

  /*  D20: (B^16 + 1)^2 = B^32 + 2*B^16 + 1 = [1, 0,..,0, 2, 0,..,0, 1].  */
  {
    arbint_limb_t * am17;
    arbint_limb_t * ep;
    size_t j;

    am17 = (arbint_limb_t *) calloc(17u, sizeof(arbint_limb_t));
    CHECK(am17 != NULL);
    am17[0] = (arbint_limb_t) 1u;
    am17[16] = (arbint_limb_t) 1u;
    set_mag_limbs(a, 1, am17, 17u);
    set_mag_limbs(b, 1, am17, 17u);
    CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);

    ep = (arbint_limb_t *) calloc(33u, sizeof(arbint_limb_t));
    CHECK(ep != NULL);
    ep[0] = (arbint_limb_t) 1u;
    for (j = 1u; j < 16u; ++j)
      ep[j] = (arbint_limb_t) 0u;
    ep[16] = (arbint_limb_t) 2u;
    for (j = 17u; j < 32u; ++j)
      ep[j] = (arbint_limb_t) 0u;
    ep[32] = (arbint_limb_t) 1u;
    CHECK(limbs_equal(r, 1, ep, 33u));

    free(ep);
    free(am17);
  }

  /*  D21: commutativity at all threshold-boundary sizes.  */
  {
    static const size_t sizes[] = {31u, 32u, 33u, 63u, 64u, 65u};
    size_t si;

    for (si = 0u; si < sizeof(sizes) / sizeof(sizes[0]) && g_failures == 0;
         ++si) {
      size_t n = sizes[si];
      fill_pattern(a, n, 1, 3, &rng);
      fill_pattern(b, n, 1, 3, &rng);
      CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
      CHECK_EQ_I(arbint_mul(t, b, a), ARBINT_OK);
      CHECK(arbint_eq(r, t));
    }
  }

  arbint_clear(t);
  arbint_clear(r);
  arbint_clear(b);
  arbint_clear(a);
}

/* ===== F: additional arbint_mul public API tests ===== */

static void test_mul_public_extra(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r, t;
  xrng_t rng;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t, &ctx), ARBINT_OK);
  xrng_seed(&rng, 0xf1f2f3f4ull);

  /*  F14: Karatsuba + aliasing (above threshold).  */
  {
    fill_pattern(a, 40u, 1, 3, &rng);
    fill_pattern(b, 40u, 1, 3, &rng);

    /*  rop == a  */
    CHECK_EQ_I(arbint_set(t, a), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t, t, b), ARBINT_OK);
    CHECK(arbint_eq(r, t));

    /*  rop == b  */
    CHECK_EQ_I(arbint_set(t, b), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t, a, t), ARBINT_OK);
    CHECK(arbint_eq(r, t));

    /*  rop == a == b (squaring alias, above threshold)  */
    CHECK_EQ_I(arbint_set(t, a), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t, t, t), ARBINT_OK);
    CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
    CHECK(arbint_eq(r, t));
  }

  /*  F15: 1 * 1 = 1 (unit multiply).  */
  CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(b, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  check_u32_value(r, 1u);

  /*  F16: 1 * x = x (multiplicative identity, multi-limb).  */
  fill_pattern(a, 50u, 1, 3, &rng);
  CHECK_EQ_I(arbint_set_u32(b, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK(arbint_eq(r, a));
  CHECK_EQ_I(arbint_mul(r, b, a), ARBINT_OK);
  CHECK(arbint_eq(r, a));

  /*  F17: (-1) * x = -x  */
  fill_pattern(a, 50u, 1, 3, &rng);
  CHECK_EQ_I(arbint_set_i32(b, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(t, a), ARBINT_OK);
  CHECK(arbint_eq(r, t));

  /*  F18: two 32-limb operands (at Karatsuba threshold).  */
  fill_pattern(a, 32u, 1, 3, &rng);
  fill_pattern(b, 32u, 1, 3, &rng);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK(!arbint_is_zero(r));
  /*  Verify via distributive law: (a+1)*b = a*b + b  */
  {
    arbint_t ap1;
    CHECK_EQ_I(arbint_init(ap1, &ctx), ARBINT_OK);
    CHECK_EQ_I(arbint_add_u32(ap1, a, 1u), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(ap1, ap1, b), ARBINT_OK);
    CHECK_EQ_I(arbint_add(t, r, b), ARBINT_OK);
    CHECK(arbint_eq(ap1, t));
    arbint_clear(ap1);
  }

  /*  F19-F21: large operands (100, 500, 1000 limbs), commutativity.  */
  {
    static const size_t sizes[] = {100u, 500u, 1000u};
    size_t si;

    for (si = 0u; si < sizeof(sizes) / sizeof(sizes[0]) && g_failures == 0;
         ++si) {
      size_t n = sizes[si];
      fill_pattern(a, n, 1, 3, &rng);
      fill_pattern(b, n, 1, 3, &rng);
      CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
      CHECK_EQ_I(arbint_mul(t, b, a), ARBINT_OK);
      CHECK(arbint_eq(r, t));
    }
  }

  arbint_clear(t);
  arbint_clear(r);
  arbint_clear(b);
  arbint_clear(a);
}

/* ===== G: additional sqr tests ===== */

static void test_sqr_extra(void) {
  arbint_ctx_t ctx;
  arbint_t a, r, t;
  xrng_t rng;
  size_t i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t, &ctx), ARBINT_OK);
  xrng_seed(&rng, 0x99887766ull);

  /*  G3: sqr(0) = 0.  */
  arbint_zero(a);
  CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
  CHECK(arbint_is_zero(r));

  /*  G4: sqr(1) = 1.  */
  CHECK_EQ_I(arbint_set_u32(a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
  check_u32_value(r, 1u);

  /*  G5: sqr(-x) = sqr(x) for multi-limb x.  */
  for (i = 0u; i < 10u && g_failures == 0; ++i) {
    size_t n = 20u + xrng_range(&rng, 80u);
    fill_pattern(a, n, 1, 3, &rng);
    CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
    CHECK_EQ_I(arbint_neg(a, a), ARBINT_OK);
    CHECK_EQ_I(arbint_sqr(t, a), ARBINT_OK);
    CHECK(arbint_eq(r, t));
  }

  /*  G6: sqr(a) == mul(a, a) for various sizes.  */
  {
    static const size_t sizes[] = {1u, 2u, 16u, 31u, 32u, 33u, 64u, 100u};
    size_t si;

    for (si = 0u; si < sizeof(sizes) / sizeof(sizes[0]) && g_failures == 0;
         ++si) {
      fill_pattern(a, sizes[si], 1, 3, &rng);
      CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
      CHECK_EQ_I(arbint_mul(t, a, a), ARBINT_OK);
      CHECK(arbint_eq(r, t));
    }
  }

  /*  G8: threshold boundary sizes for sqr.  */
  {
    static const size_t sizes[] = {31u, 32u, 33u};
    size_t si;

    for (si = 0u; si < 3u && g_failures == 0; ++si) {
      fill_pattern(a, sizes[si], 1, 0, NULL);
      CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
      CHECK_EQ_I(arbint_mul(t, a, a), ARBINT_OK);
      CHECK(arbint_eq(r, t));
    }
  }

  /*  G9: sqr(B^k - 1) for various k (all-ones).  */
  {
    static const size_t sizes[] = {1u, 5u, 32u, 33u, 64u};
    size_t si;

    for (si = 0u; si < sizeof(sizes) / sizeof(sizes[0]) && g_failures == 0;
         ++si) {
      fill_pattern(a, sizes[si], 1, 0, NULL);
      CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
      CHECK_EQ_I(arbint_mul(t, a, a), ARBINT_OK);
      CHECK(arbint_eq(r, t));
    }
  }

  /*  G10: repeated squaring to build very large numbers.  */
  CHECK_EQ_I(arbint_set_u32(a, 3u), ARBINT_OK);
  for (i = 0u; i < 15u && g_failures == 0; ++i) {
    CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t, a, a), ARBINT_OK);
    CHECK(arbint_eq(r, t));
    CHECK_EQ_I(arbint_set(a, r), ARBINT_OK);
  }
  /*  3^32768: just verify positive and nonzero.  */
  CHECK(!arbint_is_zero(a));
  CHECK(arbint_signum(a) > 0);

  arbint_clear(t);
  arbint_clear(r);
  arbint_clear(a);
}

/* ===== I10: mul_i32 multiply by -1 ===== */

static void test_mul_i32_neg1(void) {
  arbint_ctx_t ctx;
  arbint_t a, r, t;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t, &ctx), ARBINT_OK);

  /*  I10: multiply by -1 is negation.  */
  CHECK_EQ_I(arbint_set_i32(a, 42), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_i32(r, a, -1), ARBINT_OK);
  check_i32_value(r, -42);

  CHECK_EQ_I(arbint_set_i32(a, -42), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_i32(r, a, -1), ARBINT_OK);
  check_i32_value(r, 42);

  CHECK_EQ_I(arbint_set_u32(a, 0u), ARBINT_OK);
  CHECK_EQ_I(arbint_mul_i32(r, a, -1), ARBINT_OK);
  CHECK(arbint_is_zero(r));

  /*  Multi-limb: -1 * a == neg(a).  */
  {
    xrng_t rng;
    xrng_seed(&rng, 0xaabbccddull);
    fill_pattern(a, 50u, 1, 3, &rng);
    CHECK_EQ_I(arbint_mul_i32(r, a, -1), ARBINT_OK);
    CHECK_EQ_I(arbint_neg(t, a), ARBINT_OK);
    CHECK(arbint_eq(r, t));
  }

  arbint_clear(t);
  arbint_clear(r);
  arbint_clear(a);
}

/* ===== J: algebraic properties ===== */

static void test_algebraic_properties(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, c;
  arbint_t t1, t2, t3, t4, t5;
  xrng_t rng;
  size_t i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t2, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t3, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t4, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t5, &ctx), ARBINT_OK);
  xrng_seed(&rng, 0xa1a2a3a4ull);

  for (i = 0u; i < 60u && g_failures == 0; ++i) {
    size_t na = 16u + xrng_range(&rng, 80u);
    size_t nb = 16u + xrng_range(&rng, 80u);
    size_t nc = 16u + xrng_range(&rng, 80u);
    int sa = (xrng_u64(&rng) & 1u) ? 1 : -1;
    int sb = (xrng_u64(&rng) & 1u) ? 1 : -1;
    int sc = (xrng_u64(&rng) & 1u) ? 1 : -1;

    fill_pattern(a, na, sa, 3, &rng);
    fill_pattern(b, nb, sb, 3, &rng);
    fill_pattern(c, nc, sc, 3, &rng);

    /*  J2: associativity: (a * b) * c == a * (b * c).  */
    CHECK_EQ_I(arbint_mul(t1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t1, t1, c), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t2, b, c), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t2, a, t2), ARBINT_OK);
    CHECK(arbint_eq(t1, t2));

    /*  J4: a * 1 == a.  */
    CHECK_EQ_I(arbint_set_u32(t3, 1u), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t4, a, t3), ARBINT_OK);
    CHECK(arbint_eq(t4, a));

    /*  J5: a * 0 == 0.  */
    arbint_zero(t3);
    CHECK_EQ_I(arbint_mul(t4, a, t3), ARBINT_OK);
    CHECK(arbint_is_zero(t4));

    /*  J6: (-a) * b == -(a * b).  */
    CHECK_EQ_I(arbint_neg(t3, a), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t4, t3, b), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t5, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_neg(t5, t5), ARBINT_OK);
    CHECK(arbint_eq(t4, t5));

    /*  J8: (a+b)(a-b) == a^2 - b^2.  */
    CHECK_EQ_I(arbint_add(t1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(t2, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t3, t1, t2), ARBINT_OK);
    CHECK_EQ_I(arbint_sqr(t4, a), ARBINT_OK);
    CHECK_EQ_I(arbint_sqr(t5, b), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(t4, t4, t5), ARBINT_OK);
    CHECK(arbint_eq(t3, t4));

    /*  J9: (a * b) / b == a (for b != 0, already guaranteed).  */
    CHECK_EQ_I(arbint_mul(t1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_q(t2, t1, b), ARBINT_OK);
    CHECK(arbint_eq(t2, a));

    /*  J10: sign rule: (-a)*(-b) == a*b.  */
    CHECK_EQ_I(arbint_neg(t1, a), ARBINT_OK);
    CHECK_EQ_I(arbint_neg(t2, b), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t3, t1, t2), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t4, a, b), ARBINT_OK);
    CHECK(arbint_eq(t3, t4));
  }

  arbint_clear(t5);
  arbint_clear(t4);
  arbint_clear(t3);
  arbint_clear(t2);
  arbint_clear(t1);
  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
}

/* ===== K: stochastic tests targeting gaps ===== */

static void test_stochastic_threshold_edges(void) {
  /*  K2: threshold edge sizes, commutativity + distributive.  */
  static const size_t edges[] = {31u, 32u, 33u, 63u, 64u, 65u};
  arbint_ctx_t ctx;
  arbint_t a, b, r, t, ap1, ab_plus_b;
  xrng_t rng;
  size_t i;
  size_t ei;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(ap1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(ab_plus_b, &ctx), ARBINT_OK);
  xrng_seed(&rng, 0xed9eed9eull);

  for (ei = 0u; ei < sizeof(edges) / sizeof(edges[0]) && g_failures == 0;
       ++ei) {
    size_t n = edges[ei];
    for (i = 0u; i < 30u && g_failures == 0; ++i) {
      int sa = (xrng_u64(&rng) & 1u) ? 1 : -1;
      int sb = (xrng_u64(&rng) & 1u) ? 1 : -1;
      int pattern_a = (int) (xrng_u64(&rng) % 5u);
      int pattern_b = (int) (xrng_u64(&rng) % 5u);

      fill_pattern(a, n, sa, pattern_a, &rng);
      fill_pattern(b, n, sb, pattern_b, &rng);

      /*  Commutativity.  */
      CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
      CHECK_EQ_I(arbint_mul(t, b, a), ARBINT_OK);
      CHECK(arbint_eq(r, t));

      /*  Distributive: (a+1)*b == a*b + b.  */
      CHECK_EQ_I(arbint_add_u32(ap1, a, 1u), ARBINT_OK);
      CHECK_EQ_I(arbint_mul(t, ap1, b), ARBINT_OK);
      CHECK_EQ_I(arbint_add(ab_plus_b, r, b), ARBINT_OK);
      CHECK(arbint_eq(t, ab_plus_b));
    }
  }

  arbint_clear(ab_plus_b);
  arbint_clear(ap1);
  arbint_clear(t);
  arbint_clear(r);
  arbint_clear(b);
  arbint_clear(a);
}

static void test_stochastic_all_max(void) {
  /*  K4: (B^n - 1)^2 for n = 1..100.  */
  arbint_ctx_t ctx;
  arbint_t a, r;
  arbint_limb_t maxl = (arbint_limb_t) ~((arbint_limb_t) 0u);
  size_t n;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  for (n = 1u; n <= 100u && g_failures == 0; ++n) {
    arbint_limb_t * ep;
    size_t rn = 2u * n;
    size_t j;

    fill_pattern(a, n, 1, 0, NULL);
    CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);

    ep = (arbint_limb_t *) calloc(rn, sizeof(arbint_limb_t));
    CHECK(ep != NULL);
    ep[0] = (arbint_limb_t) 1u;
    for (j = 1u; j < n; ++j)
      ep[j] = (arbint_limb_t) 0u;
    ep[n] = maxl - (arbint_limb_t) 1u;
    for (j = n + 1u; j < rn; ++j)
      ep[j] = maxl;
    CHECK(limbs_equal(r, 1, ep, rn));
    free(ep);
  }

  arbint_clear(r);
  arbint_clear(a);
}

static void test_stochastic_power_of_two(void) {
  /*  K5: (B^(n-1))^2 = B^(2n-2) for n = 1..100.  */
  arbint_ctx_t ctx;
  arbint_t a, r;
  size_t n;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  for (n = 1u; n <= 100u && g_failures == 0; ++n) {
    arbint_limb_t * ep;
    size_t rn = (n == 1u) ? 1u : (2u * n - 1u);

    fill_pattern(a, n, 1, 1, NULL);
    CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);

    ep = (arbint_limb_t *) calloc(rn, sizeof(arbint_limb_t));
    CHECK(ep != NULL);
    ep[rn - 1u] = (arbint_limb_t) 1u;
    CHECK(limbs_equal(r, 1, ep, rn));
    free(ep);
  }

  arbint_clear(r);
  arbint_clear(a);
}

static void test_stochastic_mul_u32_crosscheck(void) {
  /*  K7: mul_u32 vs mul cross-check for random multi-limb a.  */
  arbint_ctx_t ctx;
  arbint_t a, b, r1, r2;
  xrng_t rng;
  size_t i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r2, &ctx), ARBINT_OK);
  xrng_seed(&rng, 0xc7c7c7c7ull);

  for (i = 0u; i < 200u && g_failures == 0; ++i) {
    size_t n = 1u + xrng_range(&rng, 120u);
    uint32_t scalar = (uint32_t) xrng_u64(&rng);
    int sign = (xrng_u64(&rng) & 1u) ? 1 : -1;

    fill_pattern(a, n, sign, 3, &rng);
    CHECK_EQ_I(arbint_set_u32(b, scalar), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(r1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_mul_u32(r2, a, scalar), ARBINT_OK);
    CHECK(arbint_eq(r1, r2));
  }

  arbint_clear(r2);
  arbint_clear(r1);
  arbint_clear(b);
  arbint_clear(a);
}

static void test_stochastic_mul_i32_crosscheck(void) {
  /*  K8: mul_i32 vs mul cross-check for random multi-limb a.  */
  arbint_ctx_t ctx;
  arbint_t a, b, r1, r2;
  xrng_t rng;
  size_t i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r2, &ctx), ARBINT_OK);
  xrng_seed(&rng, 0xc8c8c8c8ull);

  for (i = 0u; i < 200u && g_failures == 0; ++i) {
    size_t n = 1u + xrng_range(&rng, 120u);
    int32_t scalar = (int32_t) xrng_u64(&rng);
    int sign = (xrng_u64(&rng) & 1u) ? 1 : -1;

    if (scalar == 0)
      scalar = 1;

    fill_pattern(a, n, sign, 3, &rng);
    CHECK_EQ_I(arbint_set_i32(b, scalar), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(r1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_mul_i32(r2, a, scalar), ARBINT_OK);
    CHECK(arbint_eq(r1, r2));
  }

  arbint_clear(r2);
  arbint_clear(r1);
  arbint_clear(b);
  arbint_clear(a);
}

/* ===== L: Karatsuba-specific edge cases ===== */

static void test_karatsuba_edges(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r, t;
  xrng_t rng;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t, &ctx), ARBINT_OK);
  xrng_seed(&rng, 0xb1b2b3b4ull);

  /*  L1: bn == 1 with an >= threshold.  */
  fill_pattern(a, 64u, 1, 3, &rng);
  fill_pattern(b, 1u, 1, 3, &rng);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(t, b, a), ARBINT_OK);
  CHECK(arbint_eq(r, t));

  /*  L2: bn == threshold, an == threshold + 1.  */
  fill_pattern(a, 33u, 1, 3, &rng);
  fill_pattern(b, 32u, 1, 3, &rng);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(t, b, a), ARBINT_OK);
  CHECK(arbint_eq(r, t));

  /*  L3: an == 2*bn exactly at threshold.  */
  fill_pattern(a, 64u, 1, 3, &rng);
  fill_pattern(b, 32u, 1, 3, &rng);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(t, b, a), ARBINT_OK);
  CHECK(arbint_eq(r, t));

  /*  L4: an == 2*bn + 1 at threshold. Triggers schoolbook fallback.  */
  fill_pattern(a, 65u, 1, 3, &rng);
  fill_pattern(b, 32u, 1, 3, &rng);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(t, b, a), ARBINT_OK);
  CHECK(arbint_eq(r, t));

  /*  L5: B^(n-1) * B^(m-1) = B^(n+m-2) at Karatsuba sizes.  */
  {
    static const size_t sizes[][2] = {{32u, 32u}, {33u, 33u}, {64u, 32u}};
    size_t si;

    for (si = 0u; si < 3u && g_failures == 0; ++si) {
      size_t na = sizes[si][0];
      size_t nb = sizes[si][1];
      size_t rn = na + nb - 1u;
      arbint_limb_t * ep;

      fill_pattern(a, na, 1, 1, NULL);
      fill_pattern(b, nb, 1, 1, NULL);
      CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);

      ep = (arbint_limb_t *) calloc(rn, sizeof(arbint_limb_t));
      CHECK(ep != NULL);
      ep[rn - 1u] = (arbint_limb_t) 1u;
      CHECK(limbs_equal(r, 1, ep, rn));
      free(ep);
    }
  }

  /*  L6: both operands all-ones at Karatsuba sizes.  */
  {
    static const size_t sizes[] = {32u, 33u, 64u, 65u};
    size_t si;
    arbint_limb_t maxl = (arbint_limb_t) ~((arbint_limb_t) 0u);

    for (si = 0u; si < sizeof(sizes) / sizeof(sizes[0]) && g_failures == 0;
         ++si) {
      size_t n = sizes[si];
      size_t rn = 2u * n;
      size_t j;
      arbint_limb_t * ep;

      fill_pattern(a, n, 1, 0, NULL);
      CHECK_EQ_I(arbint_sqr(r, a), ARBINT_OK);

      ep = (arbint_limb_t *) calloc(rn, sizeof(arbint_limb_t));
      CHECK(ep != NULL);
      ep[0] = (arbint_limb_t) 1u;
      for (j = 1u; j < n; ++j)
        ep[j] = (arbint_limb_t) 0u;
      ep[n] = maxl - (arbint_limb_t) 1u;
      for (j = n + 1u; j < rn; ++j)
        ep[j] = maxl;
      CHECK(limbs_equal(r, 1, ep, rn));
      free(ep);
    }
  }

  /*  L7: (B^n - 1) * B^(m-1) = B^(n+m-1) - B^(m-1).
      limbs [0..m-2] = 0, limbs [m-1..n+m-2] = LIMB_MAX.  */
  {
    size_t n = 33u;
    size_t m = 32u;
    size_t rn = n + m - 1u;
    size_t j;
    arbint_limb_t maxl2 = (arbint_limb_t) ~((arbint_limb_t) 0u);
    arbint_limb_t * ep;

    fill_pattern(a, n, 1, 0, NULL);
    fill_pattern(b, m, 1, 1, NULL);
    CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);

    ep = (arbint_limb_t *) calloc(rn, sizeof(arbint_limb_t));
    CHECK(ep != NULL);
    for (j = 0u; j < m - 1u; ++j)
      ep[j] = (arbint_limb_t) 0u;
    for (j = m - 1u; j < rn; ++j)
      ep[j] = maxl2;
    CHECK(limbs_equal(r, 1, ep, rn));
    free(ep);
  }

  /*  L8: a0 + a1 carry (already covered by all-ones in L6).  */

  /*  L9: middle term carry propagation with n=33 (odd).  */
  fill_pattern(a, 33u, 1, 0, NULL);
  fill_pattern(b, 33u, -1, 0, NULL);
  CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
  CHECK(arbint_signum(r) < 0);
  /*  Verify magnitude matches sqr.  */
  fill_pattern(a, 33u, 1, 0, NULL);
  CHECK_EQ_I(arbint_sqr(t, a), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(t, t), ARBINT_OK);
  CHECK(arbint_eq(r, t));

  arbint_clear(t);
  arbint_clear(r);
  arbint_clear(b);
  arbint_clear(a);
}

/* ===== E: add_shifted (tested indirectly via Karatsuba) ===== */

static void test_add_shifted_indirect(void) {
  arbint_ctx_t ctx;
  arbint_t a, b, r, t, one;
  xrng_t rng;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(one, &ctx), ARBINT_OK);
  xrng_seed(&rng, 0xe1e2e3e4ull);

  CHECK_EQ_I(arbint_set_u32(one, 1u), ARBINT_OK);

  /*  E3/E4/E5: all-ones * (all-ones - 1) at Karatsuba size.
      Verifies via a*(a-1) = a^2 - a.  */
  {
    arbint_limb_t maxl = (arbint_limb_t) ~((arbint_limb_t) 0u);
    arbint_limb_t * bm;
    size_t n = 32u;

    fill_pattern(a, n, 1, 0, NULL);
    bm = (arbint_limb_t *) calloc(n, sizeof(arbint_limb_t));
    CHECK(bm != NULL);
    {
      size_t j;
      for (j = 0u; j < n; ++j)
        bm[j] = maxl;
    }
    bm[0] = maxl - (arbint_limb_t) 1u;
    set_mag_limbs(b, 1, bm, n);

    CHECK_EQ_I(arbint_mul(r, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_sqr(t, a), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(t, t, a), ARBINT_OK);
    CHECK(arbint_eq(r, t));
    free(bm);
  }

  /*  E6: covered by all-ones tests above.  */

  arbint_clear(one);
  arbint_clear(t);
  arbint_clear(r);
  arbint_clear(b);
  arbint_clear(a);
}

/* ===== Main ===== */

int main(void) {
  test_wide_limb_via_mul();
  test_muladd_via_mul();
  test_schoolbook_patterns();
  test_karatsuba_threshold();
  test_mul_public_extra();
  test_sqr_extra();
  test_mul_i32_neg1();
  test_algebraic_properties();
  test_stochastic_threshold_edges();
  test_stochastic_all_max();
  test_stochastic_power_of_two();
  test_stochastic_mul_u32_crosscheck();
  test_stochastic_mul_i32_crosscheck();
  test_karatsuba_edges();
  test_add_shifted_indirect();
  ARBINT_TEST_FINISH("test_mul_exhaustive");
}
