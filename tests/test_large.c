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

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

ARBINT_TEST_DECLARE_FAILURES();

#define LARGE_MAX_LIMBS 1600u

/* Keep int32 reductions safely inside range for exact host-side checks. */
#define SMALL_I32_BOUND 100000000

typedef struct {
  uint64_t state;
} large_rng_t;

static void large_seed(large_rng_t * rng, uint64_t seed) {
  if (seed == 0u)
    seed = 0x9e3779b97f4a7c15ull;
  rng->state = seed;
}

static uint64_t large_u64(large_rng_t * rng) {
  uint64_t x = rng->state;
  x ^= x >> 12;
  x ^= x << 25;
  x ^= x >> 27;
  rng->state = x;
  return x * 0x2545f4914f6cdd1dull;
}

static size_t large_range(large_rng_t * rng, size_t hi_exclusive) {
  if (hi_exclusive == 0u)
    return 0u;
  return (size_t) (large_u64(rng) % (uint64_t) hi_exclusive);
}

static arbint_limb_t large_limb(large_rng_t * rng) {
#if ARBINT_LIMB_BITS == 64
  return (arbint_limb_t) large_u64(rng);
#else
  return (arbint_limb_t) (large_u64(rng) & 0xffffffffu);
#endif
}

static size_t large_min_size(size_t a, size_t b) { return (a < b) ? a : b; }

static int32_t large_i32(large_rng_t * rng) {
  uint32_t span = (uint32_t) (2 * SMALL_I32_BOUND + 1);
  uint32_t v = (uint32_t) (large_u64(rng) % span);
  return (int32_t) v - SMALL_I32_BOUND;
}

static int32_t large_i32_nonzero_not_min(large_rng_t * rng) {
  int32_t v;
  do {
    v = (int32_t) large_u64(rng);
  } while (v == 0 || v == INT32_MIN);
  return v;
}

static void large_dump_brief(const char * label, const arbint_t x) {
  int32_t i32v;
  size_t used = arbint_abs_sz(x[0]._sz);

  fprintf(stderr, "%s sign=%d used_limbs=%zu", label, arbint_signum(x), used);
  if (arbint_fits_i32(x) && arbint_get_i32(x, &i32v) == ARBINT_OK)
    fprintf(stderr, " i32=%d", i32v);
  fputc('\n', stderr);
}

static void large_check_eq(const char * op, size_t iter, const arbint_t got,
                           const arbint_t expect, const arbint_t a,
                           const arbint_t b, const arbint_t c) {
  if (arbint_eq(got, expect))
    return;

  fprintf(stderr, "FAIL large[%s] iter=%zu\n", op, iter);
  large_dump_brief("got", got);
  large_dump_brief("expect", expect);
  large_dump_brief("a", a);
  large_dump_brief("b", b);
  large_dump_brief("c", c);
  ++g_failures;
}

static void large_check_i32(const char * op, size_t iter, const arbint_t x,
                            int32_t expect, const arbint_t a,
                            const arbint_t b) {
  int32_t out = 0;
  arbint_err_t rc = arbint_get_i32(x, &out);

  if (rc == ARBINT_OK && out == expect)
    return;

  fprintf(stderr, "FAIL large[%s] iter=%zu i32_rc=%d got=%d expect=%d\n", op,
          iter, (int) rc, out, expect);
  large_dump_brief("x", x);
  large_dump_brief("a", a);
  large_dump_brief("b", b);
  ++g_failures;
}

static void large_random_big(arbint_t x, large_rng_t * rng, size_t max_limbs,
                             int allow_zero) {
  size_t n;
  size_t i;
  uint64_t bucket;
  int sign;

  if (max_limbs == 0u)
    max_limbs = 1u;
  if (max_limbs > LARGE_MAX_LIMBS)
    max_limbs = LARGE_MAX_LIMBS;

  bucket = large_u64(rng) % 100u;
  if (allow_zero && bucket < 3u) {
    arbint_zero(x);
    return;
  }

  if (bucket < 20u) {
    n = 1u + large_range(rng, large_min_size(max_limbs, 64u));
  } else if (bucket < 50u) {
    n = 1u + large_range(rng, large_min_size(max_limbs, 256u));
  } else if (bucket < 80u) {
    n = 1u + large_range(rng, large_min_size(max_limbs, 1024u));
  } else {
    size_t near = large_min_size(max_limbs, 64u);
    n = max_limbs - large_range(rng, near);
    if (n == 0u)
      n = 1u;
  }

  CHECK_EQ_I(arbint_resize(x, n), ARBINT_OK);

  switch (large_u64(rng) & 7u) {
  case 0u:
    for (i = 0u; i < n; ++i)
      ARBINT_LIMBS(x)[i] = large_limb(rng);
    break;
  case 1u: {
    arbint_limb_t maxv = (arbint_limb_t) ~((arbint_limb_t) 0u);
    for (i = 0u; i < n; ++i)
      ARBINT_LIMBS(x)[i] = maxv;
    break;
  }
  case 2u: {
    arbint_limb_t maxv = (arbint_limb_t) ~((arbint_limb_t) 0u);
    for (i = 0u; i < n; ++i)
      ARBINT_LIMBS(x)[i] = (i & 1u) ? maxv : (arbint_limb_t) 0u;
    break;
  }
  case 3u: {
    arbint_limb_t maxv = (arbint_limb_t) ~((arbint_limb_t) 0u);
    for (i = 0u; i + 1u < n; ++i)
      ARBINT_LIMBS(x)[i] = maxv;
    ARBINT_LIMBS(x)[n - 1u] = large_limb(rng);
    break;
  }
  case 4u:
    memset(ARBINT_LIMBS(x), 0, n * sizeof(arbint_limb_t));
    for (i = 0u; i < large_min_size(n, 10u); ++i) {
      size_t bit = large_range(rng, n * ARBINT_LIMB_BITS);
      size_t li = bit / ARBINT_LIMB_BITS;
      size_t bi = bit % ARBINT_LIMB_BITS;
      ARBINT_LIMBS(x)[li] |= ((arbint_limb_t) 1u) << bi;
    }
    break;
  case 5u:
    for (i = 0u; i < n; ++i) {
      uint64_t roll = large_u64(rng) & 15u;
      if (roll < 3u)
        ARBINT_LIMBS(x)[i] = (arbint_limb_t) 0u;
      else if (roll < 8u)
        ARBINT_LIMBS(x)[i] = (arbint_limb_t) ~((arbint_limb_t) 0u);
      else
        ARBINT_LIMBS(x)[i] = large_limb(rng);
    }
    break;
  case 6u:
    for (i = 0u; i < n; ++i)
      ARBINT_LIMBS(x)[i] = (arbint_limb_t) 0u;
    ARBINT_LIMBS(x)
    [n - 1u] = ((arbint_limb_t) 1u) << large_range(rng, ARBINT_LIMB_BITS);
    break;
  default:
    for (i = 0u; i < n; ++i)
      ARBINT_LIMBS(x)[i] = large_limb(rng);
    break;
  }

  if (ARBINT_LIMBS(x)[n - 1u] == (arbint_limb_t) 0u)
    ARBINT_LIMBS(x)
  [n - 1u] = ((arbint_limb_t) 1u) << large_range(rng, ARBINT_LIMB_BITS);

  sign = (large_u64(rng) & 1u) ? 1 : -1;
  x[0]._sz = (sign < 0) ? -(ptrdiff_t) n : (ptrdiff_t) n;
}

static size_t large_mul_size(size_t iter, large_rng_t * rng) {
  static const size_t edges[] = {31u,  32u,  33u,  63u,  64u,  65u,  95u,
                                 96u,  97u,  127u, 128u, 129u, 191u, 192u,
                                 193u, 255u, 256u, 257u, 383u, 384u, 385u};
  if ((iter % 3u) == 0u)
    return edges[iter % (sizeof(edges) / sizeof(edges[0]))];
  return 48u + large_range(rng, 360u);
}

static void test_large_add_sub_properties(void) {
  enum { ITERS = 900 };
  size_t i;
  large_rng_t rng;
  arbint_ctx_t ctx;
  arbint_t a, b, c;
  arbint_t t1, t2, t3, t4;

  large_seed(&rng, 0xa20f4e1c53d9b671ull);

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t2, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t3, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t4, &ctx), ARBINT_OK);

  for (i = 0u; i < ITERS && g_failures == 0; ++i) {
    int32_t s = large_i32(&rng);

    large_random_big(a, &rng, 1200u + large_range(&rng, 400u), 1);
    large_random_big(b, &rng, 1200u + large_range(&rng, 400u), 1);
    large_random_big(c, &rng, 1200u + large_range(&rng, 400u), 1);

    /* (a + b) - b == a */
    CHECK_EQ_I(arbint_add(t1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(t2, t1, b), ARBINT_OK);
    large_check_eq("add_sub_cancel", i, t2, a, a, b, c);

    /* (a - b) + b == a */
    CHECK_EQ_I(arbint_sub(t1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_add(t2, t1, b), ARBINT_OK);
    large_check_eq("sub_add_cancel", i, t2, a, a, b, c);

    /* Commutativity of addition. */
    CHECK_EQ_I(arbint_add(t1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_add(t2, b, a), ARBINT_OK);
    large_check_eq("add_comm", i, t1, t2, a, b, c);

    /* -(b - a) == (a - b). */
    CHECK_EQ_I(arbint_sub(t1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(t2, b, a), ARBINT_OK);
    CHECK_EQ_I(arbint_neg(t2, t2), ARBINT_OK);
    large_check_eq("sub_antisym", i, t1, t2, a, b, c);

    /* Associativity: (a + b) + c == a + (b + c). */
    CHECK_EQ_I(arbint_add(t1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_add(t1, t1, c), ARBINT_OK);
    CHECK_EQ_I(arbint_add(t2, b, c), ARBINT_OK);
    CHECK_EQ_I(arbint_add(t2, a, t2), ARBINT_OK);
    large_check_eq("add_assoc", i, t1, t2, a, b, c);

    /* Alias stress for add/sub round-trip. */
    CHECK_EQ_I(arbint_set(t1, a), ARBINT_OK);
    CHECK_EQ_I(arbint_add(t1, t1, b), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(t1, t1, b), ARBINT_OK);
    large_check_eq("add_alias_roundtrip", i, t1, a, a, b, c);

    CHECK_EQ_I(arbint_set(t1, a), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(t1, t1, b), ARBINT_OK);
    CHECK_EQ_I(arbint_add(t1, t1, b), ARBINT_OK);
    large_check_eq("sub_alias_roundtrip", i, t1, a, a, b, c);

    /* Reduce huge expression to exact int32 oracle: ((a + s) - a) == s */
    CHECK_EQ_I(arbint_set(t1, a), ARBINT_OK);
    CHECK_EQ_I(arbint_add_i32(t1, t1, s), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(t1, t1, a), ARBINT_OK);
    large_check_i32("add_i32_reduce", i, t1, s, a, b);

    CHECK_EQ_I(arbint_set(t2, b), ARBINT_OK);
    CHECK_EQ_I(arbint_sub_i32(t2, t2, s), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(t2, b, t2), ARBINT_OK);
    large_check_i32("sub_i32_reduce", i, t2, s, a, b);
  }

  arbint_clear(a);
  arbint_clear(b);
  arbint_clear(c);
  arbint_clear(t1);
  arbint_clear(t2);
  arbint_clear(t3);
  arbint_clear(t4);
}

static void test_large_mul_properties(void) {
  enum { ITERS = 280 };
  size_t i;
  large_rng_t rng;
  arbint_ctx_t ctx;
  arbint_t a, b, c;
  arbint_t t1, t2, t3, t4, t5, t6;

  large_seed(&rng, 0x6d7b3c51f9a402e4ull);

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t2, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t3, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t4, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t5, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t6, &ctx), ARBINT_OK);

  for (i = 0u; i < ITERS && g_failures == 0; ++i) {
    int32_t scale = large_i32_nonzero_not_min(&rng);
    int32_t keep = large_i32(&rng);

    large_random_big(a, &rng, large_mul_size(i, &rng), 1);
    large_random_big(b, &rng, large_mul_size(i + 17u, &rng), 1);
    large_random_big(c, &rng, large_mul_size(i + 41u, &rng), 1);

    /* Multiplication commutativity. */
    CHECK_EQ_I(arbint_mul(t1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t2, b, a), ARBINT_OK);
    large_check_eq("mul_comm", i, t1, t2, a, b, c);

    /* Distributive law: (a+b)*c == a*c + b*c */
    CHECK_EQ_I(arbint_add(t1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t1, t1, c), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t2, a, c), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t3, b, c), ARBINT_OK);
    CHECK_EQ_I(arbint_add(t2, t2, t3), ARBINT_OK);
    large_check_eq("mul_distrib", i, t1, t2, a, b, c);

    /* (a+b)^2 == a^2 + 2ab + b^2 */
    CHECK_EQ_I(arbint_add(t1, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_sqr(t1, t1), ARBINT_OK);
    CHECK_EQ_I(arbint_sqr(t2, a), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t3, a, b), ARBINT_OK);
    CHECK_EQ_I(arbint_add(t3, t3, t3), ARBINT_OK);
    CHECK_EQ_I(arbint_sqr(t4, b), ARBINT_OK);
    CHECK_EQ_I(arbint_add(t2, t2, t3), ARBINT_OK);
    CHECK_EQ_I(arbint_add(t2, t2, t4), ARBINT_OK);
    large_check_eq("sqr_expand", i, t1, t2, a, b, c);

    /* Alias: (a*b)/scale recovers a exactly when built as a*scale. */
    CHECK_EQ_I(arbint_set_i32(t6, scale), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t1, a, t6), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_qr_i32(t2, t3, t1, scale), ARBINT_OK);
    large_check_eq("mul_i32_tdiv_i32_q", i, t2, a, a, b, c);
    CHECK(arbint_is_zero(t3));

    /* Reduce huge expression to exact int32 oracle: ((a*scale)+keep)-(a*scale)
     */
    CHECK_EQ_I(arbint_set_i32(t6, scale), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t4, a, t6), ARBINT_OK);
    CHECK_EQ_I(arbint_set(t5, t4), ARBINT_OK);
    CHECK_EQ_I(arbint_add_i32(t5, t5, keep), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(t6, t5, t4), ARBINT_OK);
    large_check_i32("mul_reduce_i32", i, t6, keep, a, b);
  }

  arbint_clear(a);
  arbint_clear(b);
  arbint_clear(c);
  arbint_clear(t1);
  arbint_clear(t2);
  arbint_clear(t3);
  arbint_clear(t4);
  arbint_clear(t5);
  arbint_clear(t6);
}

static void test_large_tdiv_properties(void) {
  enum { ITERS = 220 };
  size_t i;
  large_rng_t rng;
  arbint_ctx_t ctx;
  arbint_t n, d, q, r;
  arbint_t t1, t2, ncopy, s;

  large_seed(&rng, 0xf24a6dc3850b91e7ull);

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t1, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(t2, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(ncopy, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(s, &ctx), ARBINT_OK);

  for (i = 0u; i < ITERS && g_failures == 0; ++i) {
    int32_t d32 = large_i32_nonzero_not_min(&rng);
    int64_t absd64;
    int32_t rem_abs;
    int32_t rem;

    large_random_big(n, &rng, 900u + large_range(&rng, 700u), 1);
    do {
      large_random_big(d, &rng, 320u + large_range(&rng, 500u), 1);
    } while (arbint_is_zero(d));

    CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);

    /* Recomposition: q*d + r == n. */
    CHECK_EQ_I(arbint_mul(t1, q, d), ARBINT_OK);
    CHECK_EQ_I(arbint_add(t1, t1, r), ARBINT_OK);
    large_check_eq("tdiv_recompose", i, t1, n, n, d, q);

    /* Remainder contract. */
    if (!arbint_is_zero(r)) {
      CHECK_EQ_I(arbint_signum(r), arbint_signum(n));
      CHECK(arbint_cmpabs(r, d) < 0);
    }

    /* q-only and r-only APIs agree with qr API. */
    CHECK_EQ_I(arbint_tdiv_q(t1, n, d), ARBINT_OK);
    large_check_eq("tdiv_q_consistency", i, t1, q, n, d, r);

    CHECK_EQ_I(arbint_tdiv_r(t1, n, d), ARBINT_OK);
    large_check_eq("tdiv_r_consistency", i, t1, r, n, d, q);

    /* Alias paths for qr. */
    CHECK_EQ_I(arbint_set(ncopy, n), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_qr(ncopy, t1, ncopy, d), ARBINT_OK);
    large_check_eq("tdiv_alias_q", i, ncopy, q, n, d, r);
    large_check_eq("tdiv_alias_q_r", i, t1, r, n, d, q);

    CHECK_EQ_I(arbint_set(ncopy, n), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_qr(t1, ncopy, ncopy, d), ARBINT_OK);
    large_check_eq("tdiv_alias_r_q", i, t1, q, n, d, r);
    large_check_eq("tdiv_alias_r", i, ncopy, r, n, d, q);

    /* Construct exact i32-division case: n = q*d32 + rem, |rem| < |d32|. */
    CHECK_EQ_I(arbint_set_i32(s, d32), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(t1, d, s), ARBINT_OK);

    absd64 = (d32 < 0) ? -(int64_t) d32 : (int64_t) d32;
    rem_abs = (int32_t) large_range(&rng, (size_t) absd64);

    if (arbint_signum(t1) < 0)
      rem = -rem_abs;
    else if (arbint_signum(t1) > 0)
      rem = rem_abs;
    else
      rem = 0;

    CHECK_EQ_I(arbint_add_i32(t1, t1, rem), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_qr_i32(t2, r, t1, d32), ARBINT_OK);
    large_check_eq("tdiv_i32_exact_q", i, t2, d, n, t1, q);
    large_check_i32("tdiv_i32_exact_r", i, r, rem, t1, d);

    /* Same idea for positive u32 divisor. */
    {
      uint32_t du = 1u + (uint32_t) large_range(&rng, (size_t) INT32_MAX);
      int32_t urem = (int32_t) large_range(&rng, (size_t) du);
      CHECK_EQ_I(arbint_set_u32(s, du), ARBINT_OK);
      CHECK_EQ_I(arbint_mul(t1, d, s), ARBINT_OK);
      if (arbint_signum(t1) < 0)
        urem = -urem;
      CHECK_EQ_I(arbint_add_i32(t1, t1, urem), ARBINT_OK);
      CHECK_EQ_I(arbint_tdiv_qr_u32(t2, r, t1, du), ARBINT_OK);
      large_check_eq("tdiv_u32_exact_q", i, t2, d, n, t1, q);
      large_check_i32("tdiv_u32_exact_r", i, r, urem, t1, d);
    }
  }

  arbint_clear(n);
  arbint_clear(d);
  arbint_clear(q);
  arbint_clear(r);
  arbint_clear(t1);
  arbint_clear(t2);
  arbint_clear(ncopy);
  arbint_clear(s);
}

int main(void) {
  test_large_add_sub_properties();
  test_large_mul_properties();
  test_large_tdiv_properties();
  ARBINT_TEST_FINISH("test_large");
}
