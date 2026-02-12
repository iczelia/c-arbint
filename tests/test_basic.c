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
#include <stdlib.h>
#include <string.h>

typedef struct {
  size_t nonfree_calls;
  size_t alloc_calls;
  size_t realloc_calls;
  size_t free_calls;
  size_t fail_at;
} counting_alloc_state_t;

ARBINT_TEST_DECLARE_FAILURES();

static void * counting_realloc(void * ud, void * ptr, size_t new_size) {
  counting_alloc_state_t * st = (counting_alloc_state_t *) ud;

  if (new_size == 0u) {
    if (ptr != NULL)
      ++st->free_calls;
    free(ptr);
    return NULL;
  }

  ++st->nonfree_calls;
  if (st->fail_at != 0u && st->nonfree_calls >= st->fail_at)
    return NULL;

  if (ptr == NULL)
    ++st->alloc_calls;
  else
    ++st->realloc_calls;

  return realloc(ptr, new_size);
}

static arbint_ctx_t make_counting_ctx(counting_alloc_state_t * st) {
  arbint_ctx_t ctx;
  arbint_alloc_t a;
  a.ud = st;
  a.realloc = counting_realloc;
  CHECK_EQ_I(arbint_ctx_init(&ctx, &a, 0u), ARBINT_OK);
  return ctx;
}

static void test_ctx_and_base(void) {
  arbint_ctx_t ctx;
  arbint_alloc_t a;
  arbint_t x;
  arbint_t y;

  CHECK_EQ_I(arbint_ctx_init_default(NULL), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_ctx_init(NULL, &a, 0u), ARBINT_EINVAL);
  a.ud = NULL;
  a.realloc = NULL;
  CHECK_EQ_I(arbint_ctx_init(&ctx, NULL, 0u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_ctx_init(&ctx, &a, 0u), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(ctx.flags, 0u);
  CHECK_EQ_I(ctx.rng_flags, 0u);

  arbint_ctx_clear(NULL);
  arbint_ctx_clear(&ctx);

  CHECK_EQ_I(arbint_init(NULL, &ctx), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);
  CHECK_EQ_I(x[0]._cap, 0u);
  CHECK_EQ_I(x[0]._sz, 0);
  CHECK(x[0]._ptr == NULL);
  CHECK(arbint_get_ctx(x) == &ctx);
  CHECK(arbint_get_ctx(NULL) == NULL);

  CHECK_EQ_I(arbint_resize(NULL, 1u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_resize(x, 2u), ARBINT_OK);
  CHECK(x[0]._cap >= 2u);

  CHECK_EQ_I(arbint_set_u32(x, 123u), ARBINT_OK);
  CHECK_EQ_I(x[0]._sz, 1);

  arbint_zero(x);
  CHECK_EQ_I(x[0]._sz, 0);
  CHECK(x[0]._cap >= 2u);

  CHECK_EQ_I(arbint_resize(x, 1u), ARBINT_OK);
  CHECK(x[0]._cap >= 2u);

  CHECK_EQ_I(arbint_resize(x, SIZE_MAX), ARBINT_EOVERFLOW);

  CHECK_EQ_I(arbint_init(y, NULL), ARBINT_OK);
  CHECK_EQ_I(arbint_resize(y, 1u), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_set_u32(x, 111u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(y, 222u), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_init(y, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(y, 222u), ARBINT_OK);

  arbint_swap(NULL, y);
  arbint_swap(y, NULL);
  arbint_swap(y, y);

  arbint_swap(x, y);
  check_u32_value(x, 222u);
  check_u32_value(y, 111u);

  arbint_clear(NULL);
  arbint_clear(x);
  CHECK_EQ_I(x[0]._cap, 0u);
  CHECK_EQ_I(x[0]._sz, 0);
  CHECK(x[0]._ptr == NULL);
  CHECK(x[0]._ctx == NULL);

  arbint_clear(y);
}

static void test_set_ctx_and_allocators(void) {
  counting_alloc_state_t st1;
  counting_alloc_state_t st2;
  counting_alloc_state_t st_fail;
  arbint_ctx_t ctx1;
  arbint_ctx_t ctx2;
  arbint_ctx_t ctx_fail;
  arbint_t x;
  arbint_t y;

  memset(&st1, 0, sizeof(st1));
  memset(&st2, 0, sizeof(st2));
  memset(&st_fail, 0, sizeof(st_fail));

  ctx1 = make_counting_ctx(&st1);
  ctx2 = make_counting_ctx(&st2);
  st_fail.fail_at = 1u;
  ctx_fail = make_counting_ctx(&st_fail);

  CHECK_EQ_I(arbint_set_ctx(NULL, &ctx1), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_init(x, &ctx1), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(x, 123456u), ARBINT_OK);
  CHECK(st1.alloc_calls >= 1u);

  CHECK_EQ_I(arbint_set_ctx(x, &ctx1), ARBINT_OK);
  CHECK_EQ_I(arbint_set_ctx(x, NULL), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_set_ctx(x, &ctx2), ARBINT_OK);
  CHECK(arbint_get_ctx(x) == &ctx2);
  check_u32_value(x, 123456u);
  CHECK(st1.free_calls >= 1u);
  CHECK(st2.alloc_calls >= 1u);

  CHECK_EQ_I(arbint_init(y, &ctx1), ARBINT_OK);
  CHECK_EQ_I(arbint_set_u32(y, 77u), ARBINT_OK);
  CHECK_EQ_I(arbint_set_ctx(y, &ctx_fail), ARBINT_ENOMEM);
  CHECK(arbint_get_ctx(y) == &ctx1);
  check_u32_value(y, 77u);

  arbint_clear(y);
  CHECK_EQ_I(arbint_init(y, &ctx1), ARBINT_OK);
  CHECK_EQ_I(arbint_set_ctx(y, NULL), ARBINT_OK);
  CHECK(arbint_get_ctx(y) == NULL);

  {
    arbint_t z;
    arbint_ctx_t bad_ctx;

    memset(&bad_ctx, 0, sizeof(bad_ctx));

    CHECK_EQ_I(arbint_init(z, &ctx1), ARBINT_OK);
    CHECK_EQ_I(arbint_set_u32(z, 5u), ARBINT_OK);
    CHECK_EQ_I(arbint_set_ctx(z, &bad_ctx), ARBINT_EINVAL);

    z[0]._ctx = NULL;
    CHECK_EQ_I(arbint_set_ctx(z, &ctx1), ARBINT_EINVAL);

    z[0]._ctx = &ctx1;
    z[0]._cap = (SIZE_MAX / sizeof(arbint_limb_t)) + 1u;
    CHECK_EQ_I(arbint_set_ctx(z, &ctx2), ARBINT_EOVERFLOW);

    z[0]._cap = 1u;
    z[0]._sz = 2;
    CHECK_EQ_I(arbint_set_ctx(z, &ctx2), ARBINT_EINVAL);

    arbint_clear(z);
  }

  arbint_clear(x);
  arbint_clear(y);
}

static void test_assignment_get_fits(void) {
  arbint_ctx_t ctx;
  arbint_t x;
  arbint_t y;
  int32_t out_i32;
  uint32_t out_u32;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(x, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(y, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_set(NULL, x), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_set(y, NULL), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_set_i32(NULL, 0), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_set_u32(NULL, 0u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_get_i32(NULL, NULL), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_get_u32(NULL, NULL), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_set_i32(x, 0), ARBINT_OK);
  check_i32_value(x, 0);

  CHECK_EQ_I(arbint_set_i32(x, 12345), ARBINT_OK);
  check_i32_value(x, 12345);

  CHECK_EQ_I(arbint_set_i32(x, -12345), ARBINT_OK);
  check_i32_value(x, -12345);

  CHECK_EQ_I(arbint_set_i32(x, INT32_MIN), ARBINT_OK);
  check_i32_value(x, INT32_MIN);

  CHECK_EQ_I(arbint_set_u32(x, 0u), ARBINT_OK);
  check_u32_value(x, 0u);

  CHECK_EQ_I(arbint_set_u32(x, UINT32_MAX), ARBINT_OK);
  check_u32_value(x, UINT32_MAX);

  CHECK_EQ_I(arbint_set(y, x), ARBINT_OK);
  check_u32_value(y, UINT32_MAX);
  CHECK_EQ_I(arbint_set(y, y), ARBINT_OK);

  CHECK_EQ_I(arbint_set_u32(x, 0u), ARBINT_OK);
  CHECK_EQ_I(arbint_set(y, x), ARBINT_OK);
  CHECK(arbint_is_zero(y));

  y[0]._sz = 1;
  y[0]._cap = 1u;
  y[0]._ptr = NULL;
  CHECK_EQ_I(arbint_set(x, y), ARBINT_EINVAL);
  y[0]._sz = 0;

  CHECK_EQ_I(arbint_set_u32(x, (uint32_t) INT32_MAX), ARBINT_OK);
  CHECK_EQ_I(arbint_add_u32(x, x, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_get_i32(x, &out_i32), ARBINT_EOVERFLOW);

  CHECK_EQ_I(arbint_set_i32(x, INT32_MIN), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(x, x, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_get_i32(x, &out_i32), ARBINT_EOVERFLOW);

  CHECK_EQ_I(arbint_set_i32(x, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_get_u32(x, &out_u32), ARBINT_ESIGN);

  CHECK_EQ_I(arbint_set_u32(x, UINT32_MAX), ARBINT_OK);
  CHECK_EQ_I(arbint_add_u32(x, x, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_get_u32(x, &out_u32), ARBINT_EOVERFLOW);

  CHECK_EQ_I(arbint_get_i32(x, NULL), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_get_u32(x, NULL), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_set_u32(x, 255u), ARBINT_OK);
  CHECK(arbint_fits_u8(x));
  CHECK(arbint_fits_u16(x));
  CHECK(arbint_fits_u32(x));
  CHECK(arbint_fits_u64(x));
  CHECK(arbint_fits_i16(x));

  CHECK_EQ_I(arbint_set_u32(x, 256u), ARBINT_OK);
  CHECK(!arbint_fits_u8(x));

  CHECK_EQ_I(arbint_set_i32(x, -128), ARBINT_OK);
  CHECK(arbint_fits_i8(x));
  CHECK(!arbint_fits_u8(x));

  CHECK_EQ_I(arbint_set_i32(x, -129), ARBINT_OK);
  CHECK(!arbint_fits_i8(x));

  {
    arbint_limb_t big[4];
    size_t n = (64u / ARBINT_LIMB_BITS) + 1u;

    memset(big, 0, sizeof(big));
    CHECK(n < 4u);
    big[n - 1u] = 1u;

    set_mag_limbs(x, 1, big, n);
    CHECK(!arbint_fits_u64(x));
    CHECK(!arbint_fits_i64(x));

    set_mag_limbs(x, -1, big, n);
    CHECK(!arbint_fits_u64(x));
    CHECK(!arbint_fits_i64(x));
  }

  CHECK(!arbint_fits_u32(NULL));
  CHECK(!arbint_fits_i32(NULL));

  arbint_clear(x);
  arbint_clear(y);
}

static void test_cmp_and_predicates(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_signum(NULL), 0);
  CHECK(arbint_is_zero(NULL));
  CHECK(!arbint_is_one(NULL));
  CHECK(!arbint_is_neg(NULL));
  CHECK(!arbint_is_odd(NULL));
  CHECK(!arbint_is_even(NULL));

  CHECK_EQ_I(arbint_cmp(NULL, NULL), 0);
  CHECK_EQ_I(arbint_cmpabs(NULL, NULL), 0);
  CHECK(arbint_eq(NULL, NULL));
  CHECK(!arbint_ne(NULL, NULL));
  CHECK_EQ_I(arbint_cmp_u32(NULL, 0u), 0);
  CHECK_EQ_I(arbint_cmp_u32(NULL, 1u), -1);
  CHECK_EQ_I(arbint_cmp_i32(NULL, 0), 0);
  CHECK_EQ_I(arbint_cmp_i32(NULL, -1), 1);

  CHECK_EQ_I(arbint_set_i32(a, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 3), ARBINT_OK);

  CHECK_EQ_I(arbint_signum(a), -1);
  CHECK(arbint_is_neg(a));
  CHECK(arbint_is_odd(a));
  CHECK(!arbint_is_even(a));
  CHECK(!arbint_is_one(a));

  CHECK_EQ_I(arbint_signum(b), 1);
  CHECK(arbint_is_odd(b));

  CHECK_EQ_I(arbint_abs(NULL, a), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_abs(b, NULL), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_neg(NULL, a), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_neg(b, NULL), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_abs(b, a), ARBINT_OK);
  check_i32_value(b, 7);
  check_i32_value(a, -7);

  CHECK_EQ_I(arbint_neg(b, a), ARBINT_OK);
  check_i32_value(b, 7);
  CHECK_EQ_I(arbint_neg(a, a), ARBINT_OK);
  check_i32_value(a, 7);

  CHECK_EQ_I(arbint_set_u32(a, 0u), ARBINT_OK);
  CHECK_EQ_I(arbint_neg(a, a), ARBINT_OK);
  CHECK(arbint_is_zero(a));

  CHECK_EQ_I(arbint_set_i32(a, -10), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -5), ARBINT_OK);

  CHECK_EQ_I(arbint_cmp(a, b), -1);
  CHECK_EQ_I(arbint_cmpabs(a, b), 1);
  CHECK(arbint_lt(a, b));
  CHECK(arbint_le(a, b));
  CHECK(arbint_ne(a, b));
  CHECK(arbint_gt(b, a));
  CHECK(arbint_ge(b, a));

  CHECK_EQ_I(arbint_cmp_i32(a, -10), 0);
  CHECK_EQ_I(arbint_cmp_i32(a, INT32_MIN), 1);
  CHECK_EQ_I(arbint_cmp_u32(a, 0u), -1);

  CHECK(arbint_eq_i32(a, -10));
  CHECK(arbint_ne_i32(a, -9));
  CHECK(arbint_lt_i32(a, -9));
  CHECK(arbint_le_i32(a, -10));
  CHECK(arbint_gt_i32(a, -11));
  CHECK(arbint_ge_i32(a, -10));

  CHECK_EQ_I(arbint_set_u32(b, 5u), ARBINT_OK);
  CHECK(arbint_eq_u32(b, 5u));
  CHECK(arbint_ne_u32(b, 6u));
  CHECK(arbint_lt_u32(b, 6u));
  CHECK(arbint_le_u32(b, 5u));
  CHECK(arbint_gt_u32(b, 4u));
  CHECK(arbint_ge_u32(b, 5u));

  arbint_clear(a);
  arbint_clear(b);
}

static void test_add_sub_small_cases(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t r;
  uint32_t out_u32;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  CHECK_EQ_I(arbint_add(NULL, a, b), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_add(r, NULL, b), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_add(r, a, NULL), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_sub(NULL, a, b), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_sub(r, NULL, b), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_sub(r, a, NULL), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_add_u32(NULL, a, 1u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_add_u32(r, NULL, 1u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_sub_u32(NULL, a, 1u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_sub_u32(r, NULL, 1u), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_add_i32(NULL, a, 1), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_add_i32(r, NULL, 1), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_sub_i32(NULL, a, 1), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_sub_i32(r, NULL, 1), ARBINT_EINVAL);

  CHECK_EQ_I(arbint_set_i32(a, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, -3), ARBINT_OK);
  CHECK_EQ_I(arbint_add(r, a, b), ARBINT_OK);
  check_i32_value(r, 4);
  CHECK_EQ_I(arbint_sub(r, a, b), ARBINT_OK);
  check_i32_value(r, 10);

  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_add(r, a, b), ARBINT_OK);
  check_i32_value(r, 5);
  CHECK_EQ_I(arbint_sub(r, a, b), ARBINT_OK);
  check_i32_value(r, -5);

  CHECK_EQ_I(arbint_set_i32(a, 123), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 77), ARBINT_OK);

  CHECK_EQ_I(arbint_add(a, a, b), ARBINT_OK);
  check_i32_value(a, 200);

  CHECK_EQ_I(arbint_set_i32(a, 123), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 77), ARBINT_OK);
  CHECK_EQ_I(arbint_add(b, a, b), ARBINT_OK);
  check_i32_value(b, 200);

  CHECK_EQ_I(arbint_set_i32(a, 123), ARBINT_OK);
  CHECK_EQ_I(arbint_add(a, a, a), ARBINT_OK);
  check_i32_value(a, 246);

  CHECK_EQ_I(arbint_set_i32(a, 123), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 77), ARBINT_OK);
  CHECK_EQ_I(arbint_sub(a, a, b), ARBINT_OK);
  check_i32_value(a, 46);

  CHECK_EQ_I(arbint_set_i32(a, 123), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 77), ARBINT_OK);
  CHECK_EQ_I(arbint_sub(b, a, b), ARBINT_OK);
  check_i32_value(b, 46);

  CHECK_EQ_I(arbint_set_i32(a, 123), ARBINT_OK);
  CHECK_EQ_I(arbint_sub(a, a, a), ARBINT_OK);
  CHECK(arbint_is_zero(a));

  CHECK_EQ_I(arbint_set_i32(a, -50), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 20), ARBINT_OK);
  CHECK_EQ_I(arbint_add(r, a, b), ARBINT_OK);
  check_i32_value(r, -30);
  CHECK_EQ_I(arbint_sub(r, a, b), ARBINT_OK);
  check_i32_value(r, -70);

  CHECK_EQ_I(arbint_set_i32(a, 20), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 50), ARBINT_OK);
  CHECK_EQ_I(arbint_sub(r, a, b), ARBINT_OK);
  check_i32_value(r, -30);

  CHECK_EQ_I(arbint_set_u32(a, UINT32_MAX), ARBINT_OK);
  CHECK_EQ_I(arbint_add_u32(a, a, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_get_u32(a, &out_u32), ARBINT_EOVERFLOW);

  CHECK_EQ_I(arbint_set_i32(a, -5), ARBINT_OK);
  CHECK_EQ_I(arbint_add_u32(a, a, 10u), ARBINT_OK);
  check_i32_value(a, 5);

  CHECK_EQ_I(arbint_set_i32(a, -15), ARBINT_OK);
  CHECK_EQ_I(arbint_add_u32(a, a, 10u), ARBINT_OK);
  check_i32_value(a, -5);

  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(a, a, 10u), ARBINT_OK);
  check_i32_value(a, -5);

  CHECK_EQ_I(arbint_set_i32(a, -5), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_u32(a, a, 10u), ARBINT_OK);
  check_i32_value(a, -15);

  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_add_i32(a, a, INT32_MIN), ARBINT_OK);
  check_i32_value(a, INT32_MIN + 1);

  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_sub_i32(a, a, INT32_MIN), ARBINT_OK);
  check_u32_value(a, 2147483649u);

  arbint_clear(a);
  arbint_clear(b);
  arbint_clear(r);
}

static void make_pow2(arbint_t x, uint32_t bits) {
  uint32_t i;
  CHECK_EQ_I(arbint_set_u32(x, 1u), ARBINT_OK);
  for (i = 0; i < bits; ++i)
    CHECK_EQ_I(arbint_add(x, x, x), ARBINT_OK);
}

static void test_many_limb_and_complex_ops(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t c;
  arbint_t d;
  arbint_t r;
  arbint_limb_t maxv;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  maxv = (arbint_limb_t) ~((arbint_limb_t) 0u);

  {
    arbint_limb_t am[4] = {maxv, maxv, maxv, maxv};
    arbint_limb_t expect[5] = {0u, 0u, 0u, 0u, 1u};

    set_mag_limbs(a, 1, am, 4u);
    CHECK_EQ_I(arbint_set_u32(b, 1u), ARBINT_OK);

    CHECK_EQ_I(arbint_add(r, a, b), ARBINT_OK);
    CHECK(limbs_equal(r, 1, expect, 5u));

    CHECK_EQ_I(arbint_set(a, r), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(a, a, b), ARBINT_OK);
    CHECK(limbs_equal(a, 1, am, 4u));
  }

  {
    arbint_limb_t am[4] = {0u, 0u, 0u, 1u};
    arbint_limb_t expect[3] = {maxv, maxv, maxv};

    set_mag_limbs(a, 1, am, 4u);
    CHECK_EQ_I(arbint_set_u32(b, 1u), ARBINT_OK);

    CHECK_EQ_I(arbint_sub(r, a, b), ARBINT_OK);
    CHECK(limbs_equal(r, 1, expect, 3u));

    CHECK_EQ_I(arbint_set(a, r), ARBINT_OK);
    CHECK_EQ_I(arbint_add(a, a, b), ARBINT_OK);
    CHECK(limbs_equal(a, 1, am, 4u));
  }

  {
    arbint_limb_t m[5] = {5u, 6u, 7u, 8u, 9u};

    set_mag_limbs(a, 1, m, 5u);
    set_mag_limbs(b, -1, m, 5u);

    CHECK_EQ_I(arbint_add(r, a, b), ARBINT_OK);
    CHECK(arbint_is_zero(r));

    CHECK_EQ_I(arbint_sub(r, a, a), ARBINT_OK);
    CHECK(arbint_is_zero(r));
  }

  {
    uint32_t bits = (uint32_t) (ARBINT_LIMB_BITS * 5u + 11u);

    make_pow2(a, bits);
    CHECK_EQ_I(arbint_add_u32(a, a, 12345u), ARBINT_OK);

    make_pow2(b, bits - 3u);
    CHECK_EQ_I(arbint_add_u32(b, b, 99u), ARBINT_OK);

    CHECK_EQ_I(arbint_set(c, a), ARBINT_OK);
    CHECK_EQ_I(arbint_set(d, b), ARBINT_OK);

    CHECK_EQ_I(arbint_add(r, c, d), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(r, r, d), ARBINT_OK);
    CHECK(arbint_eq(r, c));

    CHECK_EQ_I(arbint_sub(r, c, d), ARBINT_OK);
    CHECK_EQ_I(arbint_add(r, r, d), ARBINT_OK);
    CHECK(arbint_eq(r, c));

    CHECK_EQ_I(arbint_set(r, c), ARBINT_OK);
    CHECK_EQ_I(arbint_add(r, r, d), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(r, r, c), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(r, r, d), ARBINT_OK);
    CHECK(arbint_is_zero(r));

    CHECK_EQ_I(arbint_set(r, c), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(r, r, d), ARBINT_OK);
    CHECK_EQ_I(arbint_add(r, r, d), ARBINT_OK);
    CHECK(arbint_eq(r, c));

    CHECK_EQ_I(arbint_set(r, c), ARBINT_OK);
    CHECK_EQ_I(arbint_add(r, r, r), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(r, r, c), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(r, r, c), ARBINT_OK);
    CHECK(arbint_is_zero(r));
  }

  arbint_clear(a);
  arbint_clear(b);
  arbint_clear(c);
  arbint_clear(d);
  arbint_clear(r);
}

static void test_zero_div_mul_and_overlap_basics(void) {
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t c;
  arbint_t d;
  arbint_t q;
  arbint_t r;
  arbint_limb_t maxv;
  int32_t n = -987654321;
  int32_t dv = 12345;
  int32_t qv = n / dv;
  int32_t rv = n % dv;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  maxv = (arbint_limb_t) ~((arbint_limb_t) 0u);

  {
    arbint_limb_t am[4] = {1u, 2u, maxv, maxv};

    set_mag_limbs(a, 1, am, 4u);
    arbint_zero(b);

    CHECK_EQ_I(arbint_mul(c, a, b), ARBINT_OK);
    CHECK(arbint_is_zero(c));

    CHECK_EQ_I(arbint_set(d, a), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(d, d, b), ARBINT_OK);
    CHECK(arbint_is_zero(d));

    CHECK_EQ_I(arbint_set(d, a), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(d, b, d), ARBINT_OK);
    CHECK(arbint_is_zero(d));
  }

  CHECK_EQ_I(arbint_set_i32(a, -1234567), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 321), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(c, a, b), ARBINT_OK);
  check_i32_value(c, -396296007);

  CHECK_EQ_I(arbint_set_i32(a, -1234567), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 321), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(a, a, b), ARBINT_OK);
  check_i32_value(a, -396296007);

  CHECK_EQ_I(arbint_set_i32(a, -1234567), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, 321), ARBINT_OK);
  CHECK_EQ_I(arbint_mul(b, a, b), ARBINT_OK);
  check_i32_value(b, -396296007);

  CHECK_EQ_I(arbint_set_i32(a, n), ARBINT_OK);
  arbint_zero(b);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, a, b), ARBINT_EZERO);
  CHECK_EQ_I(arbint_tdiv_q(q, a, b), ARBINT_EZERO);
  CHECK_EQ_I(arbint_tdiv_r(r, a, b), ARBINT_EZERO);

  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, a, 0u), ARBINT_EZERO);
  CHECK_EQ_I(arbint_tdiv_q_u32(q, a, 0u), ARBINT_EZERO);
  CHECK_EQ_I(arbint_tdiv_r_u32(r, a, 0u), ARBINT_EZERO);

  CHECK_EQ_I(arbint_tdiv_qr_i32(q, r, a, 0), ARBINT_EZERO);
  CHECK_EQ_I(arbint_tdiv_q_i32(q, a, 0), ARBINT_EZERO);
  CHECK_EQ_I(arbint_tdiv_r_i32(r, a, 0), ARBINT_EZERO);

  CHECK_EQ_I(arbint_set_i32(a, n), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(b, dv), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, a, b), ARBINT_OK);
  check_i32_value(q, qv);
  check_i32_value(r, rv);

  CHECK_EQ_I(arbint_set(c, a), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr(c, r, c, b), ARBINT_OK);
  check_i32_value(c, qv);
  check_i32_value(r, rv);

  CHECK_EQ_I(arbint_set(c, a), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr(q, c, c, b), ARBINT_OK);
  check_i32_value(q, qv);
  check_i32_value(c, rv);

  arbint_clear(a);
  arbint_clear(b);
  arbint_clear(c);
  arbint_clear(d);
  arbint_clear(q);
  arbint_clear(r);
}

int main(void) {
  ARBINT_TEST_START();
  test_ctx_and_base();
  test_set_ctx_and_allocators();
  test_assignment_get_fits();
  test_cmp_and_predicates();
  test_add_sub_small_cases();
  test_many_limb_and_complex_ops();
  test_zero_div_mul_and_overlap_basics();
  ARBINT_TEST_FINISH("test_basic");
}
