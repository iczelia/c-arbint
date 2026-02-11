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

/*  Test power-of-two optimization for division operations.

    These tests verify that division by powers of two produces correct results
    for both truncated (tdiv), floor (fdiv), and ceiling (cdiv) modes.

    Truncated division semantics:
      - Quotient: round toward zero
      - Remainder: same sign as dividend

    Examples:
       7 / 4 =  1,  7 % 4 =  3
      -7 / 4 = -1, -7 % 4 = -3
       7 / -4 = -1,  7 % -4 =  3
      -7 / -4 =  1, -7 % -4 = -3  */

#include "test_framework.h"

#include <stdint.h>
#include <string.h>

ARBINT_TEST_DECLARE_FAILURES();

/* ------------------------------------------------------------------ */
/*  Helper: verify division result using multiplication identity.      */
/*  n == q * d + r  &&  |r| < |d|  &&  sign(r) == sign(n) or r == 0    */
/* ------------------------------------------------------------------ */

static void verify_tdiv_identity(arbint_t n, arbint_t d, arbint_t q,
                                 arbint_t r, arbint_t tmp, arbint_t tmp2) {
  int nsign = arbint_signum(n);
  int rsign = arbint_signum(r);

  /*  q * d + r == n  */
  CHECK_EQ_I(arbint_mul(tmp, q, d), ARBINT_OK);
  CHECK_EQ_I(arbint_add(tmp2, tmp, r), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(tmp2, n), 0);

  /*  |r| < |d|  */
  CHECK_EQ_I(arbint_cmpabs(r, d), -1);

  /*  sign(r) == sign(n) or r == 0  */
  if (rsign != 0)
    CHECK_EQ_I(rsign, nsign);
}

/* ------------------------------------------------------------------ */
/*  Test: small power-of-two divisors (u32 path).                      */
/* ------------------------------------------------------------------ */

static void test_tdiv_pow2_u32_small(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t q;
  arbint_t r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Division by 1 (2^0).  */
  CHECK_EQ_I(arbint_set_i32(n, 17), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 1u), ARBINT_OK);
  check_i32_value(q, 17);
  check_i32_value(r, 0);

  CHECK_EQ_I(arbint_set_i32(n, -17), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 1u), ARBINT_OK);
  check_i32_value(q, -17);
  check_i32_value(r, 0);

  /*  Division by 2 (2^1).  */
  CHECK_EQ_I(arbint_set_i32(n, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 2u), ARBINT_OK);
  check_i32_value(q, 3);
  check_i32_value(r, 1);

  CHECK_EQ_I(arbint_set_i32(n, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 2u), ARBINT_OK);
  check_i32_value(q, -3);
  check_i32_value(r, -1);

  /*  Division by 4 (2^2).  */
  CHECK_EQ_I(arbint_set_i32(n, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 4u), ARBINT_OK);
  check_i32_value(q, 1);
  check_i32_value(r, 3);

  CHECK_EQ_I(arbint_set_i32(n, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 4u), ARBINT_OK);
  check_i32_value(q, -1);
  check_i32_value(r, -3);

  /*  Division by 8 (2^3).  */
  CHECK_EQ_I(arbint_set_i32(n, 17), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 8u), ARBINT_OK);
  check_i32_value(q, 2);
  check_i32_value(r, 1);

  CHECK_EQ_I(arbint_set_i32(n, -17), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 8u), ARBINT_OK);
  check_i32_value(q, -2);
  check_i32_value(r, -1);

  /*  Division by 16 (2^4).  */
  CHECK_EQ_I(arbint_set_i32(n, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 16u), ARBINT_OK);
  check_i32_value(q, 6);
  check_i32_value(r, 4);

  CHECK_EQ_I(arbint_set_i32(n, -100), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 16u), ARBINT_OK);
  check_i32_value(q, -6);
  check_i32_value(r, -4);

  /*  Large power of two: 2^30.  */
  CHECK_EQ_I(arbint_set_i32(n, 1073741824 + 123), ARBINT_OK); /* 2^30 + 123 */
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 1u << 30), ARBINT_OK);
  check_i32_value(q, 1);
  check_i32_value(r, 123);

  /*  Exact division (no remainder).  */
  CHECK_EQ_I(arbint_set_i32(n, 64), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 8u), ARBINT_OK);
  check_i32_value(q, 8);
  check_i32_value(r, 0);

  CHECK_EQ_I(arbint_set_i32(n, -64), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 8u), ARBINT_OK);
  check_i32_value(q, -8);
  check_i32_value(r, 0);

  /*  Zero dividend.  */
  CHECK_EQ_I(arbint_set_i32(n, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 4u), ARBINT_OK);
  check_i32_value(q, 0);
  check_i32_value(r, 0);

  arbint_clear(n);
  arbint_clear(q);
  arbint_clear(r);
  arbint_ctx_clear(&ctx);
}

/* ------------------------------------------------------------------ */
/*  Test: power-of-two with i32 divisor (signed).                      */
/* ------------------------------------------------------------------ */

static void test_tdiv_pow2_i32_signed(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t q;
  arbint_t r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Division by -1 (-(2^0)).  */
  CHECK_EQ_I(arbint_set_i32(n, 17), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_i32(q, r, n, -1), ARBINT_OK);
  check_i32_value(q, -17);
  check_i32_value(r, 0);

  CHECK_EQ_I(arbint_set_i32(n, -17), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_i32(q, r, n, -1), ARBINT_OK);
  check_i32_value(q, 17);
  check_i32_value(r, 0);

  /*  Division by -2 (-(2^1)).  */
  CHECK_EQ_I(arbint_set_i32(n, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_i32(q, r, n, -2), ARBINT_OK);
  check_i32_value(q, -3);
  check_i32_value(r, 1);

  CHECK_EQ_I(arbint_set_i32(n, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_i32(q, r, n, -2), ARBINT_OK);
  check_i32_value(q, 3);
  check_i32_value(r, -1);

  /*  Division by -4 (-(2^2)).  */
  CHECK_EQ_I(arbint_set_i32(n, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_i32(q, r, n, -4), ARBINT_OK);
  check_i32_value(q, -1);
  check_i32_value(r, 3);

  CHECK_EQ_I(arbint_set_i32(n, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_i32(q, r, n, -4), ARBINT_OK);
  check_i32_value(q, 1);
  check_i32_value(r, -3);

  arbint_clear(n);
  arbint_clear(q);
  arbint_clear(r);
  arbint_ctx_clear(&ctx);
}

/* ------------------------------------------------------------------ */
/*  Test: arbint / arbint power-of-two divisor.                        */
/* ------------------------------------------------------------------ */

static void test_tdiv_pow2_arbint(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t d;
  arbint_t q;
  arbint_t r;
  arbint_t tmp;
  arbint_t tmp2;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(tmp, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(tmp2, &ctx), ARBINT_OK);

  /*  Small cases via arbint/arbint.  */
  CHECK_EQ_I(arbint_set_i32(n, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(d, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  check_i32_value(q, 1);
  check_i32_value(r, 3);
  verify_tdiv_identity(n, d, q, r, tmp, tmp2);

  CHECK_EQ_I(arbint_set_i32(n, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  check_i32_value(q, -1);
  check_i32_value(r, -3);
  verify_tdiv_identity(n, d, q, r, tmp, tmp2);

  CHECK_EQ_I(arbint_set_i32(d, -4), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  check_i32_value(q, 1);
  check_i32_value(r, -3);
  verify_tdiv_identity(n, d, q, r, tmp, tmp2);

  CHECK_EQ_I(arbint_set_i32(n, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  check_i32_value(q, -1);
  check_i32_value(r, 3);
  verify_tdiv_identity(n, d, q, r, tmp, tmp2);

  arbint_clear(n);
  arbint_clear(d);
  arbint_clear(q);
  arbint_clear(r);
  arbint_clear(tmp);
  arbint_clear(tmp2);
  arbint_ctx_clear(&ctx);
}

/* ------------------------------------------------------------------ */
/*  Test: large multi-limb power-of-two divisor.                       */
/* ------------------------------------------------------------------ */

static void test_tdiv_pow2_large(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t d;
  arbint_t q;
  arbint_t r;
  arbint_t tmp;
  arbint_t tmp2;
  arbint_limb_t pow2_limb;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(tmp, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(tmp2, &ctx), ARBINT_OK);

  /*  Build a large dividend via repeated squaring: 2^128 + 123.  */
  CHECK_EQ_I(arbint_set_i32(n, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 4 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 16 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 256 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 65536 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 2^32 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 2^64 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 2^128 */
  CHECK_EQ_I(arbint_add_i32(n, n, 123), ARBINT_OK);

  /*  Divisor: 2^64 (represented as limb[1] = 1, limb[0] = 0 on 64-bit).  */
  CHECK_EQ_I(arbint_set_i32(d, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 4 */
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 16 */
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 256 */
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 65536 */
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 2^32 */
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 2^64 */

  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  verify_tdiv_identity(n, d, q, r, tmp, tmp2);

  /*  q should be 2^64, r should be 123.  */
  CHECK_EQ_I(arbint_cmp(q, d), 0);
  check_i32_value(r, 123);

  /*  Test with negative dividend.  */
  CHECK_EQ_I(arbint_neg(n, n), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  verify_tdiv_identity(n, d, q, r, tmp, tmp2);

  /*  q should be -2^64, r should be -123.  */
  CHECK_EQ_I(arbint_neg(tmp, d), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(q, tmp), 0);
  check_i32_value(r, -123);

  /*  Test with single-limb power of two that's larger than dividend.  */
#if ARBINT_LIMB_BITS == 64
  pow2_limb = (arbint_limb_t) 1u << 63;
#else
  pow2_limb = (arbint_limb_t) 1u << 31;
#endif
  CHECK_EQ_I(arbint_set_i32(n, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_resize(d, 1u), ARBINT_OK);
  ARBINT_LIMBS(d)[0] = pow2_limb;
  d[0]._sz = 1;
  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  check_i32_value(q, 0);
  check_i32_value(r, 100);

  arbint_clear(n);
  arbint_clear(d);
  arbint_clear(q);
  arbint_clear(r);
  arbint_clear(tmp);
  arbint_clear(tmp2);
  arbint_ctx_clear(&ctx);
}

/* ------------------------------------------------------------------ */
/*  Test: q-only and r-only variants.                                  */
/* ------------------------------------------------------------------ */

static void test_tdiv_pow2_q_r_only(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t q;
  arbint_t r;
  arbint_t d;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);

  /*  q-only via u32.  */
  CHECK_EQ_I(arbint_set_i32(n, -17), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_q_u32(q, n, 4u), ARBINT_OK);
  check_i32_value(q, -4);

  /*  r-only via u32.  */
  CHECK_EQ_I(arbint_tdiv_r_u32(r, n, 4u), ARBINT_OK);
  check_i32_value(r, -1);

  /*  q-only via i32.  */
  CHECK_EQ_I(arbint_tdiv_q_i32(q, n, -4), ARBINT_OK);
  check_i32_value(q, 4);

  /*  r-only via i32.  */
  CHECK_EQ_I(arbint_tdiv_r_i32(r, n, -4), ARBINT_OK);
  check_i32_value(r, -1);

  /*  q-only via arbint.  */
  CHECK_EQ_I(arbint_set_i32(d, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_q(q, n, d), ARBINT_OK);
  check_i32_value(q, -4);

  /*  r-only via arbint.  */
  CHECK_EQ_I(arbint_tdiv_r(r, n, d), ARBINT_OK);
  check_i32_value(r, -1);

  arbint_clear(n);
  arbint_clear(q);
  arbint_clear(r);
  arbint_clear(d);
  arbint_ctx_clear(&ctx);
}

/* ------------------------------------------------------------------ */
/*  Test: floor division with power-of-two.                            */
/* ------------------------------------------------------------------ */

static void test_fdiv_pow2(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t q;
  arbint_t r;
  arbint_t d;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);

  /*  Floor division: rounds toward -infinity.
      For negative dividend with positive divisor, q = tdiv_q - 1.  */
  CHECK_EQ_I(arbint_set_i32(n, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_fdiv_qr_u32(q, r, n, 4u), ARBINT_OK);
  check_i32_value(q, -2); /* -7/4 = -1.75, floor = -2 */
  check_i32_value(r, 1);  /* r = -7 - (-2)*4 = -7 + 8 = 1 */

  CHECK_EQ_I(arbint_set_i32(n, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_fdiv_qr_u32(q, r, n, 4u), ARBINT_OK);
  check_i32_value(q, 1);
  check_i32_value(r, 3);

  /*  Exact division: no adjustment.  */
  CHECK_EQ_I(arbint_set_i32(n, -8), ARBINT_OK);
  CHECK_EQ_I(arbint_fdiv_qr_u32(q, r, n, 4u), ARBINT_OK);
  check_i32_value(q, -2);
  check_i32_value(r, 0);

  /*  Via arbint divisor.  */
  CHECK_EQ_I(arbint_set_i32(n, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(d, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_fdiv_qr(q, r, n, d), ARBINT_OK);
  check_i32_value(q, -2);
  check_i32_value(r, 1);

  arbint_clear(n);
  arbint_clear(q);
  arbint_clear(r);
  arbint_clear(d);
  arbint_ctx_clear(&ctx);
}

/* ------------------------------------------------------------------ */
/*  Test: ceiling division with power-of-two.                          */
/* ------------------------------------------------------------------ */

static void test_cdiv_pow2(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t q;
  arbint_t r;
  arbint_t d;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);

  /*  Ceiling division: rounds toward +infinity.
      For positive dividend with positive divisor, q = tdiv_q + 1 if r != 0. */
  CHECK_EQ_I(arbint_set_i32(n, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(d, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_cdiv_qr(q, r, n, d), ARBINT_OK);
  check_i32_value(q, 2);  /* 7/4 = 1.75, ceil = 2 */
  check_i32_value(r, -1); /* r = 7 - 2*4 = 7 - 8 = -1 */

  CHECK_EQ_I(arbint_set_i32(n, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_cdiv_qr(q, r, n, d), ARBINT_OK);
  check_i32_value(q, -1);
  check_i32_value(r, -3);

  /*  Exact division: no adjustment.  */
  CHECK_EQ_I(arbint_set_i32(n, 8), ARBINT_OK);
  CHECK_EQ_I(arbint_cdiv_qr(q, r, n, d), ARBINT_OK);
  check_i32_value(q, 2);
  check_i32_value(r, 0);

  /*  Division by 2.  */
  CHECK_EQ_I(arbint_set_i32(n, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(d, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_cdiv_qr(q, r, n, d), ARBINT_OK);
  check_i32_value(q, 4);  /* 7/2 = 3.5, ceil = 4 */
  check_i32_value(r, -1); /* r = 7 - 4*2 = -1 */

  arbint_clear(n);
  arbint_clear(q);
  arbint_clear(r);
  arbint_clear(d);
  arbint_ctx_clear(&ctx);
}

/* ------------------------------------------------------------------ */
/*  Test: aliasing (q == n, r == n, etc.).                             */
/* ------------------------------------------------------------------ */

static void test_tdiv_pow2_aliasing(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t q;
  arbint_t r;
  arbint_t d;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);

  /*  q == n aliasing.  */
  CHECK_EQ_I(arbint_set_i32(n, 17), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(n, r, n, 4u), ARBINT_OK);
  check_i32_value(n, 4);
  check_i32_value(r, 1);

  /*  r == n aliasing.  */
  CHECK_EQ_I(arbint_set_i32(n, 17), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, n, n, 4u), ARBINT_OK);
  check_i32_value(q, 4);
  check_i32_value(n, 1);

  /*  q == r aliasing (last written wins, which is r).  */
  CHECK_EQ_I(arbint_set_i32(n, 17), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, q, n, 4u), ARBINT_OK);
  check_i32_value(q, 1); /* r is written last */

  /*  Aliasing with arbint divisor.  */
  CHECK_EQ_I(arbint_set_i32(n, 17), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(d, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr(n, r, n, d), ARBINT_OK);
  check_i32_value(n, 4);
  check_i32_value(r, 1);

  CHECK_EQ_I(arbint_set_i32(n, 17), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr(q, n, n, d), ARBINT_OK);
  check_i32_value(q, 4);
  check_i32_value(n, 1);

  arbint_clear(n);
  arbint_clear(q);
  arbint_clear(r);
  arbint_clear(d);
  arbint_ctx_clear(&ctx);
}

/* ------------------------------------------------------------------ */
/*  Test: compare pow2 path to non-pow2 for consistency.               */
/* ------------------------------------------------------------------ */

static void test_tdiv_pow2_vs_generic(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t q_pow2;
  arbint_t r_pow2;
  arbint_t q_gen;
  arbint_t r_gen;
  arbint_t d;
  int i;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q_pow2, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r_pow2, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q_gen, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r_gen, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);

  /*  Build a large number: 2^100 + 12345.  */
  CHECK_EQ_I(arbint_set_i32(n, 2), ARBINT_OK);
  for (i = 0; i < 7; ++i) {
    CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK);
  }
  /* n = 2^128, shift right to get 2^100 */
  CHECK_EQ_I(arbint_shr(n, n, 28), ARBINT_OK);
  CHECK_EQ_I(arbint_add_i32(n, n, 12345), ARBINT_OK);

  /*  Test several power-of-two divisors.  */
  uint32_t pow2_divs[] = {1u,    2u,       4u,       8u,      16u,
                          32u,   64u,      128u,     256u,    512u,
                          1024u, 1u << 16, 1u << 20, 1u << 30};
  size_t num_divs = sizeof(pow2_divs) / sizeof(pow2_divs[0]);
  size_t j;

  for (j = 0u; j < num_divs; ++j) {
    uint32_t div = pow2_divs[j];

    /*  Test positive n.  */
    CHECK_EQ_I(arbint_tdiv_qr_u32(q_pow2, r_pow2, n, div), ARBINT_OK);

    /*  Verify: n == q * d + r.  */
    CHECK_EQ_I(arbint_mul_u32(q_gen, q_pow2, div), ARBINT_OK);
    CHECK_EQ_I(arbint_add(q_gen, q_gen, r_pow2), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(q_gen, n), 0);

    /*  Test negative n.  */
    CHECK_EQ_I(arbint_neg(d, n), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_qr_u32(q_pow2, r_pow2, d, div), ARBINT_OK);

    /*  Verify: -n == q * d + r.  */
    CHECK_EQ_I(arbint_mul_u32(q_gen, q_pow2, div), ARBINT_OK);
    CHECK_EQ_I(arbint_add(q_gen, q_gen, r_pow2), ARBINT_OK);
    CHECK_EQ_I(arbint_cmp(q_gen, d), 0);
  }

  arbint_clear(n);
  arbint_clear(q_pow2);
  arbint_clear(r_pow2);
  arbint_clear(q_gen);
  arbint_clear(r_gen);
  arbint_clear(d);
  arbint_ctx_clear(&ctx);
}

/* ------------------------------------------------------------------ */
/*  Test: edge case - dividend smaller than divisor.                   */
/* ------------------------------------------------------------------ */

static void test_tdiv_pow2_dividend_smaller(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t d;
  arbint_t q;
  arbint_t r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Small n, large power-of-two divisor.  */
  CHECK_EQ_I(arbint_set_i32(n, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 1u << 20), ARBINT_OK);
  check_i32_value(q, 0);
  check_i32_value(r, 7);

  CHECK_EQ_I(arbint_set_i32(n, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 1u << 20), ARBINT_OK);
  check_i32_value(q, 0);
  check_i32_value(r, -7);

  /*  Via arbint divisor.  */
  CHECK_EQ_I(arbint_set_i32(n, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(d, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* d = 2^32 */

  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  check_i32_value(q, 0);
  check_i32_value(r, 7);

  arbint_clear(n);
  arbint_clear(d);
  arbint_clear(q);
  arbint_clear(r);
  arbint_ctx_clear(&ctx);
}

/* ------------------------------------------------------------------ */
/*  Test: trailing zeros optimization.                                  */
/*                                                                      */
/*  When dividend has trailing zeros >= divisor's k (for 2^k), the      */
/*  remainder is 0 and we can use a faster quotient computation.        */
/* ------------------------------------------------------------------ */

static void test_tdiv_pow2_trailing_zeros(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t d;
  arbint_t q;
  arbint_t r;
  arbint_t tmp;
  arbint_t tmp2;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(d, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(tmp, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(tmp2, &ctx), ARBINT_OK);

  /*  n = 64 = 2^6, divide by 4 = 2^2.
      n has 6 trailing zeros >= 2, so remainder should be 0.  */
  CHECK_EQ_I(arbint_set_i32(n, 64), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 4u), ARBINT_OK);
  check_i32_value(q, 16);
  check_i32_value(r, 0);

  /*  n = -64, divide by 4.  */
  CHECK_EQ_I(arbint_set_i32(n, -64), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 4u), ARBINT_OK);
  check_i32_value(q, -16);
  check_i32_value(r, 0);

  /*  n = 1024 = 2^10, divide by 8 = 2^3.  */
  CHECK_EQ_I(arbint_set_i32(n, 1024), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 8u), ARBINT_OK);
  check_i32_value(q, 128);
  check_i32_value(r, 0);

  /*  n = 2^32, divide by 2^16.  */
  CHECK_EQ_I(arbint_set_i32(n, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 4 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 16 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 256 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 65536 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 2^32 */

  CHECK_EQ_I(arbint_set_i32(d, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 4 */
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 16 */
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 256 */
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 2^16 */

  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  verify_tdiv_identity(n, d, q, r, tmp, tmp2);
  CHECK_EQ_I(arbint_cmp(q, d), 0); /* q = 2^16 */
  check_i32_value(r, 0);

  /*  n = 2^128, divide by 2^64.
      Both are multi-limb, n has 128 trailing zeros >= 64.  */
  CHECK_EQ_I(arbint_set_i32(n, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 4 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 16 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 256 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 65536 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 2^32 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 2^64 */
  CHECK_EQ_I(arbint_sqr(n, n), ARBINT_OK); /* 2^128 */

  CHECK_EQ_I(arbint_set_i32(d, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 4 */
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 16 */
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 256 */
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 65536 */
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 2^32 */
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 2^64 */

  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  verify_tdiv_identity(n, d, q, r, tmp, tmp2);
  CHECK_EQ_I(arbint_cmp(q, d), 0); /* q = 2^64 */
  check_i32_value(r, 0);

  /*  n = 3 * 2^64, divide by 2^32.
      n has 64 trailing zeros >= 32, so remainder is 0.  */
  CHECK_EQ_I(arbint_set_i32(tmp, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(tmp, tmp), ARBINT_OK);      /* 4 */
  CHECK_EQ_I(arbint_sqr(tmp, tmp), ARBINT_OK);      /* 16 */
  CHECK_EQ_I(arbint_sqr(tmp, tmp), ARBINT_OK);      /* 256 */
  CHECK_EQ_I(arbint_sqr(tmp, tmp), ARBINT_OK);      /* 65536 */
  CHECK_EQ_I(arbint_sqr(tmp, tmp), ARBINT_OK);      /* 2^32 */
  CHECK_EQ_I(arbint_sqr(tmp, tmp), ARBINT_OK);      /* 2^64 */
  CHECK_EQ_I(arbint_mul_i32(n, tmp, 3), ARBINT_OK); /* 3 * 2^64 */

  CHECK_EQ_I(arbint_set_i32(d, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 4 */
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 16 */
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 256 */
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 65536 */
  CHECK_EQ_I(arbint_sqr(d, d), ARBINT_OK); /* 2^32 */

  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  verify_tdiv_identity(n, d, q, r, tmp, tmp2);
  check_i32_value(r, 0);

  /*  Verify quotient is 3 * 2^32.  */
  CHECK_EQ_I(arbint_mul_i32(tmp, d, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_cmp(q, tmp), 0);

  /*  Test negative: n = -3 * 2^64, divide by 2^32.  */
  CHECK_EQ_I(arbint_set_i32(tmp, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(tmp, tmp), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(tmp, tmp), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(tmp, tmp), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(tmp, tmp), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(tmp, tmp), ARBINT_OK);
  CHECK_EQ_I(arbint_sqr(tmp, tmp), ARBINT_OK);       /* 2^64 */
  CHECK_EQ_I(arbint_mul_i32(n, tmp, -3), ARBINT_OK); /* -3 * 2^64 */

  CHECK_EQ_I(arbint_tdiv_qr(q, r, n, d), ARBINT_OK);
  verify_tdiv_identity(n, d, q, r, tmp, tmp2);
  check_i32_value(r, 0);

  /*  q-only test: ensure fast path works when r is NULL.  */
  CHECK_EQ_I(arbint_set_i32(n, 256), ARBINT_OK);       /* 2^8 */
  CHECK_EQ_I(arbint_tdiv_q_u32(q, n, 16u), ARBINT_OK); /* 2^4 */
  check_i32_value(q, 16);                              /* 2^4 */

  /*  r-only test: ensure fast path works when q is NULL.  */
  CHECK_EQ_I(arbint_tdiv_r_u32(r, n, 16u), ARBINT_OK);
  check_i32_value(r, 0);

  /*  Test case where trailing zeros < k (should use general path).  */
  CHECK_EQ_I(arbint_set_i32(n, 12),
             ARBINT_OK); /* 12 = 4 * 3, has 2 trailing 0s */
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 8u),
             ARBINT_OK); /* 2^3 > 2 trailing 0s */
  check_i32_value(q, 1);
  check_i32_value(r, 4);

  arbint_clear(n);
  arbint_clear(d);
  arbint_clear(q);
  arbint_clear(r);
  arbint_clear(tmp);
  arbint_clear(tmp2);
  arbint_ctx_clear(&ctx);
}

/* ------------------------------------------------------------------ */
/*  Test: ensure non-power-of-two still works (regression check).      */
/* ------------------------------------------------------------------ */

static void test_tdiv_non_pow2_regression(void) {
  arbint_ctx_t ctx;
  arbint_t n;
  arbint_t q;
  arbint_t r;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(n, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);

  /*  Non-power-of-two divisors should still work.  */
  CHECK_EQ_I(arbint_set_i32(n, 17), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 5u), ARBINT_OK);
  check_i32_value(q, 3);
  check_i32_value(r, 2);

  CHECK_EQ_I(arbint_set_i32(n, -17), ARBINT_OK);
  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 5u), ARBINT_OK);
  check_i32_value(q, -3);
  check_i32_value(r, -2);

  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 3u), ARBINT_OK);
  check_i32_value(q, -5);
  check_i32_value(r, -2);

  CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, n, 7u), ARBINT_OK);
  check_i32_value(q, -2);
  check_i32_value(r, -3);

  arbint_clear(n);
  arbint_clear(q);
  arbint_clear(r);
  arbint_ctx_clear(&ctx);
}

int main(void) {
  test_tdiv_pow2_u32_small();
  test_tdiv_pow2_i32_signed();
  test_tdiv_pow2_arbint();
  test_tdiv_pow2_large();
  test_tdiv_pow2_q_r_only();
  test_fdiv_pow2();
  test_cdiv_pow2();
  test_tdiv_pow2_aliasing();
  test_tdiv_pow2_vs_generic();
  test_tdiv_pow2_dividend_smaller();
  test_tdiv_pow2_trailing_zeros();
  test_tdiv_non_pow2_regression();

  ARBINT_TEST_FINISH("test_div_pow2");
}
