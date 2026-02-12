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

/*  ================================================================
    Jacobi symbol tests
    ================================================================  */

static void test_jacobi_basic(arbint_ctx_t * ctx) {
  arbint_t a, n;
  int out;

  CHECK_EQ_I(arbint_init_all(ctx, a, n, (arbint_t *) NULL), ARBINT_OK);

  /*  (a/1) = 1 for all a.  */
  CHECK_EQ_I(arbint_set_i32(n, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(a, -7), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  (0/n) = 0 for n > 1.  */
  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(n, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  CHECK_EQ_I(arbint_set_i32(n, 15), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  /*  (1/n) = 1 for all odd n > 0.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(n, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(n, 15), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  arbint_clear_all(a, n, (arbint_t *) NULL);
}

/*  (2/p) rule: 1 if p == +-1 (mod 8), -1 if p == +-3 (mod 8).  */
static void test_jacobi_two_rule(arbint_ctx_t * ctx) {
  arbint_t a, n;
  int out;

  CHECK_EQ_I(arbint_init_all(ctx, a, n, (arbint_t *) NULL), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);

  /*  p = 7: 7 == 7 (mod 8) == -1 (mod 8), so (2/7) = 1.  */
  CHECK_EQ_I(arbint_set_i32(n, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  p = 17: 17 == 1 (mod 8), so (2/17) = 1.  */
  CHECK_EQ_I(arbint_set_i32(n, 17), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  p = 3: 3 == 3 (mod 8), so (2/3) = -1.  */
  CHECK_EQ_I(arbint_set_i32(n, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, -1);

  /*  p = 5: 5 == 5 (mod 8), so (2/5) = -1.  */
  CHECK_EQ_I(arbint_set_i32(n, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, -1);

  /*  p = 11: 11 == 3 (mod 8), so (2/11) = -1.  */
  CHECK_EQ_I(arbint_set_i32(n, 11), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, -1);

  arbint_clear_all(a, n, (arbint_t *) NULL);
}

/*  (-1/p) rule: 1 if p == 1 (mod 4), -1 if p == 3 (mod 4).  */
static void test_jacobi_minus_one_rule(arbint_ctx_t * ctx) {
  arbint_t a, n;
  int out;

  CHECK_EQ_I(arbint_init_all(ctx, a, n, (arbint_t *) NULL), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(a, -1), ARBINT_OK);

  /*  p = 5: 5 == 1 (mod 4), so (-1/5) = 1.  */
  CHECK_EQ_I(arbint_set_i32(n, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  p = 13: 13 == 1 (mod 4), so (-1/13) = 1.  */
  CHECK_EQ_I(arbint_set_i32(n, 13), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  p = 3: 3 == 3 (mod 4), so (-1/3) = -1.  */
  CHECK_EQ_I(arbint_set_i32(n, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, -1);

  /*  p = 7: 7 == 3 (mod 4), so (-1/7) = -1.  */
  CHECK_EQ_I(arbint_set_i32(n, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, -1);

  arbint_clear_all(a, n, (arbint_t *) NULL);
}

/*  Test Jacobi with composite moduli.  */
static void test_jacobi_composite(arbint_ctx_t * ctx) {
  arbint_t a, n;
  int out;

  CHECK_EQ_I(arbint_init_all(ctx, a, n, (arbint_t *) NULL), ARBINT_OK);

  /*  (5/21) = (5/3)(5/7) = (2/3)(5/7) = (-1)(-1) = 1.
      5 mod 3 = 2, 5 mod 7 = 5.
      (2/3) = -1 (since 3 == 3 mod 8).
      (5/7): 5 is not a QR mod 7 (5^3 = 125 = 6 mod 7 != 1), so (5/7) = -1.  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(n, 21), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  (2/15) = (2/3)(2/5) = (-1)(-1) = 1.  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(n, 15), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  (3/35) = (3/5)(3/7) = (-1)(-1) = 1.
      3^2 = 9 = 4 mod 5, not 1, so (3/5) = -1.
      3^3 = 27 = 6 mod 7 != 1, so (3/7) = -1.  */
  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(n, 35), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  (6/35): gcd(6, 35) = 1.
      (6/35) = (6/5)(6/7) = (1/5)(6/7) = 1 * (6/7).
      6 mod 7 = 6, (6/7) = (-1/7) = -1.  */
  CHECK_EQ_I(arbint_set_i32(a, 6), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, -1);

  arbint_clear_all(a, n, (arbint_t *) NULL);
}

/*  Quadratic reciprocity: (p/q)(q/p) = (-1)^((p-1)/2 * (q-1)/2).  */
static void test_jacobi_reciprocity(arbint_ctx_t * ctx) {
  arbint_t p, q;
  int pq, qp;

  CHECK_EQ_I(arbint_init_all(ctx, p, q, (arbint_t *) NULL), ARBINT_OK);

  /*  p=3, q=7: both == 3 mod 4, so (3/7)(7/3) = -1.  */
  CHECK_EQ_I(arbint_set_i32(p, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(q, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(p, q, &pq), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(q, p, &qp), ARBINT_OK);
  CHECK_EQ_I(pq * qp, -1);

  /*  p=5, q=11: 5 == 1 mod 4, 11 == 3 mod 4, so (5/11)(11/5) = 1.  */
  CHECK_EQ_I(arbint_set_i32(p, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(q, 11), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(p, q, &pq), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(q, p, &qp), ARBINT_OK);
  CHECK_EQ_I(pq * qp, 1);

  /*  p=13, q=17: both == 1 mod 4, so (13/17)(17/13) = 1.  */
  CHECK_EQ_I(arbint_set_i32(p, 13), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(q, 17), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(p, q, &pq), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(q, p, &qp), ARBINT_OK);
  CHECK_EQ_I(pq * qp, 1);

  arbint_clear_all(p, q, (arbint_t *) NULL);
}

/*  Test Jacobi errors.  */
static void test_jacobi_errors(arbint_ctx_t * ctx) {
  arbint_t a, n;
  int out;

  CHECK_EQ_I(arbint_init_all(ctx, a, n, (arbint_t *) NULL), ARBINT_OK);

  /*  NULL arguments.  */
  CHECK_EQ_I(arbint_jacobi(NULL, n, &out), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_jacobi(a, NULL, &out), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_jacobi(a, n, NULL), ARBINT_EINVAL);

  /*  n <= 0: domain error.  */
  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(n, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_EDOM);

  CHECK_EQ_I(arbint_set_i32(n, -5), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_EDOM);

  /*  n even: domain error.  */
  CHECK_EQ_I(arbint_set_i32(n, 6), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_EDOM);

  arbint_clear_all(a, n, (arbint_t *) NULL);
}

/*  Test Jacobi with large numbers.  */
static void test_jacobi_large(arbint_ctx_t * ctx) {
  arbint_t a, n, tmp;
  int out;

  CHECK_EQ_I(arbint_init_all(ctx, a, n, tmp, (arbint_t *) NULL), ARBINT_OK);

  /*  n = 3^50 (large odd number).  */
  CHECK_EQ_I(arbint_set_i32(n, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(n, n, 50), ARBINT_OK);

  /*  (1/n) = 1.  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  (2/n): n = 3^50. 3^50 mod 8 = 3^2 mod 8 = 1 (since 3^2 = 9 = 1 mod 8).
      So (2/n) = 1.  */
  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  n = 7^30 + 2 (ensure odd). If this is even, skip.  */
  CHECK_EQ_I(arbint_set_i32(n, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(n, n, 30), ARBINT_OK);
  CHECK_EQ_I(arbint_add_i32(n, n, 2), ARBINT_OK);

  /*  Check if odd before testing.  */
  {
    uint32_t low;
    CHECK_EQ_I(arbint_fdiv_r_u32(tmp, n, 2u), ARBINT_OK);
    CHECK_EQ_I(arbint_get_u32(tmp, &low), ARBINT_OK);
    if (low == 1u) {
      CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
      CHECK_EQ_I(arbint_jacobi(a, n, &out), ARBINT_OK);
      /*  Just verify it completes without error.  */
      CHECK(out == -1 || out == 0 || out == 1);
    }
  }

  arbint_clear_all(a, n, tmp, (arbint_t *) NULL);
}

/*  ================================================================
    Kronecker symbol tests
    ================================================================  */

static void test_kronecker_basic(arbint_ctx_t * ctx) {
  arbint_t a, n;
  int out;

  CHECK_EQ_I(arbint_init_all(ctx, a, n, (arbint_t *) NULL), ARBINT_OK);

  /*  (a/0) = 1 if |a| = 1, else 0.  */
  CHECK_EQ_I(arbint_set_i32(n, 0), ARBINT_OK);

  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_kronecker(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(a, -1), ARBINT_OK);
  CHECK_EQ_I(arbint_kronecker(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_kronecker(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  CHECK_EQ_I(arbint_set_i32(a, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_kronecker(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  arbint_clear_all(a, n, (arbint_t *) NULL);
}

/*  (a/-1) = -1 if a < 0.  */
static void test_kronecker_negative_n(arbint_ctx_t * ctx) {
  arbint_t a, n;
  int out;

  CHECK_EQ_I(arbint_init_all(ctx, a, n, (arbint_t *) NULL), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(n, -1), ARBINT_OK);

  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_kronecker(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(a, -5), ARBINT_OK);
  CHECK_EQ_I(arbint_kronecker(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, -1);

  /*  (a/-7) = (a/-1)(a/7).  */
  CHECK_EQ_I(arbint_set_i32(n, -7), ARBINT_OK);

  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_kronecker(a, n, &out), ARBINT_OK);
  /*  (2/-7) = (2/-1)(2/7) = 1 * 1 = 1.  */
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(a, -2), ARBINT_OK);
  CHECK_EQ_I(arbint_kronecker(a, n, &out), ARBINT_OK);
  /*  (-2/-7) = (-2/-1)(-2/7) = -1 * (-1/7)(2/7) = -1 * (-1) * 1 = 1.  */
  CHECK_EQ_I(out, 1);

  arbint_clear_all(a, n, (arbint_t *) NULL);
}

/*  (a/2) and (a/2^k).  */
static void test_kronecker_power_of_two(arbint_ctx_t * ctx) {
  arbint_t a, n;
  int out;

  CHECK_EQ_I(arbint_init_all(ctx, a, n, (arbint_t *) NULL), ARBINT_OK);

  /*  (a/2) = 0 if a even.  */
  CHECK_EQ_I(arbint_set_i32(n, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(a, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_kronecker(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  /*  (a/2) for odd a depends on a mod 8.
      (a/2) = 1 if a == +-1 (mod 8), -1 if a == +-3 (mod 8).  */
  CHECK_EQ_I(arbint_set_i32(a, 1), ARBINT_OK); /*  1 mod 8  */
  CHECK_EQ_I(arbint_kronecker(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(a, 7), ARBINT_OK); /*  7 == -1 mod 8  */
  CHECK_EQ_I(arbint_kronecker(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK); /*  3 mod 8  */
  CHECK_EQ_I(arbint_kronecker(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, -1);

  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK); /*  5 == -3 mod 8  */
  CHECK_EQ_I(arbint_kronecker(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, -1);

  /*  (a/4) = (a/2)^2 = 1 for odd a.  */
  CHECK_EQ_I(arbint_set_i32(n, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_kronecker(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(a, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_kronecker(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  (a/8) = (a/2)^3 = (a/2) for odd a.  */
  CHECK_EQ_I(arbint_set_i32(n, 8), ARBINT_OK);
  CHECK_EQ_I(arbint_set_i32(a, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_kronecker(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, -1);

  CHECK_EQ_I(arbint_set_i32(a, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_kronecker(a, n, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  arbint_clear_all(a, n, (arbint_t *) NULL);
}

/*  Kronecker reduces to Jacobi for odd n.  */
static void test_kronecker_vs_jacobi(arbint_ctx_t * ctx) {
  arbint_t a, n;
  int kron, jac;

  CHECK_EQ_I(arbint_init_all(ctx, a, n, (arbint_t *) NULL), ARBINT_OK);

  /*  Test several values where Kronecker should match Jacobi.  */
  CHECK_EQ_I(arbint_set_i32(n, 15), ARBINT_OK);

  CHECK_EQ_I(arbint_set_i32(a, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_kronecker(a, n, &kron), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &jac), ARBINT_OK);
  CHECK_EQ_I(kron, jac);

  CHECK_EQ_I(arbint_set_i32(a, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_kronecker(a, n, &kron), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &jac), ARBINT_OK);
  CHECK_EQ_I(kron, jac);

  CHECK_EQ_I(arbint_set_i32(a, -3), ARBINT_OK);
  CHECK_EQ_I(arbint_kronecker(a, n, &kron), ARBINT_OK);
  CHECK_EQ_I(arbint_jacobi(a, n, &jac), ARBINT_OK);
  CHECK_EQ_I(kron, jac);

  arbint_clear_all(a, n, (arbint_t *) NULL);
}

/*  ================================================================
    Moebius function tests
    ================================================================  */

static void test_moebius_basic(arbint_ctx_t * ctx) {
  arbint_t n;
  int out;

  CHECK_EQ_I(arbint_init(n, ctx), ARBINT_OK);

  /*  mu(1) = 1.  */
  CHECK_EQ_I(arbint_set_i32(n, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  mu(p) = -1 for prime p.  */
  CHECK_EQ_I(arbint_set_i32(n, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_OK);
  CHECK_EQ_I(out, -1);

  CHECK_EQ_I(arbint_set_i32(n, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_OK);
  CHECK_EQ_I(out, -1);

  CHECK_EQ_I(arbint_set_i32(n, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_OK);
  CHECK_EQ_I(out, -1);

  CHECK_EQ_I(arbint_set_i32(n, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_OK);
  CHECK_EQ_I(out, -1);

  /*  mu(pq) = 1 for distinct primes p, q.  */
  CHECK_EQ_I(arbint_set_i32(n, 6), ARBINT_OK); /*  2 * 3  */
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(n, 10), ARBINT_OK); /*  2 * 5  */
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  CHECK_EQ_I(arbint_set_i32(n, 15), ARBINT_OK); /*  3 * 5  */
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  /*  mu(pqr) = -1 for distinct primes p, q, r.  */
  CHECK_EQ_I(arbint_set_i32(n, 30), ARBINT_OK); /*  2 * 3 * 5  */
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_OK);
  CHECK_EQ_I(out, -1);

  CHECK_EQ_I(arbint_set_i32(n, 42), ARBINT_OK); /*  2 * 3 * 7  */
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_OK);
  CHECK_EQ_I(out, -1);

  arbint_clear(n);
}

static void test_moebius_squared_factor(arbint_ctx_t * ctx) {
  arbint_t n;
  int out;

  CHECK_EQ_I(arbint_init(n, ctx), ARBINT_OK);

  /*  mu(p^2) = 0.  */
  CHECK_EQ_I(arbint_set_i32(n, 4), ARBINT_OK); /*  2^2  */
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  CHECK_EQ_I(arbint_set_i32(n, 9), ARBINT_OK); /*  3^2  */
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  CHECK_EQ_I(arbint_set_i32(n, 25), ARBINT_OK); /*  5^2  */
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  /*  mu(p^2 * q) = 0.  */
  CHECK_EQ_I(arbint_set_i32(n, 12), ARBINT_OK); /*  4 * 3  */
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  CHECK_EQ_I(arbint_set_i32(n, 18), ARBINT_OK); /*  2 * 9  */
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  CHECK_EQ_I(arbint_set_i32(n, 50), ARBINT_OK); /*  2 * 25  */
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  /*  mu(p^3) = 0.  */
  CHECK_EQ_I(arbint_set_i32(n, 8), ARBINT_OK);
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  CHECK_EQ_I(arbint_set_i32(n, 27), ARBINT_OK);
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_OK);
  CHECK_EQ_I(out, 0);

  arbint_clear(n);
}

static void test_moebius_errors(arbint_ctx_t * ctx) {
  arbint_t n;
  int out;

  CHECK_EQ_I(arbint_init(n, ctx), ARBINT_OK);

  /*  NULL arguments.  */
  CHECK_EQ_I(arbint_moebius(NULL, n), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_moebius(&out, NULL), ARBINT_EINVAL);

  /*  n <= 0: domain error.  */
  CHECK_EQ_I(arbint_set_i32(n, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_EDOM);

  CHECK_EQ_I(arbint_set_i32(n, -5), ARBINT_OK);
  CHECK_EQ_I(arbint_moebius(&out, n), ARBINT_EDOM);

  arbint_clear(n);
}

/*  Verify sum of mu(d) for d|n equals [n=1] (Mobius inversion).  */
static void test_moebius_sum_property(arbint_ctx_t * ctx) {
  arbint_t n, d, q, r;
  int mu_d, sum;
  uint32_t test_vals[] = {1, 2, 6, 12, 30, 60, 100};
  size_t i, j;

  CHECK_EQ_I(arbint_init_all(ctx, n, d, q, r, (arbint_t *) NULL), ARBINT_OK);

  for (i = 0; i < sizeof(test_vals) / sizeof(test_vals[0]); ++i) {
    uint32_t nv = test_vals[i];
    CHECK_EQ_I(arbint_set_u32(n, nv), ARBINT_OK);
    sum = 0;

    for (j = 1; j <= nv; ++j) {
      if (nv % j != 0)
        continue;

      CHECK_EQ_I(arbint_set_u32(d, (uint32_t) j), ARBINT_OK);
      CHECK_EQ_I(arbint_moebius(&mu_d, d), ARBINT_OK);
      sum += mu_d;
    }

    /*  Sum should be 1 if n=1, else 0.  */
    CHECK_EQ_I(sum, (nv == 1u) ? 1 : 0);
  }

  arbint_clear_all(n, d, q, r, (arbint_t *) NULL);
}

/*  ================================================================
    Carmichael function tests
    ================================================================  */

static void test_carmichael_basic(arbint_ctx_t * ctx) {
  arbint_t n, result;

  CHECK_EQ_I(arbint_init_all(ctx, n, result, (arbint_t *) NULL), ARBINT_OK);

  /*  lambda(1) = 1.  */
  CHECK_EQ_I(arbint_set_i32(n, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 1u);

  /*  lambda(2) = 1.  */
  CHECK_EQ_I(arbint_set_i32(n, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 1u);

  /*  lambda(3) = 2 (prime: p-1 = 2).  */
  CHECK_EQ_I(arbint_set_i32(n, 3), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 2u);

  /*  lambda(4) = 2.  */
  CHECK_EQ_I(arbint_set_i32(n, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 2u);

  /*  lambda(5) = 4 (prime: p-1 = 4).  */
  CHECK_EQ_I(arbint_set_i32(n, 5), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 4u);

  /*  lambda(6) = lcm(lambda(2), lambda(3)) = lcm(1, 2) = 2.  */
  CHECK_EQ_I(arbint_set_i32(n, 6), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 2u);

  /*  lambda(7) = 6 (prime: p-1 = 6).  */
  CHECK_EQ_I(arbint_set_i32(n, 7), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 6u);

  /*  lambda(8) = 2 (2^3: 2^(3-2) = 2).  */
  CHECK_EQ_I(arbint_set_i32(n, 8), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 2u);

  /*  lambda(9) = 6 (3^2: 3^1 * 2 = 6).  */
  CHECK_EQ_I(arbint_set_i32(n, 9), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 6u);

  /*  lambda(10) = lcm(lambda(2), lambda(5)) = lcm(1, 4) = 4.  */
  CHECK_EQ_I(arbint_set_i32(n, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 4u);

  arbint_clear_all(n, result, (arbint_t *) NULL);
}

static void test_carmichael_powers_of_two(arbint_ctx_t * ctx) {
  arbint_t n, result;

  CHECK_EQ_I(arbint_init_all(ctx, n, result, (arbint_t *) NULL), ARBINT_OK);

  /*  lambda(2) = 1.  */
  CHECK_EQ_I(arbint_set_i32(n, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 1u);

  /*  lambda(4) = 2.  */
  CHECK_EQ_I(arbint_set_i32(n, 4), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 2u);

  /*  lambda(8) = 2.  */
  CHECK_EQ_I(arbint_set_i32(n, 8), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 2u);

  /*  lambda(16) = 4.  */
  CHECK_EQ_I(arbint_set_i32(n, 16), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 4u);

  /*  lambda(32) = 8.  */
  CHECK_EQ_I(arbint_set_i32(n, 32), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 8u);

  /*  lambda(64) = 16.  */
  CHECK_EQ_I(arbint_set_i32(n, 64), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 16u);

  arbint_clear_all(n, result, (arbint_t *) NULL);
}

static void test_carmichael_prime_powers(arbint_ctx_t * ctx) {
  arbint_t n, result;

  CHECK_EQ_I(arbint_init_all(ctx, n, result, (arbint_t *) NULL), ARBINT_OK);

  /*  lambda(3^2 = 9) = 3^1 * 2 = 6.  */
  CHECK_EQ_I(arbint_set_i32(n, 9), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 6u);

  /*  lambda(3^3 = 27) = 3^2 * 2 = 18.  */
  CHECK_EQ_I(arbint_set_i32(n, 27), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 18u);

  /*  lambda(5^2 = 25) = 5^1 * 4 = 20.  */
  CHECK_EQ_I(arbint_set_i32(n, 25), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 20u);

  /*  lambda(7^2 = 49) = 7^1 * 6 = 42.  */
  CHECK_EQ_I(arbint_set_i32(n, 49), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 42u);

  arbint_clear_all(n, result, (arbint_t *) NULL);
}

static void test_carmichael_composite(arbint_ctx_t * ctx) {
  arbint_t n, result;

  CHECK_EQ_I(arbint_init_all(ctx, n, result, (arbint_t *) NULL), ARBINT_OK);

  /*  lambda(15) = lcm(lambda(3), lambda(5)) = lcm(2, 4) = 4.  */
  CHECK_EQ_I(arbint_set_i32(n, 15), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 4u);

  /*  lambda(21) = lcm(lambda(3), lambda(7)) = lcm(2, 6) = 6.  */
  CHECK_EQ_I(arbint_set_i32(n, 21), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 6u);

  /*  lambda(24) = lcm(lambda(8), lambda(3)) = lcm(2, 2) = 2.  */
  CHECK_EQ_I(arbint_set_i32(n, 24), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 2u);

  /*  lambda(35) = lcm(lambda(5), lambda(7)) = lcm(4, 6) = 12.  */
  CHECK_EQ_I(arbint_set_i32(n, 35), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 12u);

  /*  lambda(100) = lcm(lambda(4), lambda(25)) = lcm(2, 20) = 20.  */
  CHECK_EQ_I(arbint_set_i32(n, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_OK);
  check_u32_value(result, 20u);

  arbint_clear_all(n, result, (arbint_t *) NULL);
}

static void test_carmichael_errors(arbint_ctx_t * ctx) {
  arbint_t n, result;

  CHECK_EQ_I(arbint_init_all(ctx, n, result, (arbint_t *) NULL), ARBINT_OK);

  /*  NULL arguments.  */
  CHECK_EQ_I(arbint_carmichael(NULL, n), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_carmichael(result, NULL), ARBINT_EINVAL);

  /*  n <= 0: domain error.  */
  CHECK_EQ_I(arbint_set_i32(n, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_EDOM);

  CHECK_EQ_I(arbint_set_i32(n, -5), ARBINT_OK);
  CHECK_EQ_I(arbint_carmichael(result, n), ARBINT_EDOM);

  arbint_clear_all(n, result, (arbint_t *) NULL);
}

/*  Verify lambda(n) divides totient(n) for several values.  */
static void test_carmichael_divides_totient(arbint_ctx_t * ctx) {
  arbint_t n, lambda, totient, q, r;
  uint32_t test_ns[] = {5, 7, 12, 15, 21, 35, 100};
  size_t i;

  CHECK_EQ_I(
      arbint_init_all(ctx, n, lambda, totient, q, r, (arbint_t *) NULL),
      ARBINT_OK);

  for (i = 0; i < sizeof(test_ns) / sizeof(test_ns[0]); ++i) {
    uint32_t nv = test_ns[i];
    CHECK_EQ_I(arbint_set_u32(n, nv), ARBINT_OK);
    CHECK_EQ_I(arbint_carmichael(lambda, n), ARBINT_OK);
    CHECK_EQ_I(arbint_totient(totient, n), ARBINT_OK);

    /*  Verify totient(n) is divisible by lambda(n).  */
    CHECK_EQ_I(arbint_tdiv_qr(q, r, totient, lambda), ARBINT_OK);
    CHECK(arbint_is_zero(r));
  }

  arbint_clear_all(n, lambda, totient, q, r, (arbint_t *) NULL);
}

/*  ================================================================
    Main entry point
    ================================================================  */

int main(void) {
  ARBINT_TEST_START();
  arbint_ctx_t ctx;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);

  /*  Jacobi tests.  */
  test_jacobi_basic(&ctx);
  test_jacobi_two_rule(&ctx);
  test_jacobi_minus_one_rule(&ctx);
  test_jacobi_composite(&ctx);
  test_jacobi_reciprocity(&ctx);
  test_jacobi_errors(&ctx);
  test_jacobi_large(&ctx);

  /*  Kronecker tests.  */
  test_kronecker_basic(&ctx);
  test_kronecker_negative_n(&ctx);
  test_kronecker_power_of_two(&ctx);
  test_kronecker_vs_jacobi(&ctx);

  /*  Moebius tests.  */
  test_moebius_basic(&ctx);
  test_moebius_squared_factor(&ctx);
  test_moebius_errors(&ctx);
  test_moebius_sum_property(&ctx);

  /*  Carmichael tests.  */
  test_carmichael_basic(&ctx);
  test_carmichael_powers_of_two(&ctx);
  test_carmichael_prime_powers(&ctx);
  test_carmichael_composite(&ctx);
  test_carmichael_errors(&ctx);
  test_carmichael_divides_totient(&ctx);

  arbint_ctx_clear(&ctx);

  ARBINT_TEST_FINISH("test_numtheory");
}
