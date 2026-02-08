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

#include "arbint_cmp.h"

#include <string.h>

/*  Compare absolute value of arbint with uint32_t.
    Returns -1 if |a| < b, 0 if |a| == b, +1 if |a| > b.  */
static int arbint_cmp_abs_u32(const arbint_t a, uint32_t b) {
  size_t asz;

  if (a == NULL)
    return (b == 0u) ? 0 : -1;

  asz = arbint_abs_sz(a[0]._sz);
  if (asz == 0u)
    return (b == 0u) ? 0 : -1;
  if (asz > 1u)
    return 1;
  if (ARBINT_CLIMBS(a)[0] < (arbint_limb_t) b)
    return -1;
  if (ARBINT_CLIMBS(a)[0] > (arbint_limb_t) b)
    return 1;
  return 0;
}

/*  Get sign of arbint: -1 if negative, 0 if zero, +1 if positive.  */
int arbint_signum(const arbint_t x) {
  if (x == NULL || x[0]._sz == 0)
    return 0;
  return (x[0]._sz > 0) ? 1 : -1;
}

/*  Check if arbint is zero.  */
int arbint_is_zero(const arbint_t x) { return arbint_signum(x) == 0; }

int arbint_is_one(const arbint_t x) {
  if (x == NULL || x[0]._sz != 1)
    return 0;
  return ARBINT_CLIMBS(x)[0] == (arbint_limb_t) 1u;
}

int arbint_is_neg(const arbint_t x) { return x != NULL && x[0]._sz < 0; }

int arbint_is_odd(const arbint_t x) {
  if (x == NULL || x[0]._sz == 0)
    return 0;
  return (ARBINT_CLIMBS(x)[0] & (arbint_limb_t) 1u) != 0u;
}

int arbint_is_even(const arbint_t x) {
  if (x == NULL)
    return 0;
  return !arbint_is_odd(x);
}

/*  Set rop to absolute value of op (rop = |op|).
    Handles aliasing (rop == op) correctly.  */
arbint_err_t arbint_abs(arbint_t rop, const arbint_t op) {
  arbint_err_t rc;

  if (rop == NULL || op == NULL)
    return ARBINT_EINVAL;

  if (rop != op) {
    rc = arbint_set(rop, op);
    if (rc != ARBINT_OK)
      return rc;
  }

  if (rop[0]._sz < 0)
    rop[0]._sz = -rop[0]._sz;
  return ARBINT_OK;
}

/*  Negate arbint value (rop = -op).
    Handles aliasing (rop == op) correctly.  */
arbint_err_t arbint_neg(arbint_t rop, const arbint_t op) {
  arbint_err_t rc;

  if (rop == NULL || op == NULL)
    return ARBINT_EINVAL;

  if (rop != op) {
    rc = arbint_set(rop, op);
    if (rc != ARBINT_OK)
      return rc;
  }

  if (rop[0]._sz != 0)
    rop[0]._sz = -rop[0]._sz;
  return ARBINT_OK;
}

/*  Compare absolute values of two arbints.
    Returns -1 if |a| < |b|, 0 if |a| == |b|, +1 if |a| > |b|.  */
int arbint_cmpabs(const arbint_t a, const arbint_t b) {
  size_t asz, bsz, i;

  if (a == b)
    return 0;

  if (a == NULL)
    return (b == NULL || b[0]._sz == 0) ? 0 : -1;
  if (b == NULL)
    return (a[0]._sz == 0) ? 0 : 1;

  asz = arbint_abs_sz(a[0]._sz);
  bsz = arbint_abs_sz(b[0]._sz);
  if (asz < bsz)
    return -1;
  if (asz > bsz)
    return 1;
  if (asz == 0u)
    return 0;

  for (i = asz; i != 0u; --i) {
    arbint_limb_t al = ARBINT_CLIMBS(a)[i - 1u];
    arbint_limb_t bl = ARBINT_CLIMBS(b)[i - 1u];
    if (al < bl)
      return -1;
    if (al > bl)
      return 1;
  }
  return 0;
}

/*  Compare two arbints with sign consideration.
    Returns -1 if a < b, 0 if a == b, +1 if a > b.  */
int arbint_cmp(const arbint_t a, const arbint_t b) {
  int as, bs, c;

  if (a == b)
    return 0;

  as = arbint_signum(a);
  bs = arbint_signum(b);
  if (as < bs)
    return -1;
  if (as > bs)
    return 1;
  if (as == 0)
    return 0;

  c = arbint_cmpabs(a, b);
  return (as > 0) ? c : -c;
}

int arbint_cmp_i32(const arbint_t a, int32_t b) {
  int as = arbint_signum(a);
  int bs = (b > 0) - (b < 0);
  uint32_t bmag = (b < 0) ? (uint32_t) (-(b + 1)) + 1u : (uint32_t) b;
  int c;

  if (as < bs)
    return -1;
  if (as > bs)
    return 1;
  if (as == 0)
    return 0;

  c = arbint_cmp_abs_u32(a, bmag);
  return (as > 0) ? c : -c;
}

int arbint_cmp_u32(const arbint_t a, uint32_t b) {
  int as = arbint_signum(a);
  if (as < 0)
    return -1;
  return arbint_cmp_abs_u32(a, b);
}

int arbint_eq(const arbint_t a, const arbint_t b) {
  size_t n;

  if (a == b)
    return 1;
  if (a == NULL || b == NULL)
    return 0;
  if (a[0]._sz != b[0]._sz)
    return 0;
  if (a[0]._sz == 0)
    return 1;

  n = arbint_abs_sz(a[0]._sz);
  return memcmp(ARBINT_CLIMBS(a), ARBINT_CLIMBS(b),
                n * sizeof(arbint_limb_t)) == 0;
}

int arbint_ne(const arbint_t a, const arbint_t b) { return !arbint_eq(a, b); }

int arbint_lt(const arbint_t a, const arbint_t b) {
  return arbint_cmp(a, b) < 0;
}

int arbint_le(const arbint_t a, const arbint_t b) {
  return arbint_cmp(a, b) <= 0;
}

int arbint_gt(const arbint_t a, const arbint_t b) {
  return arbint_cmp(a, b) > 0;
}

int arbint_ge(const arbint_t a, const arbint_t b) {
  return arbint_cmp(a, b) >= 0;
}

int arbint_eq_i32(const arbint_t a, int32_t b) {
  return arbint_cmp_i32(a, b) == 0;
}

int arbint_ne_i32(const arbint_t a, int32_t b) {
  return arbint_cmp_i32(a, b) != 0;
}

int arbint_lt_i32(const arbint_t a, int32_t b) {
  return arbint_cmp_i32(a, b) < 0;
}

int arbint_le_i32(const arbint_t a, int32_t b) {
  return arbint_cmp_i32(a, b) <= 0;
}

int arbint_gt_i32(const arbint_t a, int32_t b) {
  return arbint_cmp_i32(a, b) > 0;
}

int arbint_ge_i32(const arbint_t a, int32_t b) {
  return arbint_cmp_i32(a, b) >= 0;
}

int arbint_eq_u32(const arbint_t a, uint32_t b) {
  return arbint_cmp_u32(a, b) == 0;
}

int arbint_ne_u32(const arbint_t a, uint32_t b) {
  return arbint_cmp_u32(a, b) != 0;
}

int arbint_lt_u32(const arbint_t a, uint32_t b) {
  return arbint_cmp_u32(a, b) < 0;
}

int arbint_le_u32(const arbint_t a, uint32_t b) {
  return arbint_cmp_u32(a, b) <= 0;
}

int arbint_gt_u32(const arbint_t a, uint32_t b) {
  return arbint_cmp_u32(a, b) > 0;
}

int arbint_ge_u32(const arbint_t a, uint32_t b) {
  return arbint_cmp_u32(a, b) >= 0;
}
