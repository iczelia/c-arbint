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

#ifndef ARBINT_TEST_FRAMEWORK_H
#define ARBINT_TEST_FRAMEWORK_H

#include <arbint.h>

#include "../src/arbint_base.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define ARBINT_TEST_DECLARE_FAILURES() static int g_failures = 0

#define CHECK(cond)                                                           \
  do {                                                                        \
    if (!(cond)) {                                                            \
      fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);         \
      ++g_failures;                                                           \
    }                                                                         \
  } while (0)

#define CHECK_EQ_I(a, b) CHECK((a) == (b))
#define CHECK_NE_I(a, b) CHECK((a) != (b))

static inline int arbint_test_set_mag_limbs_raw(arbint_t x, int sign,
                                                const arbint_limb_t * limbs,
                                                size_t n) {
  if (x == NULL)
    return 0;
  if (!((sign == -1) || (sign == 1) || (n == 0u)))
    return 0;
  if (n == 0u) {
    arbint_zero(x);
    return 1;
  }
  if (limbs == NULL || limbs[n - 1u] == (arbint_limb_t) 0u)
    return 0;
  if (arbint_resize(x, n) != ARBINT_OK)
    return 0;

  memcpy(ARBINT_LIMBS(x), limbs, n * sizeof(arbint_limb_t));
  x[0]._sz = (sign < 0) ? -(ptrdiff_t) n : (ptrdiff_t) n;
  return 1;
}

static inline int arbint_test_limbs_equal(const arbint_t x, int sign,
                                          const arbint_limb_t * limbs,
                                          size_t n) {
  size_t used;

  if (x == NULL)
    return 0;

  used = arbint_abs_sz(x[0]._sz);
  if (n == 0u)
    return x[0]._sz == 0;

  if (used != n)
    return 0;

  if (sign > 0 && x[0]._sz <= 0)
    return 0;
  if (sign < 0 && x[0]._sz >= 0)
    return 0;

  if (x[0]._ptr == NULL)
    return 0;

  return memcmp(ARBINT_CLIMBS(x), limbs, n * sizeof(arbint_limb_t)) == 0;
}

static inline int arbint_test_i32_eq(const arbint_t x, int32_t expected) {
  int32_t out = 0;
  return arbint_get_i32(x, &out) == ARBINT_OK && out == expected;
}

static inline int arbint_test_u32_eq(const arbint_t x, uint32_t expected) {
  uint32_t out = 0;
  return arbint_get_u32(x, &out) == ARBINT_OK && out == expected;
}

static inline int arbint_test_finish_impl(const char * test_name,
                                          int failures) {
  if (failures != 0) {
    fprintf(stderr, "%s: %d failure(s)\n", test_name, failures);
    return 1;
  }

  printf("%s: all tests passed\n", test_name);
  return 0;
}

#define set_mag_limbs(x, sign, limbs, n)                                      \
  do {                                                                        \
    CHECK(arbint_test_set_mag_limbs_raw((x), (sign), (limbs), (n)));          \
  } while (0)

#define limbs_equal arbint_test_limbs_equal

#define check_i32_value(x, expected) CHECK(arbint_test_i32_eq((x), (expected)))
#define check_u32_value(x, expected) CHECK(arbint_test_u32_eq((x), (expected)))

#define ARBINT_TEST_FINISH(test_name)                                         \
  return arbint_test_finish_impl((test_name), g_failures)

#endif /* ARBINT_TEST_FRAMEWORK_H */
