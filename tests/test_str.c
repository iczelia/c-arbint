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

/*  Test suite for arbint_set_str and arbint_get_str.  */

#include "test_framework.h"

#include <stdlib.h>
#include <string.h>

ARBINT_TEST_DECLARE_FAILURES();

/*  Helper to free string allocated by arbint context.
    Must use the same allocator that arbint_get_str used.  */
static void free_str(arbint_ctx_t * ctx, char * s) {
  if (s == NULL)
    return;
  if (ctx != NULL && ctx->a.realloc != NULL)
    ctx->a.realloc(ctx->a.ud, s, 0u);
  else
    free(s);
}

/*  Test zero handling in all bases.  */
static void test_zero(arbint_ctx_t * ctx) {
  arbint_t x;
  char * str;
  int base;

  CHECK_EQ_I(arbint_init(x, ctx), ARBINT_OK);

  for (base = 2; base <= 36; ++base) {
    /*  get_str of zero should return "0".  */
    arbint_zero(x);
    str = NULL;
    CHECK_EQ_I(arbint_get_str(x, &str, base), ARBINT_OK);
    CHECK(str != NULL);
    CHECK(strcmp(str, "0") == 0);
    free_str(ctx, str);

    /*  set_str of "0" should give zero.  */
    CHECK_EQ_I(arbint_set_str(x, "0", base), ARBINT_OK);
    CHECK(arbint_is_zero(x));

    /*  set_str of empty string (after leading zeros stripped) should give zero.  */
    CHECK_EQ_I(arbint_set_str(x, "0000", base), ARBINT_OK);
    CHECK(arbint_is_zero(x));
  }

  arbint_clear(x);
}

/*  Test small values (single limb) in various bases.  */
static void test_small_values(arbint_ctx_t * ctx) {
  arbint_t x, y;
  char * str;
  int32_t test_vals[] = {1, -1, 42, -42, 127, -128, 255, 256, 1000, -1000,
                         12345, -12345, INT32_MAX, INT32_MIN};
  size_t num_vals = sizeof(test_vals) / sizeof(test_vals[0]);
  size_t i;
  int base;

  CHECK_EQ_I(arbint_init(x, ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(y, ctx), ARBINT_OK);

  for (i = 0u; i < num_vals; ++i) {
    CHECK_EQ_I(arbint_set_i32(x, test_vals[i]), ARBINT_OK);

    for (base = 2; base <= 36; ++base) {
      /*  Round-trip: get_str then set_str should preserve value.  */
      str = NULL;
      CHECK_EQ_I(arbint_get_str(x, &str, base), ARBINT_OK);
      CHECK(str != NULL);

      CHECK_EQ_I(arbint_set_str(y, str, base), ARBINT_OK);
      CHECK(arbint_eq(x, y));

      free_str(ctx, str);
    }
  }

  arbint_clear(y);
  arbint_clear(x);
}

/*  Test power-of-2 bases with known bit patterns.  */
static void test_pow2_bases(arbint_ctx_t * ctx) {
  arbint_t x, y;
  char * str;

  CHECK_EQ_I(arbint_init(x, ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(y, ctx), ARBINT_OK);

  /*  Binary (base 2).  */
  CHECK_EQ_I(arbint_set_str(x, "11111111", 2), ARBINT_OK);
  check_u32_value(x, 255u);

  CHECK_EQ_I(arbint_set_str(x, "100000000", 2), ARBINT_OK);
  check_u32_value(x, 256u);

  CHECK_EQ_I(arbint_set_i32(x, 255), ARBINT_OK);
  str = NULL;
  CHECK_EQ_I(arbint_get_str(x, &str, 2), ARBINT_OK);
  CHECK(strcmp(str, "11111111") == 0);
  free_str(ctx, str);

  /*  Octal (base 8).  */
  CHECK_EQ_I(arbint_set_str(x, "377", 8), ARBINT_OK);
  check_u32_value(x, 255u);

  CHECK_EQ_I(arbint_set_str(x, "1000", 8), ARBINT_OK);
  check_u32_value(x, 512u);

  CHECK_EQ_I(arbint_set_i32(x, 511), ARBINT_OK);
  str = NULL;
  CHECK_EQ_I(arbint_get_str(x, &str, 8), ARBINT_OK);
  CHECK(strcmp(str, "777") == 0);
  free_str(ctx, str);

  /*  Hexadecimal (base 16).  */
  CHECK_EQ_I(arbint_set_str(x, "ff", 16), ARBINT_OK);
  check_u32_value(x, 255u);

  CHECK_EQ_I(arbint_set_str(x, "FF", 16), ARBINT_OK);  /*  Uppercase input.  */
  check_u32_value(x, 255u);

  CHECK_EQ_I(arbint_set_str(x, "DeAdBeEf", 16), ARBINT_OK);
  check_u32_value(x, 0xDEADBEEFu);

  CHECK_EQ_I(arbint_set_i32(x, (int32_t) 0x7FFFFFFF), ARBINT_OK);
  str = NULL;
  CHECK_EQ_I(arbint_get_str(x, &str, 16), ARBINT_OK);
  CHECK(strcmp(str, "7fffffff") == 0);
  free_str(ctx, str);

  /*  Base 32.  */
  CHECK_EQ_I(arbint_set_str(x, "10", 32), ARBINT_OK);
  check_u32_value(x, 32u);

  CHECK_EQ_I(arbint_set_str(x, "vvvvvv", 32), ARBINT_OK);
  /*  v = 31, so this is 31 * (32^5 + 32^4 + 32^3 + 32^2 + 32 + 1) = 1073741823.  */
  check_u32_value(x, 1073741823u);

  arbint_clear(y);
  arbint_clear(x);
}

/*  Test base 10 specifically (most common use case).  */
static void test_base10(arbint_ctx_t * ctx) {
  arbint_t x, y;
  char * str;

  CHECK_EQ_I(arbint_init(x, ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(y, ctx), ARBINT_OK);

  /*  Test known values.  */
  CHECK_EQ_I(arbint_set_str(x, "12345", 10), ARBINT_OK);
  check_i32_value(x, 12345);

  CHECK_EQ_I(arbint_set_str(x, "-67890", 10), ARBINT_OK);
  check_i32_value(x, -67890);

  /*  Test output formatting.  */
  CHECK_EQ_I(arbint_set_i32(x, 1000000), ARBINT_OK);
  str = NULL;
  CHECK_EQ_I(arbint_get_str(x, &str, 10), ARBINT_OK);
  CHECK(strcmp(str, "1000000") == 0);
  free_str(ctx, str);

  CHECK_EQ_I(arbint_set_i32(x, -999999), ARBINT_OK);
  str = NULL;
  CHECK_EQ_I(arbint_get_str(x, &str, 10), ARBINT_OK);
  CHECK(strcmp(str, "-999999") == 0);
  free_str(ctx, str);

  /*  Test round-trip for INT32_MAX and INT32_MIN.  */
  CHECK_EQ_I(arbint_set_str(x, "2147483647", 10), ARBINT_OK);
  check_i32_value(x, INT32_MAX);

  CHECK_EQ_I(arbint_set_str(x, "-2147483648", 10), ARBINT_OK);
  check_i32_value(x, INT32_MIN);

  arbint_clear(y);
  arbint_clear(x);
}

/*  Test edge cases: whitespace, signs, leading zeros.  */
static void test_edge_cases(arbint_ctx_t * ctx) {
  arbint_t x;

  CHECK_EQ_I(arbint_init(x, ctx), ARBINT_OK);

  /*  Leading whitespace should be skipped.  */
  CHECK_EQ_I(arbint_set_str(x, "   123", 10), ARBINT_OK);
  check_i32_value(x, 123);

  CHECK_EQ_I(arbint_set_str(x, "\t\n  456", 10), ARBINT_OK);
  check_i32_value(x, 456);

  /*  Explicit positive sign.  */
  CHECK_EQ_I(arbint_set_str(x, "+789", 10), ARBINT_OK);
  check_i32_value(x, 789);

  /*  Leading zeros should be stripped.  */
  CHECK_EQ_I(arbint_set_str(x, "00000123", 10), ARBINT_OK);
  check_i32_value(x, 123);

  CHECK_EQ_I(arbint_set_str(x, "-00042", 10), ARBINT_OK);
  check_i32_value(x, -42);

  /*  All zeros -> zero.  */
  CHECK_EQ_I(arbint_set_str(x, "0000000", 10), ARBINT_OK);
  CHECK(arbint_is_zero(x));

  CHECK_EQ_I(arbint_set_str(x, "-0", 10), ARBINT_OK);
  CHECK(arbint_is_zero(x));

  arbint_clear(x);
}

/*  Test invalid input handling.  */
static void test_invalid_input(arbint_ctx_t * ctx) {
  arbint_t x;

  CHECK_EQ_I(arbint_init(x, ctx), ARBINT_OK);

  /*  Invalid characters for base.  */
  CHECK_EQ_I(arbint_set_str(x, "123abc", 10), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_set_str(x, "1g", 16), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_set_str(x, "9", 8), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_set_str(x, "2", 2), ARBINT_EINVAL);

  /*  NULL pointer handling.  */
  CHECK_EQ_I(arbint_set_str(NULL, "123", 10), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_set_str(x, NULL, 10), ARBINT_EINVAL);

  /*  Invalid base.  */
  CHECK_EQ_I(arbint_set_str(x, "123", 1), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_set_str(x, "123", 37), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_set_str(x, "123", 0), ARBINT_EINVAL);

  /*  Empty strings and sign-only inputs should be rejected.  */
  CHECK_EQ_I(arbint_set_str(x, "", 10), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_set_str(x, "-", 10), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_set_str(x, "+", 10), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_set_str(x, "   ", 10), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_set_str(x, "  -", 10), ARBINT_EINVAL);
  CHECK_EQ_I(arbint_set_str(x, "  +", 10), ARBINT_EINVAL);

  /*  But "0", "-0", "+0" should be accepted as zero.  */
  CHECK_EQ_I(arbint_set_str(x, "0", 10), ARBINT_OK);
  CHECK(arbint_is_zero(x));
  CHECK_EQ_I(arbint_set_str(x, "-0", 10), ARBINT_OK);
  CHECK(arbint_is_zero(x));
  CHECK_EQ_I(arbint_set_str(x, "+0", 10), ARBINT_OK);
  CHECK(arbint_is_zero(x));

  /*  get_str invalid inputs.  */
  {
    char * str = NULL;
    CHECK_EQ_I(arbint_get_str(NULL, &str, 10), ARBINT_EINVAL);
    CHECK_EQ_I(arbint_get_str(x, NULL, 10), ARBINT_EINVAL);
    CHECK_EQ_I(arbint_get_str(x, &str, 1), ARBINT_EINVAL);
    CHECK_EQ_I(arbint_get_str(x, &str, 37), ARBINT_EINVAL);
  }

  arbint_clear(x);
}

/*  Test large numbers (multi-limb).  */
static void test_large_numbers(arbint_ctx_t * ctx) {
  arbint_t x, y, z;
  char * str;

  CHECK_EQ_I(arbint_init(x, ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(y, ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(z, ctx), ARBINT_OK);

  /*  10^30 - a number requiring multiple limbs.  */
  CHECK_EQ_I(arbint_set_str(x, "1000000000000000000000000000000", 10), ARBINT_OK);
  str = NULL;
  CHECK_EQ_I(arbint_get_str(x, &str, 10), ARBINT_OK);
  CHECK(strcmp(str, "1000000000000000000000000000000") == 0);
  free_str(ctx, str);

  /*  Negative large number.  */
  CHECK_EQ_I(arbint_set_str(x, "-12345678901234567890", 10), ARBINT_OK);
  str = NULL;
  CHECK_EQ_I(arbint_get_str(x, &str, 10), ARBINT_OK);
  CHECK(strcmp(str, "-12345678901234567890") == 0);
  free_str(ctx, str);

  /*  Build a large number by repeated squaring, then round-trip.  */
  CHECK_EQ_I(arbint_set_i32(x, 2), ARBINT_OK);
  {
    int i;
    for (i = 0; i < 20; ++i) {
      CHECK_EQ_I(arbint_mul(x, x, x), ARBINT_OK);
    }
  }

  /*  x is now 2^(2^20) = 2^1048576, a very large number.
      Actually that's too large, let's use a smaller power.  */
  CHECK_EQ_I(arbint_set_i32(x, 2), ARBINT_OK);
  {
    int i;
    for (i = 0; i < 8; ++i) {
      CHECK_EQ_I(arbint_mul(x, x, x), ARBINT_OK);
    }
  }
  /*  x is now 2^256.  */

  /*  Round-trip in decimal.  */
  str = NULL;
  CHECK_EQ_I(arbint_get_str(x, &str, 10), ARBINT_OK);
  CHECK(str != NULL);
  CHECK_EQ_I(arbint_set_str(y, str, 10), ARBINT_OK);
  CHECK(arbint_eq(x, y));
  free_str(ctx, str);

  /*  Round-trip in hex.  */
  str = NULL;
  CHECK_EQ_I(arbint_get_str(x, &str, 16), ARBINT_OK);
  CHECK(str != NULL);
  /*  2^256 in hex is 1 followed by 64 zeros.  */
  CHECK(strlen(str) == 65u);
  CHECK(str[0] == '1');
  CHECK_EQ_I(arbint_set_str(y, str, 16), ARBINT_OK);
  CHECK(arbint_eq(x, y));
  free_str(ctx, str);

  /*  Round-trip in binary.  */
  str = NULL;
  CHECK_EQ_I(arbint_get_str(x, &str, 2), ARBINT_OK);
  CHECK(str != NULL);
  CHECK(strlen(str) == 257u);  /*  2^256 has 257 bits.  */
  CHECK_EQ_I(arbint_set_str(y, str, 2), ARBINT_OK);
  CHECK(arbint_eq(x, y));
  free_str(ctx, str);

  arbint_clear(z);
  arbint_clear(y);
  arbint_clear(x);
}

/*  Test all bases 2-36 with round-trip for a few values.  */
static void test_all_bases_roundtrip(arbint_ctx_t * ctx) {
  arbint_t x, y;
  char * str;
  int base;

  CHECK_EQ_I(arbint_init(x, ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(y, ctx), ARBINT_OK);

  /*  Test with a multi-limb value.  */
  CHECK_EQ_I(arbint_set_str(x, "1234567890123456789012345", 10), ARBINT_OK);

  for (base = 2; base <= 36; ++base) {
    str = NULL;
    CHECK_EQ_I(arbint_get_str(x, &str, base), ARBINT_OK);
    CHECK(str != NULL);
    CHECK(strlen(str) > 0u);

    CHECK_EQ_I(arbint_set_str(y, str, base), ARBINT_OK);
    CHECK(arbint_eq(x, y));

    free_str(ctx, str);
  }

  /*  Test with negative value.  */
  CHECK_EQ_I(arbint_set_str(x, "-9876543210987654321", 10), ARBINT_OK);

  for (base = 2; base <= 36; ++base) {
    str = NULL;
    CHECK_EQ_I(arbint_get_str(x, &str, base), ARBINT_OK);
    CHECK(str != NULL);
    CHECK(str[0] == '-');

    CHECK_EQ_I(arbint_set_str(y, str, base), ARBINT_OK);
    CHECK(arbint_eq(x, y));

    free_str(ctx, str);
  }

  arbint_clear(y);
  arbint_clear(x);
}

/*  Stochastic testing with pseudo-random values.  */
static void test_stochastic(arbint_ctx_t * ctx) {
  arbint_t x, y;
  char * str;
  int i;
  int base;
  uint32_t seed = 12345u;

  CHECK_EQ_I(arbint_init(x, ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(y, ctx), ARBINT_OK);

  /*  Generate pseudo-random values via deterministic computation.  */
  for (i = 0; i < 50; ++i) {
    /*  Build a multi-limb number via repeated operations.  */
    seed = seed * 1103515245u + 12345u;
    CHECK_EQ_I(arbint_set_u32(x, seed), ARBINT_OK);

    /*  Make it larger by multiplying and adding.  */
    {
      int j;
      int growth = (i % 6) + 1;
      for (j = 0; j < growth; ++j) {
        seed = seed * 1103515245u + 12345u;
        CHECK_EQ_I(arbint_mul_u32(x, x, (seed >> 16) | 1u), ARBINT_OK);
        seed = seed * 1103515245u + 12345u;
        CHECK_EQ_I(arbint_add_u32(x, x, seed >> 8), ARBINT_OK);
      }
    }

    /*  Make some negative.  */
    if (i % 3 == 0 && x[0]._sz > 0)
      x[0]._sz = -x[0]._sz;

    /*  Test in a few different bases.  */
    for (base = 2; base <= 36; base += 7) {
      str = NULL;
      CHECK_EQ_I(arbint_get_str(x, &str, base), ARBINT_OK);
      CHECK(str != NULL);

      CHECK_EQ_I(arbint_set_str(y, str, base), ARBINT_OK);
      if (!arbint_eq(x, y)) {
        fprintf(stderr, "FAIL: round-trip failed for base %d, i=%d\n", base, i);
        ++g_failures;
      }

      free_str(ctx, str);
    }
  }

  arbint_clear(y);
  arbint_clear(x);
}

/*  Test divide-and-conquer threshold behavior.  */
static void test_dc_threshold(arbint_ctx_t * ctx) {
  arbint_t x, y;
  char * str;
  int i;

  CHECK_EQ_I(arbint_init(x, ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(y, ctx), ARBINT_OK);

  /*  Build a number that should trigger D&C for output.
      D&C threshold is ~32 limbs, which is about 600+ decimal digits.  */
  CHECK_EQ_I(arbint_set_i32(x, 2), ARBINT_OK);
  for (i = 0; i < 10; ++i) {
    CHECK_EQ_I(arbint_mul(x, x, x), ARBINT_OK);
  }
  /*  x = 2^1024, about 309 decimal digits.  */

  /*  Go a bit larger to ensure D&C is triggered.  */
  CHECK_EQ_I(arbint_mul(x, x, x), ARBINT_OK);
  /*  x = 2^2048, about 617 decimal digits.  */

  /*  Round-trip in decimal.  */
  str = NULL;
  CHECK_EQ_I(arbint_get_str(x, &str, 10), ARBINT_OK);
  CHECK(str != NULL);
  CHECK(strlen(str) > 600u);

  CHECK_EQ_I(arbint_set_str(y, str, 10), ARBINT_OK);
  CHECK(arbint_eq(x, y));

  free_str(ctx, str);

  /*  Round-trip in hex (should use power-of-2 fast path).  */
  str = NULL;
  CHECK_EQ_I(arbint_get_str(x, &str, 16), ARBINT_OK);
  CHECK(str != NULL);
  CHECK(strlen(str) == 513u);  /*  2^2048 in hex: 512 zeros + leading 1.  */

  CHECK_EQ_I(arbint_set_str(y, str, 16), ARBINT_OK);
  CHECK(arbint_eq(x, y));

  free_str(ctx, str);

  arbint_clear(y);
  arbint_clear(x);
}

int main(void) {
  arbint_ctx_t ctx;

  ARBINT_TEST_START();

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);

  test_zero(&ctx);
  test_small_values(&ctx);
  test_pow2_bases(&ctx);
  test_base10(&ctx);
  test_edge_cases(&ctx);
  test_invalid_input(&ctx);
  test_large_numbers(&ctx);
  test_all_bases_roundtrip(&ctx);
  test_stochastic(&ctx);
  test_dc_threshold(&ctx);

  arbint_ctx_clear(&ctx);

  ARBINT_TEST_FINISH("test_str");
}
