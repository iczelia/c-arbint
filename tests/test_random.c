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

/*  Deterministic seed for reproducible tests  */
static const uint32_t test_seed[4] = {0x123, 0x234, 0x345, 0x456};

int main(void) {
  ARBINT_TEST_DECLARE_FAILURES();

  arbint_ctx_t ctx;
  arbint_t a, b, bound;
  arbint_rng_t rng, rng2;

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(bound, &ctx), ARBINT_OK);

  /* Test 1: Manual seeding */
  CHECK_EQ_I(arbint_rng_init(&rng, test_seed, sizeof(test_seed)), ARBINT_OK);

  /* Test 2: arbint_urandomb(k=0) -> 0 */
  CHECK_EQ_I(arbint_urandomb(a, &rng, 0), ARBINT_OK);
  check_i32_value(a, 0);

  /* Test 3: arbint_urandomb(k=1) -> 0 or 1 */
  CHECK_EQ_I(arbint_urandomb(a, &rng, 1), ARBINT_OK);
  CHECK(arbint_nbits(a) <= 1);

  /* Test 4: arbint_urandomb(k=32) fills correctly */
  CHECK_EQ_I(arbint_urandomb(a, &rng, 32), ARBINT_OK);
  CHECK(arbint_nbits(a) <= 32);

  /* Test 5: arbint_urandomb(k=64) fills 1 limb */
  CHECK_EQ_I(arbint_urandomb(a, &rng, 64), ARBINT_OK);
  CHECK(arbint_nbits(a) <= 64);

  /* Test 6: arbint_urandomb(k=100) with bit masking */
  CHECK_EQ_I(arbint_urandomb(a, &rng, 100), ARBINT_OK);
  CHECK(arbint_nbits(a) <= 100);

  /* Test 7: arbint_urandomb(k=128) fills 2 limbs */
  CHECK_EQ_I(arbint_urandomb(a, &rng, 128), ARBINT_OK);
  CHECK(arbint_nbits(a) <= 128);

  /* Test 8: arbint_urandomm(bound=1) -> always 0 */
  CHECK_EQ_I(arbint_set_u32(bound, 1), ARBINT_OK);
  CHECK_EQ_I(arbint_urandomm(a, &rng, bound), ARBINT_OK);
  check_i32_value(a, 0);

  /* Test 9: arbint_urandomm(bound=2) -> 0 or 1 */
  CHECK_EQ_I(arbint_set_u32(bound, 2), ARBINT_OK);
  CHECK_EQ_I(arbint_urandomm(a, &rng, bound), ARBINT_OK);
  CHECK(arbint_cmp_u32(a, 0) >= 0);
  CHECK(arbint_cmp_u32(a, 2) < 0);

  /* Test 10: arbint_urandomm(bound=10) -> [0,9] */
  CHECK_EQ_I(arbint_set_u32(bound, 10), ARBINT_OK);
  CHECK_EQ_I(arbint_urandomm(a, &rng, bound), ARBINT_OK);
  CHECK(arbint_cmp_u32(a, 0) >= 0);
  CHECK(arbint_cmp_u32(a, 10) < 0);

  /* Test 11: arbint_urandomm(bound=100) */
  CHECK_EQ_I(arbint_set_u32(bound, 100), ARBINT_OK);
  CHECK_EQ_I(arbint_urandomm(a, &rng, bound), ARBINT_OK);
  CHECK(arbint_cmp_u32(a, 0) >= 0);
  CHECK(arbint_cmp_u32(a, 100) < 0);

  /* Test 12: Error handling - bound == 0 */
  CHECK_EQ_I(arbint_set_u32(bound, 0), ARBINT_OK);
  CHECK_EQ_I(arbint_urandomm(a, &rng, bound), ARBINT_EDOM);

  /* Test 13: Error handling - bound < 0 */
  CHECK_EQ_I(arbint_set_i32(bound, -5), ARBINT_OK);
  CHECK_EQ_I(arbint_urandomm(a, &rng, bound), ARBINT_EDOM);

  /* Test 14: Error handling - NULL rng */
  CHECK_EQ_I(arbint_urandomb(a, NULL, 64), ARBINT_EINVAL);

  /* Test 15: Error handling - NULL rop */
  CHECK_EQ_I(arbint_urandomb(NULL, &rng, 64), ARBINT_EINVAL);

  /* Test 16: Error handling - invalid seed_len */
  CHECK_EQ_I(arbint_rng_init(&rng, test_seed, 3),
             ARBINT_EINVAL); /* Not multiple of 4 */

  /* Test 17: Platform entropy seeding (seed=NULL) */
  arbint_rng_clear(&rng);
  /* May fail if no platform entropy source, so don't CHECK */
  (void) arbint_rng_init(&rng, NULL, 0);

  /* Test 18: Verify different values from different seeds */
  uint32_t seed2[4] = {0x999, 0x888, 0x777, 0x666};

  CHECK_EQ_I(arbint_rng_init(&rng, test_seed, sizeof(test_seed)), ARBINT_OK);
  CHECK_EQ_I(arbint_rng_init(&rng2, seed2, sizeof(seed2)), ARBINT_OK);

  CHECK_EQ_I(arbint_urandomb(a, &rng, 128), ARBINT_OK);
  CHECK_EQ_I(arbint_urandomb(b, &rng2, 128), ARBINT_OK);
  CHECK(arbint_cmp(a, b) != 0); /* Different seeds -> different values */

  /* Test 19: Verify deterministic sequence (same seed -> same sequence) */
  arbint_rng_clear(&rng);
  arbint_rng_clear(&rng2);
  CHECK_EQ_I(arbint_rng_init(&rng, test_seed, sizeof(test_seed)), ARBINT_OK);
  CHECK_EQ_I(arbint_rng_init(&rng2, test_seed, sizeof(test_seed)), ARBINT_OK);

  CHECK_EQ_I(arbint_urandomb(a, &rng, 256), ARBINT_OK);
  CHECK_EQ_I(arbint_urandomb(b, &rng2, 256), ARBINT_OK);
  CHECK(arbint_cmp(a, b) == 0); /* Same seed -> same value */

  /* Test 20: Large bound (via repeated multiplication) */
  CHECK_EQ_I(arbint_set_u32(bound, 1000000), ARBINT_OK);
  CHECK_EQ_I(arbint_urandomm(a, &rng, bound), ARBINT_OK);
  CHECK(arbint_cmp(a, bound) < 0);

  /* Test 21: Aliasing - rop == bound */
  CHECK_EQ_I(arbint_set_u32(bound, 1000), ARBINT_OK);
  CHECK_EQ_I(arbint_urandomm(bound, &rng, bound), ARBINT_OK);
  CHECK(arbint_cmp_u32(bound, 1000) < 0);

  /* Test 22: Multiple sequential generations (check RNG state progression) */
  {
    int i;
    int all_same = 1;

    CHECK_EQ_I(arbint_set_u32(bound, 100), ARBINT_OK);
    CHECK_EQ_I(arbint_urandomm(a, &rng, bound), ARBINT_OK);

    for (i = 0; i < 10; i++) {
      CHECK_EQ_I(arbint_urandomm(b, &rng, bound), ARBINT_OK);
      if (arbint_cmp(a, b) != 0) {
        all_same = 0;
        break;
      }
    }
    /* Should produce different values (probability of all same ~= 0) */
    CHECK(all_same == 0);
  }

  /* Test 24: Distribution uniformity for small bounds (chi-squared) */
  {
    int i, bucket;
    int counts[10] = {0};
    double chi_squared, expected;

    CHECK_EQ_I(arbint_set_u32(bound, 10), ARBINT_OK);

    /* Generate 10000 samples */
    for (i = 0; i < 10000; i++) {
      int32_t val;
      CHECK_EQ_I(arbint_urandomm(a, &rng, bound), ARBINT_OK);
      CHECK_EQ_I(arbint_get_i32(a, &val), ARBINT_OK);
      CHECK(val >= 0 && val < 10);
      counts[val]++;
    }

    /* Chi-squared test: each bucket should have ~1000 samples */
    /* We expect roughly uniform distribution */
    expected = 1000.0;
    chi_squared = 0.0;
    for (bucket = 0; bucket < 10; bucket++) {
      double deviation = counts[bucket] - expected;
      chi_squared += (deviation * deviation) / expected;
    }

    /* For 9 degrees of freedom at 95% confidence: chi_squared < 16.919 */
    /* We use a more lenient threshold of 25.0 to avoid flaky tests */
    CHECK(chi_squared < 25.0);
  }

  /*  Test 25: Boundary testing - bounds near powers of 2
      (worst case rejection)  */
  {
    /* bound = 2^k - 1: best case, no rejection */
    CHECK_EQ_I(arbint_set_u32(bound, 255), ARBINT_OK);
    CHECK_EQ_I(arbint_urandomm(a, &rng, bound), ARBINT_OK);
    CHECK(arbint_cmp_u32(a, 0) >= 0);
    CHECK(arbint_cmp(a, bound) < 0);

    /* bound = 2^k + 1: worst case, ~50% rejection rate */
    CHECK_EQ_I(arbint_set_u32(bound, 257), ARBINT_OK);
    CHECK_EQ_I(arbint_urandomm(a, &rng, bound), ARBINT_OK);
    CHECK(arbint_cmp_u32(a, 0) >= 0);
    CHECK(arbint_cmp(a, bound) < 0);

    /* bound = 2^k: exact power of 2 */
    CHECK_EQ_I(arbint_set_u32(bound, 256), ARBINT_OK);
    CHECK_EQ_I(arbint_urandomm(a, &rng, bound), ARBINT_OK);
    CHECK(arbint_cmp_u32(a, 0) >= 0);
    CHECK(arbint_cmp(a, bound) < 0);
  }

  /* Test 26: Large bit counts for arbint_urandomb */
  {
    size_t bit_counts[] = {1,   7,   8,   9,   15,  16,  17,  31,
                           32,  33,  63,  64,  65,  127, 128, 129,
                           255, 256, 257, 511, 512, 513, 1024};
    size_t i;

    for (i = 0; i < sizeof(bit_counts) / sizeof(bit_counts[0]); i++) {
      CHECK_EQ_I(arbint_urandomb(a, &rng, bit_counts[i]), ARBINT_OK);
      CHECK(arbint_nbits(a) <= bit_counts[i]);
    }
  }

  /*  Test 27: Verify rejection sampling correctness,
      bound just above power of 2.  */
  {
    int i, found_zero, found_max;
    int32_t val;

    /* bound = 3: values should be 0, 1, 2 (not 0, 1, 2, 3) */
    CHECK_EQ_I(arbint_set_u32(bound, 3), ARBINT_OK);

    found_zero = 0;
    found_max = 0;

    /* Generate many samples, verify range */
    for (i = 0; i < 1000; i++) {
      CHECK_EQ_I(arbint_urandomm(a, &rng, bound), ARBINT_OK);
      CHECK_EQ_I(arbint_get_i32(a, &val), ARBINT_OK);
      CHECK(val >= 0 && val < 3);
      if (val == 0)
        found_zero = 1;
      if (val == 2)
        found_max = 1;
    }

    /* Should have seen both extremes */
    CHECK(found_zero);
    CHECK(found_max);
  }

  /* Test 28: Large multi-limb bounds */
  {
    /* Build large bound via squaring */
    CHECK_EQ_I(arbint_set_u32(a, 1000000), ARBINT_OK);
    CHECK_EQ_I(arbint_sqr(bound, a), ARBINT_OK);

    CHECK_EQ_I(arbint_urandomm(a, &rng, bound), ARBINT_OK);
    CHECK(arbint_cmp(a, bound) < 0);
    CHECK(arbint_signum(a) >= 0);
  }

  /* Test 29: Verify no modulo bias with statistical test */
  {
    int i;
    int counts[5] = {0};
    double chi_squared, expected;

    /* bound = 5: sample 5000 times */
    CHECK_EQ_I(arbint_set_u32(bound, 5), ARBINT_OK);

    for (i = 0; i < 5000; i++) {
      int32_t val;
      CHECK_EQ_I(arbint_urandomm(a, &rng, bound), ARBINT_OK);
      CHECK_EQ_I(arbint_get_i32(a, &val), ARBINT_OK);
      CHECK(val >= 0 && val < 5);
      counts[val]++;
    }

    /* Chi-squared test: each bucket should have ~1000 samples */
    expected = 1000.0;
    chi_squared = 0.0;
    for (i = 0; i < 5; i++) {
      double deviation = counts[i] - expected;
      chi_squared += (deviation * deviation) / expected;
    }

    /* For 4 degrees of freedom at 95% confidence: chi_squared < 9.488 */
    /* Use lenient threshold of 15.0 */
    CHECK(chi_squared < 15.0);
  }

  /* Test 30: arbint_urandomb edge case - very small bit counts */
  {
    CHECK_EQ_I(arbint_urandomb(a, &rng, 1), ARBINT_OK);
    CHECK(arbint_cmp_u32(a, 0) >= 0 && arbint_cmp_u32(a, 2) < 0);

    CHECK_EQ_I(arbint_urandomb(a, &rng, 2), ARBINT_OK);
    CHECK(arbint_cmp_u32(a, 0) >= 0 && arbint_cmp_u32(a, 4) < 0);
  }

  /* Test 31: Multiple RNG instances produce different sequences */
  {
    arbint_rng_t rng3;
    uint32_t seed3[4] = {0xaaa, 0xbbb, 0xccc, 0xddd};

    CHECK_EQ_I(arbint_rng_init(&rng, test_seed, sizeof(test_seed)), ARBINT_OK);
    CHECK_EQ_I(arbint_rng_init(&rng3, seed3, sizeof(seed3)), ARBINT_OK);

    CHECK_EQ_I(arbint_urandomb(a, &rng, 128), ARBINT_OK);
    CHECK_EQ_I(arbint_urandomb(b, &rng3, 128), ARBINT_OK);

    /* Different seeds should produce different values */
    CHECK(arbint_cmp(a, b) != 0);

    arbint_rng_clear(&rng3);
  }

  /* Test 23: arbint_rng_clear zeros memory */
  CHECK_EQ_I(arbint_rng_init(&rng, test_seed, sizeof(test_seed)), ARBINT_OK);
  arbint_rng_clear(&rng);
  /*  After clear, state should be zero (can't easily verify, but ensure no
      crash)  */
  CHECK_EQ_I(arbint_rng_init(&rng, test_seed, sizeof(test_seed)), ARBINT_OK);

  /* Cleanup */
  arbint_rng_clear(&rng);
  arbint_rng_clear(&rng2);
  arbint_clear(bound);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);

  ARBINT_TEST_FINISH("test_random");
}
