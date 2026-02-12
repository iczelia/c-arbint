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

#include "arbint_cache.h"

#include "arbint_ntt.h"

#include <stdlib.h>
#include <string.h>

/*  Cached prime sieve state.
    Initialized on first call needing primes > 521.
    Grows as needed; never shrinks until cache drop.  */
static uint8_t * g_prime_sieve = NULL;
static uint32_t g_prime_sieve_limit = 0;

arbint_err_t arbint_cache_prime_sieve_ensure(uint32_t limit) {
  uint8_t * new_sieve;
  size_t new_bytes;
  uint32_t p;
  uint32_t i;
  uint32_t sqrt_limit;

  if (g_prime_sieve_limit >= limit)
    return ARBINT_OK;

  new_bytes = ((size_t) limit + 1u + 7u) / 8u;
  new_sieve = (uint8_t *) realloc(g_prime_sieve, new_bytes);
  if (new_sieve == NULL)
    return ARBINT_ENOMEM;

  if (g_prime_sieve_limit > 0) {
    size_t old_bytes = ((size_t) g_prime_sieve_limit + 1u + 7u) / 8u;
    memset(new_sieve + old_bytes, 0, new_bytes - old_bytes);
  } else {
    memset(new_sieve, 0, new_bytes);
  }

  new_sieve[0] |= 0x03;

  sqrt_limit = 1u;
  while (sqrt_limit <= limit / sqrt_limit)
    ++sqrt_limit;

  for (p = 2u; p < sqrt_limit; ++p) {
    if ((new_sieve[p / 8u] & (1u << (p % 8u))) != 0)
      continue;
    for (i = p * p; i <= limit; i += p)
      new_sieve[i / 8u] |= (uint8_t) (1u << (i % 8u));
  }

  g_prime_sieve = new_sieve;
  g_prime_sieve_limit = limit;
  return ARBINT_OK;
}

int arbint_cache_prime_sieve_is_prime(uint32_t n) {
  if (n > g_prime_sieve_limit)
    return 0;
  return (g_prime_sieve[n / 8u] & (1u << (n % 8u))) == 0;
}

void arbint_cache_drop_all(void) {
  free(g_prime_sieve);
  g_prime_sieve = NULL;
  g_prime_sieve_limit = 0;

  arbint_ntt_cache_clear_generic();
#if HAS_BMI2
  arbint_ntt_cache_clear_bmi2();
#endif
#if HAS_AVX2
  arbint_ntt_cache_clear_avx2();
#endif
}

void arbint_drop_caches(void) { arbint_cache_drop_all(); }
