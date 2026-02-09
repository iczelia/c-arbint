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

#include "arbint_bitops.h"

#include "arbint.h"

/*  Count number of significant bits in x.
    Returns 0 for x == 0, otherwise position of highest set bit + 1.  */
size_t arbint_nbits(const arbint_t x) {
  size_t used;
  arbint_limb_t top_limb;
  unsigned clz_count;

  if (x == NULL)
    return 0u;

  used = arbint_abs_sz(x[0]._sz);
  if (used == 0u)
    return 0u;

  top_limb = ARBINT_CLIMBS(x)[used - 1u];
  if (top_limb == 0u)
    return 0u; /* Should never happen if normalized */

  clz_count = arbint_clz_limb(top_limb);
  return (used - 1u) * ARBINT_LIMB_BITS + (ARBINT_LIMB_BITS - clz_count);
}

/*  Fixed-point ceil(2^16 / log2(base)) for bases 3..62.
    Entry i corresponds to base (i + 3).
    Used by arbint_sizeinbase to compute ceil(nbits / log2(base))
    without floating-point dependencies.  */
static const unsigned arbint_recip_log2_16[] = {
    /* base  3 */ 41349u, /* base  4 */ 32768u,
    /* base  5 */ 28225u, /* base  6 */ 25353u,
    /* base  7 */ 23345u, /* base  8 */ 21846u,
    /* base  9 */ 20675u, /* base 10 */ 19729u,
    /* base 11 */ 18945u, /* base 12 */ 18281u,
    /* base 13 */ 17711u, /* base 14 */ 17213u,
    /* base 15 */ 16775u, /* base 16 */ 16384u,
    /* base 17 */ 16034u, /* base 18 */ 15717u,
    /* base 19 */ 15428u, /* base 20 */ 15164u,
    /* base 21 */ 14921u, /* base 22 */ 14697u,
    /* base 23 */ 14488u, /* base 24 */ 14294u,
    /* base 25 */ 14113u, /* base 26 */ 13943u,
    /* base 27 */ 13783u, /* base 28 */ 13633u,
    /* base 29 */ 13491u, /* base 30 */ 13356u,
    /* base 31 */ 13229u, /* base 32 */ 13108u,
    /* base 33 */ 12992u, /* base 34 */ 12882u,
    /* base 35 */ 12777u, /* base 36 */ 12677u,
    /* base 37 */ 12581u, /* base 38 */ 12488u,
    /* base 39 */ 12400u, /* base 40 */ 12315u,
    /* base 41 */ 12233u, /* base 42 */ 12154u,
    /* base 43 */ 12078u, /* base 44 */ 12005u,
    /* base 45 */ 11934u, /* base 46 */ 11865u,
    /* base 47 */ 11799u, /* base 48 */ 11735u,
    /* base 49 */ 11673u, /* base 50 */ 11612u,
    /* base 51 */ 11554u, /* base 52 */ 11497u,
    /* base 53 */ 11442u, /* base 54 */ 11388u,
    /* base 55 */ 11336u, /* base 56 */ 11285u,
    /* base 57 */ 11236u, /* base 58 */ 11188u,
    /* base 59 */ 11141u, /* base 60 */ 11095u,
    /* base 61 */ 11051u, /* base 62 */ 11007u};

/*  Number of digits in the given base needed to represent |x|.
    Returns 1 for x == 0.  base must be in [2, 62].
    For base == 2 this is exactly arbint_nbits(x).
    For power-of-two bases, an exact bit-count formula is used.
    For other bases, the result may overestimate by at most 2.  */
size_t arbint_sizeinbase(const arbint_t x, int base) {
  size_t nb;
  if (x == NULL || base < 2 || base > 62)
    return 1u;
  nb = arbint_nbits(x);
  if (nb == 0u)
    return 1u;
  if (base == 2)
    return nb;
  /*  Power-of-two bases: exact via bit counts.  */
  if ((base & (base - 1)) == 0) {
    unsigned bits_per_digit = 0u;
    for (int tmp = base; tmp > 1; tmp >>= 1)
      ++bits_per_digit;
    return (nb + bits_per_digit - 1u) / bits_per_digit;
  }
  /*  General case: ceil(nbits / log2(base)) via fixed-point table.
      recip = ceil(2^16 / log2(base)), so recip >= 1/log2(base) * 2^16.
      digits = floor(nb * recip / 2^16) + 1.  */
  uint64_t recip = (uint64_t) arbint_recip_log2_16[base - 3];
  return (size_t) (((uint64_t) nb * recip) >> 16u) + 1u;
}
