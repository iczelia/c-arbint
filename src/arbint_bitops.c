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
