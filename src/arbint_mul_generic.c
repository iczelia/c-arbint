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

#include "arbint_mul.h"

#include "config.h"

#include <string.h>

/*  Multiply two limbs producing full double-width result (hi:lo = x * y).
    Uses half-limb multiplication and combining when 2x-width types
    unavailable.  */
static void arbint_mul_wide_limb(arbint_limb_t x, arbint_limb_t y,
                                 arbint_limb_t * hi, arbint_limb_t * lo) {
#if ARBINT_LIMB_BITS == 32
  uint64_t p = (uint64_t) x * (uint64_t) y;
  *lo = (arbint_limb_t) p;
  *hi = (arbint_limb_t) (p >> 32);
#else
  const arbint_limb_t mask = ARBINT_HALF_MASK;
  arbint_limb_t x0 = x & mask;
  arbint_limb_t x1 = x >> ARBINT_HALF_BITS;
  arbint_limb_t y0 = y & mask;
  arbint_limb_t y1 = y >> ARBINT_HALF_BITS;

  arbint_limb_t w0 = x0 * y0;
  arbint_limb_t t = x1 * y0 + (w0 >> ARBINT_HALF_BITS);
  arbint_limb_t w1 = t & mask;
  arbint_limb_t w2 = t >> ARBINT_HALF_BITS;

  w1 = x0 * y1 + w1;

  *hi = x1 * y1 + w2 + (w1 >> ARBINT_HALF_BITS);
  *lo = (w1 << ARBINT_HALF_BITS) | (w0 & mask);
#endif /* ARBINT_LIMB_BITS */
}

/*  Compute x*y + acc + carry, returning high limb in result and low limb
    in *out. Essential primitive for schoolbook multiplication inner loop.  */
static arbint_limb_t arbint_muladd_limb(arbint_limb_t x, arbint_limb_t y,
                                        arbint_limb_t acc, arbint_limb_t carry,
                                        arbint_limb_t * out) {
  arbint_limb_t ph;
  arbint_limb_t pl;

  arbint_mul_wide_limb(x, y, &ph, &pl);

#if ARBINT_HAVE_X86_CARRY_KERNEL
  {
    arbint_x86_carry_word_t t = (arbint_x86_carry_word_t) 0;
    arbint_x86_carry_word_t o = (arbint_x86_carry_word_t) 0;
    unsigned char c1;
    unsigned char c2;

    c1 = ARBINT_X86_ADDCARRY((unsigned char) 0, (arbint_x86_carry_word_t) pl,
                             (arbint_x86_carry_word_t) acc, &t);
    c2 = ARBINT_X86_ADDCARRY((unsigned char) 0, t,
                             (arbint_x86_carry_word_t) carry, &o);
    *out = (arbint_limb_t) o;
    return ph + (arbint_limb_t) c1 + (arbint_limb_t) c2;
  }
#else
  {
    arbint_limb_t s;
    arbint_limb_t c1;
    arbint_limb_t s2;
    arbint_limb_t c2;

    s = pl + acc;
    c1 = (s < pl) ? 1u : 0u;
    s2 = s + carry;
    c2 = (s2 < s) ? 1u : 0u;

    *out = s2;
    return ph + c1 + c2;
  }
#endif /* ARBINT_HAVE_X86_CARRY_KERNEL */
}

#define ARBINT_MUL_LIMB_1_FN arbint_mul_limb_1_generic
#define ARBINT_MUL_IMPL_FN arbint_mul_impl_generic
#define ARBINT_SQR_IMPL_FN arbint_sqr_impl_generic

#include "arbint_mul_core.inc"
