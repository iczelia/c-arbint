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
#include "arbint_ntt.h"
#include "arbint_cpu.h"

#include "config.h"

#include <assert.h>
#include <string.h>

/*  Multiply two limbs producing full double-width result (hi:lo = x * y).
    Uses half-limb multiplication and combining when 2x-width types
    unavailable.  */
static void arbint_mul_wide_limb(arbint_limb_t x, arbint_limb_t y,
                                 arbint_limb_t * hi, arbint_limb_t * lo) {
  arbint_umul_limb_generic(hi, lo, x, y);
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

/*  Barrett division step for a known normalized divisor.
    Computes *q = floor((nh:nl) / d), *r = (nh:nl) mod d.
    Precondition: nh < d, d normalized (MSB set), di = reciprocal of d.  */
static inline void arbint_div3_barrett(arbint_limb_t * q, arbint_limb_t * r,
                                       arbint_limb_t nh, arbint_limb_t nl,
                                       arbint_limb_t d, arbint_limb_t di) {
  arbint_limb_t qh;
  arbint_limb_t ql;
  arbint_limb_t _r;
  arbint_limb_t mask;

  arbint_mul_wide_limb(nh, di, &qh, &ql);

  {
    arbint_limb_t lo = ql + nl;
    arbint_limb_t carry = (lo < ql) ? (arbint_limb_t) 1u : (arbint_limb_t) 0u;
    ql = lo;
    qh = qh + (nh + (arbint_limb_t) 1u) + carry;
  }

  _r = nl - qh * d;

  mask = (arbint_limb_t) 0u - (arbint_limb_t) (_r > ql);
  qh += mask;
  _r += mask & d;

  if (_r >= d) {
    _r -= d;
    qh++;
  }

  *q = qh;
  *r = _r;
}

/*  Exact division by 3, in-place, high-to-low reciprocal division.
    Precondition: x[0..n-1] must be exactly divisible by 3.
    Returns normalized result limb count.

    Uses Barrett reciprocal to replace hardware DIV instructions with
    a single wide multiply + a few ALU ops per limb.

    Compile-time constants for divisor 3:
      64-bit: d_norm = 0xC000000000000000, di = 0x5555555555555555, shift = 62
      32-bit: d_norm = 0xC0000000, di = 0x55555555, shift = 30  */
static size_t arbint_divexact3_generic(arbint_limb_t * x, size_t n) {
  arbint_limb_t rem;
  size_t i;

#if ARBINT_LIMB_BITS == 64
  static const arbint_limb_t d_norm = UINT64_C(0xC000000000000000);
  static const arbint_limb_t di = UINT64_C(0x5555555555555555);
  static const unsigned shift = 62u;
#elif ARBINT_LIMB_BITS == 32
  static const arbint_limb_t d_norm = UINT32_C(0xC0000000);
  static const arbint_limb_t di = UINT32_C(0x55555555);
  static const unsigned shift = 30u;
#else
  #error "Unsupported ARBINT_LIMB_BITS for divexact3"
#endif

  if (n == 0u)
    return 0u;

  rem = x[n - 1u] >> (ARBINT_LIMB_BITS - shift);

  for (i = n; i != 0u;) {
    arbint_limb_t nl;
    arbint_limb_t qi;

    --i;
    nl = x[i] << shift;
    if (i >= 1u)
      nl |= x[i - 1u] >> (ARBINT_LIMB_BITS - shift);

    arbint_div3_barrett(&qi, &rem, rem, nl, d_norm, di);
    x[i] = qi;
  }

  (void) rem;
  return arbint_norm_used(x, n);
}

#define ARBINT_DIVEXACT3_FN arbint_divexact3_generic
#define ARBINT_MUL_LIMB_1_FN arbint_mul_limb_1_generic
#define ARBINT_MUL_IMPL_FN arbint_mul_impl_generic
#define ARBINT_SQR_IMPL_FN arbint_sqr_impl_generic

#include "arbint_mul_core.inc"

/*  Exported wrappers for internal functions used by addmul/submul.  */

arbint_err_t arbint_mul_mag_generic(arbint_limb_t * dst, size_t * out_used,
                                    const arbint_limb_t * a, size_t an,
                                    const arbint_limb_t * b, size_t bn,
                                    const arbint_alloc_t * alloc) {
  return arbint_mul_mag_rec(dst, out_used, a, an, b, bn, alloc);
}

size_t arbint_mulacc_1_generic(arbint_limb_t * dst, size_t dst_n,
                               size_t dst_cap, const arbint_limb_t * a,
                               size_t an, arbint_limb_t b) {
  return arbint_mulacc_1(dst, dst_n, dst_cap, a, an, b);
}

size_t arbint_mulacc_generic(arbint_limb_t * dst, size_t dst_n, size_t dst_cap,
                             const arbint_limb_t * a, size_t an,
                             const arbint_limb_t * c, size_t cn) {
  return arbint_mulacc(dst, dst_n, dst_cap, a, an, c, cn);
}

const arbint_mul_kernel_table_t * arbint_mul_kernel_table_get(void) {
  static arbint_mul_kernel_table_t table;
  static int ready = 0;

  if (!ready) {
#if HAS_BMI2_ALWAYS
    table.mul_impl = arbint_mul_impl_bmi2;
    table.sqr_impl = arbint_sqr_impl_bmi2;
    table.mul_limb_1 = arbint_mul_limb_1_bmi2;
    table.mul_mag = arbint_mul_mag_bmi2;
    table.mulacc = arbint_mulacc_bmi2;
    table.mulacc_1 = arbint_mulacc_1_bmi2;
#elif HAS_BMI2
    if (arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2)) {
      table.mul_impl = arbint_mul_impl_bmi2;
      table.sqr_impl = arbint_sqr_impl_bmi2;
      table.mul_limb_1 = arbint_mul_limb_1_bmi2;
      table.mul_mag = arbint_mul_mag_bmi2;
      table.mulacc = arbint_mulacc_bmi2;
      table.mulacc_1 = arbint_mulacc_1_bmi2;
    } else {
      table.mul_impl = arbint_mul_impl_generic;
      table.sqr_impl = arbint_sqr_impl_generic;
      table.mul_limb_1 = arbint_mul_limb_1_generic;
      table.mul_mag = arbint_mul_mag_generic;
      table.mulacc = arbint_mulacc_generic;
      table.mulacc_1 = arbint_mulacc_1_generic;
    }
#else
    table.mul_impl = arbint_mul_impl_generic;
    table.sqr_impl = arbint_sqr_impl_generic;
    table.mul_limb_1 = arbint_mul_limb_1_generic;
    table.mul_mag = arbint_mul_mag_generic;
    table.mulacc = arbint_mulacc_generic;
    table.mulacc_1 = arbint_mulacc_1_generic;
#endif
    ready = 1;
  }

  return &table;
}
