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

/*  String conversion functions: arbint_set_str and arbint_get_str.

    Performance optimizations:
    1. Power-of-2 bases (2, 4, 8, 16, 32): Pure bit manipulation, O(n)
    2. Blocking: Extract/parse multiple digits per division/multiplication
    3. Divide-and-conquer: O(n log^2 n) for large numbers  */

#include "arbint_str.h"
#include "arbint.h"
#include "arbint_internal_util.h"

#include <stdlib.h>
#include <string.h>

/*  ========== Forward Declarations ==========  */

static arbint_err_t arbint_get_str_pow2(const arbint_t op, char ** out_str,
                                        unsigned bits_per_digit);
static arbint_err_t arbint_get_str_block(const arbint_t op, char ** out_str,
                                         int base);
static arbint_err_t arbint_get_str_dc(const arbint_t op, char ** out_str,
                                      int base);

static arbint_err_t arbint_set_str_pow2(arbint_t rop, const char * s,
                                        size_t len, unsigned bits_per_digit,
                                        int sign);
static arbint_err_t arbint_set_str_block(arbint_t rop, const char * s,
                                         size_t len, int base, int sign);
static arbint_err_t arbint_set_str_dc(arbint_t rop, const char * s, size_t len,
                                      int base, int sign);

/*  ========== arbint_get_str Implementation ==========  */

arbint_err_t arbint_get_str(const arbint_t op, char ** out_str, int base) {
  unsigned pow2_log;
  size_t nlimbs;
  char * buf;

  /*  Validate inputs.  */
  if (op == NULL || out_str == NULL)
    return ARBINT_EINVAL;
  if (base < 2 || base > 36)
    return ARBINT_EINVAL;

  *out_str = NULL;

  /*  Handle zero.  */
  if (op[0]._sz == 0) {
    buf = arbint_str_alloc(op, 2u);
    if (buf == NULL)
      return ARBINT_ENOMEM;
    buf[0] = '0';
    buf[1] = '\0';
    *out_str = buf;
    return ARBINT_OK;
  }

  /*  Power-of-2 base fast path.  */
  pow2_log = arbint_base_log2_pow2(base);
  if (pow2_log > 0u)
    return arbint_get_str_pow2(op, out_str, pow2_log);

  /*  Choose algorithm based on size.  */
  nlimbs = arbint_abs_sz(op[0]._sz);
  if (nlimbs >= ARBINT_STR_DC_THRESHOLD)
    return arbint_get_str_dc(op, out_str, base);
  else
    return arbint_get_str_block(op, out_str, base);
}

/*  Power-of-2 base: extract digits via bit shifts.
    No division required.  O(n) time complexity.  */
static arbint_err_t arbint_get_str_pow2(const arbint_t op, char ** out_str,
                                        unsigned bits_per_digit) {
  size_t nbits;
  size_t num_digits;
  size_t alloc_size;
  char * buf;
  char * p;
  const arbint_limb_t * limbs;
  size_t nlimbs;
  arbint_limb_t mask;
  size_t digit_idx;

  nbits = arbint_nbits(op);
  num_digits = (nbits + bits_per_digit - 1u) / bits_per_digit;
  alloc_size = num_digits + 2u; /*  sign + NUL  */

  buf = arbint_str_alloc(op, alloc_size);
  if (buf == NULL)
    return ARBINT_ENOMEM;

  p = buf;
  if (op[0]._sz < 0)
    *p++ = '-';

  limbs = ARBINT_CLIMBS(op);
  nlimbs = arbint_abs_sz(op[0]._sz);
  mask = ((arbint_limb_t) 1u << bits_per_digit) - 1u;

  /*  Extract digits from high to low.  */
  for (digit_idx = 0u; digit_idx < num_digits; ++digit_idx) {
    size_t bit_pos;
    size_t limb_idx;
    unsigned bit_offset;
    arbint_limb_t digit;

    /*  Bit position of this digit (counting from LSB).  */
    bit_pos = (num_digits - 1u - digit_idx) * bits_per_digit;
    limb_idx = bit_pos / ARBINT_LIMB_BITS;
    bit_offset = (unsigned) (bit_pos % ARBINT_LIMB_BITS);

    if (limb_idx < nlimbs) {
      digit = (limbs[limb_idx] >> bit_offset) & mask;

      /*  Handle digit spanning two limbs.  */
      if (bit_offset + bits_per_digit > ARBINT_LIMB_BITS &&
          limb_idx + 1u < nlimbs) {
        unsigned high_bits = bit_offset + bits_per_digit - ARBINT_LIMB_BITS;
        arbint_limb_t high_part = limbs[limb_idx + 1u] &
                                  (((arbint_limb_t) 1u << high_bits) - 1u);
        digit |= high_part << (ARBINT_LIMB_BITS - bit_offset);
      }
    } else {
      digit = 0u;
    }

    *p++ = arbint_digit_to_char[digit];
  }

  *p = '\0';
  *out_str = buf;
  return ARBINT_OK;
}

/*  Blocking algorithm: extract multiple digits per division.
    Uses precomputed base^digits as divisor.  */
static arbint_err_t arbint_get_str_block(const arbint_t op, char ** out_str,
                                         int base) {
  arbint_t tmp;
  arbint_err_t rc;
  size_t est_digits;
  size_t alloc_size;
  char * buf;
  char * end;
  char * p;
  const arbint_base_block_t * blk;
  arbint_limb_t divisor;
  unsigned digits_per_block;
  int is_negative;

  is_negative = (op[0]._sz < 0);

  /*  Estimate output size.  */
  est_digits = arbint_sizeinbase(op, base);
  alloc_size = est_digits + 2u; /*  sign + NUL  */

  buf = arbint_str_alloc(op, alloc_size);
  if (buf == NULL)
    return ARBINT_ENOMEM;

  /*  Work on a copy of |op|.  */
  rc = arbint_init(tmp, op[0]._ctx);
  if (rc != ARBINT_OK) {
    arbint_str_free(op, buf);
    return rc;
  }

  rc = arbint_abs(tmp, op);
  if (rc != ARBINT_OK) {
    arbint_clear(tmp);
    arbint_str_free(op, buf);
    return rc;
  }

  blk = &arbint_base_block[base];
  divisor = blk->power;
  digits_per_block = blk->digits;

  /*  Build string from end (low to high digits).  */
  end = buf + alloc_size - 1u;
  *end = '\0';
  p = end;

  while (!arbint_is_zero(tmp)) {
    arbint_limb_t rem;
    arbint_t rem_arbint;
    unsigned i;

    /*  rem = tmp % divisor, tmp = tmp / divisor.  */
    rc = arbint_init(rem_arbint, op[0]._ctx);
    if (rc != ARBINT_OK)
      goto cleanup_fail;

    /*  Use 32-bit division if divisor fits.  */
#if ARBINT_LIMB_BITS == 64
    if (divisor <= UINT32_MAX) {
      rc = arbint_tdiv_qr_u32(tmp, rem_arbint, tmp, (uint32_t) divisor);
    } else {
      /*  Need full arbint division for large divisors.  */
      arbint_t div_arbint;
      rc = arbint_init(div_arbint, op[0]._ctx);
      if (rc != ARBINT_OK) {
        arbint_clear(rem_arbint);
        goto cleanup_fail;
      }
      rc = arbint_resize(div_arbint, 1u);
      if (rc == ARBINT_OK) {
        ARBINT_LIMBS(div_arbint)[0] = divisor;
        div_arbint[0]._sz = 1;
        rc = arbint_tdiv_qr(tmp, rem_arbint, tmp, div_arbint);
      }
      arbint_clear(div_arbint);
    }
#else
    rc = arbint_tdiv_qr_u32(tmp, rem_arbint, tmp, (uint32_t) divisor);
#endif

    if (rc != ARBINT_OK) {
      arbint_clear(rem_arbint);
      goto cleanup_fail;
    }

    /*  Extract remainder as limb value.  */
    if (rem_arbint[0]._sz == 0) {
      rem = 0u;
    } else {
      rem = ARBINT_CLIMBS(rem_arbint)[0];
    }
    arbint_clear(rem_arbint);

    /*  Convert remainder to digits.
        For intermediate blocks, write exactly digits_per_block chars.
        For the last block (when tmp becomes zero), write only significant digits.  */
    if (arbint_is_zero(tmp)) {
      /*  Last block: write only significant digits.  */
      do {
        unsigned digit = (unsigned) (rem % (arbint_limb_t) base);
        rem /= (arbint_limb_t) base;
        *--p = arbint_digit_to_char[digit];
      } while (rem != 0u && p > buf);
    } else {
      /*  Intermediate block: write full block (with leading zeros).  */
      for (i = 0u; i < digits_per_block && p > buf; ++i) {
        unsigned digit = (unsigned) (rem % (arbint_limb_t) base);
        rem /= (arbint_limb_t) base;
        *--p = arbint_digit_to_char[digit];
      }
    }
  }

  /*  Handle case where input was zero.  */
  if (p == end)
    *--p = '0';

  /*  Move result to start of buffer, adding sign if needed.  */
  {
    size_t digit_len = (size_t) (end - p);
    char * dst = buf;

    if (is_negative)
      *dst++ = '-';

    if (p != dst)
      memmove(dst, p, digit_len + 1u); /*  +1 for NUL  */
  }

  arbint_clear(tmp);
  *out_str = buf;
  return ARBINT_OK;

cleanup_fail:
  arbint_clear(tmp);
  arbint_str_free(op, buf);
  return rc;
}

/*  Divide-and-conquer: split number at base^(digits/2), recurse.
    O(n log^2 n) with fast multiplication/division.  */
static arbint_err_t arbint_get_str_dc(const arbint_t op, char ** out_str,
                                      int base) {
  size_t nlimbs;
  size_t est_digits;
  size_t half_digits;
  arbint_t divisor;
  arbint_t high;
  arbint_t low;
  arbint_t base_arbint;
  arbint_t abs_op;
  char * high_str;
  char * low_str;
  size_t high_len;
  size_t low_len;
  size_t total_len;
  char * buf;
  char * p;
  size_t pad;
  int is_negative;
  arbint_err_t rc;

  nlimbs = arbint_abs_sz(op[0]._sz);

  /*  Base case: use blocking algorithm.  */
  if (nlimbs < ARBINT_STR_DC_THRESHOLD)
    return arbint_get_str_block(op, out_str, base);

  is_negative = (op[0]._sz < 0);

  /*  Estimate digit count and compute split point.  */
  est_digits = arbint_sizeinbase(op, base);
  half_digits = est_digits / 2u;

  rc = arbint_init_all(op[0]._ctx, base_arbint, divisor, high, low, abs_op,
                       (arbint_t *) NULL);
  if (rc != ARBINT_OK)
    return rc;

  /*  Compute divisor = base^half_digits.  */
  rc = arbint_set_i32(base_arbint, base);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Use binary exponentiation for base^half_digits.  */
  if (half_digits <= UINT32_MAX) {
    rc = arbint_pow_u32(divisor, base_arbint, (uint32_t) half_digits);
  } else {
    /*  Extremely large: fall back to blocking.  */
    arbint_clear_all(base_arbint, divisor, high, low, abs_op,
                     (arbint_t *) NULL);
    return arbint_get_str_block(op, out_str, base);
  }

  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Split: high = |op| / divisor, low = |op| % divisor.  */
  rc = arbint_abs(abs_op, op);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_tdiv_qr(high, low, abs_op, divisor);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Recursively convert high and low parts.  */
  high_str = NULL;
  low_str = NULL;

  rc = arbint_get_str(high, &high_str, base);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_get_str(low, &low_str, base);
  if (rc != ARBINT_OK) {
    arbint_str_free(op, high_str);
    goto cleanup;
  }

  /*  All arbints done; free them before string assembly.  */
  arbint_clear_all(base_arbint, divisor, high, low, abs_op,
                   (arbint_t *) NULL);

  /*  Combine: high_str + zero_pad + low_str.  */
  high_len = strlen(high_str);
  low_len = strlen(low_str);

  /*  Low part must be exactly half_digits chars (pad with leading zeros).  */
  if (low_len < half_digits)
    pad = half_digits - low_len;
  else
    pad = 0u;

  total_len = high_len + pad + low_len + (is_negative ? 1u : 0u) + 1u;

  buf = arbint_str_alloc(op, total_len);
  if (buf == NULL) {
    arbint_str_free(op, high_str);
    arbint_str_free(op, low_str);
    return ARBINT_ENOMEM;
  }

  p = buf;
  if (is_negative)
    *p++ = '-';

  memcpy(p, high_str, high_len);
  p += high_len;

  /*  Pad low part with leading zeros.  */
  memset(p, '0', pad);
  p += pad;

  memcpy(p, low_str, low_len + 1u); /*  Include NUL.  */

  arbint_str_free(op, high_str);
  arbint_str_free(op, low_str);

  *out_str = buf;
  return ARBINT_OK;

cleanup:
  arbint_clear_all(base_arbint, divisor, high, low, abs_op,
                   (arbint_t *) NULL);
  return rc;
}

/*  ========== arbint_set_str Implementation ==========  */

arbint_err_t arbint_set_str(arbint_t rop, const char * s, int base) {
  int sign;
  size_t len;
  unsigned pow2_log;
  const char * start;
  int saw_digit;

  /*  Validate inputs.  */
  if (rop == NULL || s == NULL)
    return ARBINT_EINVAL;
  if (base < 2 || base > 36)
    return ARBINT_EINVAL;

  /*  Skip leading whitespace.  */
  while (*s == ' ' || *s == '\t' || *s == '\n' || *s == '\r')
    ++s;

  /*  Parse optional sign.  */
  sign = 1;
  if (*s == '-') {
    sign = -1;
    ++s;
  } else if (*s == '+') {
    ++s;
  }

  /*  Skip leading zeros (these count as digits).  */
  saw_digit = 0;
  while (*s == '0') {
    saw_digit = 1;
    ++s;
  }

  /*  Handle empty string after stripping.  */
  if (*s == '\0') {
    /*  Must have seen at least one digit (zero counts).  */
    if (!saw_digit)
      return ARBINT_EINVAL;
    arbint_zero(rop);
    return ARBINT_OK;
  }

  start = s;
  len = strlen(s);

  /*  Validate all characters for the given base.  */
  {
    size_t i;
    for (i = 0u; i < len; ++i) {
      uint8_t digit = arbint_char_to_digit[(unsigned char) s[i]];
      if (digit == 0xFFu || digit >= (unsigned) base)
        return ARBINT_EINVAL;
    }
  }

  /*  Power-of-2 base fast path.  */
  pow2_log = arbint_base_log2_pow2(base);
  if (pow2_log > 0u)
    return arbint_set_str_pow2(rop, start, len, pow2_log, sign);

  /*  Choose algorithm based on string length.  */
  if (len >= ARBINT_STR_DC_DIGIT_THRESHOLD)
    return arbint_set_str_dc(rop, start, len, base, sign);
  else
    return arbint_set_str_block(rop, start, len, base, sign);
}

/*  Power-of-2 base: pack digits directly into limbs.
    O(n) time complexity.  */
static arbint_err_t arbint_set_str_pow2(arbint_t rop, const char * s,
                                        size_t len, unsigned bits_per_digit,
                                        int sign) {
  size_t total_bits;
  size_t nlimbs;
  arbint_limb_t * limbs;
  size_t i;
  arbint_err_t rc;
  size_t used;

  /*  Check for overflow in len * bits_per_digit.  */
  if (len > SIZE_MAX / bits_per_digit)
    return ARBINT_EOVERFLOW;

  total_bits = len * bits_per_digit;
  nlimbs = (total_bits + ARBINT_LIMB_BITS - 1u) / ARBINT_LIMB_BITS;

  rc = arbint_resize(rop, nlimbs);
  if (rc != ARBINT_OK)
    return rc;

  limbs = ARBINT_LIMBS(rop);
  memset(limbs, 0, nlimbs * sizeof(arbint_limb_t));

  /*  Pack digits from high (s[0]) to low (s[len-1]).  */
  for (i = 0u; i < len; ++i) {
    uint8_t digit;
    size_t bit_pos;
    size_t limb_idx;
    unsigned bit_offset;

    digit = arbint_char_to_digit[(unsigned char) s[i]];

    /*  Bit position for this digit (counting from LSB).  */
    bit_pos = (len - 1u - i) * bits_per_digit;
    limb_idx = bit_pos / ARBINT_LIMB_BITS;
    bit_offset = (unsigned) (bit_pos % ARBINT_LIMB_BITS);

    limbs[limb_idx] |= ((arbint_limb_t) digit) << bit_offset;

    /*  Handle digit spanning two limbs.  */
    if (bit_offset + bits_per_digit > ARBINT_LIMB_BITS &&
        limb_idx + 1u < nlimbs) {
      unsigned high_bits = bit_offset + bits_per_digit - ARBINT_LIMB_BITS;
      limbs[limb_idx + 1u] |= ((arbint_limb_t) digit) >>
                              (ARBINT_LIMB_BITS - bit_offset);
      (void) high_bits; /*  Suppress unused warning.  */
    }
  }

  used = arbint_norm_used(limbs, nlimbs);
  if (!arbint_set_signed_sz(rop, used, sign))
    return ARBINT_EOVERFLOW;

  return ARBINT_OK;
}

/*  Blocking algorithm: parse multiple digits per multiplication.
    rop = rop * base^digits + chunk_value.  */
static arbint_err_t arbint_set_str_block(arbint_t rop, const char * s,
                                         size_t len, int base, int sign) {
  const arbint_base_block_t * blk;
  unsigned digits_per_block;
  arbint_limb_t block_mul;
  size_t pos;
  arbint_err_t rc;

  blk = &arbint_base_block[base];
  digits_per_block = blk->digits;
  block_mul = blk->power;

  arbint_zero(rop);

  pos = 0u;
  while (pos < len) {
    size_t block_size;
    size_t remaining;
    arbint_limb_t block_val;
    arbint_limb_t multiplier;
    size_t j;

    remaining = len - pos;

    /*  First block may be smaller if len % digits_per_block != 0.  */
    if (pos == 0u && (len % digits_per_block) != 0u)
      block_size = len % digits_per_block;
    else
      block_size = (remaining >= digits_per_block) ? digits_per_block :
                                                     remaining;

    /*  Parse block to limb value.  */
    block_val = 0u;
    for (j = 0u; j < block_size; ++j) {
      uint8_t digit = arbint_char_to_digit[(unsigned char) s[pos + j]];
      block_val = block_val * (arbint_limb_t) base + digit;
    }

    /*  Compute multiplier = base^block_size.  */
    if (block_size == digits_per_block) {
      multiplier = block_mul;
    } else {
      multiplier = 1u;
      for (j = 0u; j < block_size; ++j)
        multiplier *= (arbint_limb_t) base;
    }

    /*  rop = rop * multiplier + block_val.  */
    if (rop[0]._sz != 0) {
      /*  Need to use full arbint multiplication for large multipliers.  */
#if ARBINT_LIMB_BITS == 64
      if (multiplier <= UINT32_MAX) {
        rc = arbint_mul_u32(rop, rop, (uint32_t) multiplier);
      } else {
        arbint_t mul_arbint;
        rc = arbint_init(mul_arbint, rop[0]._ctx);
        if (rc != ARBINT_OK)
          return rc;
        rc = arbint_resize(mul_arbint, 1u);
        if (rc == ARBINT_OK) {
          ARBINT_LIMBS(mul_arbint)[0] = multiplier;
          mul_arbint[0]._sz = 1;
          rc = arbint_mul(rop, rop, mul_arbint);
        }
        arbint_clear(mul_arbint);
      }
#else
      rc = arbint_mul_u32(rop, rop, (uint32_t) multiplier);
#endif
      if (rc != ARBINT_OK)
        return rc;
    }

    /*  Add block_val.  */
#if ARBINT_LIMB_BITS == 64
    if (block_val <= UINT32_MAX) {
      rc = arbint_add_u32(rop, rop, (uint32_t) block_val);
    } else {
      arbint_t val_arbint;
      rc = arbint_init(val_arbint, rop[0]._ctx);
      if (rc != ARBINT_OK)
        return rc;
      rc = arbint_resize(val_arbint, 1u);
      if (rc == ARBINT_OK) {
        ARBINT_LIMBS(val_arbint)[0] = block_val;
        val_arbint[0]._sz = 1;
        rc = arbint_add(rop, rop, val_arbint);
      }
      arbint_clear(val_arbint);
    }
#else
    rc = arbint_add_u32(rop, rop, (uint32_t) block_val);
#endif
    if (rc != ARBINT_OK)
      return rc;

    pos += block_size;
  }

  /*  Apply sign.  */
  if (sign < 0 && rop[0]._sz > 0)
    rop[0]._sz = -rop[0]._sz;

  return ARBINT_OK;
}

/*  Divide-and-conquer: parse halves, combine with base^(len/2).
    O(n log^2 n) with fast multiplication.  */
static arbint_err_t arbint_set_str_dc(arbint_t rop, const char * s, size_t len,
                                      int base, int sign) {
  size_t half;
  arbint_t high;
  arbint_t low;
  arbint_t base_power;
  arbint_t base_arbint;
  arbint_err_t rc;

  /*  Base case: use blocking algorithm.  */
  if (len < ARBINT_STR_DC_DIGIT_THRESHOLD)
    return arbint_set_str_block(rop, s, len, base, sign);

  half = len / 2u;

  rc = arbint_init_all(rop[0]._ctx, high, low, base_arbint, base_power,
                       (arbint_t *) NULL);
  if (rc != ARBINT_OK)
    return rc;

  /*  Parse s[0..half-1] -> high (unsigned).  */
  rc = arbint_set_str_dc(high, s, half, base, 1);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Parse s[half..len-1] -> low (unsigned).  */
  rc = arbint_set_str_dc(low, s + half, len - half, base, 1);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Compute base^(len - half).  */
  rc = arbint_set_i32(base_arbint, base);
  if (rc != ARBINT_OK)
    goto cleanup;

  if ((len - half) <= UINT32_MAX) {
    rc = arbint_pow_u32(base_power, base_arbint, (uint32_t) (len - half));
  } else {
    /*  Extremely large exponent: fall back to blocking.  */
    arbint_clear_all(high, low, base_arbint, base_power, (arbint_t *) NULL);
    return arbint_set_str_block(rop, s, len, base, sign);
  }

  if (rc != ARBINT_OK)
    goto cleanup;

  /*  rop = high * base_power + low.  */
  rc = arbint_mul(rop, high, base_power);
  if (rc != ARBINT_OK)
    goto cleanup;

  rc = arbint_add(rop, rop, low);
  if (rc != ARBINT_OK)
    goto cleanup;

  /*  Apply sign.  */
  if (sign < 0 && rop[0]._sz > 0)
    rop[0]._sz = -rop[0]._sz;

cleanup:
  arbint_clear_all(high, low, base_arbint, base_power, (arbint_t *) NULL);
  return rc;
}
