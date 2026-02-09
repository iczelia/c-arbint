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
#include "config.h"
#include "arbint_cpu.h"
#include "arbint.h"

#include <limits.h>
#include <string.h>

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

/* ========== Internal magnitude helpers ========== */

/*  Increment magnitude by 1.  dst may alias src.
    Caller must ensure dst has capacity for n+1 limbs.
    Returns used limb count (may be n+1). dst and src
    may alias, careful!  */
static size_t arbint_mag_inc(arbint_limb_t * dst, const arbint_limb_t * src,
                             size_t n) {
  size_t i;
  arbint_limb_t carry = 1u;

  for (i = 0u; i < n; ++i) {
    arbint_limb_t xi = src[i];
    arbint_limb_t s = xi + carry;
    dst[i] = s;
    if (s >= xi) {
      carry = 0u;
      ++i;
      break;
    }
    carry = 1u;
  }

  if (i < n && dst != src)
    memcpy(dst + i, src + i, (n - i) * sizeof(arbint_limb_t));

  if (carry != 0u) {
    dst[n] = 1u;
    return n + 1u;
  }
  return n;
}

/*  Decrement magnitude by 1.  dst may alias src.
    Precondition: the value represented by src[0..n-1] is nonzero.
    Returns normalized used limb count.  */
static size_t arbint_mag_dec(arbint_limb_t * dst, const arbint_limb_t * src,
                             size_t n) {
  size_t i;
  arbint_limb_t borrow = 1u;

  for (i = 0u; i < n; ++i) {
    arbint_limb_t xi = src[i];
    arbint_limb_t d = xi - borrow;
    dst[i] = d;
    if (xi >= borrow) {
      borrow = 0u;
      ++i;
      break;
    }
    borrow = 1u;
  }

  if (i < n && dst != src)
    memcpy(dst + i, src + i, (n - i) * sizeof(arbint_limb_t));

  return arbint_norm_used(dst, n);
}

/*  Find index of first nonzero limb.  Returns n if all zero.  */
static size_t arbint_first_nonzero(const arbint_limb_t * p, size_t n) {
  size_t i;
  for (i = 0u; i < n; ++i) {
    if (p[i] != 0u)
      return i;
  }
  return n;
}

/* ========== Scalar popcount of limb array ========== */

static size_t arbint_popcount_limbs_scalar(const arbint_limb_t * x, size_t n) {
  size_t count = 0u;
  size_t i;
  for (i = 0u; i < n; ++i)
    count += arbint_popcount_limb(x[i]);
  return count;
}

/* ========== Dispatch for popcount on limb arrays ========== */

typedef size_t (*arbint_popcount_limbs_fn_t)(const arbint_limb_t *, size_t);

static arbint_popcount_limbs_fn_t arbint_select_popcount_limbs(void) {
#if HAS_POPCNT_ALWAYS
  return arbint__popcount_limbs_popcnt;
#elif HAS_POPCNT
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_POPCNT)
             ? arbint__popcount_limbs_popcnt
             : arbint_popcount_limbs_scalar;
#else
  return arbint_popcount_limbs_scalar;
#endif /* HAS_POPCNT_ALWAYS */
}

static size_t arbint_popcount_limbs(const arbint_limb_t * x, size_t n) {
  static arbint_popcount_limbs_fn_t impl = NULL;

  if (impl == NULL)
    impl = arbint_select_popcount_limbs();

  return impl(x, n);
}

/* ========== Scalar hamming distance of limb arrays ========== */

static size_t arbint_hamming_limbs_scalar(const arbint_limb_t * a, size_t an,
                                          const arbint_limb_t * b, size_t bn) {
  size_t count = 0u;
  size_t i;
  size_t min_n = (an < bn) ? an : bn;

  for (i = 0u; i < min_n; ++i)
    count += arbint_popcount_limb(a[i] ^ b[i]);
  for (i = min_n; i < an; ++i)
    count += arbint_popcount_limb(a[i]);
  for (i = min_n; i < bn; ++i)
    count += arbint_popcount_limb(b[i]);
  return count;
}

/* ========== Dispatch for hamming distance on limb arrays ========== */

typedef size_t (*arbint_hamming_limbs_fn_t)(const arbint_limb_t *, size_t,
                                            const arbint_limb_t *, size_t);

static arbint_hamming_limbs_fn_t arbint_select_hamming_limbs(void) {
#if HAS_POPCNT_ALWAYS
  return arbint__hamming_limbs_popcnt;
#elif HAS_POPCNT
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_POPCNT)
             ? arbint__hamming_limbs_popcnt
             : arbint_hamming_limbs_scalar;
#else
  return arbint_hamming_limbs_scalar;
#endif /* HAS_POPCNT_ALWAYS */
}

static size_t arbint_hamming_limbs(const arbint_limb_t * a, size_t an,
                                   const arbint_limb_t * b, size_t bn) {
  static arbint_hamming_limbs_fn_t impl = NULL;

  if (impl == NULL)
    impl = arbint_select_hamming_limbs();

  return impl(a, an, b, bn);
}

/* ========== Public API: testbit ========== */

/*  Test bit at bit_index in two's-complement representation of x.
    For non-negative x, this is simply the bit in the magnitude.
    For negative x, TC(-x) = ~(x-1): find first nonzero limb to
    determine which region the target falls in.  */
arbint_err_t arbint_testbit(const arbint_t x, size_t bit_index, int * out) {
  size_t limb_idx;
  unsigned bit_idx;
  size_t used;
  const arbint_limb_t * p;

  if (x == NULL || out == NULL)
    return ARBINT_EINVAL;

  used = arbint_abs_sz(x[0]._sz);
  limb_idx = bit_index / ARBINT_LIMB_BITS;
  bit_idx = (unsigned) (bit_index % ARBINT_LIMB_BITS);

  if (x[0]._sz >= 0) {
    if (limb_idx >= used) {
      *out = 0;
      return ARBINT_OK;
    }
    *out = (int) ((ARBINT_CLIMBS(x)[limb_idx] >> bit_idx) & 1u);
    return ARBINT_OK;
  }

  /*  Negative: TC(-x) = ~(x-1).  Sign-extends as all-ones.  */
  if (limb_idx >= used) {
    *out = 1;
    return ARBINT_OK;
  }

  p = ARBINT_CLIMBS(x);
  /*  Find first nonzero limb, scanning at most to limb_idx.  */
  size_t k;
  arbint_limb_t tc;

  for (k = 0u; k <= limb_idx && k < used; ++k) {
    if (p[k] != 0u)
      break;
  }

  if (k > limb_idx) {
    /*  limb_idx is in trailing-zero region: (x-1) borrows through,
        giving 0xFFFF..., so TC = ~0xFFFF... = 0.  */
    tc = 0u;
  } else if (k == limb_idx) {
    /*  First nonzero limb is at limb_idx: TC = ~(p[k] - 1).  */
    tc = ~(p[k] - 1u);
  } else {
    /*  k < limb_idx: borrow absorbed, TC = ~p[limb_idx].  */
    tc = ~p[limb_idx];
  }
  *out = (int) ((tc >> bit_idx) & 1u);
  return ARBINT_OK;
}

/* ========== Public API: setbit ========== */

/*  Set bit at bit_index in two's-complement representation of x.
    For non-negative x: straightforward magnitude bit set.
    For negative x: setting a TC bit to 1 reduces magnitude.
    Uses algebraic approach: result = -(((x-1) & ~(1<<p)) + 1).  */
arbint_err_t arbint_setbit(arbint_t x, size_t bit_index) {
  size_t limb_idx;
  unsigned bit_idx;
  size_t used;
  arbint_err_t rc;

  if (x == NULL)
    return ARBINT_EINVAL;

  limb_idx = bit_index / ARBINT_LIMB_BITS;
  bit_idx = (unsigned) (bit_index % ARBINT_LIMB_BITS);
  used = arbint_abs_sz(x[0]._sz);

  if (x[0]._sz >= 0) {
    /*  Non-negative: set bit in magnitude.  */
    size_t need = limb_idx + 1u;
    if (need == 0u)
      return ARBINT_EOVERFLOW;
    if (x[0]._cap < need) {
      rc = arbint_resize(x, need);
      if (rc != ARBINT_OK)
        return rc;
    }
    if (need > used)
      memset(ARBINT_LIMBS(x) + used, 0, (need - used) * sizeof(arbint_limb_t));
    ARBINT_LIMBS(x)[limb_idx] |= (arbint_limb_t) 1u << bit_idx;
    size_t new_used = arbint_norm_used(ARBINT_LIMBS(x), need);
    if (!arbint_set_signed_sz(x, new_used, (new_used == 0u) ? 0 : 1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  /*  Negative: setting bit p in TC(-x) = ~(x-1).
      If bit is already 1 in TC, no-op.
      Result magnitude = ((x-1) & ~(1<<p)) + 1.  */
  int cur;
  rc = arbint_testbit(x, bit_index, &cur);
  if (rc != ARBINT_OK)
    return rc;
  if (cur)
    return ARBINT_OK;

  /*  Bit is 0 in TC -- need to set it.  Setting a 0-bit to 1 in TC of a
      negative number means clearing a bit in (x-1), then adding 1.
      The result magnitude can only decrease or stay the same.  */
  const arbint_limb_t * p = ARBINT_CLIMBS(x);
  arbint_limb_t * rp;
  size_t k;
  size_t i;
  size_t new_used;

  /*  Ensure capacity (result <= used limbs).  */
  if (x[0]._cap < used) {
    rc = arbint_resize(x, used);
    if (rc != ARBINT_OK)
      return rc;
  }

  rp = ARBINT_LIMBS(x);
  p = ARBINT_CLIMBS(x);

  /*  Compute (x-1) in place, apply & ~(1<<p), then +1.
      Step 1: compute x-1 (borrow propagation).  */
  k = arbint_first_nonzero(p, used);
  /*  k < used always since x is nonzero.  */

  /*  Copy limbs and apply decrement, then clear the target bit.  */
  for (i = 0u; i < k; ++i)
    rp[i] = ~(arbint_limb_t) 0u; /*  0 - 1 = 0xFFFF... */
  rp[k] = p[k] - 1u;
  for (i = k + 1u; i < used; ++i)
    rp[i] = p[i];

  /*  Clear the target bit in (x-1).  */
  if (limb_idx < used)
    rp[limb_idx] &= ~((arbint_limb_t) 1u << bit_idx);
  /*  If limb_idx >= used, the bit was in sign-extension (all-ones in TC),
      but we already checked it was 0, so this cannot happen for setbit.  */

  /*  Add 1 to get final magnitude.  */
  new_used = arbint_mag_inc(rp, rp, used);
  /*  Normalize: magnitude could shrink.  */
  new_used = arbint_norm_used(rp, new_used);
  if (new_used == 0u)
    arbint_zero(x);
  else if (!arbint_set_signed_sz(x, new_used, -1))
    return ARBINT_EOVERFLOW;
  return ARBINT_OK;
}

/* ========== Public API: clrbit ========== */

/*  Clear bit at bit_index in two's-complement representation of x.
    For non-negative x: straightforward magnitude bit clear.
    For negative x: clearing a TC bit to 0 increases magnitude.
    Uses algebraic approach: result = -(((x-1) | (1<<p)) + 1).  */
arbint_err_t arbint_clrbit(arbint_t x, size_t bit_index) {
  size_t limb_idx;
  unsigned bit_idx;
  size_t used;
  arbint_err_t rc;

  if (x == NULL)
    return ARBINT_EINVAL;

  limb_idx = bit_index / ARBINT_LIMB_BITS;
  bit_idx = (unsigned) (bit_index % ARBINT_LIMB_BITS);
  used = arbint_abs_sz(x[0]._sz);

  if (x[0]._sz >= 0) {
    if (limb_idx >= used)
      return ARBINT_OK;
    ARBINT_LIMBS(x)[limb_idx] &= ~((arbint_limb_t) 1u << bit_idx);
    size_t new_used = arbint_norm_used(ARBINT_LIMBS(x), used);
    if (!arbint_set_signed_sz(x, new_used, (new_used == 0u) ? 0 : 1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  /*  Negative: clearing bit p in TC(-x) = ~(x-1).
      If bit is already 0 in TC, no-op.
      Result magnitude = ((x-1) | (1<<p)) + 1.  */
  int cur;
  rc = arbint_testbit(x, bit_index, &cur);
  if (rc != ARBINT_OK)
    return rc;
  if (!cur)
    return ARBINT_OK;

  /*  Bit is 1 in TC -- need to clear it.  Setting the bit in (x-1)
      then adding 1 gives the new magnitude.  Result may be larger.  */
  const arbint_limb_t * p;
  arbint_limb_t * rp;
  size_t k;
  size_t i;
  size_t need;
  size_t new_used;

  /*  Result may need limb_idx+2 limbs (setting a bit beyond current
      magnitude of (x-1), then +1 can carry).  */
  need = (limb_idx + 2u > used + 1u) ? limb_idx + 2u : used + 1u;
  if (x[0]._cap < need) {
    rc = arbint_resize(x, need);
    if (rc != ARBINT_OK)
      return rc;
  }

  rp = ARBINT_LIMBS(x);
  p = ARBINT_CLIMBS(x);

  /*  Compute (x-1).  */
  k = arbint_first_nonzero(p, used);

  for (i = 0u; i < k; ++i)
    rp[i] = ~(arbint_limb_t) 0u;
  rp[k] = p[k] - 1u;
  for (i = k + 1u; i < used; ++i)
    rp[i] = p[i];
  /*  Zero any extension beyond used.  */
  for (i = used; i < need; ++i)
    rp[i] = 0u;

  /*  Set the target bit in (x-1).  */
  rp[limb_idx] |= (arbint_limb_t) 1u << bit_idx;

  /*  Normalize (x-1) | bit, then add 1.  */
  {
    size_t xm1_used = arbint_norm_used(rp, need);
    if (xm1_used == 0u) {
      /*  (x-1) | bit == 0 should not happen since we set a bit.  */
      arbint_zero(x);
      return ARBINT_OK;
    }
    new_used = arbint_mag_inc(rp, rp, xm1_used);
    new_used = arbint_norm_used(rp, new_used);
  }
  if (new_used == 0u)
    arbint_zero(x);
  else if (!arbint_set_signed_sz(x, new_used, -1))
    return ARBINT_EOVERFLOW;
  return ARBINT_OK;
}

/* ========== Public API: ctz ========== */

/*  Count trailing zeros in |x|.  Returns EDOM if x == 0.  */
arbint_err_t arbint_ctz(const arbint_t x, size_t * out) {
  size_t used;
  size_t i;
  const arbint_limb_t * p;

  if (x == NULL || out == NULL)
    return ARBINT_EINVAL;

  used = arbint_abs_sz(x[0]._sz);
  if (used == 0u)
    return ARBINT_EDOM;

  p = ARBINT_CLIMBS(x);
  for (i = 0u; i < used; ++i) {
    if (p[i] != 0u) {
      *out = i * ARBINT_LIMB_BITS + arbint_ctz_limb(p[i]);
      return ARBINT_OK;
    }
  }

  return ARBINT_EDOM;
}

/* ========== Public API: clz ========== */

/*  Count leading zeros within canonical width W(x).
    W(x) is the smallest positive multiple of LIMB_BITS that can represent
    x in two's complement.  For x < 0, the TC representation has the MSB
    set, so clz is 0.  Returns EDOM if x == 0.  */
arbint_err_t arbint_clz(const arbint_t x, size_t * out) {
  size_t used;

  if (x == NULL || out == NULL)
    return ARBINT_EINVAL;

  used = arbint_abs_sz(x[0]._sz);
  if (used == 0u)
    return ARBINT_EDOM;

  if (x[0]._sz < 0) {
    *out = 0u;
    return ARBINT_OK;
  }

  *out = (unsigned) arbint_clz_limb(ARBINT_CLIMBS(x)[used - 1u]);
  return ARBINT_OK;
}

/* ========== Public API: popcount ========== */

/*  Count set bits in the low W(x) bits of x.
    W(x) = used * LIMB_BITS (canonical two's-complement width).
    For positive x, popcount of the magnitude limbs.
    For negative -x, popcount(TC) = popcount(~(x-1)):
      limbs before first nonzero: contribute 0.
      first nonzero limb k: LIMB_BITS - popcount(x[k] - 1).
      remaining limbs: LIMB_BITS - popcount(x[i]).  */
arbint_err_t arbint_popcount(const arbint_t x, size_t * out) {
  size_t used;

  if (x == NULL || out == NULL)
    return ARBINT_EINVAL;

  used = arbint_abs_sz(x[0]._sz);
  if (used == 0u) {
    *out = 0u;
    return ARBINT_OK;
  }

  if (x[0]._sz > 0) {
    *out = arbint_popcount_limbs(ARBINT_CLIMBS(x), used);
    return ARBINT_OK;
  }

  /*  Negative: popcount(~(x-1)) without TC buffer.  */
  const arbint_limb_t * p = ARBINT_CLIMBS(x);
  size_t k;
  size_t count;
  size_t i;

  k = arbint_first_nonzero(p, used);
  /*  k < used for normalized nonzero.
      Limbs 0..k-1: TC is 0 (contribute 0 bits).
      Limb k: TC = ~(p[k]-1).  */
  count = (size_t) ARBINT_LIMB_BITS - (size_t) arbint_popcount_limb(p[k] - 1u);
  for (i = k + 1u; i < used; ++i)
    count += (size_t) ARBINT_LIMB_BITS - (size_t) arbint_popcount_limb(p[i]);
  *out = count;
  return ARBINT_OK;
}

/* ========== Public API: hammingdist ========== */

/*  Hamming distance between TC representations of a and b within
    max(W(a), W(b)) bits.

    Both non-negative: direct magnitude hamming (hardware POPCNT dispatch).
    Both negative (-a, -b): hamming = popcount((a-1) ^ (b-1)).
    Mixed (+a, -b): hamming = max_n * LIMB_BITS - popcount(a ^ (b-1)).  */
arbint_err_t arbint_hammingdist(const arbint_t a, const arbint_t b,
                                size_t * out) {
  size_t an;
  size_t bn;
  int a_neg;
  int b_neg;

  if (a == NULL || b == NULL || out == NULL)
    return ARBINT_EINVAL;

  an = arbint_abs_sz(a[0]._sz);
  bn = arbint_abs_sz(b[0]._sz);
  a_neg = (a[0]._sz < 0);
  b_neg = (b[0]._sz < 0);

  /*  Both non-negative: direct magnitude hamming.  */
  if (!a_neg && !b_neg) {
    *out = arbint_hamming_limbs(ARBINT_CLIMBS(a), an, ARBINT_CLIMBS(b), bn);
    return ARBINT_OK;
  }

  size_t max_n = (an > bn) ? an : bn;
  const arbint_limb_t * ap = (an > 0u) ? ARBINT_CLIMBS(a) : NULL;
  const arbint_limb_t * bp = (bn > 0u) ? ARBINT_CLIMBS(b) : NULL;
  size_t count = 0u;
  size_t i;

  if (a_neg && b_neg) {
    /*  Both negative: popcount((a-1) ^ (b-1)).
        Phase 1: both borrows active for min_n limbs.
        Phase 2: shorter operand zero-extends, (shorter-1) = 0,
        so XOR = (longer-1)[i]; continue longer borrow only.  */
    size_t min_n = (an < bn) ? an : bn;
    const arbint_limb_t * longer = (an > bn) ? ap : bp;
    arbint_limb_t ba = 1u;
    arbint_limb_t bb = 1u;
    arbint_limb_t * longer_borrow = (an > bn) ? &ba : &bb;

    for (i = 0u; i < min_n; ++i) {
      arbint_limb_t ad = ap[i] - ba;
      ba = (ap[i] < ba) ? 1u : 0u;
      arbint_limb_t bd = bp[i] - bb;
      bb = (bp[i] < bb) ? 1u : 0u;
      count += arbint_popcount_limb(ad ^ bd);
    }
    for (i = min_n; i < max_n; ++i) {
      arbint_limb_t ld = longer[i] - *longer_borrow;
      *longer_borrow = (longer[i] < *longer_borrow) ? 1u : 0u;
      count += arbint_popcount_limb(ld);
    }
    *out = count;
    return ARBINT_OK;
  }

  /*  Mixed sign: ensure a is positive, b is negative.  */
  if (a_neg) {
    const arbint_limb_t * tmp_p = ap;
    size_t tmp_n = an;
    ap = bp;
    an = bn;
    bp = tmp_p;
    bn = tmp_n;
  }

  /*  (+a) vs (-b): hamming = max_n * LIMB_BITS - popcount(a ^ (b-1)).
      Phase 1 (i < min(an,bn)): both operands present.
      Phase 2 (bn <= i < an): (b-1) zero-extends to 0, XOR = a[i].
      Phase 2 (an <= i < bn): a zero-extends to 0, XOR = (b-1)[i],
      continue borrow.  */
  size_t min_n = (an < bn) ? an : bn;
  arbint_limb_t bb = 1u;
  size_t xor_pop = 0u;

  for (i = 0u; i < min_n; ++i) {
    arbint_limb_t bd = bp[i] - bb;
    bb = (bp[i] < bb) ? 1u : 0u;
    xor_pop += arbint_popcount_limb(ap[i] ^ bd);
  }
  if (an >= bn) {
    /*  Beyond bn: (b-1) = 0, XOR = a[i].  */
    for (i = min_n; i < max_n; ++i)
      xor_pop += arbint_popcount_limb(ap[i]);
  } else {
    /*  Beyond an: a = 0, XOR = (b-1)[i].  */
    for (i = min_n; i < max_n; ++i) {
      arbint_limb_t bd = bp[i] - bb;
      bb = (bp[i] < bb) ? 1u : 0u;
      xor_pop += arbint_popcount_limb(bd);
    }
  }
  *out = max_n * (size_t) ARBINT_LIMB_BITS - xor_pop;
  return ARBINT_OK;
}

/* ========== Public API: not ========== */

/*  Two's-complement bitwise NOT: not(a) = -(a+1) for a >= 0,
    not(-m) = m - 1 for m > 0.  Single-pass magnitude inc/dec.  */
arbint_err_t arbint_not(arbint_t rop, const arbint_t a) {
  size_t an;
  arbint_err_t rc;
  const arbint_limb_t * ap;
  arbint_limb_t * rp;

  if (rop == NULL || a == NULL)
    return ARBINT_EINVAL;

  an = arbint_abs_sz(a[0]._sz);

  if (a[0]._sz >= 0) {
    /*  not(a) = -(a+1).  Result has at most an+1 limbs.  */
    size_t need = an + 1u;
    size_t used;

    if (need == 0u)
      return ARBINT_EOVERFLOW;
    if (rop[0]._cap < need) {
      rc = arbint_resize(rop, need);
      if (rc != ARBINT_OK)
        return rc;
    }
    rp = ARBINT_LIMBS(rop);
    ap = ARBINT_CLIMBS(a);

    if (an == 0u) {
      /*  not(0) = -1.  */
      rp[0] = 1u;
      rop[0]._sz = -1;
      return ARBINT_OK;
    }

    used = arbint_mag_inc(rp, ap, an);
    if (!arbint_set_signed_sz(rop, used, -1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  /*  a < 0: not(-m) = m - 1.  */
  if (an == 1u && ARBINT_CLIMBS(a)[0] == 1u) {
    /*  not(-1) = 0.  */
    arbint_zero(rop);
    return ARBINT_OK;
  }

  if (rop[0]._cap < an) {
    rc = arbint_resize(rop, an);
    if (rc != ARBINT_OK)
      return rc;
  }
  rp = ARBINT_LIMBS(rop);
  ap = ARBINT_CLIMBS(a);

  size_t used = arbint_mag_dec(rp, ap, an);
  if (used == 0u) {
    arbint_zero(rop);
    return ARBINT_OK;
  }
  if (!arbint_set_signed_sz(rop, used, 1))
    return ARBINT_EOVERFLOW;
  return ARBINT_OK;
}

/* ========== Public API: and ========== */

/*  Two's-complement bitwise AND with 4-case sign decomposition.
    For negative -x (x > 0): TC(-x) = ~(x-1).

    (+a) & (+b) = a & b                                 (positive)
    (+a) & (-b) = a & ~(b-1)                            (positive)
    (-a) & (+b) = b & ~(a-1)                            (positive)
    (-a) & (-b) = -(((a-1) | (b-1)) + 1)                (negative)  */
arbint_err_t arbint_and(arbint_t rop, const arbint_t a, const arbint_t b) {
  size_t an;
  size_t bn;
  arbint_err_t rc;
  arbint_limb_t * rp;
  const arbint_limb_t * ap;
  const arbint_limb_t * bp;
  size_t i;

  if (rop == NULL || a == NULL || b == NULL)
    return ARBINT_EINVAL;

  /*  Fast paths.  */
  if (a[0]._sz == 0 || b[0]._sz == 0) {
    arbint_zero(rop);
    return ARBINT_OK;
  }
  if (a == b)
    return arbint_set(rop, a);

  an = arbint_abs_sz(a[0]._sz);
  bn = arbint_abs_sz(b[0]._sz);

  if (a[0]._sz > 0 && b[0]._sz > 0) {
    /*  (+,+): result = a & b, bounded by min(an,bn) limbs.  */
    size_t min_n = (an < bn) ? an : bn;
    size_t used;

    if (rop[0]._cap < min_n) {
      rc = arbint_resize(rop, min_n);
      if (rc != ARBINT_OK)
        return rc;
    }
    rp = ARBINT_LIMBS(rop);
    ap = ARBINT_CLIMBS(a);
    bp = ARBINT_CLIMBS(b);

    for (i = 0u; i < min_n; ++i)
      rp[i] = ap[i] & bp[i];

    used = arbint_norm_used(rp, min_n);
    if (used == 0u)
      arbint_zero(rop);
    else if (!arbint_set_signed_sz(rop, used, 1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  if (a[0]._sz > 0 && b[0]._sz < 0) {
    /*  (+a) & (-b): result = a & ~(b-1), positive, bounded by an limbs.
        Phase 1 (i < min(an,bn)): borrow-subtract b, mask with a.
        Phase 2 (bn <= i < an): (b-1) zero-extends to 0, ~0 = all-ones,
        a[i] & all-ones = a[i], so tail is a copy.  */
    size_t phase1 = (an < bn) ? an : bn;
    size_t used;
    arbint_limb_t borrow = 1u;

    if (rop[0]._cap < an) {
      rc = arbint_resize(rop, an);
      if (rc != ARBINT_OK)
        return rc;
    }
    rp = ARBINT_LIMBS(rop);
    ap = ARBINT_CLIMBS(a);
    bp = ARBINT_CLIMBS(b);

    for (i = 0u; i < phase1; ++i) {
      arbint_limb_t bd = bp[i] - borrow;
      borrow = (bp[i] < borrow) ? 1u : 0u;
      rp[i] = ap[i] & ~bd;
    }

    /*  Beyond bn: ~(b-1) = all-ones, so result = a[i].  */
    if (rp != ap && phase1 < an)
      memcpy(rp + phase1, ap + phase1,
             (an - phase1) * sizeof(arbint_limb_t));

    used = arbint_norm_used(rp, an);
    if (used == 0u)
      arbint_zero(rop);
    else if (!arbint_set_signed_sz(rop, used, 1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  if (a[0]._sz < 0 && b[0]._sz > 0) {
    /*  (-a) & (+b): symmetric to (+b) & (-a).  */
    return arbint_and(rop, b, a);
  }

  /*  (-a) & (-b): result = -(((a-1) | (b-1)) + 1), negative.
      Phase 1 (i < min_n): both borrows active.
      Phase 2 (min_n <= i < max_n): shorter operand zero-extends,
      (shorter-1) = 0, so result = (longer-1)[i] | 0 = (longer-1)[i].
      Continue the borrow of the longer operand only.
      Intermediate needs max(an,bn) limbs, +1 for possible carry.  */
  size_t min_n = (an < bn) ? an : bn;
  size_t max_n = (an > bn) ? an : bn;
  const arbint_limb_t * longer = (an > bn) ? ARBINT_CLIMBS(a)
                                           : ARBINT_CLIMBS(b);
  size_t need = max_n + 1u;
  size_t used;
  arbint_limb_t ba = 1u;
  arbint_limb_t bb = 1u;
  arbint_limb_t * longer_borrow;

  if (rop[0]._cap < need) {
    rc = arbint_resize(rop, need);
    if (rc != ARBINT_OK)
      return rc;
  }
  rp = ARBINT_LIMBS(rop);
  ap = ARBINT_CLIMBS(a);
  bp = ARBINT_CLIMBS(b);
  longer = (an > bn) ? ap : bp;

  for (i = 0u; i < min_n; ++i) {
    arbint_limb_t ad = ap[i] - ba;
    ba = (ap[i] < ba) ? 1u : 0u;
    arbint_limb_t bd = bp[i] - bb;
    bb = (bp[i] < bb) ? 1u : 0u;
    rp[i] = ad | bd;
  }

  /*  Tail: only the longer operand contributes.  */
  longer_borrow = (an > bn) ? &ba : &bb;
  for (i = min_n; i < max_n; ++i) {
    arbint_limb_t ld = longer[i] - *longer_borrow;
    *longer_borrow = (longer[i] < *longer_borrow) ? 1u : 0u;
    rp[i] = ld;
  }

  used = arbint_mag_inc(rp, rp, max_n);
  used = arbint_norm_used(rp, used);
  if (used == 0u)
    arbint_zero(rop);
  else if (!arbint_set_signed_sz(rop, used, -1))
    return ARBINT_EOVERFLOW;
  return ARBINT_OK;
}

/* ========== Public API: or ========== */

/*  Two's-complement bitwise OR with 4-case sign decomposition.

    (+a) | (+b) = a | b                                (positive)
    (+a) | (-b) = -(((b-1) & ~a) + 1)                  (negative)
    (-a) | (+b) = -(((a-1) & ~b) + 1)                  (negative)
    (-a) | (-b) = -(((a-1) & (b-1)) + 1)               (negative)  */
arbint_err_t arbint_or(arbint_t rop, const arbint_t a, const arbint_t b) {
  size_t an;
  size_t bn;
  arbint_err_t rc;
  arbint_limb_t * rp;
  const arbint_limb_t * ap;
  const arbint_limb_t * bp;
  size_t i;

  if (rop == NULL || a == NULL || b == NULL)
    return ARBINT_EINVAL;

  /*  Fast paths.  */
  if (a[0]._sz == 0)
    return arbint_set(rop, b);
  if (b[0]._sz == 0)
    return arbint_set(rop, a);
  if (a == b)
    return arbint_set(rop, a);

  an = arbint_abs_sz(a[0]._sz);
  bn = arbint_abs_sz(b[0]._sz);

  if (a[0]._sz > 0 && b[0]._sz > 0) {
    /*  (+,+): result = a | b, max(an,bn) limbs.
        Beyond the shorter operand, x | 0 = x, so the tail is a copy.  */
    size_t min_n = (an < bn) ? an : bn;
    size_t max_n = (an > bn) ? an : bn;
    const arbint_limb_t * longer = (an > bn) ? ARBINT_CLIMBS(a)
                                              : ARBINT_CLIMBS(b);
    size_t used;

    if (rop[0]._cap < max_n) {
      rc = arbint_resize(rop, max_n);
      if (rc != ARBINT_OK)
        return rc;
    }
    rp = ARBINT_LIMBS(rop);
    ap = ARBINT_CLIMBS(a);
    bp = ARBINT_CLIMBS(b);
    longer = (an > bn) ? ap : bp;

    for (i = 0u; i < min_n; ++i)
      rp[i] = ap[i] | bp[i];

    if (rp != longer && min_n < max_n)
      memcpy(rp + min_n, longer + min_n,
             (max_n - min_n) * sizeof(arbint_limb_t));

    used = arbint_norm_used(rp, max_n);
    if (used == 0u)
      arbint_zero(rop);
    else if (!arbint_set_signed_sz(rop, used, 1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  if (a[0]._sz > 0 && b[0]._sz < 0) {
    /*  (+a) | (-b): result = -(((b-1) & ~a) + 1), negative.
        Phase 1 (i < min(an,bn)): both operands present, borrow-subtract b.
        Phase 2 (bn <= i < an): (b-1) is zero-extended, bd = 0, result = 0.
        Phase 2 (an <= i < bn): ~a is all-ones, result = (b-1)[i].
        Bounded by max(an,bn) limbs + 1 for carry.  */
    size_t min_n = (an < bn) ? an : bn;
    size_t max_n = (an > bn) ? an : bn;
    size_t need = max_n + 1u;
    size_t used;
    arbint_limb_t borrow = 1u;

    if (rop[0]._cap < need) {
      rc = arbint_resize(rop, need);
      if (rc != ARBINT_OK)
        return rc;
    }
    rp = ARBINT_LIMBS(rop);
    ap = ARBINT_CLIMBS(a);
    bp = ARBINT_CLIMBS(b);

    for (i = 0u; i < min_n; ++i) {
      arbint_limb_t bd = bp[i] - borrow;
      borrow = (bp[i] < borrow) ? 1u : 0u;
      rp[i] = bd & ~ap[i];
    }

    if (an >= bn) {
      /*  Beyond bn: (b-1) zero-extends to 0, so bd & ~a[i] = 0.  */
      memset(rp + min_n, 0, (max_n - min_n) * sizeof(arbint_limb_t));
    } else {
      /*  Beyond an: ~a is all-ones, result = (b-1)[i].
          Continue borrow propagation through remaining b limbs.  */
      for (i = min_n; i < max_n; ++i) {
        arbint_limb_t bd = bp[i] - borrow;
        borrow = (bp[i] < borrow) ? 1u : 0u;
        rp[i] = bd;
      }
    }

    used = arbint_mag_inc(rp, rp, max_n);
    used = arbint_norm_used(rp, used);
    if (used == 0u)
      arbint_zero(rop);
    else if (!arbint_set_signed_sz(rop, used, -1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  if (a[0]._sz < 0 && b[0]._sz > 0) {
    /*  (-a) | (+b): symmetric.  */
    return arbint_or(rop, b, a);
  }

  /*  (-a) | (-b): result = -(((a-1) & (b-1)) + 1), negative.
      Bounded by min(an,bn) limbs (AND shrinks) + 1 for carry.  */
  size_t min_n = (an < bn) ? an : bn;
  size_t need = min_n + 1u;
  size_t used;
  arbint_limb_t ba = 1u;
  arbint_limb_t bb = 1u;

  if (rop[0]._cap < need) {
    rc = arbint_resize(rop, need);
    if (rc != ARBINT_OK)
      return rc;
  }
  rp = ARBINT_LIMBS(rop);
  ap = ARBINT_CLIMBS(a);
  bp = ARBINT_CLIMBS(b);

  for (i = 0u; i < min_n; ++i) {
    arbint_limb_t ai = ap[i];
    arbint_limb_t bi = bp[i];
    arbint_limb_t ad = ai - ba;
    ba = (ai < ba) ? 1u : 0u;
    arbint_limb_t bd = bi - bb;
    bb = (bi < bb) ? 1u : 0u;
    rp[i] = ad & bd;
  }

  used = arbint_mag_inc(rp, rp, min_n);
  used = arbint_norm_used(rp, used);
  if (used == 0u)
    arbint_zero(rop);
  else if (!arbint_set_signed_sz(rop, used, -1))
    return ARBINT_EOVERFLOW;
  return ARBINT_OK;
}

/* ========== Public API: xor ========== */

/*  Two's-complement bitwise XOR with 4-case sign decomposition.

    (+a) ^ (+b) = a ^ b                                 (positive)
    (+a) ^ (-b) = -((a ^ (b-1)) + 1)                   (negative)
    (-a) ^ (+b) = -(((a-1) ^ b) + 1)                   (negative)
    (-a) ^ (-b) = (a-1) ^ (b-1)                         (positive)  */
arbint_err_t arbint_xor(arbint_t rop, const arbint_t a, const arbint_t b) {
  size_t an;
  size_t bn;
  arbint_err_t rc;
  arbint_limb_t * rp;
  const arbint_limb_t * ap;
  const arbint_limb_t * bp;
  size_t i;

  if (rop == NULL || a == NULL || b == NULL)
    return ARBINT_EINVAL;

  /*  Fast paths.  */
  if (a[0]._sz == 0)
    return arbint_set(rop, b);
  if (b[0]._sz == 0)
    return arbint_set(rop, a);
  if (a == b) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  an = arbint_abs_sz(a[0]._sz);
  bn = arbint_abs_sz(b[0]._sz);

  if (a[0]._sz > 0 && b[0]._sz > 0) {
    /*  (+,+): result = a ^ b, max(an,bn) limbs.
        Beyond the shorter operand, x ^ 0 = x, so the tail is a copy.  */
    size_t min_n = (an < bn) ? an : bn;
    size_t max_n = (an > bn) ? an : bn;
    const arbint_limb_t * longer;
    size_t used;

    if (rop[0]._cap < max_n) {
      rc = arbint_resize(rop, max_n);
      if (rc != ARBINT_OK)
        return rc;
    }
    rp = ARBINT_LIMBS(rop);
    ap = ARBINT_CLIMBS(a);
    bp = ARBINT_CLIMBS(b);
    longer = (an > bn) ? ap : bp;

    for (i = 0u; i < min_n; ++i)
      rp[i] = ap[i] ^ bp[i];

    if (rp != longer && min_n < max_n)
      memcpy(rp + min_n, longer + min_n,
             (max_n - min_n) * sizeof(arbint_limb_t));

    used = arbint_norm_used(rp, max_n);
    if (used == 0u)
      arbint_zero(rop);
    else if (!arbint_set_signed_sz(rop, used, 1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  if (a[0]._sz > 0 && b[0]._sz < 0) {
    /*  (+a) ^ (-b): result = -((a ^ (b-1)) + 1), negative.
        Phase 1 (i < min(an,bn)): both operands present.
        Phase 2 (bn <= i < an): (b-1) zero-extends, a[i] ^ 0 = a[i], copy.
        Phase 2 (an <= i < bn): a zero-extends, 0 ^ (b-1)[i] = (b-1)[i],
        continue borrow.
        Max(an,bn) limbs + 1 for carry.  */
    size_t min_n = (an < bn) ? an : bn;
    size_t max_n = (an > bn) ? an : bn;
    size_t need = max_n + 1u;
    size_t used;
    arbint_limb_t borrow = 1u;

    if (rop[0]._cap < need) {
      rc = arbint_resize(rop, need);
      if (rc != ARBINT_OK)
        return rc;
    }
    rp = ARBINT_LIMBS(rop);
    ap = ARBINT_CLIMBS(a);
    bp = ARBINT_CLIMBS(b);

    for (i = 0u; i < min_n; ++i) {
      arbint_limb_t bd = bp[i] - borrow;
      borrow = (bp[i] < borrow) ? 1u : 0u;
      rp[i] = ap[i] ^ bd;
    }

    if (an >= bn) {
      /*  Beyond bn: (b-1) zero-extends, a[i] ^ 0 = a[i].  */
      if (rp != ap && min_n < max_n)
        memcpy(rp + min_n, ap + min_n,
               (max_n - min_n) * sizeof(arbint_limb_t));
    } else {
      /*  Beyond an: a zero-extends, 0 ^ (b-1)[i] = (b-1)[i].  */
      for (i = min_n; i < max_n; ++i) {
        arbint_limb_t bd = bp[i] - borrow;
        borrow = (bp[i] < borrow) ? 1u : 0u;
        rp[i] = bd;
      }
    }

    used = arbint_mag_inc(rp, rp, max_n);
    used = arbint_norm_used(rp, used);
    if (used == 0u)
      arbint_zero(rop);
    else if (!arbint_set_signed_sz(rop, used, -1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  if (a[0]._sz < 0 && b[0]._sz > 0) {
    /*  (-a) ^ (+b): symmetric.  */
    return arbint_xor(rop, b, a);
  }

  /*  (-a) ^ (-b): result = (a-1) ^ (b-1), positive.
      Phase 1 (i < min_n): both borrows active.
      Phase 2 (min_n <= i < max_n): shorter zero-extends,
      (longer-1)[i] ^ 0 = (longer-1)[i], continue longer borrow.
      Max(an,bn) limbs.  */
  size_t min_n = (an < bn) ? an : bn;
  size_t max_n = (an > bn) ? an : bn;
  const arbint_limb_t * longer;
  size_t used;
  arbint_limb_t ba = 1u;
  arbint_limb_t bb = 1u;
  arbint_limb_t * longer_borrow;

  if (rop[0]._cap < max_n) {
    rc = arbint_resize(rop, max_n);
    if (rc != ARBINT_OK)
      return rc;
  }
  rp = ARBINT_LIMBS(rop);
  ap = ARBINT_CLIMBS(a);
  bp = ARBINT_CLIMBS(b);
  longer = (an > bn) ? ap : bp;

  for (i = 0u; i < min_n; ++i) {
    arbint_limb_t ad = ap[i] - ba;
    ba = (ap[i] < ba) ? 1u : 0u;
    {
      arbint_limb_t bd = bp[i] - bb;
      bb = (bp[i] < bb) ? 1u : 0u;
      rp[i] = ad ^ bd;
    }
  }

  /*  Tail: only the longer operand contributes.  */
  longer_borrow = (an > bn) ? &ba : &bb;
  for (i = min_n; i < max_n; ++i) {
    arbint_limb_t ld = longer[i] - *longer_borrow;
    *longer_borrow = (longer[i] < *longer_borrow) ? 1u : 0u;
    rp[i] = ld;
  }

  used = arbint_norm_used(rp, max_n);
  if (used == 0u)
    arbint_zero(rop);
  else if (!arbint_set_signed_sz(rop, used, 1))
    return ARBINT_EOVERFLOW;
  return ARBINT_OK;
}

/* ========== Public API: immediate u32/i32 variants ========== */

/*  Sign-extend an int32_t to a full arbint_limb_t (TC representation).
    On 64-bit limbs, (uint32_t)b zero-extends; we need sign extension so
    that ~bmask has zeroes (not ones) in the upper 32 bits.  */
#if ARBINT_LIMB_BITS == 64
  #define ARBINT_I32_TO_LIMB(b) ((arbint_limb_t)(int64_t)(int32_t)(b))
#else
  #define ARBINT_I32_TO_LIMB(b) ((arbint_limb_t)(uint32_t)(b))
#endif

/*  arbint_and_u32: rop = a & b, where b is a uint32_t (single limb).
    (+a) & b: result = a[0] & b, single limb.
    (-a) & b: TC(-a)[0] & b = ~(a[0]-1) & b  (if a[0] != 0)
                            = 0              (if a[0] == 0, borrow
                                              propagates).  */
arbint_err_t arbint_and_u32(arbint_t rop, const arbint_t a, uint32_t b) {
  arbint_err_t rc;

  if (rop == NULL || a == NULL)
    return ARBINT_EINVAL;

  if (b == 0u || a[0]._sz == 0) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  if (rop[0]._cap < 1u) {
    rc = arbint_resize(rop, 1u);
    if (rc != ARBINT_OK)
      return rc;
  }

  if (a[0]._sz > 0) {
    arbint_limb_t r = ARBINT_CLIMBS(a)[0] & (arbint_limb_t) b;
    if (r == 0u) {
      arbint_zero(rop);
      return ARBINT_OK;
    }
    ARBINT_LIMBS(rop)[0] = r;
    if (!arbint_set_signed_sz(rop, 1u, 1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  /*  Negative a: TC(-a)[0] = ~(a[0]-1) if a[0] != 0, else 0 (borrow).  */
  arbint_limb_t a0 = ARBINT_CLIMBS(a)[0];
  arbint_limb_t tc0 = (a0 != 0u) ? ~(a0 - 1u) : 0u;
  arbint_limb_t r = tc0 & (arbint_limb_t) b;
  if (r == 0u) {
    arbint_zero(rop);
    return ARBINT_OK;
  }
  ARBINT_LIMBS(rop)[0] = r;
  if (!arbint_set_signed_sz(rop, 1u, 1))
    return ARBINT_EOVERFLOW;
  return ARBINT_OK;
}

/*  arbint_and_i32: rop = a & b, where b is an int32_t.
    b >= 0: delegate to and_u32.
    b < 0: TC(b) = ...FFFF_FFFF_xxxx, so AND preserves upper limbs of a.
    (+a) & (-b): mask low limb by (uint32_t)b, keep upper limbs.
    (-a) & (-b): = -(((a-1) | ~TC(b)) + 1).
                 ~TC(b)[0] = ~(uint32_t)b, upper = 0.
                 So: (a-1)[0] | ~(uint32_t)b, upper = (a-1)[i]. Then +1.  */
arbint_err_t arbint_and_i32(arbint_t rop, const arbint_t a, int32_t b) {
  size_t an;
  arbint_err_t rc;
  arbint_limb_t bmask;

  if (rop == NULL || a == NULL)
    return ARBINT_EINVAL;

  if (b >= 0)
    return arbint_and_u32(rop, a, (uint32_t) b);

  /*  b < 0: TC(b) has all-ones above 32 bits.  */
  if (a[0]._sz == 0) {
    arbint_zero(rop);
    return ARBINT_OK;
  }

  bmask = ARBINT_I32_TO_LIMB(b);
  an = arbint_abs_sz(a[0]._sz);

  if (a[0]._sz > 0) {
    /*  (+a) & (-b): mask low limb by bmask, keep upper limbs.  */
    size_t used;

    if (rop[0]._cap < an) {
      rc = arbint_resize(rop, an);
      if (rc != ARBINT_OK)
        return rc;
    }

    const arbint_limb_t * ap = ARBINT_CLIMBS(a);
    arbint_limb_t * rp = ARBINT_LIMBS(rop);
    size_t i;

    rp[0] = ap[0] & bmask;
    for (i = 1u; i < an; ++i)
      rp[i] = ap[i];

    used = arbint_norm_used(rp, an);
    if (used == 0u)
      arbint_zero(rop);
    else if (!arbint_set_signed_sz(rop, used, 1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  /*  (-a) & (-b): result = -(((a-1) | ~TC(b)) + 1).
      ~TC(b) has low limb = ~bmask, upper limbs = 0.
      So: result[0] = (a-1)[0] | ~bmask, result[i>0] = (a-1)[i].
      Then +1 for final magnitude.  */
  const arbint_limb_t * ap;
  arbint_limb_t * rp;
  size_t need = an + 1u;
  size_t k;
  size_t i;
  size_t used;

  if (rop[0]._cap < need) {
    rc = arbint_resize(rop, need);
    if (rc != ARBINT_OK)
      return rc;
  }

  ap = ARBINT_CLIMBS(a);
  rp = ARBINT_LIMBS(rop);

  /*  Compute (a-1).  */
  k = arbint_first_nonzero(ap, an);
  for (i = 0u; i < k; ++i)
    rp[i] = ~(arbint_limb_t) 0u;
  rp[k] = ap[k] - 1u;
  for (i = k + 1u; i < an; ++i)
    rp[i] = ap[i];

  /*  OR with ~bmask in low limb only.  */
  rp[0] |= ~bmask;

  /*  +1 for magnitude.  */
  used = arbint_mag_inc(rp, rp, an);
  used = arbint_norm_used(rp, used);
  if (used == 0u)
    arbint_zero(rop);
  else if (!arbint_set_signed_sz(rop, used, -1))
    return ARBINT_EOVERFLOW;
  return ARBINT_OK;
}

/*  arbint_or_u32: rop = a | b, where b is a uint32_t.
    (+a) | b: copy a, OR b into low limb.
    (-a) | b: result = -(((a-1) & ~b_ext) + 1).
              ~b_ext: limb 0 = ~b, upper = all-ones.
              (a-1) & ~b_ext: limb 0 = (a-1)[0] & ~b; upper = (a-1)[i].  */
arbint_err_t arbint_or_u32(arbint_t rop, const arbint_t a, uint32_t b) {
  size_t an;
  arbint_err_t rc;

  if (rop == NULL || a == NULL)
    return ARBINT_EINVAL;

  if (a[0]._sz == 0)
    return arbint_set_u32(rop, b);

  if (b == 0u)
    return arbint_set(rop, a);

  an = arbint_abs_sz(a[0]._sz);

  if (a[0]._sz > 0) {
    /*  (+a) | b: copy a, OR b into low limb.  */
    size_t used;

    if (rop[0]._cap < an) {
      rc = arbint_resize(rop, an);
      if (rc != ARBINT_OK)
        return rc;
    }

    const arbint_limb_t * ap = ARBINT_CLIMBS(a);
    arbint_limb_t * rp = ARBINT_LIMBS(rop);
    size_t i;

    rp[0] = ap[0] | (arbint_limb_t) b;
    for (i = 1u; i < an; ++i)
      rp[i] = ap[i];

    used = arbint_norm_used(rp, an);
    if (used == 0u)
      arbint_zero(rop);
    else if (!arbint_set_signed_sz(rop, used, 1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  /*  (-a) | b: result = -(((a-1) & ~b_ext) + 1).
      (a-1)[0] & ~(arbint_limb_t)b; upper limbs = (a-1)[i].  */
  const arbint_limb_t * ap;
  arbint_limb_t * rp;
  size_t need = an + 1u;
  size_t k;
  size_t i;
  size_t used;

  if (rop[0]._cap < need) {
    rc = arbint_resize(rop, need);
    if (rc != ARBINT_OK)
      return rc;
  }

  ap = ARBINT_CLIMBS(a);
  rp = ARBINT_LIMBS(rop);

  /*  Compute (a-1).  */
  k = arbint_first_nonzero(ap, an);
  for (i = 0u; i < k; ++i)
    rp[i] = ~(arbint_limb_t) 0u;
  rp[k] = ap[k] - 1u;
  for (i = k + 1u; i < an; ++i)
    rp[i] = ap[i];

  /*  AND with ~b in low limb only (upper limbs: & all-ones = identity).  */
  rp[0] &= ~(arbint_limb_t) b;

  /*  +1 for magnitude.  */
  used = arbint_mag_inc(rp, rp, an);
  used = arbint_norm_used(rp, used);
  if (used == 0u)
    arbint_zero(rop);
  else if (!arbint_set_signed_sz(rop, used, -1))
    return ARBINT_EOVERFLOW;
  return ARBINT_OK;
}

/*  arbint_or_i32: rop = a | b, where b is an int32_t.
    b >= 0: delegate to or_u32.
    b < 0: TC(b) = ...FFFF_xxxx.  Result is always negative.
    (+a) | (-b): = -(((b_mag-1) & ~a) + 1).
                 b_mag-1 = ~bmask.  Limb 0 = ~bmask & ~a[0]; upper = 0.
    (-a) | (-b): = -(((a-1) & (b_mag-1)) + 1).
                 b_mag-1 = ~bmask.  Limb 0 = (a-1)[0] & ~bmask; upper = 0.  */
arbint_err_t arbint_or_i32(arbint_t rop, const arbint_t a, int32_t b) {
  size_t an;
  arbint_err_t rc;
  arbint_limb_t bmask;

  if (rop == NULL || a == NULL)
    return ARBINT_EINVAL;

  if (b >= 0)
    return arbint_or_u32(rop, a, (uint32_t) b);

  /*  b < 0: TC(b) = ...FFFF_xxxx.  Result is always negative.  */
  bmask = ARBINT_I32_TO_LIMB(b);
  an = arbint_abs_sz(a[0]._sz);

  if (a[0]._sz >= 0) {
    /*  (+a) | (-b): = -(((b_mag-1) & ~a) + 1).
        b_mag = (uint32_t)(-b), b_mag-1 = ~bmask.
        Limb 0 = ~bmask & ~a[0]; upper limbs = 0 & ~a[i] = 0.  */
    arbint_limb_t a0 = (an > 0u) ? ARBINT_CLIMBS(a)[0] : 0u;
    arbint_limb_t r0 = ~bmask & ~a0;
    size_t used;

    if (rop[0]._cap < 2u) {
      rc = arbint_resize(rop, 2u);
      if (rc != ARBINT_OK)
        return rc;
    }

    ARBINT_LIMBS(rop)[0] = r0;
    used = arbint_mag_inc(ARBINT_LIMBS(rop), ARBINT_LIMBS(rop), 1u);
    used = arbint_norm_used(ARBINT_LIMBS(rop), used);
    if (used == 0u)
      arbint_zero(rop);
    else if (!arbint_set_signed_sz(rop, used, -1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  /*  (-a) | (-b): = -(((a-1) & (b_mag-1)) + 1).
      b_mag-1 = ~bmask.  Only limb 0 matters (upper limbs of ~bmask = 0).
      result = -(((a-1)[0] & ~bmask) + 1).  */
  const arbint_limb_t * ap = ARBINT_CLIMBS(a);
  size_t k;
  arbint_limb_t am1_0;
  arbint_limb_t r0;
  size_t used;

  if (rop[0]._cap < 2u) {
    rc = arbint_resize(rop, 2u);
    if (rc != ARBINT_OK)
      return rc;
  }

  /*  Compute (a-1)[0].  */
  k = arbint_first_nonzero(ap, an);
  if (k == 0u)
    am1_0 = ap[0] - 1u;
  else
    am1_0 = ~(arbint_limb_t) 0u;

  r0 = am1_0 & ~bmask;
  ARBINT_LIMBS(rop)[0] = r0;
  used = arbint_mag_inc(ARBINT_LIMBS(rop), ARBINT_LIMBS(rop), 1u);
  used = arbint_norm_used(ARBINT_LIMBS(rop), used);
  if (used == 0u)
    arbint_zero(rop);
  else if (!arbint_set_signed_sz(rop, used, -1))
    return ARBINT_EOVERFLOW;
  return ARBINT_OK;
}

/*  arbint_xor_u32: rop = a ^ b, where b is a uint32_t.
    (+a) ^ b: copy a, XOR b into low limb.
    (-a) ^ b: = -(((a-1) ^ b_ext) + 1).
              b_ext = b in limb 0, 0 above.
              Limb 0 = (a-1)[0] ^ b; upper = (a-1)[i].  */
arbint_err_t arbint_xor_u32(arbint_t rop, const arbint_t a, uint32_t b) {
  size_t an;
  arbint_err_t rc;

  if (rop == NULL || a == NULL)
    return ARBINT_EINVAL;

  if (a[0]._sz == 0)
    return arbint_set_u32(rop, b);

  if (b == 0u)
    return arbint_set(rop, a);

  an = arbint_abs_sz(a[0]._sz);

  if (a[0]._sz > 0) {
    /*  (+a) ^ b: copy a, XOR b into low limb.  */
    size_t used;

    if (rop[0]._cap < an) {
      rc = arbint_resize(rop, an);
      if (rc != ARBINT_OK)
        return rc;
    }

    const arbint_limb_t * ap = ARBINT_CLIMBS(a);
    arbint_limb_t * rp = ARBINT_LIMBS(rop);
    size_t i;

    rp[0] = ap[0] ^ (arbint_limb_t) b;
    for (i = 1u; i < an; ++i)
      rp[i] = ap[i];

    used = arbint_norm_used(rp, an);
    if (used == 0u)
      arbint_zero(rop);
    else if (!arbint_set_signed_sz(rop, used, 1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  /*  (-a) ^ b: result = -(((a-1) ^ b_ext) + 1).
      Limb 0 = (a-1)[0] ^ b; upper = (a-1)[i]. Then +1.  */
  const arbint_limb_t * ap;
  arbint_limb_t * rp;
  size_t need = an + 1u;
  size_t k;
  size_t i;
  size_t used;

  if (rop[0]._cap < need) {
    rc = arbint_resize(rop, need);
    if (rc != ARBINT_OK)
      return rc;
  }

  ap = ARBINT_CLIMBS(a);
  rp = ARBINT_LIMBS(rop);

  /*  Compute (a-1).  */
  k = arbint_first_nonzero(ap, an);
  for (i = 0u; i < k; ++i)
    rp[i] = ~(arbint_limb_t) 0u;
  rp[k] = ap[k] - 1u;
  for (i = k + 1u; i < an; ++i)
    rp[i] = ap[i];

  /*  XOR b into low limb.  */
  rp[0] ^= (arbint_limb_t) b;

  /*  +1 for magnitude.  */
  used = arbint_mag_inc(rp, rp, an);
  used = arbint_norm_used(rp, used);
  if (used == 0u)
    arbint_zero(rop);
  else if (!arbint_set_signed_sz(rop, used, -1))
    return ARBINT_EOVERFLOW;
  return ARBINT_OK;
}

/*  arbint_xor_i32: rop = a ^ b, where b is an int32_t.
    b >= 0: delegate to xor_u32.
    b < 0: TC(b) = ...FFFF_xxxx.  XOR with all-ones above complements.
    (+a) ^ (-b): = -((a ^ ~bmask_ext) + 1).
                 Limb 0 = a[0] ^ ~bmask; upper = ~a[i].
    (-a) ^ (-b): = (a-1) ^ ~bmask_ext, positive.
                 Limb 0 = (a-1)[0] ^ ~bmask; upper = ~(a-1)[i].  */
arbint_err_t arbint_xor_i32(arbint_t rop, const arbint_t a, int32_t b) {
  size_t an;
  arbint_err_t rc;
  arbint_limb_t bmask;

  if (rop == NULL || a == NULL)
    return ARBINT_EINVAL;

  if (b >= 0)
    return arbint_xor_u32(rop, a, (uint32_t) b);

  /*  b < 0: TC(b) = ...FFFF_xxxx.  */
  bmask = ARBINT_I32_TO_LIMB(b);
  an = arbint_abs_sz(a[0]._sz);

  if (a[0]._sz == 0)
    return arbint_set_i32(rop, b);

  if (a[0]._sz > 0) {
    /*  (+a) ^ (-b): result = -(a ^ (b_mag-1) + 1).
        (b_mag-1) occupies only limb 0 = ~bmask; upper = 0.
        So: limb 0 = a[0] ^ ~bmask, upper = a[i].  */
    size_t need = an + 1u;
    size_t used;

    if (rop[0]._cap < need) {
      rc = arbint_resize(rop, need);
      if (rc != ARBINT_OK)
        return rc;
    }

    const arbint_limb_t * ap = ARBINT_CLIMBS(a);
    arbint_limb_t * rp = ARBINT_LIMBS(rop);
    size_t i;

    rp[0] = ap[0] ^ ~bmask;
    for (i = 1u; i < an; ++i)
      rp[i] = ap[i];

    used = arbint_mag_inc(rp, rp, an);
    used = arbint_norm_used(rp, used);
    if (used == 0u)
      arbint_zero(rop);
    else if (!arbint_set_signed_sz(rop, used, -1))
      return ARBINT_EOVERFLOW;
    return ARBINT_OK;
  }

  /*  (-a) ^ (-b): result = (a-1) ^ (b_mag-1), positive.
      (b_mag-1) occupies only limb 0 = ~bmask; upper = 0.
      So: limb 0 = (a-1)[0] ^ ~bmask, upper = (a-1)[i].  */
  const arbint_limb_t * ap;
  arbint_limb_t * rp;
  size_t k;
  size_t i;
  size_t used;

  if (rop[0]._cap < an) {
    rc = arbint_resize(rop, an);
    if (rc != ARBINT_OK)
      return rc;
  }

  ap = ARBINT_CLIMBS(a);
  rp = ARBINT_LIMBS(rop);

  /*  Compute (a-1).  */
  k = arbint_first_nonzero(ap, an);
  for (i = 0u; i < k; ++i)
    rp[i] = ~(arbint_limb_t) 0u;
  rp[k] = ap[k] - 1u;
  for (i = k + 1u; i < an; ++i)
    rp[i] = ap[i];

  /*  XOR (b_mag-1) into limb 0 only; upper limbs unchanged.  */
  rp[0] ^= ~bmask;

  used = arbint_norm_used(rp, an);
  if (used == 0u)
    arbint_zero(rop);
  else if (!arbint_set_signed_sz(rop, used, 1))
    return ARBINT_EOVERFLOW;
  return ARBINT_OK;
}
