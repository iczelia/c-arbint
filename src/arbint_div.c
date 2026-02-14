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

#include "arbint_div.h"
#include "arbint_div_newton.h"

#include "config.h"

#include "arbint_cpu.h"
#include "arbint_internal_util.h"

#if ARBINT_COMPILER_MSVC
  #include <malloc.h>
#else
  #include <alloca.h>
#endif

#include <limits.h>
#include <string.h>

#ifndef ARBINT_TDIV_STACK_REM_LIMBS_MAX
#define ARBINT_TDIV_STACK_REM_LIMBS_MAX 8u
#endif

/*  Power-of-two truncated division: q = n / 2^k, r = n % 2^k.
    Truncated division semantics:
      - Quotient: magnitude right-shift, sign preserved
      - Remainder: low k bits of magnitude, sign of dividend
    Examples:
       7 / 4 =  1,  7 % 4 =  3
      -7 / 4 = -1, -7 % 4 = -3  */

static arbint_err_t arbint_tdiv_qr_pow2(arbint_t q, arbint_t r,
                                        const arbint_t n, unsigned k,
                                        int dsign) {
  size_t nn;
  const arbint_limb_t * np;
  size_t n_ctz;
  int nsign;
  int qsign;
  arbint_err_t rc;

  if ((q == NULL && r == NULL) || n == NULL)
    return ARBINT_EINVAL;

  nsign = arbint_signum(n);

  /*  n == 0: q = 0, r = 0.  */
  if (nsign == 0) {
    if (q != NULL)
      arbint_zero(q);
    if (r != NULL)
      arbint_zero(r);
    return ARBINT_OK;
  }

  nn = arbint_abs_sz(n[0]._sz);
  np = ARBINT_CLIMBS(n);

  /*  Quotient sign: same sign if signs match, opposite if differ.  */
  qsign = (nsign == dsign) ? 1 : -1;

  /*  Optimization: count trailing zeros in dividend.
      If n has >= k trailing zero bits, the remainder is 0.  */
  n_ctz = arbint_mag_ctz_or_size_max(np, nn);

  if (n_ctz >= (size_t) k) {
    /*  Dividend is divisible by 2^k: remainder is 0.  */
    if (r != NULL)
      arbint_zero(r);

    /*  Quotient: n >> k.
        Since n is divisible by 2^k, we can skip the low zero limbs.  */
    if (q != NULL) {
      size_t limb_skip = (size_t) (k / ARBINT_LIMB_BITS);
      unsigned bit_shift = (unsigned) (k % ARBINT_LIMB_BITS);
      size_t src_limbs = nn - limb_skip;
      size_t q_used;

      if (src_limbs == 0u) {
        arbint_zero(q);
      } else {
        arbint_limb_t * qp;
        size_t i;

        rc = arbint_resize(q, src_limbs);
        if (rc != ARBINT_OK)
          return rc;

        qp = ARBINT_LIMBS(q);

        if (bit_shift == 0u) {
          /*  Pure limb shift: just copy.  */
          for (i = 0u; i < src_limbs; ++i)
            qp[i] = np[limb_skip + i];
        } else {
          /*  Shift within limbs.  */
          for (i = 0u; i < src_limbs - 1u; ++i)
            qp[i] = (np[limb_skip + i] >> bit_shift) |
                    (np[limb_skip + i + 1u] << (ARBINT_LIMB_BITS - bit_shift));
          qp[src_limbs - 1u] = np[nn - 1u] >> bit_shift;
        }

        q_used = arbint_norm_used(qp, src_limbs);
        if (q_used == 0u) {
          arbint_zero(q);
        } else {
          if (!arbint_set_signed_sz(q, q_used, qsign))
            return ARBINT_EOVERFLOW;
        }
      }
    }
    return ARBINT_OK;
  }

  /*  General case: n has fewer than k trailing zeros.

      Aliasing considerations:
      - If q == n: arbint_shr handles this internally (shifts in place)
      - If r == n: writing r destroys n, so compute q first
      - If q == r: per aliasing contract, last write wins (r); compute q first

      Strategy: compute q first whenever r might overwrite something q needs,
      or when q == r (so r is written last).  */

  if (q != NULL && (r == n || q == r)) {
    /*  Compute q first:
        - r == n: prevents r from destroying n before we read it
        - q == r: ensures r (written later) wins per aliasing contract  */
    rc = arbint_shr(q, n, (uint32_t) k);
    if (rc != ARBINT_OK)
      return rc;
    if (!arbint_is_zero(q) && nsign != qsign)
      q[0]._sz = -q[0]._sz;
  }

  /*  Compute remainder: low k bits of magnitude with n's sign.  */
  if (r != NULL) {
    size_t limb_idx = (size_t) (k / ARBINT_LIMB_BITS);
    unsigned bit_idx = (unsigned) (k % ARBINT_LIMB_BITS);

    if (limb_idx >= nn) {
      /*  k >= total bits: remainder is n itself (magnitude).  */
      rc = arbint_set(r, n);
      if (rc != ARBINT_OK)
        return rc;
      /*  Sign of remainder = sign of dividend (already set by arbint_set).  */
    } else {
      /*  Extract low k bits into r.  */
      size_t i;
      size_t r_limbs = limb_idx + (bit_idx != 0u ? 1u : 0u);
      size_t used;

      if (r_limbs == 0u) {
        /*  k == 0: remainder is 0 (division by 1).  */
        arbint_zero(r);
      } else {
        rc = arbint_resize(r, r_limbs);
        if (rc != ARBINT_OK)
          return rc;

        arbint_limb_t * rp = ARBINT_LIMBS(r);

        /*  Copy full limbs.  */
        for (i = 0u; i < limb_idx && i < nn; ++i)
          rp[i] = np[i];

        /*  Mask partial top limb if needed.  */
        if (bit_idx != 0u && limb_idx < nn) {
          arbint_limb_t mask = ((arbint_limb_t) 1u << bit_idx) - 1u;
          rp[limb_idx] = np[limb_idx] & mask;
        }

        used = arbint_norm_used(rp, r_limbs);
        if (used == 0u) {
          arbint_zero(r);
        } else {
          /*  Remainder sign = dividend sign.  */
          if (!arbint_set_signed_sz(r, used, nsign))
            return ARBINT_EOVERFLOW;
        }
      }
    }
  }

  /*  Compute quotient if we haven't already.
      Skip if: r == n (computed above) or q == r (computed above).  */
  if (q != NULL && r != n && q != r) {
    rc = arbint_shr(q, n, (uint32_t) k);
    if (rc != ARBINT_OK)
      return rc;

    /*  arbint_shr preserves n's sign, but we need qsign.
        If nsign != qsign, flip the sign.  */
    if (!arbint_is_zero(q) && nsign != qsign)
      q[0]._sz = -q[0]._sz;
  }

  return ARBINT_OK;
}

#define ARBINT_DIV_PREPARE_FN arbint_div_prepare_barrett_limb_core
#define ARBINT_DIV_OMIT_STEP 1
#include "arbint_div_barrett_core.inc"

arbint_limb_t arbint_div_prepare_barrett_limb(arbint_limb_t d_norm) {
  return arbint_div_prepare_barrett_limb_core(d_norm);
}

/*  Function pointer types for runtime dispatch.  */

typedef arbint_err_t (*arbint_div_qr_u32_impl_fn_t)(arbint_t q, arbint_t r,
                                                    const arbint_t n,
                                                    uint32_t dmag, int dsign);

typedef arbint_err_t (*arbint_div_mag_single_limb_fn_t)(
    const arbint_limb_t * np, size_t nn, arbint_limb_t d_limb,
    arbint_limb_t * qp, size_t * q_used, arbint_limb_t * rem_out);

typedef arbint_err_t (*arbint_mod_u32_barrett_fn_t)(arbint_t x,
                                                    arbint_limb_t d_norm,
                                                    arbint_limb_t di,
                                                    unsigned shift);

typedef arbint_err_t (*arbint_tdiv_q_3_fn_t)(arbint_t q, const arbint_t n);
typedef arbint_err_t (*arbint_div_mag_two_limb_fn_t)(
    const arbint_limb_t * np, size_t nn, const arbint_limb_t * dp,
    arbint_limb_t * qp, arbint_limb_t * rp);

/*  Runtime dispatch selectors.  */

static arbint_div_qr_u32_impl_fn_t arbint_select_div_qr_u32_impl(void) {
#if HAS_BMI2_ALWAYS
  return arbint_div_qr_u32_bmi2_impl;
#elif HAS_BMI2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2)
             ? arbint_div_qr_u32_bmi2_impl
             : arbint_div_qr_u32_generic_impl;
#else
  return arbint_div_qr_u32_generic_impl;
#endif
}

static arbint_div_mag_single_limb_fn_t
arbint_select_div_mag_single_limb(void) {
#if HAS_BMI2_ALWAYS
  return arbint_div_mag_single_limb_bmi2;
#elif HAS_BMI2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2)
             ? arbint_div_mag_single_limb_bmi2
             : arbint_div_mag_single_limb_generic;
#else
  return arbint_div_mag_single_limb_generic;
#endif
}

static arbint_mod_u32_barrett_fn_t arbint_select_mod_u32_barrett(void) {
#if HAS_BMI2_ALWAYS
  return arbint_mod_u32_barrett_bmi2;
#elif HAS_BMI2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2)
             ? arbint_mod_u32_barrett_bmi2
             : arbint_mod_u32_barrett_generic;
#else
  return arbint_mod_u32_barrett_generic;
#endif
}

static arbint_tdiv_q_3_fn_t arbint_select_tdiv_q_3(void) {
#if HAS_BMI2_ALWAYS
  return arbint_tdiv_q_3_bmi2;
#elif HAS_BMI2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2)
             ? arbint_tdiv_q_3_bmi2
             : arbint_tdiv_q_3_generic;
#else
  return arbint_tdiv_q_3_generic;
#endif
}

static arbint_div_mag_two_limb_fn_t arbint_select_div_mag_two_limb(void) {
#if HAS_BMI2_ALWAYS
  return arbint_div_mag_two_limb_bmi2;
#elif HAS_BMI2
  return arbint_cpu_has_feature(ARBINT_CPU_FEATURE_BMI2)
             ? arbint_div_mag_two_limb_bmi2
             : NULL;
#else
  return NULL;
#endif
}

/*  Dispatched Barrett reduction.  */

arbint_err_t arbint_mod_u32_barrett(arbint_t x, arbint_limb_t d_norm,
                                    arbint_limb_t di, unsigned shift) {
  static arbint_mod_u32_barrett_fn_t impl = NULL;

  ARBINT_LAZY_INIT(impl, arbint_select_mod_u32_barrett);

  return impl(x, d_norm, di, shift);
}

arbint_err_t arbint_tdiv_q_3_dispatch(arbint_t q, const arbint_t n) {
  static arbint_tdiv_q_3_fn_t impl = NULL;

  ARBINT_LAZY_INIT(impl, arbint_select_tdiv_q_3);

  return impl(q, n);
}

/*  Dispatched single-limb remainder (no quotient).  */

arbint_err_t arbint_mod_mag_single_limb(const arbint_limb_t * np, size_t nn,
                                        arbint_limb_t d_limb,
                                        arbint_limb_t * rem_out) {
  static arbint_div_mag_single_limb_fn_t impl = NULL;

  ARBINT_LAZY_INIT(impl, arbint_select_div_mag_single_limb);

  return impl(np, nn, d_limb, NULL, NULL, rem_out);
}

/*  Dispatched u32 division.  */

arbint_err_t arbint_div_qr_u32_dispatch(arbint_t q, arbint_t r,
                                        const arbint_t n, uint32_t dmag,
                                        int dsign) {
  static arbint_div_qr_u32_impl_fn_t impl = NULL;

  /*  Power-of-two fast path: use shift/mask instead of reciprocal.  */
  if (arbint_u32_is_pow2(dmag))
    return arbint_tdiv_qr_pow2(q, r, n, (unsigned) arbint_ctz_limb(dmag),
                               dsign);

  ARBINT_LAZY_INIT(impl, arbint_select_div_qr_u32_impl);

  return impl(q, r, n, dmag, dsign);
}

/*  Magnitude view and signed assignment.  */

arbint_err_t arbint_get_mag_view(const arbint_t x, const arbint_limb_t ** xp,
                                 size_t * xn, int * sign) {
  size_t used;

  if (x == NULL || xp == NULL || xn == NULL || sign == NULL)
    return ARBINT_EINVAL;

  used = arbint_abs_sz(x[0]._sz);
  if (used == 0u) {
    *xp = NULL;
    *xn = 0u;
    *sign = 0;
    return ARBINT_OK;
  }

  if (x[0]._ptr == NULL || used > x[0]._cap)
    return ARBINT_EINVAL;

  *xp = ARBINT_CLIMBS(x);
  *xn = arbint_norm_used(*xp, used);
  *sign = (x[0]._sz < 0) ? -1 : 1;
  if (*xn == 0u)
    *sign = 0;
  return ARBINT_OK;
}

static arbint_err_t arbint_set_mag_signed(arbint_t x,
                                          const arbint_limb_t * mag,
                                          size_t used, int sign) {
  arbint_err_t rc;

  if (x == NULL)
    return ARBINT_OK;

  if (used == 0u) {
    arbint_zero(x);
    return ARBINT_OK;
  }

  rc = arbint_resize(x, used);
  if (rc != ARBINT_OK)
    return rc;

  memcpy(ARBINT_LIMBS(x), mag, used * sizeof(arbint_limb_t));
  if (!arbint_set_signed_sz(x, used, sign))
    return ARBINT_EOVERFLOW;

  return ARBINT_OK;
}

/*  Core truncated magnitude division.  */

static arbint_err_t
arbint_tdiv_qr_mag_impl(arbint_t q, arbint_t r, const arbint_limb_t * np,
                        size_t nn, int nsign, const arbint_limb_t * dp,
                        size_t dn, int dsign, const arbint_alloc_t * alloc) {
  size_t qcap;
  size_t rcap;
  arbint_limb_t * qmag;
  arbint_limb_t * rmag;
  int rmag_needs_free;
  size_t q_used = 0u;
  size_t r_used;
  int qsign;
  int rsign;
  arbint_err_t rc;

  if (dn == 0u)
    return ARBINT_EZERO;
  if ((q == NULL && r == NULL) || alloc == NULL)
    return ARBINT_EINVAL;

  if (nn == 0u) {
    if (q != NULL)
      arbint_zero(q);
    if (r != NULL)
      arbint_zero(r);
    return ARBINT_OK;
  }

  if (nn < dn) {
    if (q != NULL)
      arbint_zero(q);
    if (r != NULL)
      return arbint_set_mag_signed(r, np, nn, nsign);
    return ARBINT_OK;
  }

  qcap = nn;
  if (dn > SIZE_MAX - 1u)
    return ARBINT_EOVERFLOW;
  rcap = dn + 1u;

  qmag = NULL;
  if (q != NULL) {
    qmag = arbint_alloc_limbs(alloc, qcap);
    if (qmag == NULL)
      return ARBINT_ENOMEM;
  }

  rmag_needs_free = 0;
  if (r == NULL && rcap <= ARBINT_TDIV_STACK_REM_LIMBS_MAX) {
    rmag = (arbint_limb_t *) alloca(rcap * sizeof(arbint_limb_t));
  } else {
    rmag = arbint_alloc_limbs(alloc, rcap);
    if (rmag == NULL) {
      arbint_free_limbs(alloc, qmag);
      return ARBINT_ENOMEM;
    }
    rmag_needs_free = 1;
  }

  if (dn >= 2u) {
    if (dn == 2u) {
      static arbint_div_mag_two_limb_fn_t impl2 = NULL;
      ARBINT_LAZY_INIT(impl2, arbint_select_div_mag_two_limb);
      if (impl2 != NULL)
        rc = impl2(np, nn, dp, qmag, rmag);
      else
        rc = arbint_div_mag_knuth(np, nn, dp, dn, qmag, rmag);
    } else {
      /*  Choose between Newton-Raphson and Knuth Algorithm D.
          Newton is O(M(n)) vs Knuth's O(n*m), but has higher constant factor.
          Use Newton for large divisors where asymptotic advantage dominates.  */
      if (dn >= ARBINT_NEWTON_DIV_THRESHOLD) {
        rc = arbint_div_mag_newton(np, nn, dp, dn, qmag, rmag, alloc);
      } else {
        rc = arbint_div_mag_knuth(np, nn, dp, dn, qmag, rmag);
      }
    }
    if (rc != ARBINT_OK) {
      if (rmag_needs_free)
        arbint_free_limbs(alloc, rmag);
      arbint_free_limbs(alloc, qmag);
      return rc;
    }
    if (qmag != NULL)
      q_used = arbint_norm_used(qmag, nn - dn + 1u);
    else
      q_used = 0u;
    r_used = arbint_norm_used(rmag, dn);
  } else {
    static arbint_div_mag_single_limb_fn_t impl = NULL;
    arbint_limb_t rem_limb = 0u;

    ARBINT_LAZY_INIT(impl, arbint_select_div_mag_single_limb);

    rc = impl(np, nn, dp[0], qmag, &q_used, &rem_limb);
    if (rc != ARBINT_OK) {
      if (rmag_needs_free)
        arbint_free_limbs(alloc, rmag);
      arbint_free_limbs(alloc, qmag);
      return rc;
    }
    rmag[0] = rem_limb;
    r_used = (rem_limb != 0u) ? 1u : 0u;
  }

  qsign = (q == NULL || q_used == 0u || nsign == 0)
              ? 0
              : ((nsign == dsign) ? 1 : -1);
  rsign = (r_used == 0u || nsign == 0) ? 0 : nsign;

  if (q != NULL) {
    rc = arbint_set_mag_signed(q, qmag, q_used, qsign);
    if (rc != ARBINT_OK) {
      if (rmag_needs_free)
        arbint_free_limbs(alloc, rmag);
      arbint_free_limbs(alloc, qmag);
      return rc;
    }
  }

  if (r != NULL) {
    rc = arbint_set_mag_signed(r, rmag, r_used, rsign);
    if (rc != ARBINT_OK) {
      if (rmag_needs_free)
        arbint_free_limbs(alloc, rmag);
      arbint_free_limbs(alloc, qmag);
      return rc;
    }
  }

  if (rmag_needs_free)
    arbint_free_limbs(alloc, rmag);
  arbint_free_limbs(alloc, qmag);
  return ARBINT_OK;
}

/*  Truncated division implementation (arbint / arbint).  */

arbint_err_t arbint_tdiv_qr_impl(arbint_t q, arbint_t r, const arbint_t n,
                                 const arbint_t d) {
  const arbint_limb_t * np;
  const arbint_limb_t * dp;
  size_t nn;
  size_t dn;
  int nsign;
  int dsign;
  const arbint_alloc_t * alloc;
  arbint_err_t rc;
  size_t pow2_bit;

  if ((q == NULL && r == NULL) || n == NULL || d == NULL)
    return ARBINT_EINVAL;

  rc = arbint_get_mag_view(n, &np, &nn, &nsign);
  if (rc != ARBINT_OK)
    return rc;

  rc = arbint_get_mag_view(d, &dp, &dn, &dsign);
  if (rc != ARBINT_OK)
    return rc;
  if (dn == 0u)
    return ARBINT_EZERO;

  /*  Power-of-two fast path: divisor is 2^k for some k.  */
  pow2_bit = arbint_mag_pow2_bit_or_size_max(dp, dn);
  if (pow2_bit != SIZE_MAX) {
    /*  k may exceed 32 bits for large divisors (2^64, 2^128, etc.).
        arbint_tdiv_qr_pow2 takes unsigned k, which may truncate.
        For k > UINT_MAX, we need to handle differently: the quotient
        is 0 unless n has more than k bits.

        For practical purposes, k fits in unsigned for divisors up to
        2^(UINT_MAX), which is astronomically large. We cast carefully.  */
    if (pow2_bit <= (size_t) UINT_MAX) {
      return arbint_tdiv_qr_pow2(q, r, n, (unsigned) pow2_bit, dsign);
    } else {
      return ARBINT_EOVERFLOW; /* Too large to handle! */
    }
  }

  alloc = arbint_pick_alloc4(q, r, n, d);
  if (alloc == NULL)
    return ARBINT_EINVAL;

  return arbint_tdiv_qr_mag_impl(q, r, np, nn, nsign, dp, dn, dsign, alloc);
}

/*  Divisibility testing: arbint_divisible_u32.
    Tests if n is divisible by d without computing the full quotient.
    Uses Barrett reduction for O(n) time with no division in the
    inner loop.  */

arbint_err_t arbint_divisible_u32(const arbint_t n, uint32_t d, int * out) {
  size_t nn;
  const arbint_limb_t * np;
  uint64_t rem;
  uint64_t m;
  size_t i;

  if (n == NULL || out == NULL)
    return ARBINT_EINVAL;
  if (d == 0u)
    return ARBINT_EZERO;

  *out = 0;

  /* Zero is divisible by any nonzero number. */
  if (n[0]._sz == 0) {
    *out = 1;
    return ARBINT_OK;
  }

  nn = arbint_abs_sz(n[0]._sz);
  np = ARBINT_CLIMBS(n);

  /* Power of 2: O(1) - check if low bits are zero. */
  if (arbint_u32_is_pow2(d)) {
    uint32_t mask = d - 1u;
#if ARBINT_LIMB_BITS == 64
    *out = ((uint32_t) np[0] & mask) == 0u;
#else
    *out = (np[0] & mask) == 0u;
#endif
    return ARBINT_OK;
  }

  /* Barrett reduction: precompute reciprocal m = floor((2^64 - 1) / d).
     For each 32-bit chunk, compute rem = (rem * 2^32 + chunk) mod d
     using q_approx = floor(rem * m / 2^64), then rem -= q_approx * d.
     The approximation may be off by 1, so we apply a single correction.

     Note: We use UINT64_MAX / d instead of the canonical 2^64 / d. This
     underestimates by at most 1 when d divides 2^64 exactly (i.e., d is
     a power of 2), but those cases are already handled by the fast path
     above, so this is safe.  */
  m = UINT64_MAX / (uint64_t) d;
  rem = 0u;

  /* Process limbs from high to low, 32-bit chunks at a time. */
  for (i = nn; i != 0u; --i) {
    arbint_limb_t limb = np[i - 1u];
#if ARBINT_LIMB_BITS == 64
    /* High 32 bits. */
    {
      uint64_t x = (rem << 32) | (uint32_t) (limb >> 32);
      /* Compute (x * m) >> 64 using half-word arithmetic.
         x and m are both 64-bit, so we decompose:
         x = x_hi * 2^32 + x_lo, m = m_hi * 2^32 + m_lo
         (x * m) >> 64 = x_hi * m_hi + ((x_hi * m_lo + x_lo * m_hi) >> 32)
         The x_lo * m_lo term contributes at most 1 to the high word. */
      uint64_t x_lo = x & 0xFFFFFFFFu;
      uint64_t x_hi = x >> 32;
      uint64_t m_lo = m & 0xFFFFFFFFu;
      uint64_t m_hi = m >> 32;
      uint64_t mid = x_hi * m_lo + x_lo * m_hi;
      uint64_t q_hi = x_hi * m_hi + (mid >> 32);
      rem = x - q_hi * (uint64_t) d;
      if (rem >= (uint64_t) d)
        rem -= (uint64_t) d;
    }
    /* Low 32 bits. */
    {
      uint64_t x = (rem << 32) | (uint32_t) limb;
      uint64_t x_lo = x & 0xFFFFFFFFu;
      uint64_t x_hi = x >> 32;
      uint64_t m_lo = m & 0xFFFFFFFFu;
      uint64_t m_hi = m >> 32;
      uint64_t mid = x_hi * m_lo + x_lo * m_hi;
      uint64_t q_hi = x_hi * m_hi + (mid >> 32);
      rem = x - q_hi * (uint64_t) d;
      if (rem >= (uint64_t) d)
        rem -= (uint64_t) d;
    }
#else
    /* 32-bit limb: single chunk. */
    {
      uint64_t x = (rem << 32) | (uint64_t) limb;
      uint64_t x_lo = x & 0xFFFFFFFFu;
      uint64_t x_hi = x >> 32;
      uint64_t m_lo = m & 0xFFFFFFFFu;
      uint64_t m_hi = m >> 32;
      uint64_t mid = x_hi * m_lo + x_lo * m_hi;
      uint64_t q_hi = x_hi * m_hi + (mid >> 32);
      rem = x - q_hi * (uint64_t) d;
      if (rem >= (uint64_t) d)
        rem -= (uint64_t) d;
    }
#endif
  }

  *out = (rem == 0u);
  return ARBINT_OK;
}

/*  Divisibility testing: arbint_divisible (multi-limb divisor).
    Tests if n is divisible by d.
    Key optimizations:
      1. Power of 2: O(1) - check trailing zeros
      2. Single-limb divisor: delegate to arbint_divisible_u32
      3. General: compute remainder via tdiv and check if zero  */

arbint_err_t arbint_divisible(const arbint_t n, const arbint_t d, int * out) {
  size_t nn;
  size_t dn;
  const arbint_limb_t * np;
  const arbint_limb_t * dp;
  int nsign;
  int dsign;
  arbint_err_t rc;
  size_t pow2_bit;

  if (n == NULL || d == NULL || out == NULL)
    return ARBINT_EINVAL;

  *out = 0;

  rc = arbint_get_mag_view(d, &dp, &dn, &dsign);
  if (rc != ARBINT_OK)
    return rc;
  if (dn == 0u)
    return ARBINT_EZERO;

  /* Zero is divisible by any nonzero number. */
  if (n[0]._sz == 0) {
    *out = 1;
    return ARBINT_OK;
  }

  rc = arbint_get_mag_view(n, &np, &nn, &nsign);
  if (rc != ARBINT_OK)
    return rc;

  /* Quick check: if |n| < |d|, n is not divisible (and n != 0). */
  if (nn < dn) {
    *out = 0;
    return ARBINT_OK;
  }

  /* Power of 2: check if n has enough trailing zeros. */
  pow2_bit = arbint_mag_pow2_bit_or_size_max(dp, dn);
  if (pow2_bit != SIZE_MAX) {
    /* d = 2^pow2_bit. n is divisible iff n has >= pow2_bit trailing zeros. */
    *out = (arbint_mag_ctz_or_size_max(np, nn) >= pow2_bit);
    return ARBINT_OK;
  }

  /* Single-limb divisor: use optimized u32/u64 path. */
  if (dn == 1u) {
#if ARBINT_LIMB_BITS == 64
    /* For 64-bit limbs, d may exceed uint32_t. */
    arbint_limb_t d_limb = dp[0];
    if (d_limb <= (arbint_limb_t) UINT32_MAX) {
      return arbint_divisible_u32(n, (uint32_t) d_limb, out);
    }
    /* 64-bit single-limb divisor: use full remainder computation.
       Fall through to the general multi-limb path below. */
#else
    /* 32-bit limbs: delegate to u32 path. */
    return arbint_divisible_u32(n, (uint32_t) dp[0], out);
#endif
  }

  /* Multi-limb divisor (or large single-limb on 64-bit):
     compute remainder and check if zero. */
  {
    arbint_t rem;
    rc = arbint_init(rem, n[0]._ctx);
    if (rc != ARBINT_OK)
      return rc;

    rc = arbint_tdiv_r(rem, n, d);
    if (rc != ARBINT_OK) {
      arbint_clear(rem);
      return rc;
    }

    *out = arbint_is_zero(rem);
    arbint_clear(rem);
    return ARBINT_OK;
  }
}

/*  Generalized division with rounding mode.  */

arbint_err_t arbint_div_qr_mode_impl(arbint_t q, arbint_t r, const arbint_t n,
                                     const arbint_t d,
                                     arbint_div_mode_t mode) {
  arbint_t q_tmp;
  arbint_t r_tmp;
  arbint_err_t rc;
  int nsign;
  int dsign;
  int need_adjust;

  if (mode == ARBINT_DIV_TRUNC)
    return arbint_tdiv_qr_impl(q, r, n, d);

  nsign = arbint_signum(n);
  dsign = arbint_signum(d);

  if (dsign == 0)
    return ARBINT_EZERO;

  if (nsign == 0) {
    if (q != NULL)
      arbint_zero(q);
    if (r != NULL)
      arbint_zero(r);
    return ARBINT_OK;
  }

  if (q == NULL) {
    rc = arbint_init(q_tmp, n[0]._ctx);
    if (rc != ARBINT_OK)
      return rc;
  } else {
    q_tmp[0] = q[0];
  }

  if (r == NULL) {
    rc = arbint_init(r_tmp, n[0]._ctx);
    if (rc != ARBINT_OK) {
      if (q == NULL)
        arbint_clear(q_tmp);
      return rc;
    }
  } else {
    r_tmp[0] = r[0];
  }

  rc = arbint_tdiv_qr_impl(q_tmp, r_tmp, n, d);
  if (rc != ARBINT_OK) {
    if (q == NULL)
      arbint_clear(q_tmp);
    if (r == NULL)
      arbint_clear(r_tmp);
    return rc;
  }

  if (arbint_is_zero(r_tmp)) {
    need_adjust = 0;
  } else if (mode == ARBINT_DIV_FLOOR) {
    need_adjust = (nsign != dsign);
  } else {
    need_adjust = (nsign == dsign);
  }

  if (need_adjust) {
    if (mode == ARBINT_DIV_FLOOR) {
      rc = arbint_sub_i32(q_tmp, q_tmp, 1);
      if (rc == ARBINT_OK)
        rc = arbint_add(r_tmp, r_tmp, d);
    } else {
      rc = arbint_add_i32(q_tmp, q_tmp, 1);
      if (rc == ARBINT_OK)
        rc = arbint_sub(r_tmp, r_tmp, d);
    }
    if (rc != ARBINT_OK) {
      if (q == NULL)
        arbint_clear(q_tmp);
      if (r == NULL)
        arbint_clear(r_tmp);
      return rc;
    }
  }

  if (q == NULL) {
    arbint_clear(q_tmp);
  } else {
    q[0] = q_tmp[0];
  }

  if (r == NULL) {
    arbint_clear(r_tmp);
  } else {
    r[0] = r_tmp[0];
  }

  return ARBINT_OK;
}
