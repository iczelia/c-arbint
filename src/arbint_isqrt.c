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

/*  Karatsuba square root implementation.

    Based on "Karatsuba Square Root" by Paul Zimmermann (INRIA).
    Reference: https://hal.inria.fr/inria-00072854

    The algorithm computes floor(sqrt(N)) for an n-limb integer N in O(M(n))
    time, where M(n) is the time for n-limb multiplication.

    Key insight: Instead of O(log n) Newton iterations each costing O(M(n)),
    we recursively split the problem:
      1. Compute sqrt of high half (recursive call on n/2 limbs)
      2. One division to refine with low half
      3. One squaring to verify/adjust

    The recurrence T(n) = T(n/2) + O(M(n)) solves to T(n) = O(M(n)).  */

#include "arbint_isqrt.h"
#include "arbint.h"
#include "arbint_div.h"
#include "arbint_mul.h"
#include "arbint_shift.h"

#include <stdio.h>
#include <string.h>

/*  Table for 8-bit 1/sqrt approximation.
    Entry i approximates floor(256 / sqrt((256+i)/256)) - 256 for i in 0..127.  */
const unsigned char arbint_invsqrt_tab[128] = {
    /*  sqrt(1/0x100)..sqrt(1/0x107)  */
    0xff, 0xfd, 0xfb, 0xf9, 0xf7, 0xf5, 0xf3, 0xf2,
    /*  sqrt(1/0x108)..sqrt(1/0x10f)  */
    0xf0, 0xee, 0xec, 0xea, 0xe9, 0xe7, 0xe5, 0xe4,
    /*  sqrt(1/0x110)..sqrt(1/0x117)  */
    0xe2, 0xe0, 0xdf, 0xdd, 0xdb, 0xda, 0xd8, 0xd7,
    /*  sqrt(1/0x118)..sqrt(1/0x11f)  */
    0xd5, 0xd4, 0xd2, 0xd1, 0xcf, 0xce, 0xcc, 0xcb,
    /*  sqrt(1/0x120)..sqrt(1/0x127)  */
    0xc9, 0xc8, 0xc6, 0xc5, 0xc4, 0xc2, 0xc1, 0xc0,
    /*  sqrt(1/0x128)..sqrt(1/0x12f)  */
    0xbe, 0xbd, 0xbc, 0xba, 0xb9, 0xb8, 0xb7, 0xb5,
    /*  sqrt(1/0x130)..sqrt(1/0x137)  */
    0xb4, 0xb3, 0xb2, 0xb0, 0xaf, 0xae, 0xad, 0xac,
    /*  sqrt(1/0x138)..sqrt(1/0x13f)  */
    0xaa, 0xa9, 0xa8, 0xa7, 0xa6, 0xa5, 0xa4, 0xa3,
    /*  sqrt(1/0x140)..sqrt(1/0x147)  */
    0xa2, 0xa0, 0x9f, 0x9e, 0x9d, 0x9c, 0x9b, 0x9a,
    /*  sqrt(1/0x148)..sqrt(1/0x14f)  */
    0x99, 0x98, 0x97, 0x96, 0x95, 0x94, 0x93, 0x92,
    /*  sqrt(1/0x150)..sqrt(1/0x157)  */
    0x91, 0x90, 0x8f, 0x8e, 0x8d, 0x8c, 0x8c, 0x8b,
    /*  sqrt(1/0x158)..sqrt(1/0x15f)  */
    0x8a, 0x89, 0x88, 0x87, 0x86, 0x85, 0x84, 0x83,
    /*  sqrt(1/0x160)..sqrt(1/0x167)  */
    0x83, 0x82, 0x81, 0x80, 0x7f, 0x7e, 0x7e, 0x7d,
    /*  sqrt(1/0x168)..sqrt(1/0x16f)  */
    0x7c, 0x7b, 0x7a, 0x79, 0x79, 0x78, 0x77, 0x76,
    /*  sqrt(1/0x170)..sqrt(1/0x177)  */
    0x76, 0x75, 0x74, 0x73, 0x72, 0x72, 0x71, 0x70,
    /*  sqrt(1/0x178)..sqrt(1/0x17f)  */
    0x6f, 0x6f, 0x6e, 0x6d, 0x6d, 0x6c, 0x6b, 0x6a
};

/*  Single-limb square root using Newton iteration.

    Uses classic Newton iteration: x_{k+1} = (x_k + a / x_k) / 2
    Converges quadratically from above when x_0 >= sqrt(a).

    For a single limb, we use a simpler approach that avoids the
    tricky 1/sqrt conversion and denormalization issues.  */
arbint_limb_t arbint_isqrt_1(arbint_limb_t a0) {
  arbint_limb_t x, x_prev;

  if (a0 == 0)
    return 0;

  if (a0 == 1)
    return 1;

  /*  Initial guess: x = 1 << ((bits(a0) + 1) / 2).
      This ensures x >= sqrt(a0).  */
  {
    unsigned bits = ARBINT_LIMB_BITS - arbint_clz_limb(a0);
    x = (arbint_limb_t) 1 << ((bits + 1) / 2);
  }

  /*  Newton iteration: x = (x + a0/x) / 2.
      Converges in O(log(bits)) iterations.  */
  do {
    x_prev = x;
    x = (x + a0 / x) / 2;
  } while (x < x_prev);

  /*  x_prev is floor(sqrt(a0)).  */
  return x_prev;
}

/*  Two-limb square root with remainder.
    Input: {np, 2} with np[1] >= 2^(LIMB_BITS-2) (normalized).
    Output: sqrt in *sp, low remainder limb in *rp.
    Returns high remainder limb (0 or 1).  */
static arbint_limb_t arbint_isqrtrem_2(arbint_limb_t * sp, arbint_limb_t * rp,
                                       const arbint_limb_t * np) {
  arbint_limb_t np1 = np[1];
  arbint_limb_t np0 = np[0];
  arbint_limb_t s0, r0, q, u, q2;
  int cc;

  /*  Compute sqrt of high limb.  */
  s0 = arbint_isqrt_1(np1);

  /*  r0 = np1 - s0^2.  */
  r0 = np1 - s0 * s0;

  /*  Combine remainder with low half: (r0 << HALF_BITS) + (np0 >> HALF_BITS).
      Then divide by s0 to get next digit.  */
#if ARBINT_LIMB_BITS == 64
  {
    /*  Form 2-limb remainder and divide by s0.  */
    arbint_limb_t r_hi = (r0 << 31) | (np0 >> 33);
    q = r_hi / s0;
    if (q > 0xFFFFFFFFULL)
      q = 0xFFFFFFFFULL;
    u = r_hi - q * s0;

    /*  Combine: s = (s0 << 32) | q.  */
    s0 = (s0 << 32) | q;

    /*  Compute remainder and adjust.  */
    cc = (int) (u >> 31);
    r0 = ((u << 33) & 0xFFFFFFFFFFFFFFFFULL) |
         (np0 & ((1ULL << 33) - 1));

    q2 = q * q;
    cc -= (r0 < q2) ? 1 : 0;
    r0 -= q2;

    if (cc < 0) {
      r0 += s0;
      cc += (r0 < s0) ? 1 : 0;
      --s0;
      r0 += s0;
      cc += (r0 < s0) ? 1 : 0;
    }
  }
#else
  {
    /*  32-bit version.  */
    arbint_limb_t r_hi = (r0 << 15) | (np0 >> 17);
    q = r_hi / s0;
    if (q > 0xFFFFu)
      q = 0xFFFFu;
    u = r_hi - q * s0;

    s0 = (s0 << 16) | q;

    cc = (int) (u >> 15);
    r0 = ((u << 17) & 0xFFFFFFFFu) | (np0 & ((1u << 17) - 1));

    q2 = q * q;
    cc -= (r0 < q2) ? 1 : 0;
    r0 -= q2;

    if (cc < 0) {
      r0 += s0;
      cc += (r0 < s0) ? 1 : 0;
      --s0;
      r0 += s0;
      cc += (r0 < s0) ? 1 : 0;
    }
  }
#endif

  *sp = s0;
  *rp = r0;
  return (arbint_limb_t) cc;
}

/*  Two-limb square root wrapper (without exposing remainder).  */
arbint_limb_t arbint_isqrt_2(arbint_limb_t * sp, const arbint_limb_t * np) {
  arbint_limb_t rp0;
  return arbint_isqrtrem_2(sp, &rp0, np);
}

typedef size_t (*arbint_mulacc_1_fn_t)(arbint_limb_t * dst, size_t dst_n,
                                       size_t dst_cap,
                                       const arbint_limb_t * a, size_t an,
                                       arbint_limb_t b);

/*  Shared mulacc_1 dispatch for correction steps.
    We preserve rp[n] because the original local helper only touched rp[0..n-1]
    and returned carry separately.  */
static arbint_limb_t arbint_isqrt_addmul_1(arbint_limb_t * rp,
                                           const arbint_limb_t * ap, size_t n,
                                           arbint_limb_t b) {
  static arbint_mulacc_1_fn_t mulacc_1_impl = NULL;
  arbint_limb_t saved_hi;
  arbint_limb_t carry_out;

  if (mulacc_1_impl == NULL)
    mulacc_1_impl = arbint_mul_kernel_table_get()->mulacc_1;

  saved_hi = rp[n];
  rp[n] = 0u;
  (void) mulacc_1_impl(rp, n, n + 1u, ap, n, b);
  carry_out = rp[n];
  rp[n] = saved_hi;
  return carry_out;
}

/*  Compute floor(sqrt(a)) with integer Newton iteration.
    Shared by public arbint_isqrt small-input path and DC base cases.  */
arbint_err_t arbint_isqrt_newton_mag(arbint_t rop, const arbint_t a) {
  arbint_t x, t;
  arbint_err_t rc;
  size_t nbits;

  rc = arbint_init_all(a[0]._ctx, x, t, (arbint_t *) NULL);
  if (rc != ARBINT_OK)
    return rc;

  nbits = arbint_nbits(a);
  rc = arbint_set_i32(x, 1);
  if (rc != ARBINT_OK)
    goto cleanup;
  rc = arbint_shl(x, x, (uint32_t) ((nbits + 1u) / 2u));
  if (rc != ARBINT_OK)
    goto cleanup;

  for (;;) {
    rc = arbint_tdiv_q(t, a, x);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_add(t, t, x);
    if (rc != ARBINT_OK)
      goto cleanup;
    rc = arbint_shr(t, t, 1u);
    if (rc != ARBINT_OK)
      goto cleanup;
    if (arbint_cmp(t, x) >= 0)
      break;
    rc = arbint_set(x, t);
    if (rc != ARBINT_OK)
      goto cleanup;
  }

  rc = arbint_set(rop, x);

cleanup:
  arbint_clear_all(x, t, (arbint_t *) NULL);
  return rc;
}

/*  Base-case helper for small 2n-limb operands (n==2).
    Computes sqrt and remainder using local Newton iteration.  */
static int arbint_isqrt_dc_base_2(arbint_limb_t * sp, arbint_limb_t * np,
                                  const arbint_alloc_t * alloc) {
  arbint_ctx_t local_ctx;
  arbint_t a, s, sq;
  arbint_err_t rc;
  size_t an;
  size_t sn;
  arbint_limb_t q = 0;

  if (alloc != NULL && alloc->realloc != NULL)
    rc = arbint_ctx_init(&local_ctx, alloc, 0u);
  else
    rc = arbint_ctx_init_default(&local_ctx);
  if (rc != ARBINT_OK)
    return -1;

  rc = arbint_init_all(&local_ctx, a, s, sq, (arbint_t *) NULL);
  if (rc != ARBINT_OK) {
    arbint_ctx_clear(&local_ctx);
    return -1;
  }

  an = arbint_norm_used(np, 4u);
  if (an == 0u) {
    sp[0] = 0u;
    sp[1] = 0u;
    np[0] = 0u;
    np[1] = 0u;
    np[2] = 0u;
    np[3] = 0u;
    arbint_clear_all(a, s, sq, (arbint_t *) NULL);
    arbint_ctx_clear(&local_ctx);
    return 0;
  }

  rc = arbint_resize(a, an);
  if (rc != ARBINT_OK)
    goto fail;
  memcpy(ARBINT_LIMBS(a), np, an * sizeof(arbint_limb_t));
  a[0]._sz = (ptrdiff_t) an;

  rc = arbint_isqrt_newton_mag(s, a);
  if (rc != ARBINT_OK)
    goto fail;

  sp[0] = 0u;
  sp[1] = 0u;
  sn = arbint_abs_sz(s[0]._sz);
  if (sn > 2u)
    goto fail;
  if (sn != 0u)
    memcpy(sp, ARBINT_CLIMBS(s), sn * sizeof(arbint_limb_t));

  rc = arbint_sqr(sq, s);
  if (rc != ARBINT_OK)
    goto fail;
  rc = arbint_sub(a, a, sq);
  if (rc != ARBINT_OK)
    goto fail;

  np[0] = 0u;
  np[1] = 0u;
  np[2] = 0u;
  np[3] = 0u;
  an = arbint_abs_sz(a[0]._sz);
  if (an > 4u)
    goto fail;
  if (an != 0u) {
    const arbint_limb_t * ap = ARBINT_CLIMBS(a);
    if (an > 0u)
      np[0] = ap[0];
    if (an > 1u)
      np[1] = ap[1];
    if (an > 2u)
      q = ap[2];
  }

  arbint_clear_all(a, s, sq, (arbint_t *) NULL);
  arbint_ctx_clear(&local_ctx);
  return (int) q;

fail:
  arbint_clear_all(a, s, sq, (arbint_t *) NULL);
  arbint_ctx_clear(&local_ctx);
  return -1;
}

/*  Karatsuba divide-and-conquer square root.
    Writes floor(sqrt({np,2n})) into {sp,n}.
    Overwrites {np,2n} with remainder (low n limbs in np[0..n-1], high carry in
    return value). Returns 0 for perfect square, 1 otherwise, or -1 on failure.
    Requires scratch capacity >= floor(n/2)+1 limbs for n>1.  */
int arbint_isqrt_dc(arbint_limb_t * sp, arbint_limb_t * np, size_t n,
                    arbint_limb_t * scratch, const arbint_alloc_t * alloc) {
  size_t l, h;
  arbint_limb_t q;
  int c, b;
  size_t sq_used, i;

  if (sp == NULL || np == NULL || n == 0u)
    return -1;
  if (n > 1u && scratch == NULL)
    return -1;

  if (n == 1u) {
    q = arbint_isqrtrem_2(sp, np, np);
    np[1] = q;
    return (int) q;
  }
  if (n == 2u)
    return arbint_isqrt_dc_base_2(sp, np, alloc);

  l = n / 2u;
  h = n - l;

  if (h == 1u) {
    q = arbint_isqrtrem_2(sp + l, np + 2u * l, np + 2u * l);
  } else {
    int qi = arbint_isqrt_dc(sp + l, np + 2u * l, h, scratch, alloc);
    if (qi < 0)
      return -1;
    q = (arbint_limb_t) qi;
  }

  if (q != 0u)
    arbint_limb_sub_n(np + 2u * l, np + 2u * l, sp + l, h);

  if (arbint_div_mag_knuth(np + l, n, sp + l, h, scratch, np + l) != ARBINT_OK)
    return -1;
  q += scratch[l];
  c = (int) (scratch[0] & 1u);

  arbint_rshift_limbs_inplace(scratch, l, 1u);
  memcpy(sp, scratch, l * sizeof(arbint_limb_t));
  sp[l - 1u] |= (q << (ARBINT_LIMB_BITS - 1u));
  q >>= 1u;

  if (c != 0)
    c = (int) arbint_limb_add_n(np + l, np + l, sp + l, h);

  if (arbint_mul_mag_generic(np + n, &sq_used, sp, l, sp, l, alloc) != ARBINT_OK)
    return -1;
  for (i = sq_used; i < 2u * l; ++i)
    np[n + i] = 0u;

  b = (int) q + (int) arbint_limb_sub_n(np, np, np + n, 2u * l);
  if (l == h)
    c -= b;
  else
    c -=
        (int) arbint_limb_sub_1(np + 2u * l, np + 2u * l, 1u, (arbint_limb_t) b);

  if (c < 0) {
    q = arbint_limb_add_1(sp + l, sp + l, h, q);
    c += (int) arbint_isqrt_addmul_1(np, sp, n, 2u) + 2 * (int) q;
    c -= (int) arbint_limb_sub_1(np, np, n, 1u);
    q -= arbint_limb_sub_1(sp, sp, n, 1u);
  }
  return c;
}
