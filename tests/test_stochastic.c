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

#include <stdint.h>
#include <string.h>

ARBINT_TEST_DECLARE_FAILURES();

#define STOCH_MAX_LIMBS 512u
#define STOCH_MAX_INPUT_LIMBS 180u

typedef struct {
  uint64_t state;
} stoch_rng_t;

typedef struct {
  int sign; /* -1, 0, +1 */
  size_t n; /* used limbs */
  arbint_limb_t limbs[STOCH_MAX_LIMBS];
} refint_t;

static void stoch_seed(stoch_rng_t * rng, uint64_t seed) {
  if (seed == 0u)
    seed = 0x9e3779b97f4a7c15ull;
  rng->state = seed;
}

static uint64_t stoch_u64(stoch_rng_t * rng) {
  uint64_t x = rng->state;
  x ^= x >> 12;
  x ^= x << 25;
  x ^= x >> 27;
  rng->state = x;
  return x * 0x2545f4914f6cdd1dull;
}

static size_t stoch_range(stoch_rng_t * rng, size_t hi_exclusive) {
  if (hi_exclusive == 0u)
    return 0u;
  return (size_t) (stoch_u64(rng) % (uint64_t) hi_exclusive);
}

static arbint_limb_t stoch_limb(stoch_rng_t * rng) {
#if ARBINT_LIMB_BITS == 64
  return (arbint_limb_t) stoch_u64(rng);
#else
  return (arbint_limb_t) (stoch_u64(rng) & 0xffffffffu);
#endif /* ARBINT_LIMB_BITS */
}

static void ref_zero(refint_t * x) {
  x->sign = 0;
  x->n = 0u;
}

static void ref_copy(refint_t * dst, const refint_t * src) {
  dst->sign = src->sign;
  dst->n = src->n;
  if (src->n != 0u)
    memcpy(dst->limbs, src->limbs, src->n * sizeof(arbint_limb_t));
}

static void ref_norm(refint_t * x) {
  while (x->n != 0u && x->limbs[x->n - 1u] == (arbint_limb_t) 0u)
    --x->n;
  if (x->n == 0u)
    x->sign = 0;
  else if (x->sign == 0)
    x->sign = 1;
}

static int ref_cmp_mag_limbs(const arbint_limb_t * a, size_t an,
                             const arbint_limb_t * b, size_t bn) {
  size_t i;
  if (an > bn)
    return 1;
  if (an < bn)
    return -1;
  i = an;
  while (i != 0u) {
    --i;
    if (a[i] > b[i])
      return 1;
    if (a[i] < b[i])
      return -1;
  }
  return 0;
}

static int ref_cmp_abs(const refint_t * a, const refint_t * b) {
  return ref_cmp_mag_limbs(a->limbs, a->n, b->limbs, b->n);
}

static int ref_eq(const refint_t * a, const refint_t * b) {
  if (a->sign != b->sign || a->n != b->n)
    return 0;
  if (a->n == 0u)
    return 1;
  return memcmp(a->limbs, b->limbs, a->n * sizeof(arbint_limb_t)) == 0;
}

static int ref_add_mag(arbint_limb_t * out, size_t * out_n,
                       const arbint_limb_t * a, size_t an,
                       const arbint_limb_t * b, size_t bn) {
  size_t i;
  size_t n = (an > bn) ? an : bn;
  arbint_limb_t carry = 0u;

  if (n > STOCH_MAX_LIMBS)
    return 0;

  for (i = 0u; i < n; ++i) {
    arbint_limb_t ai = (i < an) ? a[i] : (arbint_limb_t) 0u;
    arbint_limb_t bi = (i < bn) ? b[i] : (arbint_limb_t) 0u;
    arbint_limb_t s1 = ai + bi;
    arbint_limb_t c1 = (s1 < ai) ? 1u : 0u;
    arbint_limb_t s2 = s1 + carry;
    arbint_limb_t c2 = (s2 < s1) ? 1u : 0u;
    out[i] = s2;
    carry = (arbint_limb_t) (c1 + c2);
  }

  if (carry != (arbint_limb_t) 0u) {
    if (n >= STOCH_MAX_LIMBS)
      return 0;
    out[n++] = carry;
  }

  *out_n = n;
  while (*out_n != 0u && out[*out_n - 1u] == (arbint_limb_t) 0u)
    --(*out_n);
  return 1;
}

static void ref_sub_mag(arbint_limb_t * out, size_t * out_n,
                        const arbint_limb_t * a, size_t an,
                        const arbint_limb_t * b, size_t bn) {
  size_t i;
  arbint_limb_t borrow = 0u;

  for (i = 0u; i < an; ++i) {
    arbint_limb_t ai = a[i];
    arbint_limb_t bi = (i < bn) ? b[i] : (arbint_limb_t) 0u;
    arbint_limb_t t = ai - bi;
    arbint_limb_t b1 = (ai < bi) ? 1u : 0u;
    arbint_limb_t d = t - borrow;
    arbint_limb_t b2 = (t < borrow) ? 1u : 0u;
    out[i] = d;
    borrow = (arbint_limb_t) (b1 + b2);
  }

  *out_n = an;
  while (*out_n != 0u && out[*out_n - 1u] == (arbint_limb_t) 0u)
    --(*out_n);
}

#if ARBINT_LIMB_BITS == 64
static void ref_mul64_wide(uint64_t a, uint64_t b, uint64_t * hi,
                           uint64_t * lo) {
  uint64_t a0 = (uint32_t) a;
  uint64_t a1 = a >> 32u;
  uint64_t b0 = (uint32_t) b;
  uint64_t b1 = b >> 32u;
  uint64_t p00 = a0 * b0;
  uint64_t p01 = a0 * b1;
  uint64_t p10 = a1 * b0;
  uint64_t p11 = a1 * b1;
  uint64_t t = (p00 >> 32u) + (uint32_t) p01 + (uint32_t) p10;

  *lo = (p00 & 0xffffffffull) | (t << 32u);
  *hi = p11 + (p01 >> 32u) + (p10 >> 32u) + (t >> 32u);
}

static int ref_mul_add_limb(arbint_limb_t a, arbint_limb_t b, arbint_limb_t in,
                            arbint_limb_t * carry_io, arbint_limb_t * out) {
  uint64_t hi;
  uint64_t lo;
  uint64_t prev;
  uint64_t carry = (uint64_t) *carry_io;

  ref_mul64_wide((uint64_t) a, (uint64_t) b, &hi, &lo);

  prev = lo;
  lo += (uint64_t) in;
  if (lo < prev) {
    if (hi == UINT64_MAX)
      return 0;
    ++hi;
  }

  prev = lo;
  lo += carry;
  if (lo < prev) {
    if (hi == UINT64_MAX)
      return 0;
    ++hi;
  }

  *out = (arbint_limb_t) lo;
  *carry_io = (arbint_limb_t) hi;
  return 1;
}
#elif ARBINT_LIMB_BITS == 32
static int ref_mul_add_limb(arbint_limb_t a, arbint_limb_t b, arbint_limb_t in,
                            arbint_limb_t * carry_io, arbint_limb_t * out) {
  uint64_t acc =
      (uint64_t) a * (uint64_t) b + (uint64_t) in + (uint64_t) *carry_io;
  *out = (arbint_limb_t) acc;
  *carry_io = (arbint_limb_t) (acc >> 32u);
  return 1;
}
#endif /* ARBINT_LIMB_BITS */

static int ref_add(refint_t * out, const refint_t * a, const refint_t * b) {
  int cmp;

  if (a->sign == 0) {
    ref_copy(out, b);
    return 1;
  }
  if (b->sign == 0) {
    ref_copy(out, a);
    return 1;
  }

  if (a->sign == b->sign) {
    if (!ref_add_mag(out->limbs, &out->n, a->limbs, a->n, b->limbs, b->n))
      return 0;
    out->sign = (out->n == 0u) ? 0 : a->sign;
    return 1;
  }

  cmp = ref_cmp_abs(a, b);
  if (cmp == 0) {
    ref_zero(out);
    return 1;
  }

  if (cmp > 0) {
    ref_sub_mag(out->limbs, &out->n, a->limbs, a->n, b->limbs, b->n);
    out->sign = (out->n == 0u) ? 0 : a->sign;
  } else {
    ref_sub_mag(out->limbs, &out->n, b->limbs, b->n, a->limbs, a->n);
    out->sign = (out->n == 0u) ? 0 : b->sign;
  }

  return 1;
}

static int ref_sub(refint_t * out, const refint_t * a, const refint_t * b) {
  refint_t nb;
  ref_copy(&nb, b);
  nb.sign = -nb.sign;
  return ref_add(out, a, &nb);
}

static int ref_mul(refint_t * out, const refint_t * a, const refint_t * b) {
  size_t i;
  size_t j;
  size_t n;

  if (a->sign == 0 || b->sign == 0) {
    ref_zero(out);
    return 1;
  }

  n = a->n + b->n;
  if (n > STOCH_MAX_LIMBS)
    return 0;

  memset(out->limbs, 0, n * sizeof(arbint_limb_t));
  for (i = 0u; i < a->n; ++i) {
    arbint_limb_t carry = (arbint_limb_t) 0u;
    for (j = 0u; j < b->n; ++j) {
      size_t k = i + j;
      if (!ref_mul_add_limb(a->limbs[i], b->limbs[j], out->limbs[k], &carry,
                            &out->limbs[k]))
        return 0;
    }

    {
      size_t k = i + b->n;
      while (carry != 0u) {
        if (k >= n)
          return 0;
        {
          arbint_limb_t prev = out->limbs[k];
          out->limbs[k] = (arbint_limb_t) (out->limbs[k] + carry);
          carry =
              (out->limbs[k] < prev) ? (arbint_limb_t) 1u : (arbint_limb_t) 0u;
        }
        ++k;
      }
    }
  }

  out->n = n;
  out->sign = (a->sign == b->sign) ? 1 : -1;
  ref_norm(out);
  return 1;
}

static size_t ref_nbits_mag(const arbint_limb_t * x, size_t n) {
  size_t bits;
  arbint_limb_t w;

  if (n == 0u)
    return 0u;

  bits = (n - 1u) * ARBINT_LIMB_BITS;
  w = x[n - 1u];
  while (w != (arbint_limb_t) 0u) {
    ++bits;
    w >>= 1u;
  }
  return bits;
}

static int ref_mag_get_bit(const arbint_limb_t * x, size_t bit) {
  size_t li = bit / ARBINT_LIMB_BITS;
  size_t bi = bit % ARBINT_LIMB_BITS;
  return (int) ((x[li] >> bi) & (arbint_limb_t) 1u);
}

static int ref_mag_shl1(arbint_limb_t * x, size_t * n, size_t cap) {
  size_t i;
  arbint_limb_t carry = 0u;

  for (i = 0u; i < *n; ++i) {
    arbint_limb_t w = x[i];
    arbint_limb_t nc = w >> (ARBINT_LIMB_BITS - 1u);
    x[i] = (w << 1u) | carry;
    carry = nc;
  }

  if (carry != (arbint_limb_t) 0u) {
    if (*n >= cap)
      return 0;
    x[*n] = carry;
    ++(*n);
  }
  return 1;
}

static int ref_div_mag(refint_t * q, refint_t * r, const refint_t * n,
                       const refint_t * d) {
  size_t bits;
  size_t b;

  ref_zero(q);
  ref_zero(r);

  if (d->n == 0u)
    return 0;
  if (n->n == 0u)
    return 1;

  memset(q->limbs, 0, sizeof(q->limbs));
  memset(r->limbs, 0, sizeof(r->limbs));

  bits = ref_nbits_mag(n->limbs, n->n);
  for (b = bits; b != 0u; --b) {
    size_t qi;
    size_t qb;
    if (!ref_mag_shl1(r->limbs, &r->n, STOCH_MAX_LIMBS))
      return 0;
    if (ref_mag_get_bit(n->limbs, b - 1u)) {
      if (r->n == 0u) {
        r->n = 1u;
        r->limbs[0] = (arbint_limb_t) 1u;
      } else {
        r->limbs[0] |= (arbint_limb_t) 1u;
      }
    }
    if (ref_cmp_mag_limbs(r->limbs, r->n, d->limbs, d->n) >= 0) {
      ref_sub_mag(r->limbs, &r->n, r->limbs, r->n, d->limbs, d->n);
      qi = (b - 1u) / ARBINT_LIMB_BITS;
      qb = (b - 1u) % ARBINT_LIMB_BITS;
      if (qi >= STOCH_MAX_LIMBS)
        return 0;
      q->limbs[qi] |= ((arbint_limb_t) 1u) << qb;
      if (q->n < qi + 1u)
        q->n = qi + 1u;
    }
  }

  q->sign = (q->n == 0u) ? 0 : 1;
  r->sign = (r->n == 0u) ? 0 : 1;
  ref_norm(q);
  ref_norm(r);
  return 1;
}

static int ref_tdiv(refint_t * q, refint_t * r, const refint_t * n,
                    const refint_t * d) {
  refint_t an;
  refint_t ad;
  refint_t qmag;
  refint_t rmag;

  if (d->sign == 0)
    return 0;

  if (n->sign == 0) {
    ref_zero(q);
    ref_zero(r);
    return 1;
  }

  ref_copy(&an, n);
  ref_copy(&ad, d);
  an.sign = (an.n == 0u) ? 0 : 1;
  ad.sign = (ad.n == 0u) ? 0 : 1;

  if (ref_cmp_abs(&an, &ad) < 0) {
    ref_zero(q);
    ref_copy(r, n);
    return 1;
  }

  if (!ref_div_mag(&qmag, &rmag, &an, &ad))
    return 0;

  ref_copy(q, &qmag);
  ref_copy(r, &rmag);

  if (q->n != 0u)
    q->sign = (n->sign == d->sign) ? 1 : -1;
  else
    q->sign = 0;

  if (r->n != 0u)
    r->sign = n->sign;
  else
    r->sign = 0;

  return 1;
}

static void ref_from_u32(refint_t * out, uint32_t v) {
  if (v == 0u) {
    ref_zero(out);
    return;
  }
  out->sign = 1;
  out->n = 1u;
  out->limbs[0] = (arbint_limb_t) v;
}

static void ref_from_i32(refint_t * out, int32_t v) {
  if (v == 0) {
    ref_zero(out);
    return;
  }
  out->sign = (v < 0) ? -1 : 1;
  out->n = 1u;
  if (v < 0)
    out->limbs[0] = (arbint_limb_t) ((uint32_t) (-(v + 1)) + 1u);
  else
    out->limbs[0] = (arbint_limb_t) (uint32_t) v;
}

static arbint_err_t ref_to_arbint(arbint_t x, const refint_t * in) {
  arbint_err_t rc;
  if (in->n == 0u) {
    arbint_zero(x);
    return ARBINT_OK;
  }

  rc = arbint_resize(x, in->n);
  if (rc != ARBINT_OK)
    return rc;
  memcpy(ARBINT_LIMBS(x), in->limbs, in->n * sizeof(arbint_limb_t));
  x[0]._sz = (in->sign < 0) ? -(ptrdiff_t) in->n : (ptrdiff_t) in->n;
  return ARBINT_OK;
}

static int ref_from_arbint(refint_t * out, const arbint_t x) {
  size_t used;

  if (x == NULL)
    return 0;
  used = arbint_abs_sz(x[0]._sz);
  if (used == 0u) {
    if (x[0]._sz != 0)
      return 0;
    ref_zero(out);
    return 1;
  }
  if (used > STOCH_MAX_LIMBS || x[0]._ptr == NULL || x[0]._cap < used)
    return 0;
  if (ARBINT_CLIMBS(x)[used - 1u] == (arbint_limb_t) 0u)
    return 0;

  out->n = used;
  out->sign = (x[0]._sz < 0) ? -1 : 1;
  memcpy(out->limbs, ARBINT_CLIMBS(x), used * sizeof(arbint_limb_t));
  return 1;
}

static void ref_dump(const char * label, const refint_t * x) {
  size_t i;
  fprintf(stderr, "%s: sign=%d n=%zu", label, x->sign, x->n);
  if (x->n != 0u) {
    fprintf(stderr, " limbs=");
    for (i = x->n; i != 0u; --i) {
#if ARBINT_LIMB_BITS == 64
      fprintf(stderr, "%016llx", (unsigned long long) x->limbs[i - 1u]);
#else
      fprintf(stderr, "%08x", (unsigned int) x->limbs[i - 1u]);
#endif /* ARBINT_LIMB_BITS */
      if (i > 1u)
        fputc('_', stderr);
    }
  }
  fputc('\n', stderr);
}

static void check_ref_match(const char * op, size_t iter, const refint_t * a,
                            const refint_t * b, const refint_t * expect,
                            const arbint_t got) {
  enum { MAX_VERBOSE_STOCH_FAILURES = 12 };
  static int g_stoch_suppressed = 0;
  refint_t actual;

  if (ref_from_arbint(&actual, got) && ref_eq(&actual, expect))
    return;

  ++g_failures;
  if (g_failures > MAX_VERBOSE_STOCH_FAILURES) {
    if (!g_stoch_suppressed) {
      fprintf(stderr,
              "stochastic: more mismatches suppressed after %d failures\n",
              MAX_VERBOSE_STOCH_FAILURES);
      g_stoch_suppressed = 1;
    }
    return;
  }

  fprintf(stderr, "FAIL stochastic[%s] iter=%zu\n", op, iter);
  if (!ref_from_arbint(&actual, got))
    fprintf(stderr, "got: INVALID-NONNORMALIZED\n");
  else
    ref_dump("got", &actual);
  ref_dump("a", a);
  ref_dump("b", b);
  ref_dump("expect", expect);
}

static void ref_random(refint_t * out, stoch_rng_t * rng, size_t max_limbs,
                       int allow_zero) {
  size_t i;
  size_t n;
  uint64_t bucket;
  uint64_t pattern;
  arbint_limb_t maxv = (arbint_limb_t) ~((arbint_limb_t) 0u);

  if (max_limbs == 0u)
    max_limbs = 1u;
  if (max_limbs > STOCH_MAX_INPUT_LIMBS)
    max_limbs = STOCH_MAX_INPUT_LIMBS;

  bucket = stoch_u64(rng) % 100u;
  if (allow_zero && bucket < 8u)
    n = 0u;
  else if (bucket < 30u)
    n = 1u + stoch_range(rng, 4u);
  else if (bucket < 70u)
    n = 1u + stoch_range(rng, (max_limbs < 24u) ? max_limbs : 24u);
  else
    n = 1u + stoch_range(rng, max_limbs);

  if (n == 0u) {
    ref_zero(out);
    return;
  }

  out->n = n;
  memset(out->limbs, 0, n * sizeof(arbint_limb_t));
  pattern = stoch_u64(rng) % 7u;
  switch (pattern) {
  case 0u:
    for (i = 0u; i < n; ++i)
      out->limbs[i] = stoch_limb(rng);
    break;
  case 1u:
    for (i = 0u; i < n; ++i)
      out->limbs[i] = maxv;
    break;
  case 2u:
    for (i = 0u; i < n; ++i)
      out->limbs[i] = (i & 1u) ? maxv : (arbint_limb_t) 0u;
    break;
  case 3u: {
    size_t bit = stoch_range(rng, n * ARBINT_LIMB_BITS);
    size_t li = bit / ARBINT_LIMB_BITS;
    size_t bi = bit % ARBINT_LIMB_BITS;
    out->limbs[li] = ((arbint_limb_t) 1u) << bi;
    if (bit != 0u && (stoch_u64(rng) & 1u)) {
      size_t j;
      for (j = 0u; j < li; ++j)
        out->limbs[j] = maxv;
    }
    break;
  }
  case 4u:
    for (i = 0u; i < n; ++i) {
      out->limbs[i] = (i < n / 2u) ? maxv : stoch_limb(rng);
    }
    break;
  case 5u:
    for (i = 0u; i < n; ++i)
      out->limbs[i] = (arbint_limb_t) ((stoch_u64(rng) & 3u) ? 0u : maxv);
    break;
  default:
    for (i = 0u; i < n; ++i) {
      uint64_t roll = stoch_u64(rng);
      if ((roll & 7u) == 0u)
        out->limbs[i] = maxv;
      else if ((roll & 7u) == 1u)
        out->limbs[i] = (arbint_limb_t) 1u;
      else
        out->limbs[i] = stoch_limb(rng);
    }
    break;
  }

  if (out->limbs[n - 1u] == (arbint_limb_t) 0u)
    out->limbs[n - 1u] = (arbint_limb_t) 1u
                         << stoch_range(rng, ARBINT_LIMB_BITS);

  out->sign = (stoch_u64(rng) & 1u) ? 1 : -1;
  ref_norm(out);
}

static void test_stochastic_arithmetic(void) {
  enum {
    ADD_SUB_ITERS = 3200,
    MUL_ITERS = 900,
    TDIV_ITERS = 1400,
    TDIV_U32_ITERS = 900,
    TDIV_I32_ITERS = 900
  };
  const uint64_t seed = 0xcafef00d1234abcdull;
  stoch_rng_t rng;
  arbint_ctx_t ctx;
  arbint_t a;
  arbint_t b;
  arbint_t c;
  arbint_t q;
  arbint_t r;
  arbint_t ncopy;
  refint_t ra;
  refint_t rb;
  refint_t rc;
  refint_t expect;
  refint_t expect_q;
  refint_t expect_r;
  size_t i;

  stoch_seed(&rng, seed);

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(q, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(r, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(ncopy, &ctx), ARBINT_OK);

  for (i = 0u; i < ADD_SUB_ITERS && g_failures == 0; ++i) {
    size_t lim = (i % 19u == 0u)  ? STOCH_MAX_INPUT_LIMBS
                 : (i % 5u == 0u) ? 48u
                                  : 4u + stoch_range(&rng, 28u);
    uint32_t u32v = (uint32_t) stoch_u64(&rng);
    int32_t i32v = (int32_t) stoch_u64(&rng);

    ref_random(&ra, &rng, lim, 1);
    ref_random(&rb, &rng, lim, 1);

    CHECK(ref_add(&expect, &ra, &rb));
    CHECK(ref_sub(&expect_q, &ra, &rb));

    CHECK_EQ_I(ref_to_arbint(a, &ra), ARBINT_OK);
    CHECK_EQ_I(ref_to_arbint(b, &rb), ARBINT_OK);

    CHECK_EQ_I(arbint_add(c, a, b), ARBINT_OK);
    check_ref_match("add", i, &ra, &rb, &expect, c);

    CHECK_EQ_I(ref_to_arbint(a, &ra), ARBINT_OK);
    CHECK_EQ_I(ref_to_arbint(b, &rb), ARBINT_OK);
    CHECK_EQ_I(arbint_add(a, a, b), ARBINT_OK);
    check_ref_match("add_alias_a", i, &ra, &rb, &expect, a);

    CHECK_EQ_I(ref_to_arbint(a, &ra), ARBINT_OK);
    CHECK_EQ_I(ref_to_arbint(b, &rb), ARBINT_OK);
    CHECK_EQ_I(arbint_add(b, a, b), ARBINT_OK);
    check_ref_match("add_alias_b", i, &ra, &rb, &expect, b);

    CHECK_EQ_I(ref_to_arbint(a, &ra), ARBINT_OK);
    CHECK_EQ_I(ref_to_arbint(b, &rb), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(c, a, b), ARBINT_OK);
    check_ref_match("sub", i, &ra, &rb, &expect_q, c);

    CHECK_EQ_I(ref_to_arbint(a, &ra), ARBINT_OK);
    CHECK_EQ_I(ref_to_arbint(b, &rb), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(a, a, b), ARBINT_OK);
    check_ref_match("sub_alias_a", i, &ra, &rb, &expect_q, a);

    CHECK_EQ_I(ref_to_arbint(a, &ra), ARBINT_OK);
    CHECK_EQ_I(ref_to_arbint(b, &rb), ARBINT_OK);
    CHECK_EQ_I(arbint_sub(b, a, b), ARBINT_OK);
    check_ref_match("sub_alias_b", i, &ra, &rb, &expect_q, b);

    ref_from_u32(&rc, u32v);
    CHECK(ref_add(&expect, &ra, &rc));
    CHECK(ref_sub(&expect_q, &ra, &rc));

    CHECK_EQ_I(ref_to_arbint(a, &ra), ARBINT_OK);
    CHECK_EQ_I(arbint_add_u32(c, a, u32v), ARBINT_OK);
    check_ref_match("add_u32", i, &ra, &rc, &expect, c);

    CHECK_EQ_I(ref_to_arbint(a, &ra), ARBINT_OK);
    CHECK_EQ_I(arbint_sub_u32(c, a, u32v), ARBINT_OK);
    check_ref_match("sub_u32", i, &ra, &rc, &expect_q, c);

    ref_from_i32(&rc, i32v);
    CHECK(ref_add(&expect, &ra, &rc));
    CHECK(ref_sub(&expect_q, &ra, &rc));

    CHECK_EQ_I(ref_to_arbint(a, &ra), ARBINT_OK);
    CHECK_EQ_I(arbint_add_i32(c, a, i32v), ARBINT_OK);
    check_ref_match("add_i32", i, &ra, &rc, &expect, c);

    CHECK_EQ_I(ref_to_arbint(a, &ra), ARBINT_OK);
    CHECK_EQ_I(arbint_sub_i32(c, a, i32v), ARBINT_OK);
    check_ref_match("sub_i32", i, &ra, &rc, &expect_q, c);
  }

  for (i = 0u; i < MUL_ITERS && g_failures == 0; ++i) {
    size_t lim_a = (i % 13u == 0u) ? STOCH_MAX_INPUT_LIMBS / 2u
                                   : 2u + stoch_range(&rng, 40u);
    size_t lim_b = (i % 11u == 0u) ? STOCH_MAX_INPUT_LIMBS / 2u
                                   : 2u + stoch_range(&rng, 40u);

    ref_random(&ra, &rng, lim_a, 1);
    ref_random(&rb, &rng, lim_b, 1);
    CHECK(ref_mul(&expect, &ra, &rb));

    CHECK_EQ_I(ref_to_arbint(a, &ra), ARBINT_OK);
    CHECK_EQ_I(ref_to_arbint(b, &rb), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(c, a, b), ARBINT_OK);
    check_ref_match("mul", i, &ra, &rb, &expect, c);

    CHECK_EQ_I(ref_to_arbint(a, &ra), ARBINT_OK);
    CHECK_EQ_I(ref_to_arbint(b, &rb), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(a, a, b), ARBINT_OK);
    check_ref_match("mul_alias_a", i, &ra, &rb, &expect, a);

    CHECK_EQ_I(ref_to_arbint(a, &ra), ARBINT_OK);
    CHECK_EQ_I(ref_to_arbint(b, &rb), ARBINT_OK);
    CHECK_EQ_I(arbint_mul(b, a, b), ARBINT_OK);
    check_ref_match("mul_alias_b", i, &ra, &rb, &expect, b);

    CHECK(ref_mul(&expect_q, &ra, &ra));
    CHECK_EQ_I(ref_to_arbint(a, &ra), ARBINT_OK);
    CHECK_EQ_I(arbint_sqr(c, a), ARBINT_OK);
    check_ref_match("sqr", i, &ra, &ra, &expect_q, c);

    CHECK_EQ_I(ref_to_arbint(a, &ra), ARBINT_OK);
    CHECK_EQ_I(arbint_sqr(a, a), ARBINT_OK);
    check_ref_match("sqr_alias", i, &ra, &ra, &expect_q, a);
  }

  for (i = 0u; i < TDIV_ITERS; ++i) {
    size_t lim_n =
        (i % 17u == 0u) ? STOCH_MAX_INPUT_LIMBS : 2u + stoch_range(&rng, 54u);
    size_t lim_d = (i % 23u == 0u) ? STOCH_MAX_INPUT_LIMBS / 2u
                                   : 1u + stoch_range(&rng, 36u);

    ref_random(&ra, &rng, lim_n, 1);
    do {
      ref_random(&rb, &rng, lim_d, 1);
    } while (rb.sign == 0);

    CHECK(ref_tdiv(&expect_q, &expect_r, &ra, &rb));

    CHECK_EQ_I(ref_to_arbint(a, &ra), ARBINT_OK);
    CHECK_EQ_I(ref_to_arbint(b, &rb), ARBINT_OK);

    CHECK_EQ_I(arbint_tdiv_qr(q, r, a, b), ARBINT_OK);
    check_ref_match("tdiv_qr_q", i, &ra, &rb, &expect_q, q);
    check_ref_match("tdiv_qr_r", i, &ra, &rb, &expect_r, r);

    CHECK_EQ_I(arbint_tdiv_q(c, a, b), ARBINT_OK);
    check_ref_match("tdiv_q", i, &ra, &rb, &expect_q, c);
    CHECK_EQ_I(arbint_tdiv_r(c, a, b), ARBINT_OK);
    check_ref_match("tdiv_r", i, &ra, &rb, &expect_r, c);

    if (!arbint_is_zero(r)) {
      CHECK_EQ_I(arbint_signum(r), arbint_signum(a));
      CHECK(arbint_cmpabs(r, b) < 0);
    }

    CHECK_EQ_I(arbint_set(ncopy, a), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_qr(ncopy, r, ncopy, b), ARBINT_OK);
    check_ref_match("tdiv_alias_q", i, &ra, &rb, &expect_q, ncopy);
    check_ref_match("tdiv_alias_q_r", i, &ra, &rb, &expect_r, r);

    CHECK_EQ_I(arbint_set(ncopy, a), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_qr(q, ncopy, ncopy, b), ARBINT_OK);
    check_ref_match("tdiv_alias_r_q", i, &ra, &rb, &expect_q, q);
    check_ref_match("tdiv_alias_r", i, &ra, &rb, &expect_r, ncopy);
  }

  for (i = 0u; i < TDIV_U32_ITERS; ++i) {
    uint32_t d = (uint32_t) stoch_u64(&rng) | 1u;

    ref_random(&ra, &rng, STOCH_MAX_INPUT_LIMBS, 1);
    ref_from_u32(&rb, d);
    CHECK(ref_tdiv(&expect_q, &expect_r, &ra, &rb));

    CHECK_EQ_I(ref_to_arbint(a, &ra), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_qr_u32(q, r, a, d), ARBINT_OK);
    check_ref_match("tdiv_qr_u32_q", i, &ra, &rb, &expect_q, q);
    check_ref_match("tdiv_qr_u32_r", i, &ra, &rb, &expect_r, r);
  }

  for (i = 0u; i < TDIV_I32_ITERS; ++i) {
    int32_t d;
    do {
      d = (int32_t) stoch_u64(&rng);
    } while (d == 0);

    ref_random(&ra, &rng, STOCH_MAX_INPUT_LIMBS, 1);
    ref_from_i32(&rb, d);
    CHECK(ref_tdiv(&expect_q, &expect_r, &ra, &rb));

    CHECK_EQ_I(ref_to_arbint(a, &ra), ARBINT_OK);
    CHECK_EQ_I(arbint_tdiv_qr_i32(q, r, a, d), ARBINT_OK);
    check_ref_match("tdiv_qr_i32_q", i, &ra, &rb, &expect_q, q);
    check_ref_match("tdiv_qr_i32_r", i, &ra, &rb, &expect_r, r);
  }

  arbint_clear(a);
  arbint_clear(b);
  arbint_clear(c);
  arbint_clear(q);
  arbint_clear(r);
  arbint_clear(ncopy);
}

int main(void) {
  ARBINT_TEST_START();
  test_stochastic_arithmetic();
  ARBINT_TEST_FINISH("test_stochastic");
}
