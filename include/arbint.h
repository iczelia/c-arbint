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

/*  The following features are currently unimplemented:
      - arbint_set_str, arbint_get_str
      - arbint_moebius, arbint_totient, arbint_carmichael,
      - arbint_xgcd, arbint_legendre, arbint_jacobi,
      - arbint_kronecker, arbint_removefactor_u32,
      - arbint_is_power, arbint_isprime,
      - arbint_nextprime, arbint_prevprime, arbint_inv_mod,
      - arbint_inv_mod_u32  */

#ifndef ARBINT_H
#define ARBINT_H

#include <stddef.h>
#include <stdint.h>

/* ---------------- Symbol visibility ---------------- */
#if defined(ARBINT_STATIC)
  #define ARBINT_API
#elif defined(_WIN32) || defined(__CYGWIN__)
  #if defined(ARBINT_BUILD_DLL)
    #if defined(__GNUC__)
      #define ARBINT_API __attribute__((dllexport))
    #else
      #define ARBINT_API __declspec(dllexport)
    #endif /* defined(__GNUC__) */
  #elif defined(ARBINT_USE_DLL)
    #if defined(__GNUC__)
      #define ARBINT_API __attribute__((dllimport))
    #else
      #define ARBINT_API __declspec(dllimport)
    #endif /* defined(__GNUC__) */
  #else
    #define ARBINT_API
  #endif /* defined(ARBINT_BUILD_DLL) */
#else
  #if defined(__GNUC__) || defined(__clang__)
    #define ARBINT_API __attribute__((visibility("default")))
  #else
    #define ARBINT_API
  #endif /* defined(__GNUC__) || defined(__clang__) */
#endif   /* defined(ARBINT_STATIC) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ---------------- Errors ---------------- */
typedef enum arbint_err {
  ARBINT_OK = 0,
  ARBINT_EINVAL,
  ARBINT_ENOMEM,
  ARBINT_EOVERFLOW,
  ARBINT_ESIGN,
  ARBINT_EDOM, /* domain error (e.g., sqrt of negative if disallowed) */
  ARBINT_EZERO /* division by zero */
} arbint_err_t;

/* ---------------- Context / allocator ---------------- */
/*  Reallocator contract: behaves like realloc(ptr, new_size).
    - ptr == NULL allocates.
    - new_size == 0 frees and may return NULL.
    - ud is the allocator user cookie from arbint_alloc_t.  */
typedef void * (*arbint_realloc_fn)(void * ud, void * ptr, size_t new_size);

typedef struct arbint_alloc {
  void * ud;
  arbint_realloc_fn realloc;
} arbint_alloc_t;

typedef struct arbint_ctx {
  arbint_alloc_t a;
  uint32_t flags;
  uint32_t rng_flags;
} arbint_ctx_t;

/* Initialize context with libc allocator and default flags. */
ARBINT_API arbint_err_t arbint_ctx_init_default(arbint_ctx_t * ctx);
/* Initialize context with caller allocator and explicit flags. */
ARBINT_API arbint_err_t arbint_ctx_init(arbint_ctx_t * ctx,
                                        const arbint_alloc_t * a,
                                        uint32_t flags);
/* Release context-owned transient state and drop global runtime caches. */
ARBINT_API void arbint_ctx_clear(arbint_ctx_t * ctx);
/* Drop global runtime caches (idempotent). */
ARBINT_API void arbint_drop_caches(void);

/* ---------------- Core type ---------------- */
/*  Internal normalization invariant:
    - _sz == 0 represents zero.
    - abs(_sz) is the number of used limbs.
    - sign(_sz) is the sign of the integer.
    - limbs are stored little-endian in _ptr.  */
typedef struct {
  size_t _cap;         /* allocated limbs in _ptr */
  ptrdiff_t _sz;       /* used limbs; sign is sign of _sz */
  void * _ptr;         /* little-endian limbs */
  arbint_ctx_t * _ctx; /* for allocations; optional (can be NULL if no further
                          allocations needed) */
} _arbint_struct;

typedef _arbint_struct arbint_t[1];

/* ---------------- Lifecycle / memory ---------------- */
/* Initialize x to numeric zero and bind it to ctx for future allocations. */
ARBINT_API arbint_err_t arbint_init(arbint_t x, arbint_ctx_t * ctx);
/*  Initialize a NULL-terminated list of arbint_t values with the same
    context. If any initialization fails, all previously initialized values
    are cleared. Usage: arbint_init_all(ctx, a, b, c, (arbint_t *) NULL).  */
ARBINT_API arbint_err_t arbint_init_all(arbint_ctx_t * ctx, arbint_t a, ...);
/* Release storage owned by x and reset it to an empty state. */
ARBINT_API void arbint_clear(arbint_t x);
/*  Clear a NULL-terminated list of arbint_t values.
    Usage: arbint_clear_all(a, b, c, (arbint_t *) NULL).  */
ARBINT_API void arbint_clear_all(arbint_t a, ...);
/* Return the context currently associated with x (may be NULL). */
ARBINT_API arbint_ctx_t * arbint_get_ctx(const arbint_t x);
/*  Rebind x to ctx. If x owns allocated storage, implementation may migrate
    it so future resize/clear calls remain allocator-safe.  */
ARBINT_API arbint_err_t arbint_set_ctx(arbint_t x, arbint_ctx_t * ctx);

/* Set x to exact zero without releasing reserved capacity. */
ARBINT_API void arbint_zero(arbint_t x);
/* Ensure capacity of x is at least new_cap limbs (value-preserving). */
ARBINT_API arbint_err_t arbint_resize(arbint_t x, size_t new_cap);
/* Swap complete object state (value, capacity, and context pointer). */
ARBINT_API void arbint_swap(arbint_t a, arbint_t b);

/* ---------------- Assignment (set) ---------------- */
ARBINT_API arbint_err_t arbint_set(arbint_t rop, const arbint_t op);

ARBINT_API arbint_err_t arbint_set_i32(arbint_t rop, int32_t v);
ARBINT_API arbint_err_t arbint_set_u32(arbint_t rop, uint32_t v);

/* Parse optional sign + digits in base 2..36 and assign the parsed value. */
ARBINT_API arbint_err_t arbint_set_str(arbint_t rop, const char * s, int base);
/*  Format op in base 2..36 into a newly allocated NUL-terminated string.
    Allocation is performed through op's context allocator and ownership
    transfers to the caller.  */
ARBINT_API arbint_err_t arbint_get_str(const arbint_t op, char ** out_str,
                                       int base);

/* ---------------- Extraction (get) ---------------- */
/* Convert to i32, EOVERFLOW if doesn't fit */
ARBINT_API arbint_err_t arbint_get_i32(const arbint_t op, int32_t * out);
/* Convert to u32. ESIGN if negative; EOVERFLOW if not fit */
ARBINT_API arbint_err_t arbint_get_u32(const arbint_t op, uint32_t * out);

/* fits queries (no allocation) */
ARBINT_API int arbint_fits_u8(const arbint_t x);
ARBINT_API int arbint_fits_u16(const arbint_t x);
ARBINT_API int arbint_fits_u32(const arbint_t x);
ARBINT_API int arbint_fits_u64(const arbint_t x);
ARBINT_API int arbint_fits_i8(const arbint_t x);
ARBINT_API int arbint_fits_i16(const arbint_t x);
ARBINT_API int arbint_fits_i32(const arbint_t x);
ARBINT_API int arbint_fits_i64(const arbint_t x);

/* ---------------- Sign / magnitude ---------------- */
ARBINT_API int arbint_signum(const arbint_t x); /* -1,0,+1 */
ARBINT_API int arbint_is_zero(const arbint_t x);
ARBINT_API int arbint_is_one(const arbint_t x);
ARBINT_API int arbint_is_neg(const arbint_t x);
ARBINT_API int arbint_is_odd(const arbint_t x);
ARBINT_API int arbint_is_even(const arbint_t x);

ARBINT_API arbint_err_t arbint_abs(arbint_t rop, const arbint_t op);
ARBINT_API arbint_err_t arbint_neg(arbint_t rop, const arbint_t op);

/* ---------------- Comparison ---------------- */
ARBINT_API int arbint_cmp(const arbint_t a, const arbint_t b); /* -1/0/+1 */
ARBINT_API int arbint_cmpabs(const arbint_t a, const arbint_t b);

ARBINT_API int arbint_cmp_i32(const arbint_t a, int32_t b);
ARBINT_API int arbint_cmp_u32(const arbint_t a, uint32_t b);

/* predicate helpers */
ARBINT_API int arbint_eq(const arbint_t a, const arbint_t b);
ARBINT_API int arbint_ne(const arbint_t a, const arbint_t b);
ARBINT_API int arbint_lt(const arbint_t a, const arbint_t b);
ARBINT_API int arbint_le(const arbint_t a, const arbint_t b);
ARBINT_API int arbint_gt(const arbint_t a, const arbint_t b);
ARBINT_API int arbint_ge(const arbint_t a, const arbint_t b);

ARBINT_API int arbint_eq_i32(const arbint_t a, int32_t b);
ARBINT_API int arbint_ne_i32(const arbint_t a, int32_t b);
ARBINT_API int arbint_lt_i32(const arbint_t a, int32_t b);
ARBINT_API int arbint_le_i32(const arbint_t a, int32_t b);
ARBINT_API int arbint_gt_i32(const arbint_t a, int32_t b);
ARBINT_API int arbint_ge_i32(const arbint_t a, int32_t b);
ARBINT_API int arbint_eq_u32(const arbint_t a, uint32_t b);
ARBINT_API int arbint_ne_u32(const arbint_t a, uint32_t b);
ARBINT_API int arbint_lt_u32(const arbint_t a, uint32_t b);
ARBINT_API int arbint_le_u32(const arbint_t a, uint32_t b);
ARBINT_API int arbint_gt_u32(const arbint_t a, uint32_t b);
ARBINT_API int arbint_ge_u32(const arbint_t a, uint32_t b);

/* ---------------- Basic arithmetic ---------------- */
/*  Aliasing contract:
    - Any arbint_t argument may alias any input argument.
    - Functions may internally reuse storage from aliased inputs.
    - For functions with multiple outputs (e.g., *_qr), output aliases are
      accepted; because only one object exists, the final value in that object
      is whichever output is written last by that routine.  */
ARBINT_API arbint_err_t arbint_add(arbint_t rop, const arbint_t a,
                                   const arbint_t b);
ARBINT_API arbint_err_t arbint_sub(arbint_t rop, const arbint_t a,
                                   const arbint_t b);
ARBINT_API arbint_err_t arbint_mul(arbint_t rop, const arbint_t a,
                                   const arbint_t b);
ARBINT_API arbint_err_t arbint_sqr(arbint_t rop, const arbint_t a);

/* small-int operand variants (destination is always arbint_t) */
ARBINT_API arbint_err_t arbint_add_i32(arbint_t rop, const arbint_t a,
                                       int32_t b);
ARBINT_API arbint_err_t arbint_add_u32(arbint_t rop, const arbint_t a,
                                       uint32_t b);
ARBINT_API arbint_err_t arbint_sub_i32(arbint_t rop, const arbint_t a,
                                       int32_t b);
ARBINT_API arbint_err_t arbint_sub_u32(arbint_t rop, const arbint_t a,
                                       uint32_t b);
ARBINT_API arbint_err_t arbint_mul_i32(arbint_t rop, const arbint_t a,
                                       int32_t b);
ARBINT_API arbint_err_t arbint_mul_u32(arbint_t rop, const arbint_t a,
                                       uint32_t b);

/* fused forms: a += b*c, a -= b*c */
ARBINT_API arbint_err_t arbint_addmul(arbint_t a, const arbint_t b,
                                      const arbint_t c);
ARBINT_API arbint_err_t arbint_submul(arbint_t a, const arbint_t b,
                                      const arbint_t c);
ARBINT_API arbint_err_t arbint_addmul_u32(arbint_t a, const arbint_t b,
                                          uint32_t c);
ARBINT_API arbint_err_t arbint_submul_u32(arbint_t a, const arbint_t b,
                                          uint32_t c);
ARBINT_API arbint_err_t arbint_addmul_i32(arbint_t a, const arbint_t b,
                                          int32_t c);
ARBINT_API arbint_err_t arbint_submul_i32(arbint_t a, const arbint_t b,
                                          int32_t c);

/* ---------------- Division, remainder, rounding modes ---------------- */
/*  Truncated quotient/remainder: q = trunc(n/d), r = n - q*d (|r| < |d|,
    sign(r)=sign(n) or r=0).  */
ARBINT_API arbint_err_t arbint_tdiv_qr(arbint_t q, arbint_t r,
                                       const arbint_t n, const arbint_t d);
ARBINT_API arbint_err_t arbint_tdiv_q(arbint_t q, const arbint_t n,
                                      const arbint_t d);
ARBINT_API arbint_err_t arbint_tdiv_r(arbint_t r, const arbint_t n,
                                      const arbint_t d);

/*  Floor division: q = floor(n/d), r = n - q*d (0 <= r < |d| if d>0;
    adjust consistently for d<0).  */
ARBINT_API arbint_err_t arbint_fdiv_qr(arbint_t q, arbint_t r,
                                       const arbint_t n, const arbint_t d);
ARBINT_API arbint_err_t arbint_fdiv_q(arbint_t q, const arbint_t n,
                                      const arbint_t d);
ARBINT_API arbint_err_t arbint_fdiv_r(arbint_t r, const arbint_t n,
                                      const arbint_t d);

/* Ceil division: q = ceil(n/d), r = n - q*d */
ARBINT_API arbint_err_t arbint_cdiv_qr(arbint_t q, arbint_t r,
                                       const arbint_t n, const arbint_t d);
ARBINT_API arbint_err_t arbint_cdiv_q(arbint_t q, const arbint_t n,
                                      const arbint_t d);
ARBINT_API arbint_err_t arbint_cdiv_r(arbint_t r, const arbint_t n,
                                      const arbint_t d);

/* small-int divisors/moduli variants (d == 0 returns ARBINT_EZERO) */
ARBINT_API arbint_err_t arbint_tdiv_qr_u32(arbint_t q, arbint_t r,
                                           const arbint_t n, uint32_t d);
ARBINT_API arbint_err_t arbint_tdiv_q_u32(arbint_t q, const arbint_t n,
                                          uint32_t d);
ARBINT_API arbint_err_t arbint_tdiv_r_u32(arbint_t r, const arbint_t n,
                                          uint32_t d);
ARBINT_API arbint_err_t arbint_tdiv_qr_i32(arbint_t q, arbint_t r,
                                           const arbint_t n, int32_t d);
ARBINT_API arbint_err_t arbint_tdiv_q_i32(arbint_t q, const arbint_t n,
                                          int32_t d);
ARBINT_API arbint_err_t arbint_tdiv_r_i32(arbint_t r, const arbint_t n,
                                          int32_t d);
ARBINT_API arbint_err_t arbint_fdiv_qr_u32(arbint_t q, arbint_t r,
                                           const arbint_t n, uint32_t d);
ARBINT_API arbint_err_t arbint_fdiv_q_u32(arbint_t q, const arbint_t n,
                                          uint32_t d);
ARBINT_API arbint_err_t arbint_fdiv_r_u32(arbint_t r, const arbint_t n,
                                          uint32_t d);

/* Divisibility: on success writes out = 1 iff d divides n exactly. */
ARBINT_API arbint_err_t arbint_divisible(const arbint_t n, const arbint_t d,
                                         int * out);
ARBINT_API arbint_err_t arbint_divisible_u32(const arbint_t n, uint32_t d,
                                             int * out);

/* ---------------- Powers ---------------- */
ARBINT_API arbint_err_t arbint_pow_u32(arbint_t rop, const arbint_t base,
                                       uint32_t exp);

/*  Modular exponentiation: rop = base^exp mod mod.
    Truncated division semantics. exp >= 0 and mod != 0.  */
ARBINT_API arbint_err_t arbint_pow_tmod(arbint_t rop, const arbint_t base,
                                        const arbint_t exp,
                                        const arbint_t mod);
ARBINT_API arbint_err_t arbint_pow_u32_tmod(arbint_t rop, const arbint_t base,
                                            uint32_t exp, const arbint_t mod);

/*  Modular exponentiation with u32 modulus: rop = base^exp mod mod.
    Uses precomputed Barrett reduction for O(1) reciprocal overhead.
    Truncated division semantics.  */
ARBINT_API arbint_err_t arbint_pow_u32u32_tmod(arbint_t rop,
                                               const arbint_t base,
                                               uint32_t exp, uint32_t mod);

/* ---------------- Shifts ---------------- */
ARBINT_API arbint_err_t arbint_shl(arbint_t rop, const arbint_t a,
                                   uint32_t k); /* a * 2^k */
ARBINT_API arbint_err_t arbint_shr(arbint_t rop, const arbint_t a,
                                   uint32_t k); /* trunc toward 0 */

/* ---------------- Bit operations / queries ---------------- */
/*  Digits in base; for x=0 returns 1.  */
ARBINT_API size_t arbint_sizeinbase(const arbint_t x, int base);
/*  Number of significant bits; 0 for x=0.  */
ARBINT_API size_t arbint_nbits(const arbint_t x);

/*  Bit semantics use two's-complement integers.
    - testbit/or/and/xor/not/setbit/clrbit behave as if over infinite
      sign-extended two's-complement representations, then normalize.
    - not(a) is therefore equivalent to (-a - 1).
    For APIs that must return finite counts, define canonical width W(x):
    the smallest positive multiple of the implementation limb width that can
    represent x in two's-complement form.  */
ARBINT_API arbint_err_t arbint_testbit(const arbint_t x, size_t bit_index,
                                       int * out); /* out = 0/1 */
ARBINT_API arbint_err_t arbint_setbit(arbint_t x, size_t bit_index);
ARBINT_API arbint_err_t arbint_clrbit(arbint_t x, size_t bit_index);
/*  Count trailing zeros; EDOM if x=0.  */
ARBINT_API arbint_err_t arbint_ctz(const arbint_t x, size_t * out);
/*  Count leading zero bits within canonical width W(x); for x<0 this is 0.
    Returns EDOM if x=0.  */
ARBINT_API arbint_err_t arbint_clz(const arbint_t x, size_t * out);
/*  Count set bits in the low W(x) bits of x.  */
ARBINT_API arbint_err_t arbint_popcount(const arbint_t x, size_t * out);

/*  Hamming distance between low max(W(a), W(b)) bits of a and b.  */
ARBINT_API arbint_err_t arbint_hammingdist(const arbint_t a, const arbint_t b,
                                           size_t * out);

ARBINT_API arbint_err_t arbint_or(arbint_t rop, const arbint_t a,
                                  const arbint_t b);
ARBINT_API arbint_err_t arbint_and(arbint_t rop, const arbint_t a,
                                   const arbint_t b);
ARBINT_API arbint_err_t arbint_xor(arbint_t rop, const arbint_t a,
                                   const arbint_t b);
/*  Two's-complement bitwise not.  */
ARBINT_API arbint_err_t arbint_not(arbint_t rop, const arbint_t a);

/*  Immediate-operand bitwise operations.  */
ARBINT_API arbint_err_t arbint_and_u32(arbint_t rop, const arbint_t a,
                                       uint32_t b);
ARBINT_API arbint_err_t arbint_and_i32(arbint_t rop, const arbint_t a,
                                       int32_t b);
ARBINT_API arbint_err_t arbint_or_u32(arbint_t rop, const arbint_t a,
                                      uint32_t b);
ARBINT_API arbint_err_t arbint_or_i32(arbint_t rop, const arbint_t a,
                                      int32_t b);
ARBINT_API arbint_err_t arbint_xor_u32(arbint_t rop, const arbint_t a,
                                       uint32_t b);
ARBINT_API arbint_err_t arbint_xor_i32(arbint_t rop, const arbint_t a,
                                       int32_t b);

/* ---------------- Number theory ---------------- */
ARBINT_API arbint_err_t arbint_gcd(arbint_t g, const arbint_t a,
                                   const arbint_t b);
ARBINT_API arbint_err_t arbint_lcm(arbint_t l, const arbint_t a,
                                   const arbint_t b);
ARBINT_API arbint_err_t arbint_gcd_u32(arbint_t g, const arbint_t a,
                                       uint32_t b);
ARBINT_API arbint_err_t arbint_lcm_u32(arbint_t l, const arbint_t a,
                                       uint32_t b);
/*  out=0 if n has squared prime factor; else out=(-1)^(number of distinct
 * prime factors).  */
ARBINT_API arbint_err_t arbint_moebius(int * out, const arbint_t n);
/*  Count of positive integers <= n that are coprime to n.  */
ARBINT_API arbint_err_t arbint_totient(arbint_t rop, const arbint_t n);
/*  lcm of lambda(p_i^{k_i}) for prime factorization n = prod p_i^{k_i};
 * lambda(p^k) = p^{k-1}*(p-1) for odd p or k<=2; else lambda(2^k) = 2^{k-2}.
 */
ARBINT_API arbint_err_t arbint_carmichael(arbint_t rop, const arbint_t n);

/* extended gcd: g = ax + by */
ARBINT_API arbint_err_t arbint_xgcd(arbint_t g, arbint_t x, arbint_t y,
                                    const arbint_t a, const arbint_t b);

/* Jacobi/Legendre/Kronecker symbols return in *out: -1,0,+1 */
ARBINT_API arbint_err_t arbint_jacobi(const arbint_t a, const arbint_t n,
                                      int * out);
ARBINT_API arbint_err_t arbint_legendre(const arbint_t a, const arbint_t p,
                                        int * out);
ARBINT_API arbint_err_t arbint_kronecker(const arbint_t a, const arbint_t n,
                                         int * out);

/* removefactor: write a = p^k * rest; return k and rest */
ARBINT_API arbint_err_t arbint_removefactor_u32(arbint_t rest,
                                                const arbint_t a, uint32_t p,
                                                uint32_t * k);

/* factorial, binomial, fibonacci, lucas */
ARBINT_API arbint_err_t arbint_fac_u32(arbint_t rop, uint32_t n);
ARBINT_API arbint_err_t arbint_bin_u32u32(arbint_t rop, uint32_t n,
                                          uint32_t k);
ARBINT_API arbint_err_t arbint_fib_u32(arbint_t rop, uint32_t n);
ARBINT_API arbint_err_t arbint_lucas_u32(arbint_t rop, uint32_t n);

/* ---------------- Roots / perfect powers ---------------- */
/*  floor(sqrt(a)) for a>=0.  */
ARBINT_API arbint_err_t arbint_isqrt(arbint_t rop, const arbint_t a);
/*  out=1 if perfect square.  */
ARBINT_API arbint_err_t arbint_is_square(const arbint_t a, int * out);
/*  out=1 iff exists integers b, k>=2 such that a = b^k.  */
ARBINT_API arbint_err_t arbint_is_power(const arbint_t a, int * out);
/*  floor(a^(1/k)) for a>=0; EDOM if k=0.  */
ARBINT_API arbint_err_t arbint_root(arbint_t rop, const arbint_t a,
                                    uint32_t k);

/* ---------------- Primality and neighboring primes ---------------- */
/*  out=0/1; reps controls accuracy.  */
ARBINT_API arbint_err_t arbint_isprime(const arbint_t n, int reps, int * out);
ARBINT_API arbint_err_t arbint_nextprime(arbint_t rop, const arbint_t n);
ARBINT_API arbint_err_t arbint_prevprime(arbint_t rop, const arbint_t n);

/* ---------------- Modular inverses ---------------- */
/*  rop = a^{-1} mod mod; EDOM if gcd(a,mod)!=1.  */
ARBINT_API arbint_err_t arbint_inv_mod(arbint_t rop, const arbint_t a,
                                       const arbint_t mod);
/*  rop = a^{-1} mod mod; EDOM if gcd(a,mod)!=1.  */
ARBINT_API arbint_err_t arbint_inv_mod_u32(arbint_t rop, const arbint_t a,
                                           uint32_t mod);

/* ---------------- Import + Export ---------------- */
/* Import magnitude from little-endian byte buffer and apply sign. */
ARBINT_API arbint_err_t arbint_import(arbint_t rop, const void * buf,
                                      size_t nbytes);
/* Export magnitude to a newly allocated little-endian byte buffer. */
ARBINT_API arbint_err_t arbint_export(const arbint_t op, void ** out_buf,
                                      size_t * out_nbytes);

/* ---------------- Hashing ---------------- */
/*  Hash canonical integer value into hash_len bytes at out_hash.
    Intended for hash tables only; algorithm and output are intentionally
    unstable across library versions and may also vary across process runs.  */
ARBINT_API arbint_err_t arbint_hash_slow(const arbint_t op, uint8_t * out_hash,
                                         size_t hash_len);
/*  Hash canonical integer value into a 32-bit unsigned integer value.
    Intended for hash tables only; algorithm and output are intentionally
    unstable across library versions and may also vary across process runs.  */
ARBINT_API arbint_err_t arbint_hash_fast(const arbint_t op,
                                         uint32_t * out_hash);

/* ---------------- Random Number Generation ---------------- */

/*  Mersenne Twister MT19937 PRNG state.
    Users should generally treat this as an opaque type and only manipulate
    via arbint_rng_init, arbint_rng_clear, arbint_urandomb, arbint_urandomm.
    State size: 2504 bytes (624 x 32-bit state + 1 x 32-bit index).  */
typedef struct arbint_rng {
  uint32_t mt[624]; /* MT19937 state array */
  unsigned int mti; /* Current index in state array */
} arbint_rng_t;

/*  Initialize RNG state.
    If seed is NULL, attempts to seed from platform entropy source:
      - Windows: BCryptGenRandom or CryptGenRandom
      - macOS/BSD: arc4random_buf
      - Linux/POSIX: /dev/urandom
    If seed is non-NULL, uses provided bytes as seed. seed_len must be
    a positive multiple of 4 (sizeof(uint32_t)).
    Returns:
      ARBINT_OK on success
      ARBINT_EINVAL if rng is NULL, or if seed is NULL and no platform
                    entropy source is available, or if seed_len is invalid  */
ARBINT_API arbint_err_t arbint_rng_init(arbint_rng_t * rng, const void * seed,
                                        size_t seed_len);

/*  Clear RNG state (zeros memory for security).  */
ARBINT_API void arbint_rng_clear(arbint_rng_t * rng);

/*  Generate uniform random integer in [0, 2^k).
    Fills rop with k random bits from MT19937 PRNG.
    Returns:
      ARBINT_OK on success
      ARBINT_EINVAL if rop or rng is NULL
      ARBINT_ENOMEM if allocation fails  */
ARBINT_API arbint_err_t arbint_urandomb(arbint_t rop, arbint_rng_t * rng,
                                        size_t k);

/*  Generate uniform random integer in [0, bound).
    Uses rejection sampling to ensure true uniform distribution with no
    modulo bias. Expected iterations: ~1.5 (worst case ~2).
    Returns:
      ARBINT_OK on success
      ARBINT_EINVAL if rop, rng, or bound is NULL
      ARBINT_EDOM if bound <= 0, or if rejection sampling fails after 256 tries
      ARBINT_ENOMEM if allocation fails  */
ARBINT_API arbint_err_t arbint_urandomm(arbint_t rop, arbint_rng_t * rng,
                                        const arbint_t bound);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* ARBINT_H */
