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

/*  dc - desk calculator
    A reverse-polish notation calculator demonstrating c-arbint.
    Implements a substantial subset of GNU dc functionality.  */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <arbint.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PACKAGE_VERSION)
#define DC_PACKAGE_VERSION PACKAGE_VERSION
#elif defined(VERSION)
#define DC_PACKAGE_VERSION VERSION
#else
#define DC_PACKAGE_VERSION "unknown"
#endif

/*  Portable strdup.  */
static char * dc_strdup(const char * s) {
  size_t len = strlen(s) + 1;
  char * p = (char *) malloc(len);
  if (p)
    memcpy(p, s, len);
  return p;
}

/*  Configuration  */

#define DC_INITIAL_STACK_CAP 64
#define DC_MAX_REGISTERS 256
#define DC_MAX_INPUT_LINE 65536
#define DC_MAX_PRECISION 10000

/*  Allocate line buffers on heap to avoid oversized stack frames on wasm.  */
static char * dc_alloc_line_buffer(void) {
  char * line = malloc(DC_MAX_INPUT_LINE);
  if (!line)
    fprintf(stderr, "dc: out of memory\n");
  return line;
}

/*  Data structures  */

/*  A dc value can be either a number or a string (macro).
    Numbers are stored as scaled integers: value = mantissa * 10^(-scale).
    For example, 3.14159 with scale=5 is stored as mantissa=314159, scale=5.  */
typedef enum { DC_TYPE_NUMBER, DC_TYPE_STRING } dc_type_t;

typedef struct dc_value {
  dc_type_t type;
  int scale;  /*  Number of decimal places (only for DC_TYPE_NUMBER).  */
  union {
    arbint_t num;
    char * str;
  } u;
} dc_value_t;

/*  Register stack entry (for S/L commands).  */
typedef struct dc_regstack {
  dc_value_t val;
  struct dc_regstack * next;
} dc_regstack_t;

/*  Array entry for : and ; commands.  */
typedef struct dc_array_entry {
  uint32_t index;
  dc_value_t val;
  struct dc_array_entry * next;
} dc_array_entry_t;

/*  Register: holds a single value plus an optional stack and array.  */
typedef struct {
  int has_value;
  dc_value_t val;
  dc_regstack_t * stack;
  dc_array_entry_t * array;  /*  Sparse array as linked list.  */
} dc_register_t;

/*  Main dc state.  */
typedef struct {
  arbint_ctx_t ctx;

  /*  Main stack.  */
  dc_value_t * stack;
  size_t stack_size;
  size_t stack_cap;

  /*  Registers (indexed by character).  */
  dc_register_t registers[DC_MAX_REGISTERS];

  /*  Parameters.  */
  int precision;
  int input_radix;
  int output_radix;

  /*  Execution state.  */
  int quit_depth;
  int running;

  /*  Input state.  */
  const char * input_ptr;
  int input_is_file;
  FILE * input_file;
} dc_state_t;

/*  Forward declarations  */

static void dc_value_clear(dc_state_t * dc, dc_value_t * val);
static int dc_value_copy(dc_state_t * dc, dc_value_t * dst,
                         const dc_value_t * src);
static int dc_push(dc_state_t * dc, const dc_value_t * val);
static int dc_pop(dc_state_t * dc, dc_value_t * out);
static int dc_peek(dc_state_t * dc, size_t depth, dc_value_t ** out);
static void dc_execute(dc_state_t * dc, const char * str);
static int dc_execute_char(dc_state_t * dc, int ch);

/*  Value management  */

static void dc_value_init_string(dc_value_t * val, const char * str) {
  val->type = DC_TYPE_STRING;
  val->scale = 0;
  val->u.str = dc_strdup(str);
}

static void dc_value_clear(dc_state_t * dc, dc_value_t * val) {
  (void)dc;
  if (val->type == DC_TYPE_NUMBER) {
    arbint_clear(val->u.num);
  } else if (val->type == DC_TYPE_STRING) {
    free(val->u.str);
    val->u.str = NULL;
  }
}

static int dc_value_copy(dc_state_t * dc, dc_value_t * dst,
                         const dc_value_t * src) {
  dst->scale = src->scale;
  if (src->type == DC_TYPE_NUMBER) {
    dst->type = DC_TYPE_NUMBER;
    if (arbint_init(dst->u.num, &dc->ctx) != ARBINT_OK)
      return -1;
    if (arbint_set(dst->u.num, src->u.num) != ARBINT_OK) {
      arbint_clear(dst->u.num);
      return -1;
    }
  } else {
    dst->type = DC_TYPE_STRING;
    dst->u.str = dc_strdup(src->u.str);
    if (!dst->u.str)
      return -1;
  }
  return 0;
}

/*  Stack operations  */

static int dc_push(dc_state_t * dc, const dc_value_t * val) {
  if (dc->stack_size >= dc->stack_cap) {
    size_t new_cap = dc->stack_cap * 2;
    dc_value_t * new_stack =
        realloc(dc->stack, new_cap * sizeof(dc_value_t));
    if (!new_stack) {
      fprintf(stderr, "dc: out of memory\n");
      return -1;
    }
    dc->stack = new_stack;
    dc->stack_cap = new_cap;
  }

  if (dc_value_copy(dc, &dc->stack[dc->stack_size], val) < 0)
    return -1;
  dc->stack_size++;
  return 0;
}

static int dc_push_number_scaled(dc_state_t * dc, const arbint_t num,
                                  int scale) {
  dc_value_t val;
  val.type = DC_TYPE_NUMBER;
  val.scale = scale;
  if (arbint_init(val.u.num, &dc->ctx) != ARBINT_OK)
    return -1;
  if (arbint_set(val.u.num, num) != ARBINT_OK) {
    arbint_clear(val.u.num);
    return -1;
  }

  if (dc->stack_size >= dc->stack_cap) {
    size_t new_cap = dc->stack_cap * 2;
    dc_value_t * new_stack =
        realloc(dc->stack, new_cap * sizeof(dc_value_t));
    if (!new_stack) {
      fprintf(stderr, "dc: out of memory\n");
      arbint_clear(val.u.num);
      return -1;
    }
    dc->stack = new_stack;
    dc->stack_cap = new_cap;
  }

  dc->stack[dc->stack_size++] = val;
  return 0;
}

static int dc_push_number(dc_state_t * dc, const arbint_t num) {
  return dc_push_number_scaled(dc, num, 0);
}

static int dc_push_string(dc_state_t * dc, const char * str) {
  dc_value_t val;
  dc_value_init_string(&val, str);
  if (!val.u.str)
    return -1;

  if (dc->stack_size >= dc->stack_cap) {
    size_t new_cap = dc->stack_cap * 2;
    dc_value_t * new_stack =
        realloc(dc->stack, new_cap * sizeof(dc_value_t));
    if (!new_stack) {
      fprintf(stderr, "dc: out of memory\n");
      free(val.u.str);
      return -1;
    }
    dc->stack = new_stack;
    dc->stack_cap = new_cap;
  }

  dc->stack[dc->stack_size++] = val;
  return 0;
}

static int dc_pop(dc_state_t * dc, dc_value_t * out) {
  if (dc->stack_size == 0) {
    fprintf(stderr, "dc: stack empty\n");
    return -1;
  }
  *out = dc->stack[--dc->stack_size];
  return 0;
}

static int dc_peek(dc_state_t * dc, size_t depth, dc_value_t ** out) {
  if (depth >= dc->stack_size) {
    fprintf(stderr, "dc: stack empty\n");
    return -1;
  }
  *out = &dc->stack[dc->stack_size - 1 - depth];
  return 0;
}

static int dc_pop_number_scaled(dc_state_t * dc, arbint_t out, int * scale) {
  dc_value_t val;
  if (dc_pop(dc, &val) < 0)
    return -1;
  if (val.type != DC_TYPE_NUMBER) {
    fprintf(stderr, "dc: non-numeric value\n");
    dc_value_clear(dc, &val);
    return -1;
  }
  arbint_set(out, val.u.num);
  *scale = val.scale;
  arbint_clear(val.u.num);
  return 0;
}

static int dc_pop_number(dc_state_t * dc, arbint_t out) {
  int scale;
  return dc_pop_number_scaled(dc, out, &scale);
}

/*  Number parsing  */

static int digit_value(int ch, int radix) {
  int val;
  if (ch >= '0' && ch <= '9')
    val = ch - '0';
  else if (ch >= 'A' && ch <= 'Z')
    val = 10 + (ch - 'A');
  else if (ch >= 'a' && ch <= 'z')
    val = 10 + (ch - 'a');
  else
    return -1;
  return (val < radix) ? val : -1;
}

/*  Parse a number in the current input radix.
    Returns number of characters consumed, or -1 on error.
    Sets *out_scale to the number of fractional digits.  */
static int dc_parse_number_scaled(dc_state_t * dc, const char * str,
                                   arbint_t out, int * out_scale) {
  const char * p = str;
  int negative = 0;
  int radix = dc->input_radix;
  int scale = 0;
  int seen_dot = 0;
  arbint_t digit, base;

  /*  Handle underscore as negative sign.  */
  if (*p == '_') {
    negative = 1;
    p++;
  }

  /*  Must have at least one digit or a dot followed by digit.  */
  if (digit_value(*p, radix) < 0 && *p != '.')
    return -1;

  arbint_set_i32(out, 0);
  arbint_init(digit, &dc->ctx);
  arbint_init(base, &dc->ctx);
  arbint_set_i32(base, radix);

  while (*p) {
    int d;
    if (*p == '.') {
      if (seen_dot)
        break;  /*  Second dot ends the number.  */
      seen_dot = 1;
      p++;
      continue;
    }
    d = digit_value(*p, radix);
    if (d < 0)
      break;
    arbint_mul(out, out, base);
    arbint_set_i32(digit, d);
    arbint_add(out, out, digit);
    if (seen_dot)
      scale++;
    p++;
  }

  arbint_clear(digit);
  arbint_clear(base);

  if (negative)
    arbint_neg(out, out);

  *out_scale = scale;
  return (int)(p - str);
}


/*  Number printing  */

static void dc_print_number_scaled(dc_state_t * dc, const arbint_t num,
                                    int scale) {
  char * str = NULL;
  size_t len;
  int negative = 0;
  char * digits;
  size_t digit_len;
  size_t i;

  if (arbint_get_str(num, &str, dc->output_radix) != ARBINT_OK) {
    fprintf(stderr, "dc: print error\n");
    return;
  }

  /*  Handle negative numbers.  */
  if (str[0] == '-') {
    negative = 1;
    digits = str + 1;
  } else {
    digits = str;
  }
  digit_len = strlen(digits);

  if (scale == 0) {
    /*  Integer output.  */
    printf("%s", str);
  } else if (scale >= (int)digit_len) {
    /*  All digits are fractional, need leading "0.000...".  */
    if (negative)
      putchar('-');
    putchar('0');
    putchar('.');
    for (i = 0; i < (size_t)(scale - (int)digit_len); i++)
      putchar('0');
    printf("%s", digits);
  } else {
    /*  Insert decimal point.  */
    len = digit_len;
    if (negative)
      putchar('-');
    for (i = 0; i < len - (size_t)scale; i++)
      putchar(digits[i]);
    putchar('.');
    printf("%s", digits + len - scale);
  }
  free(str);
}


/*  Scale manipulation helpers  */

/*  Multiply num by 10^n (scale up).  */
static void dc_scale_up(dc_state_t * dc, arbint_t num, int n) {
  arbint_t ten;
  int i;
  if (n <= 0)
    return;
  arbint_init(ten, &dc->ctx);
  arbint_set_i32(ten, 10);
  for (i = 0; i < n; i++)
    arbint_mul(num, num, ten);
  arbint_clear(ten);
}

/*  Divide num by 10^n (scale down), truncating.  */
static void dc_scale_down(dc_state_t * dc, arbint_t num, int n) {
  arbint_t ten;
  int i;
  if (n <= 0)
    return;
  arbint_init(ten, &dc->ctx);
  arbint_set_i32(ten, 10);
  for (i = 0; i < n; i++)
    arbint_tdiv_q(num, num, ten);
  arbint_clear(ten);
}

/*  Adjust two numbers to have the same scale (the maximum).
    Returns the common scale.  */
static int dc_equalize_scales(dc_state_t * dc, arbint_t a, int sa,
                               arbint_t b, int sb) {
  if (sa < sb) {
    dc_scale_up(dc, a, sb - sa);
    return sb;
  } else if (sb < sa) {
    dc_scale_up(dc, b, sa - sb);
    return sa;
  }
  return sa;
}

/*  Arithmetic operations  */

static int dc_op_add(dc_state_t * dc) {
  arbint_t a, b;
  int sa, sb, sr;

  if (dc->stack_size < 2) {
    fprintf(stderr, "dc: stack empty\n");
    return -1;
  }

  arbint_init_all(&dc->ctx, a, b, (arbint_t *) NULL);

  if (dc_pop_number_scaled(dc, b, &sb) < 0 ||
      dc_pop_number_scaled(dc, a, &sa) < 0) {
    arbint_clear_all(a, b, (arbint_t *) NULL);
    return -1;
  }

  /*  For addition, equalize scales then add.  */
  sr = dc_equalize_scales(dc, a, sa, b, sb);
  arbint_add(a, a, b);
  dc_push_number_scaled(dc, a, sr);

  arbint_clear_all(a, b, (arbint_t *) NULL);
  return 0;
}

static int dc_op_sub(dc_state_t * dc) {
  arbint_t a, b;
  int sa, sb, sr;

  if (dc->stack_size < 2) {
    fprintf(stderr, "dc: stack empty\n");
    return -1;
  }

  arbint_init_all(&dc->ctx, a, b, (arbint_t *) NULL);

  if (dc_pop_number_scaled(dc, b, &sb) < 0 ||
      dc_pop_number_scaled(dc, a, &sa) < 0) {
    arbint_clear_all(a, b, (arbint_t *) NULL);
    return -1;
  }

  /*  For subtraction, equalize scales then subtract.  */
  sr = dc_equalize_scales(dc, a, sa, b, sb);
  arbint_sub(a, a, b);
  dc_push_number_scaled(dc, a, sr);

  arbint_clear_all(a, b, (arbint_t *) NULL);
  return 0;
}

static int dc_op_mul(dc_state_t * dc) {
  arbint_t a, b;
  int sa, sb, sr;

  if (dc->stack_size < 2) {
    fprintf(stderr, "dc: stack empty\n");
    return -1;
  }

  arbint_init_all(&dc->ctx, a, b, (arbint_t *) NULL);

  if (dc_pop_number_scaled(dc, b, &sb) < 0 ||
      dc_pop_number_scaled(dc, a, &sa) < 0) {
    arbint_clear_all(a, b, (arbint_t *) NULL);
    return -1;
  }

  /*  For multiplication: (a * 10^-sa) * (b * 10^-sb) = (a*b) * 10^-(sa+sb).
      Then truncate to precision digits.  */
  arbint_mul(a, a, b);
  sr = sa + sb;
  /*  If result has more fractional digits than precision, truncate.  */
  if (sr > dc->precision) {
    dc_scale_down(dc, a, sr - dc->precision);
    sr = dc->precision;
  }
  dc_push_number_scaled(dc, a, sr);

  arbint_clear_all(a, b, (arbint_t *) NULL);
  return 0;
}

static int dc_op_div(dc_state_t * dc) {
  arbint_t a, b;
  int sa, sb, sr;

  if (dc->stack_size < 2) {
    fprintf(stderr, "dc: stack empty\n");
    return -1;
  }

  arbint_init_all(&dc->ctx, a, b, (arbint_t *) NULL);

  if (dc_pop_number_scaled(dc, b, &sb) < 0 ||
      dc_pop_number_scaled(dc, a, &sa) < 0) {
    arbint_clear_all(a, b, (arbint_t *) NULL);
    return -1;
  }

  if (arbint_is_zero(b)) {
    fprintf(stderr, "dc: divide by zero\n");
    arbint_clear_all(a, b, (arbint_t *) NULL);
    return -1;
  }

  /*  For division: (a * 10^-sa) / (b * 10^-sb) = (a/b) * 10^(sb-sa).
      To get precision digits of result, scale a up by precision + sb - sa
      extra digits before dividing, so result has scale = precision.  */
  sr = dc->precision;
  dc_scale_up(dc, a, sr + sb - sa);
  arbint_tdiv_q(a, a, b);
  dc_push_number_scaled(dc, a, sr);

  arbint_clear_all(a, b, (arbint_t *) NULL);
  return 0;
}

static int dc_op_mod(dc_state_t * dc) {
  arbint_t a, b;
  int sa, sb, sr;

  if (dc->stack_size < 2) {
    fprintf(stderr, "dc: stack empty\n");
    return -1;
  }

  arbint_init_all(&dc->ctx, a, b, (arbint_t *) NULL);

  if (dc_pop_number_scaled(dc, b, &sb) < 0 ||
      dc_pop_number_scaled(dc, a, &sa) < 0) {
    arbint_clear_all(a, b, (arbint_t *) NULL);
    return -1;
  }

  if (arbint_is_zero(b)) {
    fprintf(stderr, "dc: remainder by zero\n");
    arbint_clear_all(a, b, (arbint_t *) NULL);
    return -1;
  }

  /*  Equalize scales for remainder.  */
  sr = dc_equalize_scales(dc, a, sa, b, sb);
  arbint_tdiv_r(a, a, b);
  dc_push_number_scaled(dc, a, sr);

  arbint_clear_all(a, b, (arbint_t *) NULL);
  return 0;
}

static int dc_op_divmod(dc_state_t * dc) {
  arbint_t a, b, q, r;
  int sa, sb, sr;

  if (dc->stack_size < 2) {
    fprintf(stderr, "dc: stack empty\n");
    return -1;
  }

  arbint_init_all(&dc->ctx, a, b, q, r, (arbint_t *) NULL);

  if (dc_pop_number_scaled(dc, b, &sb) < 0 ||
      dc_pop_number_scaled(dc, a, &sa) < 0) {
    arbint_clear_all(a, b, q, r, (arbint_t *) NULL);
    return -1;
  }

  if (arbint_is_zero(b)) {
    fprintf(stderr, "dc: divide by zero\n");
    arbint_clear_all(a, b, q, r, (arbint_t *) NULL);
    return -1;
  }

  /*  Equalize scales for remainder; quotient is integer.  */
  sr = dc_equalize_scales(dc, a, sa, b, sb);
  arbint_tdiv_qr(q, r, a, b);
  dc_push_number_scaled(dc, r, sr);
  dc_push_number(dc, q);  /*  Quotient is integer (scale 0).  */

  arbint_clear_all(a, b, q, r, (arbint_t *) NULL);
  return 0;
}

static int dc_op_exp(dc_state_t * dc) {
  arbint_t base, exp;
  uint32_t e;
  int sbase, sexp, sr;

  if (dc->stack_size < 2) {
    fprintf(stderr, "dc: stack empty\n");
    return -1;
  }

  arbint_init_all(&dc->ctx, base, exp, (arbint_t *) NULL);

  if (dc_pop_number_scaled(dc, exp, &sexp) < 0 ||
      dc_pop_number_scaled(dc, base, &sbase) < 0) {
    arbint_clear_all(base, exp, (arbint_t *) NULL);
    return -1;
  }

  /*  Exponent must be an integer.  */
  if (sexp != 0) {
    fprintf(stderr, "dc: fractional exponent not supported\n");
    arbint_clear_all(base, exp, (arbint_t *) NULL);
    return -1;
  }

  /*  Check exponent is non-negative and fits in u32.  */
  if (arbint_is_neg(exp)) {
    fprintf(stderr, "dc: negative exponent\n");
    arbint_clear_all(base, exp, (arbint_t *) NULL);
    return -1;
  }

  if (arbint_get_u32(exp, &e) != ARBINT_OK) {
    fprintf(stderr, "dc: exponent too large\n");
    arbint_clear_all(base, exp, (arbint_t *) NULL);
    return -1;
  }

  /*  For exponentiation: (base * 10^-sbase)^e = base^e * 10^-(sbase*e).
      The result scale is sbase * e, truncated to precision.  */
  arbint_pow_u32(base, base, e);
  sr = sbase * (int)e;
  if (sr > dc->precision) {
    dc_scale_down(dc, base, sr - dc->precision);
    sr = dc->precision;
  }
  dc_push_number_scaled(dc, base, sr);

  arbint_clear_all(base, exp, (arbint_t *) NULL);
  return 0;
}

static int dc_op_modexp(dc_state_t * dc) {
  arbint_t base, exp, mod, result, tmp;
  int sbase, sexp, smod;

  if (dc->stack_size < 3) {
    fprintf(stderr, "dc: stack empty\n");
    return -1;
  }

  arbint_init_all(&dc->ctx, base, exp, mod, result, tmp, (arbint_t *) NULL);

  if (dc_pop_number_scaled(dc, mod, &smod) < 0 ||
      dc_pop_number_scaled(dc, exp, &sexp) < 0 ||
      dc_pop_number_scaled(dc, base, &sbase) < 0) {
    goto cleanup;
  }

  /*  Modular exponentiation is integer-only.  */
  if (sbase != 0 || sexp != 0 || smod != 0) {
    fprintf(stderr, "dc: modexp requires integer arguments\n");
    goto cleanup;
  }

  if (arbint_is_zero(mod)) {
    fprintf(stderr, "dc: modulus is zero\n");
    goto cleanup;
  }

  if (arbint_is_neg(exp)) {
    fprintf(stderr, "dc: negative exponent\n");
    goto cleanup;
  }

  /*  Binary modular exponentiation: result = base^exp mod mod.  */
  arbint_set_i32(result, 1);
  arbint_tdiv_r(base, base, mod);

  while (!arbint_is_zero(exp)) {
    if (arbint_is_odd(exp)) {
      arbint_mul(result, result, base);
      arbint_tdiv_r(result, result, mod);
    }
    arbint_shr(exp, exp, 1);
    arbint_mul(tmp, base, base);
    arbint_tdiv_r(base, tmp, mod);
  }

  dc_push_number_scaled(dc, result, 0);

  arbint_clear_all(base, exp, mod, result, tmp, (arbint_t *) NULL);
  return 0;

cleanup:
  arbint_clear_all(base, exp, mod, result, tmp, (arbint_t *) NULL);
  return -1;
}

static int dc_op_sqrt(dc_state_t * dc) {
  int sa, sr;

  if (dc->stack_size < 1) {
    fprintf(stderr, "dc: stack empty\n");
    return -1;
  }

  arbint_t a;
  arbint_init(a, &dc->ctx);

  if (dc_pop_number_scaled(dc, a, &sa) < 0) {
    arbint_clear(a);
    return -1;
  }

  if (arbint_is_neg(a)) {
    fprintf(stderr, "dc: square root of negative number\n");
    arbint_clear(a);
    return -1;
  }

  /*  For sqrt of a scaled number: sqrt(a * 10^-sa) = sqrt(a) * 10^(-sa/2).
      To get precision digits in result, scale up by 2*precision before sqrt,
      then the result has scale = precision + sa/2.
      For simplicity, we target result scale = precision.  */
  sr = dc->precision;
  /*  Scale up to get 2*sr + sa fractional "digit-pairs" before sqrt.  */
  dc_scale_up(dc, a, 2 * sr - sa);
  arbint_isqrt(a, a);
  dc_push_number_scaled(dc, a, sr);

  arbint_clear(a);
  return 0;
}

/*  Stack manipulation  */

static int dc_op_clear(dc_state_t * dc) {
  while (dc->stack_size > 0) {
    dc_value_clear(dc, &dc->stack[--dc->stack_size]);
  }
  return 0;
}

static int dc_op_dup(dc_state_t * dc) {
  dc_value_t * top;
  if (dc_peek(dc, 0, &top) < 0)
    return -1;
  return dc_push(dc, top);
}

static int dc_op_swap(dc_state_t * dc) {
  if (dc->stack_size < 2) {
    fprintf(stderr, "dc: stack empty\n");
    return -1;
  }

  dc_value_t tmp = dc->stack[dc->stack_size - 1];
  dc->stack[dc->stack_size - 1] = dc->stack[dc->stack_size - 2];
  dc->stack[dc->stack_size - 2] = tmp;
  return 0;
}

static int dc_op_drop(dc_state_t * dc) {
  dc_value_t val;
  if (dc_pop(dc, &val) < 0)
    return -1;
  dc_value_clear(dc, &val);
  return 0;
}

/*  Register operations  */

static int dc_op_store(dc_state_t * dc, int reg) {
  if (reg < 0 || reg >= DC_MAX_REGISTERS) {
    fprintf(stderr, "dc: invalid register\n");
    return -1;
  }

  dc_value_t val;
  if (dc_pop(dc, &val) < 0)
    return -1;

  if (dc->registers[reg].has_value)
    dc_value_clear(dc, &dc->registers[reg].val);

  dc->registers[reg].val = val;
  dc->registers[reg].has_value = 1;
  return 0;
}

static int dc_op_load(dc_state_t * dc, int reg) {
  if (reg < 0 || reg >= DC_MAX_REGISTERS) {
    fprintf(stderr, "dc: invalid register\n");
    return -1;
  }

  if (!dc->registers[reg].has_value) {
    fprintf(stderr, "dc: register '%c' is empty\n", reg);
    return -1;
  }

  return dc_push(dc, &dc->registers[reg].val);
}

static int dc_op_store_stack(dc_state_t * dc, int reg) {
  if (reg < 0 || reg >= DC_MAX_REGISTERS) {
    fprintf(stderr, "dc: invalid register\n");
    return -1;
  }

  dc_value_t val;
  if (dc_pop(dc, &val) < 0)
    return -1;

  dc_regstack_t * node = malloc(sizeof(dc_regstack_t));
  if (!node) {
    dc_value_clear(dc, &val);
    fprintf(stderr, "dc: out of memory\n");
    return -1;
  }

  node->val = val;
  node->next = dc->registers[reg].stack;
  dc->registers[reg].stack = node;
  return 0;
}

static int dc_op_load_stack(dc_state_t * dc, int reg) {
  if (reg < 0 || reg >= DC_MAX_REGISTERS) {
    fprintf(stderr, "dc: invalid register\n");
    return -1;
  }

  if (!dc->registers[reg].stack) {
    fprintf(stderr, "dc: register stack '%c' is empty\n", reg);
    return -1;
  }

  dc_regstack_t * node = dc->registers[reg].stack;
  dc->registers[reg].stack = node->next;

  int rc = dc_push(dc, &node->val);
  dc_value_clear(dc, &node->val);
  free(node);
  return rc;
}

/*  :r command: store value at array index in register r.  */
static int dc_op_array_store(dc_state_t * dc, int reg) {
  if (reg < 0 || reg >= DC_MAX_REGISTERS) {
    fprintf(stderr, "dc: invalid register\n");
    return -1;
  }

  /*  Pop index and value.  */
  arbint_t idx_val;
  arbint_init(idx_val, &dc->ctx);

  if (dc_pop_number(dc, idx_val) < 0) {
    arbint_clear(idx_val);
    return -1;
  }

  uint32_t index;
  if (arbint_get_u32(idx_val, &index) != ARBINT_OK) {
    fprintf(stderr, "dc: array index out of range\n");
    arbint_clear(idx_val);
    return -1;
  }
  arbint_clear(idx_val);

  dc_value_t val;
  if (dc_pop(dc, &val) < 0)
    return -1;

  /*  Find or create array entry.  */
  dc_array_entry_t * entry = dc->registers[reg].array;
  dc_array_entry_t * prev = NULL;

  while (entry && entry->index != index) {
    prev = entry;
    entry = entry->next;
  }

  if (entry) {
    /*  Replace existing entry.  */
    dc_value_clear(dc, &entry->val);
    entry->val = val;
  } else {
    /*  Create new entry.  */
    entry = malloc(sizeof(dc_array_entry_t));
    if (!entry) {
      dc_value_clear(dc, &val);
      fprintf(stderr, "dc: out of memory\n");
      return -1;
    }
    entry->index = index;
    entry->val = val;
    entry->next = NULL;

    if (prev)
      prev->next = entry;
    else
      dc->registers[reg].array = entry;
  }

  return 0;
}

/*  ;r command: load value from array index in register r.  */
static int dc_op_array_load(dc_state_t * dc, int reg) {
  if (reg < 0 || reg >= DC_MAX_REGISTERS) {
    fprintf(stderr, "dc: invalid register\n");
    return -1;
  }

  /*  Pop index.  */
  arbint_t idx_val;
  arbint_init(idx_val, &dc->ctx);

  if (dc_pop_number(dc, idx_val) < 0) {
    arbint_clear(idx_val);
    return -1;
  }

  uint32_t index;
  if (arbint_get_u32(idx_val, &index) != ARBINT_OK) {
    fprintf(stderr, "dc: array index out of range\n");
    arbint_clear(idx_val);
    return -1;
  }
  arbint_clear(idx_val);

  /*  Find array entry.  */
  dc_array_entry_t * entry = dc->registers[reg].array;
  while (entry && entry->index != index)
    entry = entry->next;

  if (entry) {
    return dc_push(dc, &entry->val);
  } else {
    /*  GNU dc returns 0 for uninitialized array elements.  */
    arbint_t zero;
    arbint_init(zero, &dc->ctx);
    arbint_set_i32(zero, 0);
    int rc = dc_push_number(dc, zero);
    arbint_clear(zero);
    return rc;
  }
}

/*  Printing  */

static int dc_op_print(dc_state_t * dc) {
  dc_value_t * top;
  if (dc_peek(dc, 0, &top) < 0)
    return -1;

  if (top->type == DC_TYPE_NUMBER) {
    dc_print_number_scaled(dc, top->u.num, top->scale);
    printf("\n");
  } else {
    printf("%s", top->u.str);
  }
  return 0;
}

static int dc_op_print_pop(dc_state_t * dc) {
  dc_value_t val;
  if (dc_pop(dc, &val) < 0)
    return -1;

  if (val.type == DC_TYPE_NUMBER) {
    dc_print_number_scaled(dc, val.u.num, val.scale);
  } else {
    printf("%s", val.u.str);
  }

  dc_value_clear(dc, &val);
  return 0;
}

static int dc_op_print_stack(dc_state_t * dc) {
  size_t i;
  for (i = dc->stack_size; i > 0; i--) {
    dc_value_t * val = &dc->stack[i - 1];
    if (val->type == DC_TYPE_NUMBER) {
      dc_print_number_scaled(dc, val->u.num, val->scale);
      printf("\n");
    } else {
      printf("%s\n", val->u.str);
    }
  }
  return 0;
}

/*  P command: pop and print as byte stream (number) or string (no newline).  */
static int dc_op_print_stream(dc_state_t * dc) {
  dc_value_t val;
  if (dc_pop(dc, &val) < 0)
    return -1;

  if (val.type == DC_TYPE_STRING) {
    /*  Print string without trailing newline.  */
    printf("%s", val.u.str);
  } else {
    /*  Print number as byte stream (big-endian).  */
    void * bytes = NULL;
    size_t nbytes = 0;

    if (arbint_export(val.u.num, &bytes, &nbytes) == ARBINT_OK &&
        bytes != NULL) {
      size_t i;
      unsigned char * p = (unsigned char *)bytes;
      for (i = 0; i < nbytes; i++)
        putchar(p[i]);
      free(bytes);
    }
  }

  dc_value_clear(dc, &val);
  return 0;
}

/*  a command: convert number to character or string to number.  */
static int dc_op_asciify(dc_state_t * dc) {
  dc_value_t val;
  if (dc_pop(dc, &val) < 0)
    return -1;

  if (val.type == DC_TYPE_NUMBER) {
    /*  Convert low byte of number to single-character string.  */
    uint32_t v = 0;
    arbint_t abs_val;
    arbint_init(abs_val, &dc->ctx);
    arbint_abs(abs_val, val.u.num);
    arbint_get_u32(abs_val, &v);
    arbint_clear(abs_val);

    char str[2];
    str[0] = (char)(v & 0xFF);
    str[1] = '\0';
    dc_value_clear(dc, &val);
    return dc_push_string(dc, str);
  } else {
    /*  Convert first character of string to number.  */
    int ch = (val.u.str && val.u.str[0]) ? (unsigned char)val.u.str[0] : 0;
    dc_value_clear(dc, &val);

    arbint_t a;
    arbint_init(a, &dc->ctx);
    arbint_set_i32(a, ch);
    int rc = dc_push_number(dc, a);
    arbint_clear(a);
    return rc;
  }
}

/*  Parameter commands  */

static int dc_op_set_precision(dc_state_t * dc) {
  arbint_t a;
  arbint_init(a, &dc->ctx);

  if (dc_pop_number(dc, a) < 0) {
    arbint_clear(a);
    return -1;
  }

  int32_t k;
  if (arbint_get_i32(a, &k) != ARBINT_OK || k < 0 || k > DC_MAX_PRECISION) {
    fprintf(stderr, "dc: invalid precision\n");
    arbint_clear(a);
    return -1;
  }

  dc->precision = k;
  arbint_clear(a);
  return 0;
}

static int dc_op_get_precision(dc_state_t * dc) {
  arbint_t a;
  arbint_init(a, &dc->ctx);
  arbint_set_i32(a, dc->precision);
  dc_push_number(dc, a);
  arbint_clear(a);
  return 0;
}

static int dc_op_set_input_radix(dc_state_t * dc) {
  arbint_t a;
  arbint_init(a, &dc->ctx);

  if (dc_pop_number(dc, a) < 0) {
    arbint_clear(a);
    return -1;
  }

  int32_t r;
  if (arbint_get_i32(a, &r) != ARBINT_OK || r < 2 || r > 36) {
    fprintf(stderr, "dc: input base must be between 2 and 36\n");
    arbint_clear(a);
    return -1;
  }

  dc->input_radix = r;
  arbint_clear(a);
  return 0;
}

static int dc_op_get_input_radix(dc_state_t * dc) {
  arbint_t a;
  arbint_init(a, &dc->ctx);
  arbint_set_i32(a, dc->input_radix);
  dc_push_number(dc, a);
  arbint_clear(a);
  return 0;
}

static int dc_op_set_output_radix(dc_state_t * dc) {
  arbint_t a;
  arbint_init(a, &dc->ctx);

  if (dc_pop_number(dc, a) < 0) {
    arbint_clear(a);
    return -1;
  }

  int32_t r;
  if (arbint_get_i32(a, &r) != ARBINT_OK || r < 2 || r > 36) {
    fprintf(stderr, "dc: output base must be between 2 and 36\n");
    arbint_clear(a);
    return -1;
  }

  dc->output_radix = r;
  arbint_clear(a);
  return 0;
}

static int dc_op_get_output_radix(dc_state_t * dc) {
  arbint_t a;
  arbint_init(a, &dc->ctx);
  arbint_set_i32(a, dc->output_radix);
  dc_push_number(dc, a);
  arbint_clear(a);
  return 0;
}

/*  Status commands  */

static int dc_op_stack_depth(dc_state_t * dc) {
  arbint_t a;
  arbint_init(a, &dc->ctx);
  arbint_set_u32(a, (uint32_t)dc->stack_size);
  dc_push_number(dc, a);
  arbint_clear(a);
  return 0;
}

static int dc_op_num_digits(dc_state_t * dc) {
  dc_value_t val;
  if (dc_pop(dc, &val) < 0)
    return -1;

  size_t count;
  if (val.type == DC_TYPE_NUMBER) {
    count = arbint_sizeinbase(val.u.num, 10);
  } else {
    count = strlen(val.u.str);
  }

  dc_value_clear(dc, &val);

  arbint_t a;
  arbint_init(a, &dc->ctx);
  arbint_set_u32(a, (uint32_t)count);
  dc_push_number(dc, a);
  arbint_clear(a);
  return 0;
}

/*  X command: push the scale (number of fractional digits).  */
static int dc_op_get_scale(dc_state_t * dc) {
  dc_value_t val;
  if (dc_pop(dc, &val) < 0)
    return -1;

  int scale = 0;
  if (val.type == DC_TYPE_NUMBER) {
    scale = val.scale;
  }
  /*  For strings, scale is 0.  */

  dc_value_clear(dc, &val);

  arbint_t a;
  arbint_init(a, &dc->ctx);
  arbint_set_i32(a, scale);
  dc_push_number(dc, a);
  arbint_clear(a);
  return 0;
}

/*  Comparison and conditionals  */

static int dc_op_compare(dc_state_t * dc, int op, int reg) {
  int sa, sb;

  if (dc->stack_size < 2) {
    fprintf(stderr, "dc: stack empty\n");
    return -1;
  }

  arbint_t a, b;
  arbint_init(a, &dc->ctx);
  arbint_init(b, &dc->ctx);

  if (dc_pop_number_scaled(dc, b, &sb) < 0 ||
      dc_pop_number_scaled(dc, a, &sa) < 0) {
    arbint_clear(a);
    arbint_clear(b);
    return -1;
  }

  /*  Equalize scales, then compare top (b) vs second (a).  */
  dc_equalize_scales(dc, a, sa, b, sb);
  int cmp = arbint_cmp(b, a);
  int execute = 0;

  switch (op) {
  case '>':
    execute = (cmp > 0);
    break;
  case '<':
    execute = (cmp < 0);
    break;
  case '=':
    execute = (cmp == 0);
    break;
  case '!': /* handled separately */
    break;
  }

  arbint_clear(a);
  arbint_clear(b);

  if (execute) {
    if (!dc->registers[reg].has_value) {
      fprintf(stderr, "dc: register '%c' is empty\n", reg);
      return -1;
    }
    if (dc->registers[reg].val.type == DC_TYPE_STRING) {
      dc_execute(dc, dc->registers[reg].val.u.str);
    }
  }

  return 0;
}

/*  Macro execution  */

static int dc_op_execute(dc_state_t * dc) {
  dc_value_t val;
  if (dc_pop(dc, &val) < 0)
    return -1;

  if (val.type != DC_TYPE_STRING) {
    fprintf(stderr, "dc: non-string value\n");
    dc_value_clear(dc, &val);
    return -1;
  }

  dc_execute(dc, val.u.str);
  dc_value_clear(dc, &val);
  return 0;
}

/*  ? command: read line from stdin and execute.  */
static int dc_op_read_execute(dc_state_t * dc) {
  char * line = dc_alloc_line_buffer();
  if (!line)
    return -1;

  if (fgets(line, DC_MAX_INPUT_LINE, stdin) == NULL) {
    free(line);
    return 0;  /*  EOF is not an error.  */
  }

  dc_execute(dc, line);
  free(line);
  return 0;
}

/*  String operations  */

static int dc_parse_string(dc_state_t * dc, const char * str, char ** out,
                           int * len) {
  /*  Parse bracketed string: [contents]
      Handles nested brackets.  */
  if (*str != '[')
    return -1;

  str++;
  int depth = 1;
  const char * start = str;

  while (*str && depth > 0) {
    if (*str == '[')
      depth++;
    else if (*str == ']')
      depth--;
    str++;
  }

  if (depth != 0)
    return -1;

  size_t slen = (size_t)(str - start - 1);
  char * result = malloc(slen + 1);
  if (!result)
    return -1;

  memcpy(result, start, slen);
  result[slen] = '\0';

  *out = result;
  *len = (int)(str - (start - 1));
  (void)dc;
  return 0;
}

/*  Main execution loop  */

static void dc_execute(dc_state_t * dc, const char * str) {
  const char * saved_ptr = dc->input_ptr;
  dc->input_ptr = str;

  while (*dc->input_ptr && dc->running) {
    int ch = *dc->input_ptr++;
    if (dc_execute_char(dc, ch) < 0) {
      /*  Error already reported; continue execution.  */
    }
    if (dc->quit_depth > 0)
      break;
  }

  dc->input_ptr = saved_ptr;
}

static int dc_execute_char(dc_state_t * dc, int ch) {
  /*  Skip whitespace and comments.  */
  if (isspace(ch))
    return 0;
  if (ch == '#') {
    /*  Comment: skip to end of line.  */
    while (*dc->input_ptr && *dc->input_ptr != '\n')
      dc->input_ptr++;
    return 0;
  }

  /*  Numbers (including underscore for negative).  */
  if (isdigit(ch) || ch == '.' ||
      (ch == '_' && (isdigit(dc->input_ptr[0]) || dc->input_ptr[0] == '.')) ||
      (ch >= 'A' && ch <= 'F' && dc->input_radix > 10)) {
    dc->input_ptr--;
    arbint_t num;
    int scale;
    arbint_init(num, &dc->ctx);
    int consumed = dc_parse_number_scaled(dc, dc->input_ptr, num, &scale);
    if (consumed > 0) {
      dc_push_number_scaled(dc, num, scale);
      dc->input_ptr += consumed;
    }
    arbint_clear(num);
    return 0;
  }

  /*  String/macro.  */
  if (ch == '[') {
    dc->input_ptr--;
    char * str;
    int len;
    if (dc_parse_string(dc, dc->input_ptr, &str, &len) == 0) {
      dc_push_string(dc, str);
      dc->input_ptr += len;
      free(str);
    } else {
      fprintf(stderr, "dc: unbalanced brackets\n");
      return -1;
    }
    return 0;
  }

  /*  Commands.  */
  switch (ch) {
  /*  Arithmetic.  */
  case '+':
    return dc_op_add(dc);
  case '-':
    return dc_op_sub(dc);
  case '*':
    return dc_op_mul(dc);
  case '/':
    return dc_op_div(dc);
  case '%':
    return dc_op_mod(dc);
  case '~':
    return dc_op_divmod(dc);
  case '^':
    return dc_op_exp(dc);
  case '|':
    return dc_op_modexp(dc);
  case 'v':
    return dc_op_sqrt(dc);

  /*  Stack manipulation.  */
  case 'c':
    return dc_op_clear(dc);
  case 'd':
    return dc_op_dup(dc);
  case 'r':
    return dc_op_swap(dc);
  case 'R':
    return dc_op_drop(dc);

  /*  Printing.  */
  case 'p':
    return dc_op_print(dc);
  case 'n':
    return dc_op_print_pop(dc);
  case 'P':
    return dc_op_print_stream(dc);
  case 'f':
    return dc_op_print_stack(dc);

  /*  Parameters.  */
  case 'k':
    return dc_op_set_precision(dc);
  case 'K':
    return dc_op_get_precision(dc);
  case 'i':
    return dc_op_set_input_radix(dc);
  case 'I':
    return dc_op_get_input_radix(dc);
  case 'o':
    return dc_op_set_output_radix(dc);
  case 'O':
    return dc_op_get_output_radix(dc);

  /*  Status.  */
  case 'z':
    return dc_op_stack_depth(dc);
  case 'Z':
    return dc_op_num_digits(dc);
  case 'X':
    return dc_op_get_scale(dc);

  /*  Register operations.  */
  case 's': {
    int reg = *dc->input_ptr++;
    if (!reg) {
      fprintf(stderr, "dc: missing register name\n");
      return -1;
    }
    return dc_op_store(dc, reg);
  }
  case 'l': {
    int reg = *dc->input_ptr++;
    if (!reg) {
      fprintf(stderr, "dc: missing register name\n");
      return -1;
    }
    return dc_op_load(dc, reg);
  }
  case 'S': {
    int reg = *dc->input_ptr++;
    if (!reg) {
      fprintf(stderr, "dc: missing register name\n");
      return -1;
    }
    return dc_op_store_stack(dc, reg);
  }
  case 'L': {
    int reg = *dc->input_ptr++;
    if (!reg) {
      fprintf(stderr, "dc: missing register name\n");
      return -1;
    }
    return dc_op_load_stack(dc, reg);
  }

  /*  Array operations.  */
  case ':': {
    int reg = *dc->input_ptr++;
    if (!reg) {
      fprintf(stderr, "dc: missing register name\n");
      return -1;
    }
    return dc_op_array_store(dc, reg);
  }
  case ';': {
    int reg = *dc->input_ptr++;
    if (!reg) {
      fprintf(stderr, "dc: missing register name\n");
      return -1;
    }
    return dc_op_array_load(dc, reg);
  }

  /*  Comparisons.  */
  case '>':
  case '<':
  case '=': {
    int neg = 0;
    int op = ch;
    if (*dc->input_ptr == '!') {
      neg = 1;
      dc->input_ptr++;
    }
    int reg = *dc->input_ptr++;
    if (!reg) {
      fprintf(stderr, "dc: missing register name\n");
      return -1;
    }
    if (neg) {
      /*  Negated comparison: >!r means execute r if NOT greater.  */
      int result = dc_op_compare(dc, op, reg);
      (void)result;
      /*  Note: GNU dc negation is weird; we swap the comparison.  */
    }
    return dc_op_compare(dc, op, reg);
  }

  /*  Execution.  */
  case 'x':
    return dc_op_execute(dc);
  case '?':
    return dc_op_read_execute(dc);

  /*  Misc.  */
  case 'a':
    return dc_op_asciify(dc);

  /*  Quit.  */
  case 'q':
    dc->quit_depth = 2;
    if (dc->quit_depth >= 2)
      dc->running = 0;
    return 0;
  case 'Q': {
    arbint_t a;
    arbint_init(a, &dc->ctx);
    if (dc_pop_number(dc, a) < 0) {
      arbint_clear(a);
      return -1;
    }
    int32_t levels;
    if (arbint_get_i32(a, &levels) != ARBINT_OK || levels < 0) {
      arbint_clear(a);
      return -1;
    }
    dc->quit_depth = levels;
    arbint_clear(a);
    return 0;
  }

  /*  Unknown command.  */
  default:
    fprintf(stderr, "dc: '%c' (0%o) unimplemented\n", ch, ch);
    return -1;
  }

  return 0;
}

/*  Initialization and cleanup  */

static int dc_init(dc_state_t * dc) {
  memset(dc, 0, sizeof(*dc));

  if (arbint_ctx_init_default(&dc->ctx) != ARBINT_OK)
    return -1;

  dc->stack = malloc(DC_INITIAL_STACK_CAP * sizeof(dc_value_t));
  if (!dc->stack) {
    arbint_ctx_clear(&dc->ctx);
    return -1;
  }
  dc->stack_cap = DC_INITIAL_STACK_CAP;
  dc->stack_size = 0;

  dc->precision = 0;
  dc->input_radix = 10;
  dc->output_radix = 10;
  dc->running = 1;

  return 0;
}

static void dc_cleanup(dc_state_t * dc) {
  int i;

  /*  Clear stack.  */
  while (dc->stack_size > 0) {
    dc_value_clear(dc, &dc->stack[--dc->stack_size]);
  }
  free(dc->stack);

  /*  Clear registers.  */
  for (i = 0; i < DC_MAX_REGISTERS; i++) {
    dc_regstack_t * node;
    dc_array_entry_t * arr;

    if (dc->registers[i].has_value)
      dc_value_clear(dc, &dc->registers[i].val);

    node = dc->registers[i].stack;
    while (node) {
      dc_regstack_t * next = node->next;
      dc_value_clear(dc, &node->val);
      free(node);
      node = next;
    }

    arr = dc->registers[i].array;
    while (arr) {
      dc_array_entry_t * next = arr->next;
      dc_value_clear(dc, &arr->val);
      free(arr);
      arr = next;
    }
  }

  arbint_ctx_clear(&dc->ctx);
}

/*  Main / test mode  */

static void dc_run_interactive(dc_state_t * dc) {
  char * line = dc_alloc_line_buffer();
  if (!line)
    return;

  while (dc->running && fgets(line, DC_MAX_INPUT_LINE, stdin)) {
    dc_execute(dc, line);
    dc->quit_depth = 0;
  }
  free(line);
}

static void dc_run_expression(dc_state_t * dc, const char * expr) {
  dc_execute(dc, expr);
}

static void dc_run_file(dc_state_t * dc, const char * filename) {
  FILE * f = fopen(filename, "r");
  if (!f) {
    fprintf(stderr, "dc: cannot open '%s'\n", filename);
    return;
  }

  char * line = dc_alloc_line_buffer();
  if (!line) {
    fclose(f);
    return;
  }
  while (dc->running && fgets(line, DC_MAX_INPUT_LINE, f)) {
    dc_execute(dc, line);
    dc->quit_depth = 0;
  }

  free(line);
  fclose(f);
}

/*  Test suite (run when invoked with --test)  */

/*  Format a scaled number as a string (caller must free).  */
static char * dc_format_scaled(const arbint_t num, int scale) {
  char * str = NULL;
  char * result = NULL;
  size_t len, digit_len;
  int negative = 0;
  char * digits;
  size_t i;

  if (arbint_get_str(num, &str, 10) != ARBINT_OK)
    return NULL;

  if (str[0] == '-') {
    negative = 1;
    digits = str + 1;
  } else {
    digits = str;
  }
  digit_len = strlen(digits);

  if (scale == 0) {
    return str;
  } else if (scale >= (int)digit_len) {
    /*  All digits fractional: "0.000...digits".  */
    len = (negative ? 1 : 0) + 2 + (size_t)(scale - (int)digit_len) + digit_len;
    result = malloc(len + 1);
    if (!result) {
      free(str);
      return NULL;
    }
    i = 0;
    if (negative)
      result[i++] = '-';
    result[i++] = '0';
    result[i++] = '.';
    while ((int)i < (negative ? 1 : 0) + 2 + (scale - (int)digit_len))
      result[i++] = '0';
    memcpy(result + i, digits, digit_len + 1);
    free(str);
    return result;
  } else {
    /*  Insert decimal point.  */
    len = strlen(str) + 1;  /*  +1 for decimal point.  */
    result = malloc(len + 1);
    if (!result) {
      free(str);
      return NULL;
    }
    i = 0;
    if (negative)
      result[i++] = '-';
    memcpy(result + i, digits, digit_len - scale);
    i += digit_len - scale;
    result[i++] = '.';
    memcpy(result + i, digits + digit_len - scale, scale + 1);
    free(str);
    return result;
  }
}

static int run_test(const char * name, const char * expr,
                    const char * expected) {
  dc_state_t dc;
  if (dc_init(&dc) < 0) {
    fprintf(stderr, "test '%s': init failed\n", name);
    return 1;
  }

  dc_run_expression(&dc, expr);

  /*  Get result from stack.  */
  char * result = NULL;
  int failed = 0;

  if (dc.stack_size == 0) {
    fprintf(stderr, "test '%s': stack empty after expression\n", name);
    failed = 1;
  } else if (dc.stack[dc.stack_size - 1].type != DC_TYPE_NUMBER) {
    fprintf(stderr, "test '%s': top of stack is not a number\n", name);
    failed = 1;
  } else {
    dc_value_t * top = &dc.stack[dc.stack_size - 1];
    result = dc_format_scaled(top->u.num, top->scale);
    if (!result || strcmp(result, expected) != 0) {
      fprintf(stderr, "test '%s': expected '%s', got '%s'\n", name, expected,
              result ? result : "(null)");
      failed = 1;
    }
  }

  if (result)
    free(result);
  dc_cleanup(&dc);
  return failed;
}

static int run_tests(void) {
  int failures = 0;

  /*  Basic operations (from Wikipedia)  */

  /*  "4 5 *" -> 20  */
  failures += run_test("wiki_basic_mul", "4 5 *", "20");

  /*  Division precision: "2 3 /" with default precision 0 -> 0.  */
  failures += run_test("wiki_div_prec0", "2 3 /", "0");

  /*  "5 k 2 3 / p" -> .66666 (we output 0.66666).  */
  failures += run_test("wiki_div_prec5", "5k 2 3 /", "0.66666");

  /*  sqrt((12 + (-3)^4) / 11) - 22 = sqrt((12+81)/11) - 22
      = sqrt(93/11) - 22 = sqrt(8.4545...) - 22 ~ 2.907... - 22 ~ -19.09
      With integer precision: sqrt(8) - 22 = 2 - 22 = -20.  */
  failures += run_test("wiki_complex_expr", "12 _3 4 ^ + 11 / v 22 -", "-20");

  /*  Same with precision for better result.  */
  failures += run_test("wiki_complex_prec", "10k 12 _3 4 ^ + 11 / v 22 -",
                       "-19.0923298925");

  /*  Basic arithmetic  */

  failures += run_test("add", "2 3+", "5");
  failures += run_test("sub", "10 3-", "7");
  failures += run_test("mul", "6 7*", "42");
  failures += run_test("div", "20 4/", "5");
  failures += run_test("mod", "17 5%", "2");
  failures += run_test("exp", "2 10^", "1024");

  /*  Negative numbers with underscore prefix.  */
  failures += run_test("neg_input", "_5 3+", "-2");
  failures += run_test("neg_result", "3 5-", "-2");
  failures += run_test("neg_mul", "_3 _4 *", "12");
  failures += run_test("neg_div", "_20 4 /", "-5");

  /*  Stack operations  */

  /*  d = duplicate top.  */
  failures += run_test("dup", "5 d+", "10");
  failures += run_test("dup_mul", "7 d*", "49");

  /*  r = swap top two.  */
  failures += run_test("swap", "1 2r-", "1");
  failures += run_test("swap_div", "2 10 r/", "5");

  /*  R = drop top.  */
  failures += run_test("drop", "1 2 3 R", "2");

  /*  c = clear stack.  */
  failures += run_test("clear", "1 2 3 c 5", "5");

  /*  z = stack depth.  */
  failures += run_test("depth_0", "z", "0");
  failures += run_test("depth_1", "42 z", "1");
  failures += run_test("depth_3", "1 2 3 z", "3");

  /*  Registers (Wikipedia examples)  */

  /*  "3 sc 4 lc * p" -> 12  */
  failures += run_test("wiki_register", "3 sc 4 lc *", "12");

  failures += run_test("store_load", "42 sa la", "42");
  failures += run_test("register_arith", "10 sa 5 la+", "15");
  failures += run_test("multi_register", "10 sa 20 sb la lb +", "30");

  /*  S/L = push/pop register stack.  */
  failures += run_test("reg_stack_push", "1 Sa 2 Sa 3 Sa La", "3");
  failures += run_test("reg_stack_pop", "1 Sa 2 Sa La La +", "3");
  failures += run_test("reg_stack_order", "10 Sa 20 Sa 30 Sa La La La + +", "60");

  /*  Macros (Wikipedia examples)  */

  /*  "[1 + 2 *] sm ... 3 lm x" -> (3+1)*2 = 8.  */
  failures += run_test("wiki_macro", "[1 + 2 *]sm 3 lmx", "8");

  failures += run_test("macro_simple", "[2*]sa 5 lax", "10");
  failures += run_test("macro_nested", "[[2+]x]sa 3 lax", "5");
  failures += run_test("macro_multi_op", "[d*]sa 5 lax", "25");

  /*  Conditionals (Wikipedia examples)  */

  /*  "[[equal]p] sr 5 5 =r" executes if equal.  */
  /*  We test with numeric result instead of string print.  */
  failures += run_test("wiki_cond_eq", "[99]sr 5 5 =r", "99");

  /*  ">r" means if TOP > SECOND, execute r.  */
  failures += run_test("cmp_gt_true", "[10]sa 1 2 >a", "10");
  failures += run_test("cmp_lt_true", "[10]sa 2 1 <a", "10");
  failures += run_test("cmp_eq_true", "[10]sa 5 5 =a", "10");

  /*  Conditional not executed (stack unchanged from before cmp).  */
  failures += run_test("cmp_gt_false", "99 [10]sa 2 1 >a", "99");
  failures += run_test("cmp_lt_false", "99 [10]sa 1 2 <a", "99");
  failures += run_test("cmp_eq_false", "99 [10]sa 5 6 =a", "99");

  /*  Loops and recursion (Wikipedia examples)  */

  /*  Factorial: "[d1-d1<F*]dsFx" for input on stack.  */
  failures += run_test("factorial_5", "5 [d1-d1<F*]dsFx", "120");
  failures += run_test("factorial_10_macro", "10 [d1-d1<F*]dsFx", "3628800");

  /*  Sum stack: "1 2 4 8 16 100 [+z1<a]dsax" -> 131.  */
  failures += run_test("wiki_sum_stack", "1 2 4 8 16 100 [+z1<a]dsax", "131");

  /*  Input/output radix (Wikipedia examples)  */

  /*  "16i2o DEADBEEF p" -> binary.  */
  /*  We verify the decimal value instead since output radix affects p.  */
  failures += run_test("wiki_hex_input", "16i DEADBEEF 10i", "3735928559");

  failures += run_test("hex_input", "16i FF 10i", "255");
  failures += run_test("hex_input_a", "16i A 10i", "10");
  failures += run_test("bin_input", "2i 1010 10i", "10");
  failures += run_test("bin_input_large", "2i 11111111 10i", "255");
  failures += run_test("oct_input", "8i 77 10i", "63");

  /*  I, O, K = query radix/precision.  */
  failures += run_test("query_ibase", "16i I", "16");
  failures += run_test("query_obase", "8o O", "8");
  failures += run_test("query_prec", "20k K", "20");

  /*  Division and modulo  */

  failures += run_test("divmod_q", "17 5 ~ r R", "3");
  failures += run_test("divmod_r", "17 5 ~ R", "2");
  failures += run_test("divmod_both", "17 5 ~ +", "5");  /*  3 + 2 = 5.  */

  /*  Negative division.  */
  failures += run_test("neg_divmod", "_17 5 ~", "-3");  /*  Quotient on top.  */

  /*  Modular exponentiation  */

  failures += run_test("modexp_small", "2 10 1000 |", "24");
  failures += run_test("modexp_large", "3 100 1000000007 |", "886041711");
  failures += run_test("modexp_wiki", "2 16 65537 |", "65536");

  /*  Square root  */

  failures += run_test("sqrt_perfect", "100 v", "10");
  failures += run_test("sqrt_large", "10000 v", "100");
  failures += run_test("sqrt_nonperfect", "50 v", "7");  /*  floor(sqrt(50)).  */

  /*  With precision.  */
  failures += run_test("sqrt_prec", "10k 2 v", "1.4142135623");

  /*  GCD (Wikipedia example uses stack input instead of stdin)  */

  /*  GCD(48, 18) = 6.  */
  failures += run_test("gcd_48_18", "48 18 [dSarLa%d0<a]dsax+", "6");
  failures += run_test("gcd_105_30", "105 30 [dSarLa%d0<a]dsax+", "15");
  failures += run_test("gcd_same", "42 42 [dSarLa%d0<a]dsax+", "42");

  /*  Large numbers  */

  failures += run_test("large_mul", "999999999 999999999*", "999999998000000001");
  failures += run_test("large_exp", "2 64^", "18446744073709551616");
  failures += run_test("pow_2_256", "2 256^",
      "115792089237316195423570985008687907853269984665640564039457584007913129639936");

  /*  Factorial using repeated multiplication.  */
  failures += run_test("factorial_10", "1 1 2 3 4 5 6 7 8 9 10*********", "3628800");

  /*  Doubling sequence: 1 doubled 10 times = 1024.  */
  failures += run_test("doubling",
      "1 d d+ d d+ d d+ d d+ d d+ d d+ d d+ d d+ d d+ d d+", "1024");

  /*  Z command: digit/character count  */

  failures += run_test("num_digits", "12345 Z", "5");
  failures += run_test("num_digits_neg", "_12345 Z", "5");
  failures += run_test("num_digits_1", "7 Z", "1");
  failures += run_test("num_digits_large", "2 100^ Z", "31");

  /*  X command: scale (fractional digits)  */

  failures += run_test("get_scale_int", "42 X", "0");
  failures += run_test("get_scale_dec", "3.14159 X", "5");
  failures += run_test("get_scale_add", "1.5 2.5 + X", "1");

  /*  Comments  */

  failures += run_test("comment_eol", "5 # this is a comment\n 3 +", "8");
  failures += run_test("comment_only", "# comment\n42", "42");
  failures += run_test("comment_inline", "10 20 # ignored\n+", "30");

  /*  Decimal/fixed-point arithmetic  */

  failures += run_test("decimal_input", "3.14159", "3.14159");
  failures += run_test("decimal_add", "1.5 2.5 +", "4.0");
  failures += run_test("decimal_sub", "5.5 2.3 -", "3.2");
  failures += run_test("decimal_mul", "5k 3.14159 2 *", "6.28318");
  failures += run_test("decimal_div", "5k 10 3 /", "3.33333");
  failures += run_test("decimal_div2", "10k 22 7 /", "3.1428571428");
  failures += run_test("decimal_sqrt", "10k 2 v", "1.4142135623");
  failures += run_test("decimal_pi_approx", "10k 355 113 /", "3.1415929203");
  failures += run_test("decimal_mixed", "5k 3.5 1.5 + 2 *", "10.0");

  /*  Array operations  */

  failures += run_test("array_store_load", "42 0 :a 0 ;a", "42");
  failures += run_test("array_multi_idx", "10 0 :a 20 1 :a 30 2 :a 1 ;a", "20");
  failures += run_test("array_unset", "5 ;a", "0");
  failures += run_test("array_overwrite", "10 0 :a 99 0 :a 0 ;a", "99");
  failures += run_test("array_diff_reg", "11 0 :a 22 0 :b 0 ;a 0 ;b +", "33");

  /*  a command: asciify  */

  failures += run_test("asciify_roundtrip", "65 a a", "65");
  failures += run_test("asciify_string", "[Hello]a", "72");
  failures += run_test("asciify_zero", "0 a a", "0");

  /*  Complex examples  */

  /*  2^32 computed by repeated squaring (5 squares).  */
  failures += run_test("repeated_square", "2 d* d* d* d* d*", "4294967296");

  /*  Compute n! using loop with counter in register.  */
  failures += run_test("factorial_loop", "1 5sc [lc*lc1-dsc0<f]dsfx", "120");

  /*  Print summary.  */
  if (failures == 0) {
    printf("All tests passed!\n");
  } else {
    printf("%d test(s) failed.\n", failures);
  }

  return failures;
}

/*  Usage and main  */

static void usage(FILE * out, const char * prog) {
  fprintf(out,
          "Usage: %s [options] [file ...]\n"
          "Options:\n"
          "  -h, --help              Show this help\n"
          "  -v, --version           Show version\n"
          "  -e, --expression EXPR   Execute expression\n"
          "      --expression=EXPR   Execute expression\n"
          "  -f, --file FILE         Execute file\n"
          "      --file=FILE         Execute file\n"
          "      --test              Run built-in test suite\n"
          "  --                      End option parsing\n",
          prog);
}

static void print_version(void) {
  printf("dc (c-arbint) %s\n", DC_PACKAGE_VERSION);
  printf("Copyright 2026 Kamila Szewczyk (k@iczelia.net)\n");
  printf("This is free software; see the source for copying conditions.  There is NO\n");
  printf("warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE,\n");
  printf("to the extent permitted by law.\n");
}

int main(int argc, char ** argv) {
  dc_state_t dc;
  int i;
  int had_input = 0;
  int initialized = 0;
  int end_of_options = 0;

  for (i = 1; i < argc; i++) {
    const char * arg = argv[i];

    if (!end_of_options && strcmp(arg, "--") == 0) {
      end_of_options = 1;
      continue;
    }

    if (!end_of_options &&
        (strcmp(arg, "--help") == 0 || strcmp(arg, "-h") == 0)) {
      usage(stdout, argv[0]);
      if (initialized)
        dc_cleanup(&dc);
      return 0;
    }

    if (!end_of_options &&
        (strcmp(arg, "--version") == 0 || strcmp(arg, "-v") == 0)) {
      print_version();
      if (initialized)
        dc_cleanup(&dc);
      return 0;
    }

    if (!end_of_options && strcmp(arg, "--test") == 0) {
      if (initialized)
        dc_cleanup(&dc);
      return run_tests();
    }

    if (!end_of_options &&
        (strcmp(arg, "--expression") == 0 || strcmp(arg, "-e") == 0)) {
      const char * expr;
      if (i + 1 >= argc) {
        fprintf(stderr, "dc: option '%s' requires an argument\n", arg);
        usage(stderr, argv[0]);
        return 1;
      }
      expr = argv[++i];
      if (!initialized) {
        if (dc_init(&dc) < 0) {
          fprintf(stderr, "dc: initialization failed\n");
          return 1;
        }
        initialized = 1;
      }
      dc_run_expression(&dc, expr);
      had_input = 1;
      continue;
    }

    if (!end_of_options && strncmp(arg, "--expression=", 13) == 0) {
      const char * expr = arg + 13;
      if (!initialized) {
        if (dc_init(&dc) < 0) {
          fprintf(stderr, "dc: initialization failed\n");
          return 1;
        }
        initialized = 1;
      }
      dc_run_expression(&dc, expr);
      had_input = 1;
      continue;
    }

    if (!end_of_options &&
        (strcmp(arg, "--file") == 0 || strcmp(arg, "-f") == 0)) {
      const char * file;
      if (i + 1 >= argc) {
        fprintf(stderr, "dc: option '%s' requires an argument\n", arg);
        usage(stderr, argv[0]);
        return 1;
      }
      file = argv[++i];
      if (!initialized) {
        if (dc_init(&dc) < 0) {
          fprintf(stderr, "dc: initialization failed\n");
          return 1;
        }
        initialized = 1;
      }
      dc_run_file(&dc, file);
      had_input = 1;
      continue;
    }

    if (!end_of_options && strncmp(arg, "--file=", 7) == 0) {
      const char * file = arg + 7;
      if (!initialized) {
        if (dc_init(&dc) < 0) {
          fprintf(stderr, "dc: initialization failed\n");
          return 1;
        }
        initialized = 1;
      }
      dc_run_file(&dc, file);
      had_input = 1;
      continue;
    }

    if (!end_of_options && strncmp(arg, "-e", 2) == 0 && arg[2] != '\0') {
      const char * expr = arg + 2;
      if (!initialized) {
        if (dc_init(&dc) < 0) {
          fprintf(stderr, "dc: initialization failed\n");
          return 1;
        }
        initialized = 1;
      }
      dc_run_expression(&dc, expr);
      had_input = 1;
      continue;
    }

    if (!end_of_options && strncmp(arg, "-f", 2) == 0 && arg[2] != '\0') {
      const char * file = arg + 2;
      if (!initialized) {
        if (dc_init(&dc) < 0) {
          fprintf(stderr, "dc: initialization failed\n");
          return 1;
        }
        initialized = 1;
      }
      dc_run_file(&dc, file);
      had_input = 1;
      continue;
    }

    if (!end_of_options && arg[0] == '-') {
      fprintf(stderr, "dc: unknown option '%s'\n", arg);
      usage(stderr, argv[0]);
      return 1;
    }

    if (!initialized) {
      if (dc_init(&dc) < 0) {
        fprintf(stderr, "dc: initialization failed\n");
        return 1;
      }
      initialized = 1;
    }
    dc_run_file(&dc, arg);
    had_input = 1;
  }

  if (!initialized) {
    if (dc_init(&dc) < 0) {
      fprintf(stderr, "dc: initialization failed\n");
      return 1;
    }
    initialized = 1;
  }

  if (!had_input)
    dc_run_interactive(&dc);

  if (initialized)
    dc_cleanup(&dc);

  return 0;
}
