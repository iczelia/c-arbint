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

ARBINT_TEST_DECLARE_FAILURES();

static void build_large_operand(arbint_t x, uint32_t shift, uint32_t tail) {
  CHECK_EQ_I(arbint_set_u32(x, 1u), ARBINT_OK);
  CHECK_EQ_I(arbint_shl(x, x, shift), ARBINT_OK);
  CHECK_EQ_I(arbint_add_u32(x, x, tail), ARBINT_OK);
}

int main(void) {
  ARBINT_TEST_START();
  arbint_ctx_t ctx;
  arbint_t a, b, c;
  int out = 0;

  arbint_drop_caches();
  arbint_drop_caches();

  CHECK_EQ_I(arbint_ctx_init_default(&ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(a, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(b, &ctx), ARBINT_OK);
  CHECK_EQ_I(arbint_init(c, &ctx), ARBINT_OK);

  /*  Prime-sieve cache build, drop, and rebuild.  */
  CHECK_EQ_I(arbint_set_u32(a, 2u), ARBINT_OK);
  CHECK_EQ_I(arbint_pow_u32(a, a, 541u), ARBINT_OK);
  CHECK_EQ_I(arbint_is_power(a, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  arbint_drop_caches();
  CHECK_EQ_I(arbint_is_power(a, &out), ARBINT_OK);
  CHECK_EQ_I(out, 1);

  arbint_drop_caches();
  arbint_drop_caches();

  /*  NTT cache build and drop via explicit API.  */
  build_large_operand(a, 70000u, 12345u);
  build_large_operand(b, 69936u, 67890u);
  CHECK_EQ_I(arbint_mul(c, a, b), ARBINT_OK);
  CHECK(!arbint_is_zero(c));

  arbint_drop_caches();

  CHECK_EQ_I(arbint_mul(c, a, b), ARBINT_OK);
  CHECK(!arbint_is_zero(c));

  /*  Context clear path also triggers cache drop and remains reusable.  */
  arbint_ctx_clear(&ctx);
  CHECK_EQ_I(arbint_mul(c, a, b), ARBINT_OK);
  CHECK(!arbint_is_zero(c));

  arbint_clear(c);
  arbint_clear(b);
  arbint_clear(a);
  arbint_ctx_clear(&ctx);

  ARBINT_TEST_FINISH("test_cache");
}
