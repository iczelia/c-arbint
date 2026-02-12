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

#ifndef ARBINT_STR_H
#define ARBINT_STR_H

#include "arbint_base.h"

#include <stddef.h>
#include <stdint.h>

/*  ========== Algorithm Thresholds ==========  */

/*  Limb count threshold for switching from blocking to divide-and-conquer.
    D&C has higher constant factors but better asymptotic complexity.  */
#define ARBINT_STR_DC_THRESHOLD 32u

/*  Digit count threshold for set_str D&C (approx 32 limbs * 19 digits).  */
#define ARBINT_STR_DC_DIGIT_THRESHOLD 600u

/*  ========== Character Lookup Tables ==========  */

/*  Map ASCII character to digit value 0-35.
    Returns 0xFF for invalid characters.
    Accepts both uppercase and lowercase letters.  */
static const uint8_t arbint_char_to_digit[256] = {
    /*  0x00-0x0F  */
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    /*  0x10-0x1F  */
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    /*  0x20-0x2F: space and punctuation  */
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    /*  0x30-0x3F: '0'-'9' = 0-9, then punctuation  */
    0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u,
    8u, 9u, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    /*  0x40-0x4F: '@', 'A'-'O' = 10-24  */
    0xFFu, 10u, 11u, 12u, 13u, 14u, 15u, 16u,
    17u, 18u, 19u, 20u, 21u, 22u, 23u, 24u,
    /*  0x50-0x5F: 'P'-'Z' = 25-35, then '[\\]^_'  */
    25u, 26u, 27u, 28u, 29u, 30u, 31u, 32u,
    33u, 34u, 35u, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    /*  0x60-0x6F: '`', 'a'-'o' = 10-24  */
    0xFFu, 10u, 11u, 12u, 13u, 14u, 15u, 16u,
    17u, 18u, 19u, 20u, 21u, 22u, 23u, 24u,
    /*  0x70-0x7F: 'p'-'z' = 25-35, then '{|}~' and DEL  */
    25u, 26u, 27u, 28u, 29u, 30u, 31u, 32u,
    33u, 34u, 35u, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    /*  0x80-0xFF: high bytes, all invalid  */
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu,
    0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu, 0xFFu
};

/*  Map digit value 0-35 to lowercase ASCII character.  */
static const char arbint_digit_to_char[36] = {
    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
    'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
    'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
    'u', 'v', 'w', 'x', 'y', 'z'
};

/*  ========== Base Blocking Parameters ==========  */

/*  For each base 2-36, precompute:
    - digits: max digits extractable per single-limb division
    - power: base^digits (fits in a limb)

    Key insight: dividing a bigint by base^digits extracts `digits` output
    digits at once, reducing the number of expensive divisions by a factor
    of `digits`.  */

typedef struct {
  uint8_t digits;        /*  Max digits per limb division  */
  arbint_limb_t power;   /*  base^digits (fits in limb)  */
} arbint_base_block_t;

#if ARBINT_LIMB_BITS == 64

/*  64-bit limb blocking table.  Index = base (0-1 unused).
    For base b: digits = max k such that b^k < 2^64, power = b^digits.  */
static const arbint_base_block_t arbint_base_block[37] = {
    {0u, 0u},                                      /*  base 0: unused  */
    {0u, 0u},                                      /*  base 1: unused  */
    {63u, UINT64_C(9223372036854775808)},          /*  base 2: 2^63  */
    {40u, UINT64_C(12157665459056928801)},         /*  base 3: 3^40  */
    {31u, UINT64_C(4611686018427387904)},          /*  base 4: 4^31  */
    {27u, UINT64_C(7450580596923828125)},          /*  base 5: 5^27  */
    {24u, UINT64_C(4738381338321616896)},          /*  base 6: 6^24  */
    {22u, UINT64_C(3909821048582988049)},          /*  base 7: 7^22  */
    {21u, UINT64_C(9223372036854775808)},          /*  base 8: 8^21  */
    {20u, UINT64_C(12157665459056928801)},         /*  base 9: 9^20  */
    {19u, UINT64_C(10000000000000000000)},         /*  base 10: 10^19  */
    {18u, UINT64_C(5559917313492231481)},          /*  base 11: 11^18  */
    {17u, UINT64_C(2218611106740436992)},          /*  base 12: 12^17  */
    {17u, UINT64_C(8650415919381337933)},          /*  base 13: 13^17  */
    {16u, UINT64_C(2177953337809371136)},          /*  base 14: 14^16  */
    {16u, UINT64_C(6568408355712890625)},          /*  base 15: 15^16  */
    {15u, UINT64_C(1152921504606846976)},          /*  base 16: 16^15  */
    {15u, UINT64_C(2862423051509815793)},          /*  base 17: 17^15  */
    {15u, UINT64_C(6746640616477458432)},          /*  base 18: 18^15  */
    {15u, UINT64_C(15181127029874798299)},         /*  base 19: 19^15  */
    {14u, UINT64_C(1638400000000000000)},          /*  base 20: 20^14  */
    {14u, UINT64_C(3243919932521508681)},          /*  base 21: 21^14  */
    {14u, UINT64_C(6221821273427820544)},          /*  base 22: 22^14  */
    {14u, UINT64_C(11592836324538749809)},         /*  base 23: 23^14  */
    {13u, UINT64_C(876488338465357824)},           /*  base 24: 24^13  */
    {13u, UINT64_C(1490116119384765625)},          /*  base 25: 25^13  */
    {13u, UINT64_C(2481152873203736576)},          /*  base 26: 26^13  */
    {13u, UINT64_C(4052555153018976267)},          /*  base 27: 27^13  */
    {13u, UINT64_C(6502111422497947648)},          /*  base 28: 28^13  */
    {13u, UINT64_C(10260628712958602189)},         /*  base 29: 29^13  */
    {13u, UINT64_C(15943230000000000000)},         /*  base 30: 30^13  */
    {12u, UINT64_C(787662783788549761)},           /*  base 31: 31^12  */
    {12u, UINT64_C(1152921504606846976)},          /*  base 32: 32^12  */
    {12u, UINT64_C(1667889514952984961)},          /*  base 33: 33^12  */
    {12u, UINT64_C(2386420683693101056)},          /*  base 34: 34^12  */
    {12u, UINT64_C(3379220508056640625)},          /*  base 35: 35^12  */
    {12u, UINT64_C(4738381338321616896)}           /*  base 36: 36^12  */
};

#elif ARBINT_LIMB_BITS == 32

/*  32-bit limb blocking table.  */
static const arbint_base_block_t arbint_base_block[37] = {
    {0u, 0u},                          /*  base 0: unused  */
    {0u, 0u},                          /*  base 1: unused  */
    {31u, UINT32_C(0x80000000)},       /*  base 2: 2^31  */
    {20u, UINT32_C(3486784401)},       /*  base 3: 3^20  */
    {15u, UINT32_C(0x40000000)},       /*  base 4: 4^15  */
    {13u, UINT32_C(1220703125)},       /*  base 5: 5^13  */
    {12u, UINT32_C(2176782336)},       /*  base 6: 6^12  */
    {11u, UINT32_C(1977326743)},       /*  base 7: 7^11  */
    {10u, UINT32_C(0x40000000)},       /*  base 8: 8^10  */
    {10u, UINT32_C(3486784401)},       /*  base 9: 9^10  */
    {9u, UINT32_C(1000000000)},        /*  base 10: 10^9  */
    {9u, UINT32_C(2357947691)},        /*  base 11: 11^9  */
    {8u, UINT32_C(429981696)},         /*  base 12: 12^8  */
    {8u, UINT32_C(815730721)},         /*  base 13: 13^8  */
    {8u, UINT32_C(1475789056)},        /*  base 14: 14^8  */
    {8u, UINT32_C(2562890625)},        /*  base 15: 15^8  */
    {7u, UINT32_C(0x10000000)},        /*  base 16: 16^7  */
    {7u, UINT32_C(410338673)},         /*  base 17: 17^7  */
    {7u, UINT32_C(612220032)},         /*  base 18: 18^7  */
    {7u, UINT32_C(893871739)},         /*  base 19: 19^7  */
    {7u, UINT32_C(1280000000)},        /*  base 20: 20^7  */
    {7u, UINT32_C(1801088541)},        /*  base 21: 21^7  */
    {7u, UINT32_C(2494357888)},        /*  base 22: 22^7  */
    {7u, UINT32_C(3404825447)},        /*  base 23: 23^7  */
    {6u, UINT32_C(191102976)},         /*  base 24: 24^6  */
    {6u, UINT32_C(244140625)},         /*  base 25: 25^6  */
    {6u, UINT32_C(308915776)},         /*  base 26: 26^6  */
    {6u, UINT32_C(387420489)},         /*  base 27: 27^6  */
    {6u, UINT32_C(481890304)},         /*  base 28: 28^6  */
    {6u, UINT32_C(594823321)},         /*  base 29: 29^6  */
    {6u, UINT32_C(729000000)},         /*  base 30: 30^6  */
    {6u, UINT32_C(887503681)},         /*  base 31: 31^6  */
    {6u, UINT32_C(0x40000000)},        /*  base 32: 32^6  */
    {6u, UINT32_C(1291467969)},        /*  base 33: 33^6  */
    {6u, UINT32_C(1544804416)},        /*  base 34: 34^6  */
    {6u, UINT32_C(1838265625)},        /*  base 35: 35^6  */
    {6u, UINT32_C(2176782336)}         /*  base 36: 36^6  */
};

#endif /* ARBINT_LIMB_BITS */

/*  ========== Helper Functions ==========  */

/*  Return log2(base) if base is a power of 2, else 0.
    Valid for base in [2, 32].  */
static inline unsigned arbint_base_log2_pow2(int base) {
  unsigned log2;

  if (base < 2 || base > 32)
    return 0u;
  if ((base & (base - 1)) != 0)
    return 0u;

  log2 = 0u;
  while ((1 << log2) < base)
    ++log2;
  return log2;
}

/*  Allocate a string buffer using the arbint's context allocator.
    Falls back to default allocator if context is NULL.  */
static inline char * arbint_str_alloc(const arbint_t x, size_t n) {
  if (x != NULL && x[0]._ctx != NULL && x[0]._ctx->a.realloc != NULL)
    return (char *) x[0]._ctx->a.realloc(x[0]._ctx->a.ud, NULL, n);
  return (char *) arbint_alloc(NULL, NULL, n);
}

/*  Free a string buffer using the arbint's context allocator.  */
static inline void arbint_str_free(const arbint_t x, char * s) {
  if (s == NULL)
    return;
  if (x != NULL && x[0]._ctx != NULL && x[0]._ctx->a.realloc != NULL)
    x[0]._ctx->a.realloc(x[0]._ctx->a.ud, s, 0u);
  else
    arbint_alloc(NULL, s, 0u);
}

#endif /* ARBINT_STR_H */
