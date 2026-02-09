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

#include "arbint_hash.h"

#include "config.h"

#include <arm_neon.h>

/*  SHA-256 round constants.  */
static const uint32_t arbint_sha256_K[64] = {
    0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu,
    0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u, 0xd807aa98u, 0x12835b01u,
    0x243185beu, 0x550c7dc3u, 0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u,
    0xc19bf174u, 0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu,
    0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau, 0x983e5152u,
    0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u,
    0x06ca6351u, 0x14292967u, 0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu,
    0x53380d13u, 0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
    0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u, 0xd192e819u,
    0xd6990624u, 0xf40e3585u, 0x106aa070u, 0x19a4c116u, 0x1e376c08u,
    0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu,
    0x682e6ff3u, 0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
    0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u};

/*  Perform four SHA-256 rounds on ABEF/CDGH with message schedule w
    and round constants from K[i..i+3].

    ARM SHA-2 crypto instructions:
      vsha256hq_u32(ABEF, CDGH, WK)  - 4 rounds updating ABEF
      vsha256h2q_u32(CDGH, ABEF, WK) - 4 rounds updating CDGH
    Both must use the SAME original ABEF for the second call.  */
#define SHA256_4ROUNDS(abef, cdgh, w, i)                                      \
  do {                                                                         \
    uint32x4_t wk_ = vaddq_u32((w), vld1q_u32(&arbint_sha256_K[(i)]));        \
    uint32x4_t abef_prev_ = (abef);                                           \
    (abef) = vsha256hq_u32((abef), (cdgh), wk_);                              \
    (cdgh) = vsha256h2q_u32((cdgh), abef_prev_, wk_);                         \
  } while (0)

/*  SHA-256 single-block compression using ARM Crypto Extensions.

    The ARM SHA-2 instructions operate on two uint32x4_t registers:
      ABEF = [A, B, E, F]
      CDGH = [C, D, G, H]

    vsha256su0q_u32 and vsha256su1q_u32 perform message schedule
    expansion (sigma0 and sigma1 respectively).  */
void arbint_sha256_compress_arm(uint32_t state[8], const uint8_t block[64]) {
  uint32x4_t abef;
  uint32x4_t cdgh;
  uint32x4_t abef_save;
  uint32x4_t cdgh_save;
  uint32x4_t w0;
  uint32x4_t w1;
  uint32x4_t w2;
  uint32x4_t w3;

  /*  Load state.  Input layout: state[0..3] = A B C D, state[4..7] = E F G H.
      SHA-2 instructions expect ABEF and CDGH packing.  */
  {
    uint32x4_t s0 = vld1q_u32(state);
    uint32x4_t s1 = vld1q_u32(state + 4);

    /*  s0 = [A, B, C, D], s1 = [E, F, G, H]
        Need: abef = [A, B, E, F], cdgh = [C, D, G, H]  */
    abef = vcombine_u32(vget_low_u32(s0), vget_low_u32(s1));
    cdgh = vcombine_u32(vget_high_u32(s0), vget_high_u32(s1));
  }

  abef_save = abef;
  cdgh_save = cdgh;

  /*  Load message block as big-endian 32-bit words.  */
  w0 = vreinterpretq_u32_u8(vrev32q_u8(vld1q_u8(block)));
  w1 = vreinterpretq_u32_u8(vrev32q_u8(vld1q_u8(block + 16)));
  w2 = vreinterpretq_u32_u8(vrev32q_u8(vld1q_u8(block + 32)));
  w3 = vreinterpretq_u32_u8(vrev32q_u8(vld1q_u8(block + 48)));

  /*  Rounds 0-3  */
  SHA256_4ROUNDS(abef, cdgh, w0, 0);
  /*  Rounds 4-7  */
  SHA256_4ROUNDS(abef, cdgh, w1, 4);
  /*  Rounds 8-11  */
  SHA256_4ROUNDS(abef, cdgh, w2, 8);
  /*  Rounds 12-15  */
  SHA256_4ROUNDS(abef, cdgh, w3, 12);

  /*  Rounds 16-19  */
  w0 = vsha256su1q_u32(vsha256su0q_u32(w0, w1), w2, w3);
  SHA256_4ROUNDS(abef, cdgh, w0, 16);
  /*  Rounds 20-23  */
  w1 = vsha256su1q_u32(vsha256su0q_u32(w1, w2), w3, w0);
  SHA256_4ROUNDS(abef, cdgh, w1, 20);
  /*  Rounds 24-27  */
  w2 = vsha256su1q_u32(vsha256su0q_u32(w2, w3), w0, w1);
  SHA256_4ROUNDS(abef, cdgh, w2, 24);
  /*  Rounds 28-31  */
  w3 = vsha256su1q_u32(vsha256su0q_u32(w3, w0), w1, w2);
  SHA256_4ROUNDS(abef, cdgh, w3, 28);

  /*  Rounds 32-35  */
  w0 = vsha256su1q_u32(vsha256su0q_u32(w0, w1), w2, w3);
  SHA256_4ROUNDS(abef, cdgh, w0, 32);
  /*  Rounds 36-39  */
  w1 = vsha256su1q_u32(vsha256su0q_u32(w1, w2), w3, w0);
  SHA256_4ROUNDS(abef, cdgh, w1, 36);
  /*  Rounds 40-43  */
  w2 = vsha256su1q_u32(vsha256su0q_u32(w2, w3), w0, w1);
  SHA256_4ROUNDS(abef, cdgh, w2, 40);
  /*  Rounds 44-47  */
  w3 = vsha256su1q_u32(vsha256su0q_u32(w3, w0), w1, w2);
  SHA256_4ROUNDS(abef, cdgh, w3, 44);

  /*  Rounds 48-51  */
  w0 = vsha256su1q_u32(vsha256su0q_u32(w0, w1), w2, w3);
  SHA256_4ROUNDS(abef, cdgh, w0, 48);
  /*  Rounds 52-55  */
  w1 = vsha256su1q_u32(vsha256su0q_u32(w1, w2), w3, w0);
  SHA256_4ROUNDS(abef, cdgh, w1, 52);
  /*  Rounds 56-59  */
  w2 = vsha256su1q_u32(vsha256su0q_u32(w2, w3), w0, w1);
  SHA256_4ROUNDS(abef, cdgh, w2, 56);
  /*  Rounds 60-63  */
  w3 = vsha256su1q_u32(vsha256su0q_u32(w3, w0), w1, w2);
  SHA256_4ROUNDS(abef, cdgh, w3, 60);

  /*  Add saved state.  */
  abef = vaddq_u32(abef, abef_save);
  cdgh = vaddq_u32(cdgh, cdgh_save);

  /*  Unpack from ABEF/CDGH back to [A B C D] [E F G H].  */
  vst1q_u32(state, vcombine_u32(vget_low_u32(abef), vget_low_u32(cdgh)));
  vst1q_u32(state + 4, vcombine_u32(vget_high_u32(abef), vget_high_u32(cdgh)));
}
