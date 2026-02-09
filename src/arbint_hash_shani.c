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

#include <immintrin.h>

/*  SHA-256 round constants packed as __m128i pairs for SHA-NI rounds.
    Each pair holds K[4i..4i+3] in big-endian word order as the intrinsic
    expects.  */
static const uint32_t arbint_sha256_K_aligned[64] = {
  0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u,
  0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
  0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
  0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
  0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu,
  0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
  0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u,
  0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
  0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
  0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
  0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u,
  0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
  0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u,
  0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
  0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
  0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u
};

/*  SHA-256 single-block compression using Intel SHA-NI extensions.

    The SHA-NI instructions operate on two 128-bit registers that
    together hold the 8 working variables A..H:
      state0 = [A, B, E, F]  (after rearrangement)
      state1 = [C, D, G, H]

    _mm_sha256rnds2_epu32 processes two rounds at a time, using the
    low 64 bits (two 32-bit words) of a message+constant register.
    _mm_sha256msg1_epu32 and _mm_sha256msg2_epu32 perform partial
    and complete message schedule expansion.  */
void arbint_sha256_compress_shani(uint32_t state[8],
                                  const uint8_t block[64]) {
  __m128i state0, state1;
  __m128i state0_save, state1_save;
  __m128i msg0, msg1, msg2, msg3;
  __m128i tmp;
  __m128i mask;

  /*  Byte-swap shuffle mask: convert little-endian loads to big-endian
      32-bit words.  */
  mask = _mm_set_epi64x((long long) 0x0c0d0e0f08090a0bULL,
                         (long long) 0x0405060700010203ULL);

  /*  Load state into SHA-NI register layout.
      Input layout:  state[0..3] = A B C D,  state[4..7] = E F G H
      SHA-NI layout: state0 = A B E F,  state1 = C D G H  */
  tmp = _mm_loadu_si128((const __m128i *) state);
  state1 = _mm_loadu_si128((const __m128i *) (state + 4));

  tmp = _mm_shuffle_epi32(tmp, 0xB1);       /*  B A D C  */
  state1 = _mm_shuffle_epi32(state1, 0x1B);  /*  H G F E  */
  state0 = _mm_alignr_epi8(tmp, state1, 8);  /*  A B E F  */
  state1 = _mm_blend_epi16(state1, tmp, 0xF0); /*  C D G H  */

  state0_save = state0;
  state1_save = state1;

  /*  Rounds 0-3  */
  msg0 = _mm_loadu_si128((const __m128i *) block);
  msg0 = _mm_shuffle_epi8(msg0, mask);
  tmp = _mm_add_epi32(msg0, _mm_loadu_si128((const __m128i *) &arbint_sha256_K_aligned[0]));
  state1 = _mm_sha256rnds2_epu32(state1, state0, tmp);
  tmp = _mm_shuffle_epi32(tmp, 0x0E);
  state0 = _mm_sha256rnds2_epu32(state0, state1, tmp);

  /*  Rounds 4-7  */
  msg1 = _mm_loadu_si128((const __m128i *) (block + 16));
  msg1 = _mm_shuffle_epi8(msg1, mask);
  tmp = _mm_add_epi32(msg1, _mm_loadu_si128((const __m128i *) &arbint_sha256_K_aligned[4]));
  state1 = _mm_sha256rnds2_epu32(state1, state0, tmp);
  tmp = _mm_shuffle_epi32(tmp, 0x0E);
  state0 = _mm_sha256rnds2_epu32(state0, state1, tmp);
  msg0 = _mm_sha256msg1_epu32(msg0, msg1);

  /*  Rounds 8-11  */
  msg2 = _mm_loadu_si128((const __m128i *) (block + 32));
  msg2 = _mm_shuffle_epi8(msg2, mask);
  tmp = _mm_add_epi32(msg2, _mm_loadu_si128((const __m128i *) &arbint_sha256_K_aligned[8]));
  state1 = _mm_sha256rnds2_epu32(state1, state0, tmp);
  tmp = _mm_shuffle_epi32(tmp, 0x0E);
  state0 = _mm_sha256rnds2_epu32(state0, state1, tmp);
  msg1 = _mm_sha256msg1_epu32(msg1, msg2);

  /*  Rounds 12-15  */
  msg3 = _mm_loadu_si128((const __m128i *) (block + 48));
  msg3 = _mm_shuffle_epi8(msg3, mask);
  tmp = _mm_add_epi32(msg3, _mm_loadu_si128((const __m128i *) &arbint_sha256_K_aligned[12]));
  state1 = _mm_sha256rnds2_epu32(state1, state0, tmp);
  tmp = _mm_shuffle_epi32(tmp, 0x0E);
  state0 = _mm_sha256rnds2_epu32(state0, state1, tmp);
  msg0 = _mm_sha256msg2_epu32(_mm_add_epi32(msg0, _mm_alignr_epi8(msg3, msg2, 4)), msg3);
  msg2 = _mm_sha256msg1_epu32(msg2, msg3);

  /*  Rounds 16-19  */
  tmp = _mm_add_epi32(msg0, _mm_loadu_si128((const __m128i *) &arbint_sha256_K_aligned[16]));
  state1 = _mm_sha256rnds2_epu32(state1, state0, tmp);
  tmp = _mm_shuffle_epi32(tmp, 0x0E);
  state0 = _mm_sha256rnds2_epu32(state0, state1, tmp);
  msg1 = _mm_sha256msg2_epu32(_mm_add_epi32(msg1, _mm_alignr_epi8(msg0, msg3, 4)), msg0);
  msg3 = _mm_sha256msg1_epu32(msg3, msg0);

  /*  Rounds 20-23  */
  tmp = _mm_add_epi32(msg1, _mm_loadu_si128((const __m128i *) &arbint_sha256_K_aligned[20]));
  state1 = _mm_sha256rnds2_epu32(state1, state0, tmp);
  tmp = _mm_shuffle_epi32(tmp, 0x0E);
  state0 = _mm_sha256rnds2_epu32(state0, state1, tmp);
  msg2 = _mm_sha256msg2_epu32(_mm_add_epi32(msg2, _mm_alignr_epi8(msg1, msg0, 4)), msg1);
  msg0 = _mm_sha256msg1_epu32(msg0, msg1);

  /*  Rounds 24-27  */
  tmp = _mm_add_epi32(msg2, _mm_loadu_si128((const __m128i *) &arbint_sha256_K_aligned[24]));
  state1 = _mm_sha256rnds2_epu32(state1, state0, tmp);
  tmp = _mm_shuffle_epi32(tmp, 0x0E);
  state0 = _mm_sha256rnds2_epu32(state0, state1, tmp);
  msg3 = _mm_sha256msg2_epu32(_mm_add_epi32(msg3, _mm_alignr_epi8(msg2, msg1, 4)), msg2);
  msg1 = _mm_sha256msg1_epu32(msg1, msg2);

  /*  Rounds 28-31  */
  tmp = _mm_add_epi32(msg3, _mm_loadu_si128((const __m128i *) &arbint_sha256_K_aligned[28]));
  state1 = _mm_sha256rnds2_epu32(state1, state0, tmp);
  tmp = _mm_shuffle_epi32(tmp, 0x0E);
  state0 = _mm_sha256rnds2_epu32(state0, state1, tmp);
  msg0 = _mm_sha256msg2_epu32(_mm_add_epi32(msg0, _mm_alignr_epi8(msg3, msg2, 4)), msg3);
  msg2 = _mm_sha256msg1_epu32(msg2, msg3);

  /*  Rounds 32-35  */
  tmp = _mm_add_epi32(msg0, _mm_loadu_si128((const __m128i *) &arbint_sha256_K_aligned[32]));
  state1 = _mm_sha256rnds2_epu32(state1, state0, tmp);
  tmp = _mm_shuffle_epi32(tmp, 0x0E);
  state0 = _mm_sha256rnds2_epu32(state0, state1, tmp);
  msg1 = _mm_sha256msg2_epu32(_mm_add_epi32(msg1, _mm_alignr_epi8(msg0, msg3, 4)), msg0);
  msg3 = _mm_sha256msg1_epu32(msg3, msg0);

  /*  Rounds 36-39  */
  tmp = _mm_add_epi32(msg1, _mm_loadu_si128((const __m128i *) &arbint_sha256_K_aligned[36]));
  state1 = _mm_sha256rnds2_epu32(state1, state0, tmp);
  tmp = _mm_shuffle_epi32(tmp, 0x0E);
  state0 = _mm_sha256rnds2_epu32(state0, state1, tmp);
  msg2 = _mm_sha256msg2_epu32(_mm_add_epi32(msg2, _mm_alignr_epi8(msg1, msg0, 4)), msg1);
  msg0 = _mm_sha256msg1_epu32(msg0, msg1);

  /*  Rounds 40-43  */
  tmp = _mm_add_epi32(msg2, _mm_loadu_si128((const __m128i *) &arbint_sha256_K_aligned[40]));
  state1 = _mm_sha256rnds2_epu32(state1, state0, tmp);
  tmp = _mm_shuffle_epi32(tmp, 0x0E);
  state0 = _mm_sha256rnds2_epu32(state0, state1, tmp);
  msg3 = _mm_sha256msg2_epu32(_mm_add_epi32(msg3, _mm_alignr_epi8(msg2, msg1, 4)), msg2);
  msg1 = _mm_sha256msg1_epu32(msg1, msg2);

  /*  Rounds 44-47  */
  tmp = _mm_add_epi32(msg3, _mm_loadu_si128((const __m128i *) &arbint_sha256_K_aligned[44]));
  state1 = _mm_sha256rnds2_epu32(state1, state0, tmp);
  tmp = _mm_shuffle_epi32(tmp, 0x0E);
  state0 = _mm_sha256rnds2_epu32(state0, state1, tmp);
  msg0 = _mm_sha256msg2_epu32(_mm_add_epi32(msg0, _mm_alignr_epi8(msg3, msg2, 4)), msg3);
  msg2 = _mm_sha256msg1_epu32(msg2, msg3);

  /*  Rounds 48-51  */
  tmp = _mm_add_epi32(msg0, _mm_loadu_si128((const __m128i *) &arbint_sha256_K_aligned[48]));
  state1 = _mm_sha256rnds2_epu32(state1, state0, tmp);
  tmp = _mm_shuffle_epi32(tmp, 0x0E);
  state0 = _mm_sha256rnds2_epu32(state0, state1, tmp);
  msg1 = _mm_sha256msg2_epu32(_mm_add_epi32(msg1, _mm_alignr_epi8(msg0, msg3, 4)), msg0);
  msg3 = _mm_sha256msg1_epu32(msg3, msg0);

  /*  Rounds 52-55  */
  tmp = _mm_add_epi32(msg1, _mm_loadu_si128((const __m128i *) &arbint_sha256_K_aligned[52]));
  state1 = _mm_sha256rnds2_epu32(state1, state0, tmp);
  tmp = _mm_shuffle_epi32(tmp, 0x0E);
  state0 = _mm_sha256rnds2_epu32(state0, state1, tmp);
  msg2 = _mm_sha256msg2_epu32(_mm_add_epi32(msg2, _mm_alignr_epi8(msg1, msg0, 4)), msg1);

  /*  Rounds 56-59  */
  tmp = _mm_add_epi32(msg2, _mm_loadu_si128((const __m128i *) &arbint_sha256_K_aligned[56]));
  state1 = _mm_sha256rnds2_epu32(state1, state0, tmp);
  tmp = _mm_shuffle_epi32(tmp, 0x0E);
  state0 = _mm_sha256rnds2_epu32(state0, state1, tmp);
  msg3 = _mm_sha256msg2_epu32(_mm_add_epi32(msg3, _mm_alignr_epi8(msg2, msg1, 4)), msg2);

  /*  Rounds 60-63  */
  tmp = _mm_add_epi32(msg3, _mm_loadu_si128((const __m128i *) &arbint_sha256_K_aligned[60]));
  state1 = _mm_sha256rnds2_epu32(state1, state0, tmp);
  tmp = _mm_shuffle_epi32(tmp, 0x0E);
  state0 = _mm_sha256rnds2_epu32(state0, state1, tmp);

  /*  Add saved state.  */
  state0 = _mm_add_epi32(state0, state0_save);
  state1 = _mm_add_epi32(state1, state1_save);

  /*  Unpack from SHA-NI layout back to linear [A B C D] [E F G H].  */
  tmp = _mm_shuffle_epi32(state0, 0x1B);      /*  F E B A  */
  state1 = _mm_shuffle_epi32(state1, 0xB1);    /*  D C H G  */
  state0 = _mm_blend_epi16(tmp, state1, 0xF0); /*  D C B A  */
  state1 = _mm_alignr_epi8(state1, tmp, 8);    /*  H G F E  */

  _mm_storeu_si128((__m128i *) state, state0);
  _mm_storeu_si128((__m128i *) (state + 4), state1);
}
