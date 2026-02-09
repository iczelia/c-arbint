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

#ifndef ARBINT_HASH_H
#define ARBINT_HASH_H

#include "arbint_base.h"

/*  SHA-256 intermediate state.  */
typedef struct arbint_sha256_state {
  uint32_t h[8];      /*  Chain values.  */
  uint8_t buf[64];    /*  Partial block buffer.  */
  size_t buf_len;     /*  Bytes currently in buf (0..63).  */
  uint64_t total_len; /*  Total bytes fed so far.  */
} arbint_sha256_state_t;

/*  SHA-256 single-block compression (generic portable).
    state: in/out 8 x uint32_t chain values.
    block: pointer to exactly 64 bytes of message data.  */
void arbint_sha256_compress_generic(uint32_t state[8],
                                    const uint8_t block[64]);

#if HAS_SHA_NI
/*  SHA-NI accelerated single-block compression.  */
void arbint_sha256_compress_shani(uint32_t state[8], const uint8_t block[64]);
#endif /* HAS_SHA_NI */

#if HAS_ARM_SHA2
/*  ARM Crypto Extensions accelerated single-block compression.  */
void arbint_sha256_compress_arm(uint32_t state[8], const uint8_t block[64]);
#endif /* HAS_ARM_SHA2 */

/*  CRC32C over a byte buffer (generic portable, table-based).  */
uint32_t arbint_crc32c_generic(uint32_t crc, const uint8_t * data, size_t len);

#if HAS_SSE42_CRC32
/*  CRC32C using SSE4.2 _mm_crc32_* intrinsics.  */
uint32_t arbint_crc32c_sse42(uint32_t crc, const uint8_t * data, size_t len);
#endif /* HAS_SSE42_CRC32 */

#if HAS_ARM_CRC32 && HAS_ARM_PMULL
/*  CRC32C using ARM CRC32 instructions and NEON PMULL folding.  */
uint32_t arbint_crc32c_neon(uint32_t crc, const uint8_t * data, size_t len);
#endif /* HAS_ARM_CRC32 && HAS_ARM_PMULL */

#endif /* ARBINT_HASH_H */
