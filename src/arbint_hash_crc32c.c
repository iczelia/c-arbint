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

#include <nmmintrin.h>
#include <string.h>

/*  CRC32C using SSE4.2 hardware intrinsics.

    Processes data in 8-byte (64-bit only), 4-byte, and single-byte
    chunks for maximum throughput.  Uses memcpy for safe unaligned
    loads (never casts void * to wider integer types).  */
uint32_t arbint_crc32c_sse42(uint32_t crc, const uint8_t * data, size_t len) {
  size_t i = 0u;

#if ARBINT_LIMB_BITS == 64
  /*  Process 8 bytes at a time with _mm_crc32_u64.  */
  while (i + 8u <= len) {
    uint64_t val;
    memcpy(&val, data + i, 8u);
    crc = (uint32_t) _mm_crc32_u64((uint64_t) crc, val);
    i += 8u;
  }
#endif /* ARBINT_LIMB_BITS == 64 */

  /*  Process 4 bytes at a time.  */
  while (i + 4u <= len) {
    uint32_t val;
    memcpy(&val, data + i, 4u);
    crc = _mm_crc32_u32(crc, val);
    i += 4u;
  }

  /*  Process remaining bytes.  */
  while (i < len) {
    crc = _mm_crc32_u8(crc, data[i]);
    ++i;
  }

  return crc;
}
