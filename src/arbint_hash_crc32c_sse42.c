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

/*  CRC32C using SSE4.2 hardware intrinsics.  */
static uint64_t clmul(uint32_t a, uint32_t b) {
  uint64_t ret = 0, bp = b;
  for (uint32_t i = 0; i < 32; i++)
    ret ^= (a & (1 << i)) * bp;
  return ret;
}

uint32_t arbint_crc32c_sse42(uint32_t crc, const uint8_t * buf, size_t len) {
  for (; len && ((uintptr_t) buf & 7); --len)
    crc = _mm_crc32_u8(crc, *buf++);
  for (; len >= 4088; len -= 4088) {
    const uint8_t * end = buf + 1352;
    uint32_t acc1 = 0, acc2 = 0;
    do {
      crc = _mm_crc32_u64(crc, *(const uint64_t *) buf);
      acc1 = _mm_crc32_u64(acc1, *(const uint64_t *) (buf + 1360));
      acc2 = _mm_crc32_u64(acc2, *(const uint64_t *) (buf + 2720));
      buf += 8;
    } while (buf <= end);
    uint64_t vc = clmul(crc, 0x7b454cb3) ^ clmul(acc1, 0x79113270);
    buf += 2728;
    crc = _mm_crc32_u64(acc2, *(const uint64_t *) (buf - 8) ^ vc);
  }
  for (; len >= 8; buf += 8, len -= 8)
    crc = _mm_crc32_u64(crc, *(const uint64_t *) buf);
  for (; len; --len)
    crc = _mm_crc32_u8(crc, *buf++);
  return crc;
}
