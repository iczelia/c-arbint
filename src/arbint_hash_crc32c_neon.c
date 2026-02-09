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

#include <arm_acle.h>
#include <arm_neon.h>
#include <stddef.h>
#include <stdint.h>

/*  Carryless multiply low half with fused EOR (PMULL + EOR).
    Computes: result = fold_lo XOR data, where fold_lo = clmul_low(a, b).  */
static inline __attribute__((always_inline)) uint64x2_t
clmul_lo_e(uint64x2_t a, uint64x2_t b, uint64x2_t c) {
  uint64x2_t result;
  __asm__("pmull %0.1q, %2.1d, %3.1d \n\t"
          "eor   %0.16b, %0.16b, %1.16b \n\t"
          : "=&w"(result)
          : "w"(c), "w"(a), "w"(b));
  return result;
}

/*  Carryless multiply high half with fused EOR (PMULL2 + EOR).
    Computes: result = fold_hi XOR data, where fold_hi = clmul_high(a, b).  */
static inline __attribute__((always_inline)) uint64x2_t
clmul_hi_e(uint64x2_t a, uint64x2_t b, uint64x2_t c) {
  uint64x2_t result;
  __asm__("pmull2 %0.1q, %2.2d, %3.2d \n\t"
          "eor    %0.16b, %0.16b, %1.16b \n\t"
          : "=&w"(result)
          : "w"(c), "w"(a), "w"(b));
  return result;
}

/*  CRC32C using ARM CRC32 instructions and NEON PMULL for parallel
    stream folding.  Processes 64 bytes per iteration in the hot loop
    with 4 parallel CRC streams reduced via carryless multiplication.  */
uint32_t arbint_crc32c_neon(uint32_t crc, const uint8_t * buf, size_t len) {
  for (; len && ((uintptr_t) buf & 7); --len)
    crc = __crc32cb(crc, *buf++);
  if (len >= 64) {
    int64_t len_s = (int64_t) len;
    uint32_t acc1 = 0, acc2 = 0, acc3 = 0;
    static const uint64_t __attribute__((aligned(16))) kc[] = {
        0x740eef02, 0x9e4addf8,  /*  fold-by-4  */
        0xf20c0dfe, 0x493c7d27,  /*  fold-to-128 + fold-by-1  */
        0x3da6d0cb, 0xba4fc28e}; /*  Barrett reduction  */
    do {
      crc = __crc32cd(crc, *(const uint64_t *) (buf + 0));
      acc1 = __crc32cd(acc1, *(const uint64_t *) (buf + 8));
      acc2 = __crc32cd(acc2, *(const uint64_t *) (buf + 16));
      acc3 = __crc32cd(acc3, *(const uint64_t *) (buf + 24));
      crc = __crc32cd(crc, *(const uint64_t *) (buf + 32));
      acc1 = __crc32cd(acc1, *(const uint64_t *) (buf + 40));
      acc2 = __crc32cd(acc2, *(const uint64_t *) (buf + 48));
      acc3 = __crc32cd(acc3, *(const uint64_t *) (buf + 56));
      buf += 64;
      len_s -= 64;
    } while (len_s >= 128);
    uint64x2_t vc0 = {crc, 0}, vc1 = {acc1, 0};
    uint64x2_t vc2 = {acc2, 0}, vc3 = {acc3, 0};
    uint64x2_t vk = vld1q_u64(kc);
    vc0 = clmul_lo_e(vc0, vk, vreinterpretq_u64_u8(vld1q_u8(buf + 0)));
    vc0 = clmul_hi_e(vc0, vk, vc1);
    vc2 = clmul_lo_e(vc2, vk, vreinterpretq_u64_u8(vld1q_u8(buf + 16)));
    vc2 = clmul_hi_e(vc2, vk, vc3);
    vk = vld1q_u64(kc + 2);
    vc0 = clmul_lo_e(vc0, vk, vc2);
    vc0 = clmul_hi_e(vc0, vk, vreinterpretq_u64_u8(vld1q_u8(buf + 32)));
    uint64x2_t vbr = vld1q_u64(kc + 4);
    uint64x2_t vt;
    vt = vreinterpretq_u64_p128(
        vmull_p64(vgetq_lane_u64(vc0, 0), vgetq_lane_u64(vbr, 0)));
    vt = vreinterpretq_u64_p128(
        vmull_p64(vgetq_lane_u64(vt, 0), vgetq_lane_u64(vbr, 1)));
    vc0 = veorq_u64(vc0, vt);
    crc = vgetq_lane_u32(vreinterpretq_u32_u64(vc0), 1);
    crc = __crc32cd(crc, *(const uint64_t *) (buf + 40));
    crc = __crc32cd(crc, *(const uint64_t *) (buf + 48));
    crc = __crc32cd(crc, *(const uint64_t *) (buf + 56));
    buf += 64;
    len_s -= 64;
    len = (size_t) len_s;
  }
  for (; len >= 16; len -= 16) {
    crc = __crc32cd(crc, *(const uint64_t *) buf);
    crc = __crc32cd(crc, *(const uint64_t *) (buf + 8));
    buf += 16;
  }
  if (len >= 8) {
    crc = __crc32cd(crc, *(const uint64_t *) buf);
    buf += 8;
    len -= 8;
  }
  for (; len; --len)
    crc = __crc32cb(crc, *buf++);
  return crc;
}
