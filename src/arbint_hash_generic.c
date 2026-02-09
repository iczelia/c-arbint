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

/* ========== SHA-256 (FIPS 180-4) Generic Compression ========== */

/*  First 32 bits of the fractional parts of the cube roots of the
    first 64 prime numbers.  */
static const uint32_t arbint_sha256_K[64] = {
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

static inline uint32_t arbint_sha256_rotr(uint32_t x, unsigned n) {
  return (x >> n) | (x << (32u - n));
}

#define ARBINT_SHA256_CH(x, y, z)  (((x) & (y)) ^ (~(x) & (z)))
#define ARBINT_SHA256_MAJ(x, y, z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))

#define ARBINT_SHA256_BSIG0(x) \
  (arbint_sha256_rotr((x), 2) ^ arbint_sha256_rotr((x), 13) ^ arbint_sha256_rotr((x), 22))
#define ARBINT_SHA256_BSIG1(x) \
  (arbint_sha256_rotr((x), 6) ^ arbint_sha256_rotr((x), 11) ^ arbint_sha256_rotr((x), 25))
#define ARBINT_SHA256_SSIG0(x) \
  (arbint_sha256_rotr((x), 7) ^ arbint_sha256_rotr((x), 18) ^ ((x) >> 3))
#define ARBINT_SHA256_SSIG1(x) \
  (arbint_sha256_rotr((x), 17) ^ arbint_sha256_rotr((x), 19) ^ ((x) >> 10))

/*  Load a 32-bit big-endian word from a byte pointer.  */
static inline uint32_t arbint_sha256_load_be32(const uint8_t * p) {
  return ((uint32_t) p[0] << 24u)
       | ((uint32_t) p[1] << 16u)
       | ((uint32_t) p[2] << 8u)
       | ((uint32_t) p[3]);
}

void arbint_sha256_compress_generic(uint32_t state[8],
                                    const uint8_t block[64]) {
  uint32_t W[64];
  uint32_t a, b, c, d, e, f, g, h;
  unsigned i;

  /*  Message schedule: load 16 words big-endian, expand to 64.  */
  for (i = 0u; i < 16u; ++i)
    W[i] = arbint_sha256_load_be32(block + i * 4u);
  for (i = 16u; i < 64u; ++i)
    W[i] = ARBINT_SHA256_SSIG1(W[i - 2u]) + W[i - 7u]
         + ARBINT_SHA256_SSIG0(W[i - 15u]) + W[i - 16u];

  /*  Initialize working variables from current chain value.  */
  a = state[0]; b = state[1]; c = state[2]; d = state[3];
  e = state[4]; f = state[5]; g = state[6]; h = state[7];

  /*  64 rounds.  */
  for (i = 0u; i < 64u; ++i) {
    uint32_t T1 = h + ARBINT_SHA256_BSIG1(e) + ARBINT_SHA256_CH(e, f, g)
                + arbint_sha256_K[i] + W[i];
    uint32_t T2 = ARBINT_SHA256_BSIG0(a) + ARBINT_SHA256_MAJ(a, b, c);
    h = g; g = f; f = e; e = d + T1;
    d = c; c = b; b = a; a = T1 + T2;
  }

  /*  Accumulate into state.  */
  state[0] += a; state[1] += b; state[2] += c; state[3] += d;
  state[4] += e; state[5] += f; state[6] += g; state[7] += h;
}

/* ========== CRC32C (Castagnoli Polynomial) Generic ========== */

/*  Precomputed CRC32C lookup table.
    Polynomial: 0x82F63B78 (bit-reversed form of 0x1EDC6F41).  */
static const uint32_t arbint_crc32c_table[256] = {
  0x00000000u, 0xf26b8303u, 0xe13b70f7u, 0x1350f3f4u,
  0xc79a971fu, 0x35f1141cu, 0x26a1e7e8u, 0xd4ca64ebu,
  0x8ad958cfu, 0x78b2dbccu, 0x6be22838u, 0x9989ab3bu,
  0x4d43cfd0u, 0xbf284cd3u, 0xac78bf27u, 0x5e133c24u,
  0x105ec76fu, 0xe235446cu, 0xf165b798u, 0x030e349bu,
  0xd7c45070u, 0x25afd373u, 0x36ff2087u, 0xc494a384u,
  0x9a879fa0u, 0x68ec1ca3u, 0x7bbcef57u, 0x89d76c54u,
  0x5d1d08bfu, 0xaf768bbcu, 0xbc267848u, 0x4e4dfb4bu,
  0x20bd8edeu, 0xd2d60dddu, 0xc186fe29u, 0x33ed7d2au,
  0xe72719c1u, 0x154c9ac2u, 0x061c6936u, 0xf477ea35u,
  0xaa64d611u, 0x580f5512u, 0x4b5fa6e6u, 0xb93425e5u,
  0x6dfe410eu, 0x9f95c20du, 0x8cc531f9u, 0x7eaeb2fau,
  0x30e349b1u, 0xc288cab2u, 0xd1d83946u, 0x23b3ba45u,
  0xf779deaeu, 0x05125dadu, 0x1642ae59u, 0xe4292d5au,
  0xba3a117eu, 0x4851927du, 0x5b016189u, 0xa96ae28au,
  0x7da08661u, 0x8fcb0562u, 0x9c9bf696u, 0x6ef07595u,
  0x417b1dbcu, 0xb3109ebfu, 0xa0406d4bu, 0x522bee48u,
  0x86e18aa3u, 0x748a09a0u, 0x67dafa54u, 0x95b17957u,
  0xcba24573u, 0x39c9c670u, 0x2a993584u, 0xd8f2b687u,
  0x0c38d26cu, 0xfe53516fu, 0xed03a29bu, 0x1f682198u,
  0x5125dad3u, 0xa34e59d0u, 0xb01eaa24u, 0x42752927u,
  0x96bf4dccu, 0x64d4cecfu, 0x77843d3bu, 0x85efbe38u,
  0xdbfc821cu, 0x2997011fu, 0x3ac7f2ebu, 0xc8ac71e8u,
  0x1c661503u, 0xee0d9600u, 0xfd5d65f4u, 0x0f36e6f7u,
  0x61c69362u, 0x93ad1061u, 0x80fde395u, 0x72966096u,
  0xa65c047du, 0x5437877eu, 0x4767748au, 0xb50cf789u,
  0xeb1fcbadu, 0x197448aeu, 0x0a24bb5au, 0xf84f3859u,
  0x2c855cb2u, 0xdeeedfb1u, 0xcdbe2c45u, 0x3fd5af46u,
  0x7198540du, 0x83f3d70eu, 0x90a324fau, 0x62c8a7f9u,
  0xb602c312u, 0x44694011u, 0x5739b3e5u, 0xa55230e6u,
  0xfb410cc2u, 0x092a8fc1u, 0x1a7a7c35u, 0xe811ff36u,
  0x3cdb9bddu, 0xceb018deu, 0xdde0eb2au, 0x2f8b6829u,
  0x82f63b78u, 0x709db87bu, 0x63cd4b8fu, 0x91a6c88cu,
  0x456cac67u, 0xb7072f64u, 0xa457dc90u, 0x563c5f93u,
  0x082f63b7u, 0xfa44e0b4u, 0xe9141340u, 0x1b7f9043u,
  0xcfb5f4a8u, 0x3dde77abu, 0x2e8e845fu, 0xdc45075cu,
  0x9208fc17u, 0x60637f14u, 0x73338ce0u, 0x81580fe3u,
  0x55926b08u, 0xa7f9e80bu, 0xb4a91bffu, 0x469298fcu,
  0x1881a4d8u, 0xeaea27dbu, 0xf9bad42fu, 0x0bd1572cu,
  0xdf1b33c7u, 0x2d70b0c4u, 0x3e204330u, 0xcc4bc033u,
  0xc30c8ea1u, 0x316745a2u, 0x2237b656u, 0xd05c3555u,
  0x049651beu, 0xf6fdd2bdu, 0xe5ad2149u, 0x17c6a24au,
  0x49d59e6eu, 0xbbbed16du, 0xa8ee6299u, 0x5a85e19au,
  0x8e4f8571u, 0x7c240672u, 0x6f74f586u, 0x9d1f7685u,
  0xd3528dceu, 0x21390ecdu, 0x3269fd39u, 0xc0027e3au,
  0x14c81ad1u, 0xe6a399d2u, 0xf5f36a26u, 0x0798e925u,
  0x598bd501u, 0xab605602u, 0xb830a5f6u, 0x4a5b26f5u,
  0x9e91421eu, 0x6cfac11du, 0x7faa32e9u, 0x8dc1b1eau,
  0xe3317b7fu, 0x115af87cu, 0x020a0b88u, 0xf061888bu,
  0x24abec60u, 0xd6c06f63u, 0xc5909c97u, 0x37fb1f94u,
  0x69e823b0u, 0x9b83a0b3u, 0x88d35347u, 0x7ab8d044u,
  0xae72b4afu, 0x5c1937acu, 0x4f49c458u, 0xbd22475bu,
  0xf36fbc10u, 0x01043f13u, 0x1254cce7u, 0xe03f4fe4u,
  0x34f52b0fu, 0xc69ea80cu, 0xd5ce5bf8u, 0x27a5d8fbu,
  0x79b6e4dfu, 0x8bdd67dcu, 0x988d9428u, 0x6ae6172bu,
  0xbe2c73c0u, 0x4c47f0c3u, 0x5f170337u, 0xad7c8034u,
  0xa2392284u, 0x5052a187u, 0x43025273u, 0xb169d170u,
  0x65a3b59bu, 0x97c83698u, 0x8498c56cu, 0x76f3466fu,
  0x28e07a4bu, 0xda8bf948u, 0xc9db0abcu, 0x3bb089bfu,
  0xef7aed54u, 0x1d116e57u, 0x0e419da3u, 0xfc2a1ea0u,
  0xb267e5ebu, 0x400c66e8u, 0x535c951cu, 0xa137161fu,
  0x75fd72f4u, 0x8796f1f7u, 0x94c60203u, 0x66ad8100u,
  0x38bead24u, 0xcad52e27u, 0xd985ddd3u, 0x2bee5ed0u,
  0xff243a3bu, 0x0d4fb938u, 0x1e1f4accu, 0xec74c9cfu
};

uint32_t arbint_crc32c_generic(uint32_t crc, const uint8_t * data,
                                size_t len) {
  size_t i;
  for (i = 0u; i < len; ++i)
    crc = arbint_crc32c_table[(crc ^ data[i]) & 0xffu] ^ (crc >> 8u);
  return crc;
}
