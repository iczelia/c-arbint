/* arbint - portable arbitrary-precision computation library
 *
 * Copyright (C) 2026 Kamila Szewczyk (k@iczelia.net)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include "arbint_cpu.h"
#include "arbint_base.h"

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#if ARBINT_TARGET_X86_FAMILY && ARBINT_COMPILER_GNU_CLANG
  #define ARBINT_CPU_CAN_QUERY_X86_CPUID 1
#else
  #define ARBINT_CPU_CAN_QUERY_X86_CPUID 0
#endif

#define ARBINT_CPUID1_ECX_SSE3 (1u << 0)
#define ARBINT_CPUID1_ECX_PCLMULQDQ (1u << 1)
#define ARBINT_CPUID1_ECX_SSSE3 (1u << 9)
#define ARBINT_CPUID1_ECX_FMA (1u << 12)
#define ARBINT_CPUID1_ECX_SSE41 (1u << 19)
#define ARBINT_CPUID1_ECX_SSE42 (1u << 20)
#define ARBINT_CPUID1_ECX_POPCNT (1u << 23)
#define ARBINT_CPUID1_ECX_XSAVE (1u << 26)
#define ARBINT_CPUID1_ECX_OSXSAVE (1u << 27)
#define ARBINT_CPUID1_ECX_AVX (1u << 28)

#define ARBINT_CPUID1_EDX_SSE (1u << 25)
#define ARBINT_CPUID1_EDX_SSE2 (1u << 26)

#define ARBINT_CPUID7_EBX_BMI1 (1u << 3)
#define ARBINT_CPUID7_EBX_AVX2 (1u << 5)
#define ARBINT_CPUID7_EBX_BMI2 (1u << 8)
#define ARBINT_CPUID7_EBX_AVX512F (1u << 16)
#define ARBINT_CPUID7_EBX_SHA (1u << 29)
#define ARBINT_CPUID7_EBX_AVX512BW (1u << 30)
#define ARBINT_CPUID7_EBX_AVX512VL (1u << 31)

static void arbint_cpu_cpuid_count(unsigned int leaf, unsigned int subleaf,
                                   unsigned int * out_eax,
                                   unsigned int * out_ebx,
                                   unsigned int * out_ecx,
                                   unsigned int * out_edx) {
  unsigned int eax = 0u;
  unsigned int ebx = 0u;
  unsigned int ecx = 0u;
  unsigned int edx = 0u;

#if ARBINT_CPU_CAN_QUERY_X86_CPUID
  #if defined(__i386__) && defined(__PIC__)
  __asm__ volatile("xchgl %%ebx, %1\n\t"
                   "cpuid\n\t"
                   "xchgl %%ebx, %1\n\t"
                   : "=a"(eax), "=&r"(ebx), "=c"(ecx), "=d"(edx)
                   : "0"(leaf), "2"(subleaf)
                   : "cc");
  #else
  __asm__ volatile("cpuid"
                   : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx)
                   : "a"(leaf), "c"(subleaf)
                   : "cc");
  #endif
#else
  (void) leaf;
  (void) subleaf;
#endif

  if (out_eax != NULL)
    *out_eax = eax;
  if (out_ebx != NULL)
    *out_ebx = ebx;
  if (out_ecx != NULL)
    *out_ecx = ecx;
  if (out_edx != NULL)
    *out_edx = edx;
}

static unsigned int arbint_cpu_cpuid_max_leaf(void) {
  unsigned int eax = 0u;
  arbint_cpu_cpuid_count(0u, 0u, &eax, NULL, NULL, NULL);
  return eax;
}

typedef struct {
  unsigned int initialized : 1;
  unsigned int has_cpuid : 1;
  unsigned int has_xgetbv : 1;
  unsigned int sse : 1;
  unsigned int sse2 : 1;
  unsigned int sse3 : 1;
  unsigned int ssse3 : 1;
  unsigned int sse41 : 1;
  unsigned int sse42 : 1;
  unsigned int popcnt : 1;
  unsigned int crc32 : 1;
  unsigned int pclmulqdq : 1;
  unsigned int bmi1 : 1;
  unsigned int bmi2 : 1;
  unsigned int avx_hw : 1;
  unsigned int avx_os : 1;
  unsigned int avx : 1;
  unsigned int avx2 : 1;
  unsigned int avx512f : 1;
  unsigned int avx512bw : 1;
  unsigned int avx512vl : 1;
  unsigned int fma : 1;
  unsigned int sha : 1;
  uint64_t xcr0;
} arbint_cpu_caps_t;

static int arbint_cpu_can_call_cpuid(void) {
#if !ARBINT_TARGET_X86_FAMILY || !ARBINT_COMPILER_GNU_CLANG
  return 0;
#elif defined(__x86_64__) || defined(_M_X64)
  return 1;
#elif defined(__i386__) || defined(_M_IX86)
  unsigned int eflags_before;
  unsigned int eflags_after;

  __asm__ volatile("pushfl\n\t"
                   "popl %0\n\t"
                   "movl %0, %1\n\t"
                   "xorl $0x200000, %1\n\t"
                   "pushl %1\n\t"
                   "popfl\n\t"
                   "pushfl\n\t"
                   "popl %1\n\t"
                   "pushl %0\n\t"
                   "popfl\n\t"
                   : "=&r"(eflags_before), "=&r"(eflags_after)
                   :
                   : "cc");

  return ((eflags_before ^ eflags_after) & 0x200000u) != 0u;
#else
  return 0;
#endif
}

static int arbint_cpu_read_xcr0(uint64_t * out_xcr0) {
#if ARBINT_TARGET_X86_FAMILY && ARBINT_COMPILER_GNU_CLANG
  unsigned int eax;
  unsigned int edx;

  if (out_xcr0 == NULL)
    return 0;

  __asm__ volatile(".byte 0x0f, 0x01, 0xd0" : "=a"(eax), "=d"(edx) : "c"(0u));
  *out_xcr0 = ((uint64_t) edx << 32u) | (uint64_t) eax;
  return 1;
#else
  (void) out_xcr0;
  return 0;
#endif
}

static arbint_cpu_caps_t arbint_cpu_probe_caps(void) {
  arbint_cpu_caps_t caps;

  memset(&caps, 0, sizeof(caps));
  caps.initialized = 1u;

#if ARBINT_CPU_CAN_QUERY_X86_CPUID
  {
    unsigned int max_leaf;
    unsigned int eax;
    unsigned int ebx;
    unsigned int ecx;
    unsigned int edx;
    int avx_state_ok;
    int avx512_state_ok;

    if (!arbint_cpu_can_call_cpuid())
      return caps;
    caps.has_cpuid = 1u;

    max_leaf = arbint_cpu_cpuid_max_leaf();
    if (max_leaf < 1u)
      return caps;
    arbint_cpu_cpuid_count(1u, 0u, &eax, &ebx, &ecx, &edx);

    caps.sse = (edx & ARBINT_CPUID1_EDX_SSE) != 0u;
    caps.sse2 = (edx & ARBINT_CPUID1_EDX_SSE2) != 0u;
    caps.sse3 = (ecx & ARBINT_CPUID1_ECX_SSE3) != 0u;
    caps.ssse3 = (ecx & ARBINT_CPUID1_ECX_SSSE3) != 0u;
    caps.sse41 = (ecx & ARBINT_CPUID1_ECX_SSE41) != 0u;
    caps.sse42 = (ecx & ARBINT_CPUID1_ECX_SSE42) != 0u;
    caps.crc32 = caps.sse42;

    caps.popcnt = (ecx & ARBINT_CPUID1_ECX_POPCNT) != 0u;
    caps.pclmulqdq = (ecx & ARBINT_CPUID1_ECX_PCLMULQDQ) != 0u;
    caps.fma = (ecx & ARBINT_CPUID1_ECX_FMA) != 0u;
    caps.avx_hw = (ecx & ARBINT_CPUID1_ECX_AVX) != 0u;

    if ((ecx & (ARBINT_CPUID1_ECX_XSAVE | ARBINT_CPUID1_ECX_OSXSAVE)) ==
            (ARBINT_CPUID1_ECX_XSAVE | ARBINT_CPUID1_ECX_OSXSAVE) &&
        arbint_cpu_read_xcr0(&caps.xcr0)) {
      caps.has_xgetbv = 1u;
    }

    avx_state_ok = (caps.has_xgetbv && (caps.xcr0 & 0x6u) == 0x6u);
    avx512_state_ok = avx_state_ok && (caps.xcr0 & 0xe0u) == 0xe0u;

    caps.avx_os = avx_state_ok ? 1u : 0u;
    caps.avx = (caps.avx_hw && caps.avx_os) ? 1u : 0u;

    if (max_leaf >= 7u) {
      arbint_cpu_cpuid_count(7u, 0u, &eax, &ebx, &ecx, &edx);
      caps.bmi1 = (ebx & ARBINT_CPUID7_EBX_BMI1) != 0u;
      caps.bmi2 = (ebx & ARBINT_CPUID7_EBX_BMI2) != 0u;
      caps.avx2 =
          ((ebx & ARBINT_CPUID7_EBX_AVX2) != 0u && avx_state_ok) ? 1u : 0u;
      caps.avx512f =
          ((ebx & ARBINT_CPUID7_EBX_AVX512F) != 0u && avx512_state_ok) ? 1u
                                                                       : 0u;
      caps.avx512bw =
          ((ebx & ARBINT_CPUID7_EBX_AVX512BW) != 0u && avx512_state_ok) ? 1u
                                                                        : 0u;
      caps.avx512vl =
          ((ebx & ARBINT_CPUID7_EBX_AVX512VL) != 0u && avx512_state_ok) ? 1u
                                                                        : 0u;
      caps.sha = (ebx & ARBINT_CPUID7_EBX_SHA) != 0u;
    }
  }
#endif

  return caps;
}

int arbint_cpu_has_feature(arbint_cpu_feature_t feature) {
  static int ready = 0;
  static arbint_cpu_caps_t caps;

  if (!ready) {
    caps = arbint_cpu_probe_caps();
    ready = 1;
  }

  switch (feature) {
  case ARBINT_CPU_FEATURE_BMI2:
    return (int) caps.bmi2;
  case ARBINT_CPU_FEATURE_BMI1:
    return (int) caps.bmi1;
  case ARBINT_CPU_FEATURE_SSE:
    return (int) caps.sse;
  case ARBINT_CPU_FEATURE_SSE2:
    return (int) caps.sse2;
  case ARBINT_CPU_FEATURE_SSE3:
    return (int) caps.sse3;
  case ARBINT_CPU_FEATURE_SSSE3:
    return (int) caps.ssse3;
  case ARBINT_CPU_FEATURE_SSE41:
    return (int) caps.sse41;
  case ARBINT_CPU_FEATURE_SSE42:
    return (int) caps.sse42;
  case ARBINT_CPU_FEATURE_POPCNT:
    return (int) caps.popcnt;
  case ARBINT_CPU_FEATURE_CRC32:
    return (int) caps.crc32;
  case ARBINT_CPU_FEATURE_PCLMULQDQ:
    return (int) caps.pclmulqdq;
  case ARBINT_CPU_FEATURE_AVX:
    return (int) caps.avx;
  case ARBINT_CPU_FEATURE_AVX2:
    return (int) caps.avx2;
  case ARBINT_CPU_FEATURE_AVX512F:
    return (int) caps.avx512f;
  case ARBINT_CPU_FEATURE_AVX512BW:
    return (int) caps.avx512bw;
  case ARBINT_CPU_FEATURE_AVX512VL:
    return (int) caps.avx512vl;
  case ARBINT_CPU_FEATURE_FMA:
    return (int) caps.fma;
  case ARBINT_CPU_FEATURE_SHA:
    return (int) caps.sha;
  default:
    return 0;
  }
}
