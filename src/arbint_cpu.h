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

#ifndef ARBINT_CPU_H
#define ARBINT_CPU_H

typedef enum arbint_cpu_feature {
  ARBINT_CPU_FEATURE_BMI2 = 0,
  ARBINT_CPU_FEATURE_BMI1,
  ARBINT_CPU_FEATURE_SSE,
  ARBINT_CPU_FEATURE_SSE2,
  ARBINT_CPU_FEATURE_SSE3,
  ARBINT_CPU_FEATURE_SSSE3,
  ARBINT_CPU_FEATURE_SSE41,
  ARBINT_CPU_FEATURE_SSE42,
  ARBINT_CPU_FEATURE_POPCNT,
  ARBINT_CPU_FEATURE_CRC32,
  ARBINT_CPU_FEATURE_PCLMULQDQ,
  ARBINT_CPU_FEATURE_AVX,
  ARBINT_CPU_FEATURE_AVX2,
  ARBINT_CPU_FEATURE_AVX512F,
  ARBINT_CPU_FEATURE_AVX512BW,
  ARBINT_CPU_FEATURE_AVX512VL,
  ARBINT_CPU_FEATURE_FMA,
  ARBINT_CPU_FEATURE_SHA,
  ARBINT_CPU_FEATURE_ARM_CRC32,
  ARBINT_CPU_FEATURE_ARM_PMULL
} arbint_cpu_feature_t;

int arbint_cpu_has_feature(arbint_cpu_feature_t feature);

#endif /* ARBINT_CPU_H */
