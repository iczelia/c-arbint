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

#ifndef ARBINT_RANDOM_H
#define ARBINT_RANDOM_H

#include "arbint.h"
#include "arbint_base.h"

/*  MT19937 constants  */
#define ARBINT_MT_N 624
#define ARBINT_MT_M 397
#define ARBINT_MT_MATRIX_A 0x9908b0dfUL
#define ARBINT_MT_UPPER_MASK 0x80000000UL
#define ARBINT_MT_LOWER_MASK 0x7fffffffUL

/*  Platform entropy source (returns 0 on success, nonzero on failure)  */
typedef int (*arbint_entropy_fn_t)(uint8_t * dst, size_t len);

/*  Select platform entropy source  */
arbint_entropy_fn_t arbint_select_entropy_source(void);

/*  Platform-specific entropy implementations  */
#if defined(ARBINT_HAS_URANDOM)
int arbint_entropy_urandom(uint8_t * dst, size_t len);
#endif /* defined(ARBINT_HAS_URANDOM) */

#if defined(ARBINT_HAS_WINAPI_ENTROPY)
int arbint_entropy_winapi(uint8_t * dst, size_t len);
#endif /* defined(ARBINT_HAS_WINAPI_ENTROPY) */

#if defined(ARBINT_HAS_ARC4RANDOM)
int arbint_entropy_arc4(uint8_t * dst, size_t len);
#endif /* defined(ARBINT_HAS_ARC4RANDOM) */

#endif /* ARBINT_RANDOM_H */
