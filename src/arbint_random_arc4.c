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

#include "arbint_random.h"

#if defined(ARBINT_HAS_ARC4RANDOM)

  #include <stdlib.h> /* arc4random_buf */

/*  Read entropy from arc4random_buf.
    Uses system-provided CSPRNG (ChaCha20 on modern macOS/OpenBSD).
    Automatically seeded from kernel entropy.
    Thread-safety: arc4random_buf is thread-safe per POSIX.1-2024.
    Parameters:
      dst - Output buffer for random bytes
      len - Number of bytes to generate
    Returns:
      0 on success (cannot fail)
    Algorithm:
      1. Call arc4random_buf to fill buffer directly
      2. Return 0 (arc4random_buf has void return type)  */
int arbint_entropy_arc4(uint8_t * dst, size_t len) {
  arc4random_buf(dst, len);
  return 0;
}

#endif /* defined(ARBINT_HAS_ARC4RANDOM) */
