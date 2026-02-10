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

#if defined(ARBINT_HAS_GETENTROPY)

  #include <unistd.h>

/*  Read entropy via getentropy().
    This is the entropy path on Emscripten, where it routes to
    crypto.getRandomValues.
    Parameters:
      dst - Output buffer for random bytes
      len - Number of bytes to generate
    Returns:
      0 on success, nonzero on failure
    Note: getentropy() is limited to 256 bytes per call on most
    implementations, so we loop over chunks.  */
int arbint_entropy_getentropy(uint8_t * dst, size_t len) {
  while (len > 0) {
    size_t chunk = (len > 256) ? 256 : len;
    if (getentropy(dst, chunk) != 0)
      return -1;
    dst += chunk;
    len -= chunk;
  }
  return 0;
}

#endif /* defined(ARBINT_HAS_GETENTROPY) */
