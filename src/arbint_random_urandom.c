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

#if defined(ARBINT_HAS_URANDOM)

  #include <fcntl.h>
  #include <unistd.h>

/*  Read entropy from /dev/urandom.
    Opens, reads, and closes /dev/urandom in a single call, leaving no
    persistent file descriptor.  Callers should minimise the number of
    invocations (ideally one per RNG lifetime).
    Parameters:
      dst - Output buffer for random bytes
      len - Number of bytes to generate
    Returns:
      0 on success, nonzero on failure  */
int arbint_entropy_urandom(uint8_t * dst, size_t len) {
  int fd;
  ssize_t nread, total_read;
  fd = open("/dev/urandom", O_RDONLY);
  if (fd == -1)
    return -1;
  total_read = 0;
  while ((size_t) total_read < len) {
    nread = read(fd, dst + total_read, len - (size_t) total_read);
    if (nread <= 0) {
      close(fd);
      return -1;
    }
    total_read += nread;
  }
  close(fd);
  return 0;
}

#endif /* defined(ARBINT_HAS_URANDOM) */
