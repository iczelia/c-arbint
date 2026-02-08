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

#if defined(ARBINT_HAS_WINAPI_ENTROPY)

  #include <windows.h>

  #ifdef HAVE_BCRYPT_H
    #include <bcrypt.h>
    #pragma comment(lib, "bcrypt.lib")
    #define USE_BCRYPT 1
  #else
    #include <wincrypt.h>
    #pragma comment(lib, "advapi.lib")
    #define USE_BCRYPT 0
  #endif

/*  Read entropy from Windows CSPRNG.
    Uses BCryptGenRandom (Vista+) or CryptGenRandom (XP+).
    Acquires and releases all handles within a single call, leaving no
    persistent state.  Callers should minimise the number of invocations.
    Parameters:
      dst - Output buffer for random bytes
      len - Number of bytes to generate
    Returns:
      0 on success, nonzero on failure  */
int arbint_entropy_winapi(uint8_t * dst, size_t len) {
  #if USE_BCRYPT
  NTSTATUS status;
  status = BCryptGenRandom(NULL, (PUCHAR) dst, (ULONG) len,
                           BCRYPT_USE_SYSTEM_PREFERRED_RNG);
  return (status == STATUS_SUCCESS) ? 0 : -1;
  #else
  HCRYPTPROV hProvider = 0;
  BOOL success;
  success = CryptAcquireContextW(&hProvider, NULL, NULL, PROV_RSA_FULL,
                                 CRYPT_VERIFYCONTEXT | CRYPT_SILENT);
  if (!success)
    return -1;
  success = CryptGenRandom(hProvider, (DWORD) len, (BYTE *) dst);
  CryptReleaseContext(hProvider, 0);
  return success ? 0 : -1;
  #endif
}

#endif
