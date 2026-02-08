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

  #if defined(HAVE_BCRYPT_H)
    #include <bcrypt.h>
    #pragma comment(lib, "bcrypt.lib")
    #define ARBINT_WIN_ENTROPY_BCRYPT 1
  #elif defined(HAVE_WINCRYPT_H)
    #include <wincrypt.h>
    #pragma comment(lib, "advapi32.lib")
    #define ARBINT_WIN_ENTROPY_WINCRYPT 1
  #else
    /*  Fallback: RtlGenRandom (aka SystemFunction036) from advapi32.dll.
        Available since Windows XP.  Needs only windows.h for the types;
        we declare the prototype manually because ntsecapi.h may be absent.
        Loaded at runtime via GetProcAddress to avoid a link-time dependency
        on advapi32.lib, which may not be available in minimal toolchains.  */
    #define ARBINT_WIN_ENTROPY_RTLGENRANDOM 1
  #endif /* HAVE_BCRYPT_H */

/*  Read entropy from Windows CSPRNG.
    Uses BCryptGenRandom (Vista+), CryptGenRandom (XP+), or
    RtlGenRandom (XP+, no special headers required) in that order.
    Acquires and releases all handles within a single call, leaving no
    persistent state.  Callers should minimise the number of invocations.
    Parameters:
      dst - Output buffer for random bytes
      len - Number of bytes to generate
    Returns:
      0 on success, nonzero on failure  */
int arbint_entropy_winapi(uint8_t * dst, size_t len) {
  #if defined(ARBINT_WIN_ENTROPY_BCRYPT)
  NTSTATUS status;
  status = BCryptGenRandom(NULL, (PUCHAR) dst, (ULONG) len,
                           BCRYPT_USE_SYSTEM_PREFERRED_RNG);
  return (status == STATUS_SUCCESS) ? 0 : -1;
  #elif defined(ARBINT_WIN_ENTROPY_WINCRYPT)
  HCRYPTPROV hProvider = 0;
  BOOL success;
  success = CryptAcquireContextW(&hProvider, NULL, NULL, PROV_RSA_FULL,
                                 CRYPT_VERIFYCONTEXT | CRYPT_SILENT);
  if (!success)
    return -1;
  success = CryptGenRandom(hProvider, (DWORD) len, (BYTE *) dst);
  CryptReleaseContext(hProvider, 0);
  return success ? 0 : -1;
  #elif defined(ARBINT_WIN_ENTROPY_RTLGENRANDOM)
  /*  SystemFunction036 is the exported name of RtlGenRandom in
      advapi32.dll.  We load it dynamically to avoid requiring an
      import library.  */
  typedef BOOLEAN(APIENTRY * RtlGenRandom_fn)(PVOID, ULONG);
  HMODULE hLib;
  RtlGenRandom_fn pfn;
  BOOLEAN ok;

  hLib = LoadLibraryA("advapi32.dll");
  if (hLib == NULL)
    return -1;
  pfn = (RtlGenRandom_fn) (void (*)(void)) GetProcAddress(hLib,
                                                          "SystemFunction036");
  if (pfn == NULL) {
    FreeLibrary(hLib);
    return -1;
  }
  ok = pfn((PVOID) dst, (ULONG) len);
  FreeLibrary(hLib);
  return ok ? 0 : -1;
  #endif /* ARBINT_WIN_ENTROPY_BCRYPT */
}

#endif /* defined(ARBINT_HAS_WINAPI_ENTROPY) */
