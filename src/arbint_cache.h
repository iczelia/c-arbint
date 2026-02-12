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

#ifndef ARBINT_CACHE_H
#define ARBINT_CACHE_H

#include "arbint_base.h"

arbint_err_t arbint_cache_prime_sieve_ensure(uint32_t limit);
int arbint_cache_prime_sieve_is_prime(uint32_t n);

void arbint_cache_drop_all(void);

#endif /* ARBINT_CACHE_H */
