# c-arbint

c-arbint - portable arbitrary-precision integer library for C, providing signed multi-precision arithmetic with runtime optimization selection and explicit memory management.
Licensed under the terms of GNU Affero GPL version 3 or later - see COPYING.
Report issues to Kamila Szewczyk k@iczelia.net.
Project homepage: https://github.com/iczelia/c-arbint

[![CI](https://github.com/iczelia/c-arbint/actions/workflows/ci.yml/badge.svg?branch=trunk)](https://github.com/iczelia/c-arbint/actions/workflows/ci.yml)

## Requirements

- C compiler with C99 support
- autotools + libtool (`autoconf`, `automake`, `libtool`), `make` (for building)
- `nroff` (for man page generation, optional)
- `pkg-config` (for consumer integration)

## Build from git

```sh
./bootstrap
./configure
make -j"$(nproc)"
make check
```

## Install

```sh
sudo make install
```

## Configure options

- `--with-base-optimization=default|tuned`
- `--with-march=<arch>` (used with `tuned`, `native` if not specified)
- `--enable-warnings`
- `--enable-lto=auto|yes|no`
  - `auto` (default): enables `-flto` if supported by the compiler/linker
  - `yes`: requires `-flto` support and fails otherwise
  - `no`: disables LTO
- `--disable-entropy` - do not compile platform entropy sources; `arbint_rng_init`
  with a NULL seed will always fail and callers must supply their own seed buffer

## Link in your project

```sh
cc your_file.c $(pkg-config --cflags --libs c-arbint)
```

## Supported platforms

The platforms that the code has been verified to compile and pass tests on include:
- Linux: x86-64 (Ubuntu, Fedora, Alpine; gcc/clang/tcc), arm64 (Ubuntu; gcc/clang/tcc), i386 (Debian; gcc/clang/tcc), riscv64 (Ubuntu; gcc/clang), s390x (Ubuntu; gcc/clang), ppc64le (Ubuntu; clang).
- macOS: x86-64, arm64 (clang).
- Windows: x86-64, i686 (via MinGW-GCC and TCC and MSVC).

Note: The following configurations are disabled in CI due to compiler segfaults/internal errors:
- ppc64le/s390x with GCC (all LTO settings)
- riscv64 with clang