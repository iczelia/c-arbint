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

## WebAssembly / Emscripten

c-arbint builds as a wasm side module for use in browser and Node.js environments.
Both wasm32 and wasm64 (Memory64) are supported.

### Build

```sh
./bootstrap
emconfigure ./configure --enable-shared --disable-static \
  --host=wasm32-unknown-emscripten
emmake make -j"$(nproc)"
emmake make check    # requires Node.js (24+ for wasm64)
```

For wasm64, use `--host=wasm64-unknown-emscripten` and pass
`CFLAGS="-sMEMORY64" LDFLAGS="-sMEMORY64"` to configure.

### Output

The build produces `src/libarbint.wasm`, a SIDE_MODULE that can be loaded at
runtime by an emscripten MAIN_MODULE application.  Symbol visibility matches the
native shared library (`-fvisibility=hidden` + `ARBINT_API` exports).

### Link a test program

```sh
emcc -sMAIN_MODULE=2 your_file.c src/libarbint.wasm -I include -o your_file.js
node your_file.js
```

## Supported platforms

The platforms that the code has been verified to compile and pass tests on include:

Primary platforms (with LTO support):
- Linux: x86-64 (Ubuntu, Fedora, Alpine; gcc/clang/tcc), arm64 (Ubuntu; gcc/clang/tcc), i386 (Debian; gcc/clang/tcc)
- macOS: arm64 (clang)
- Windows: x86-64, i686 (MinGW-GCC, TCC, MSVC)
- WebAssembly: wasm32, wasm64 (Emscripten; with and without LTO)

Exotic architectures (no LTO, tested via Docker QEMU emulation):
- Linux: riscv64 (Debian sid, Alpine edge; gcc), s390x (Debian, Alpine; clang/gcc), ppc64le (Debian, Alpine; clang/gcc), mips64le (Debian; gcc/clang), arm32v7 (Debian, Alpine; gcc), arm32v5 (Debian; gcc)
- Linux (Debian ports): ppc64 big-endian (gcc), powerpc 32-bit (gcc), hppa/PA-RISC (gcc), m68k/Motorola 68000 (gcc), sh4/SuperH (gcc), alpha/DEC Alpha (gcc)
