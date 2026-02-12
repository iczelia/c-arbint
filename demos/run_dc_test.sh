#!/bin/sh
# run_dc_test.sh - wrapper to run dc demo in test mode
#
# This script is used by 'make check' to run dc with --test flag.
# Works for both native and emscripten builds.

set -e

# DC_BINARY is set by TESTS_ENVIRONMENT in Makefile.am
if [ -z "$DC_BINARY" ]; then
    # Fallback for manual execution
    DC_BINARY="./dc"
fi

# For emscripten builds, use the wasm test runner
if [ -n "$WASM_RUNNER" ] && [ -f "$WASM_RUNNER" ]; then
    exec "$WASM_RUNNER" "$DC_BINARY" --test
fi

# Native build - run directly
exec "$DC_BINARY" --test
