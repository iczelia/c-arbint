#!/bin/sh
# run_wasm_test.sh -- run an emscripten test binary under node.
#
# Usage: run_wasm_test.sh ./test_name
#
# The argument is a libtool wrapper script.  We derive the actual
# .js file produced by emcc (in .libs/) and execute it with node.

test_path="$1"
test_name=$(basename "$test_path")
test_dir=$(cd "$(dirname "$test_path")" && pwd)
libs_dir="$test_dir/.libs"

# emcc/libtool may name the output with or without .js extension.
js_file=""
for candidate in "$libs_dir/$test_name" "$libs_dir/${test_name}.js"; do
  if [ -f "$candidate" ]; then
    js_file="$candidate"
    break
  fi
done

if [ -z "$js_file" ]; then
  echo "run_wasm_test.sh: no JS file found for $test_name in $libs_dir" >&2
  ls -la "$libs_dir/" 2>/dev/null >&2
  exit 99
fi

# Ensure the side module is available next to the test .js file.
# The emscripten runtime loads neededDynlibs relative to the .js location.
src_lib="$test_dir/../src/.libs/libarbint.so"
if [ -f "$src_lib" ] && [ ! -f "$libs_dir/libarbint.so" ]; then
  cp -f "$src_lib" "$libs_dir/"
fi

# Run from the test source directory so that data files (pi10k.txt,
# e10k.txt) are found by tests that open them relative to cwd.
cd "$test_dir"
exec node "$js_file"
