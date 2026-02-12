#!/bin/sh
# run_wasm_test.sh -- run an emscripten test binary under node.
#
# Usage: run_wasm_test.sh ./test_name
#
# The argument is a libtool wrapper script (or the JS file itself).
# We derive the actual .js file produced by emcc and execute it with
# node.

test_path="$1"
shift
extra_args="$*"
test_name=$(basename "$test_path")
test_dir=$(cd "$(dirname "$test_path")" && pwd)
libs_dir="$test_dir/.libs"
node_bin="${NODE:-node}"

# emcc + libtool may place the JS output in several locations depending
# on the emscripten and libtool versions:
#   1. .libs/<name>      -- extensionless JS (emscripten <= 3.x)
#   2. .libs/<name>.js   -- explicit .js extension (emscripten >= 4.x)
#   3. <name>.js         -- next to the wrapper when libtool skips .libs
# As a last resort, the wrapper itself may be the JS file.
js_file=""
for candidate in \
  "$libs_dir/$test_name" \
  "$libs_dir/${test_name}.js" \
  "$test_dir/${test_name}.js" \
  "$test_path"
do
  [ -f "$candidate" ] || continue
  # Skip libtool shell wrappers (they start with "#! /bin/bash" or
  # "#!/bin/sh" followed by a comment mentioning "wrapper script").
  case "$(head -n 2 "$candidate" 2>/dev/null)" in
    *"wrapper script"*) continue ;;
  esac
  js_file="$candidate"
  break
done

if [ -z "$js_file" ]; then
  echo "run_wasm_test.sh: no JS file found for $test_name" >&2
  echo "  looked in: $libs_dir/ and $test_dir/" >&2
  echo "  libs_dir contents:" >&2
  ls -la "$libs_dir/" 2>/dev/null >&2 || echo "    (directory does not exist)" >&2
  echo "  test_dir contents matching $test_name*:" >&2
  ls -la "$test_dir/$test_name"* 2>/dev/null >&2 || echo "    (none)" >&2
  exit 99
fi

# Keep an absolute JS path before changing cwd for data files.
js_file=$(cd "$(dirname "$js_file")" && pwd)/$(basename "$js_file")

# Ensure the side module (wasm) is next to the JS file so the
# emscripten runtime can find it via neededDynlibs.
# The libtool-bypass build places libarbint.wasm in src/;
# the legacy libtool build places libarbint.so in src/.libs/.
js_dir=$(dirname "$js_file")
for lib_name in libarbint.wasm libarbint.so; do
  src_lib="$test_dir/../src/$lib_name"
  if [ ! -f "$src_lib" ]; then
    src_lib="$test_dir/../src/.libs/$lib_name"
  fi
  if [ -f "$src_lib" ] && [ ! -f "$js_dir/$lib_name" ]; then
    cp -f "$src_lib" "$js_dir/"
  fi
done

# Also ensure the .wasm companion file is next to the JS if the JS
# lives outside .libs/ (e.g. tests/test_name.js alongside tests/test_name.wasm).
wasm_in_libs="$libs_dir/${test_name}.wasm"
wasm_next_to_js="$js_dir/${test_name}.wasm"
if [ -f "$wasm_in_libs" ] && [ ! -f "$wasm_next_to_js" ]; then
  cp -f "$wasm_in_libs" "$js_dir/"
fi

# Run from the source directory when available so tests that open data files
# relative to cwd (pi10k.txt, e10k.txt) work in out-of-tree builds.
run_dir="$test_dir"
if [ -n "$srcdir" ]; then
  src_abs=$(cd "$srcdir" 2>/dev/null && pwd)
  if [ -n "$src_abs" ]; then
    run_dir="$src_abs"
  fi
fi
cd "$run_dir"

# If this build targets Memory64 but the runtime lacks support,
# skip rather than reporting a spurious test failure.
if grep -q "'index': 'i64'" "$js_file"; then
  if ! "$node_bin" -e "new WebAssembly.Memory({initial: 1n, index: 'i64'});" \
      >/dev/null 2>&1; then
    echo "run_wasm_test.sh: skipping $test_name (runtime has no wasm Memory64 support)" >&2
    exit 77
  fi
fi

exec "$node_bin" "$js_file" $extra_args
