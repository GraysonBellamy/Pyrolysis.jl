#!/usr/bin/env bash
# Build ThermaKin2Ds (Version 7) for Linux from a local, non-redistributable
# source copy. Nothing from the source tree is committed: the patched copy and
# the executable go to an ignored directory under docs/refs/.
#
#   build_thermakin.sh [SOURCE_DIR] [BUILD_DIR]
#
# The only portability change replaces the MSVC-only localtime_s/asctime_s
# start-time banner calls with their POSIX equivalents; numerics are untouched.
set -euo pipefail

repo=$(cd "$(dirname "$0")/../../../.." && pwd)
src=${1:-"$repo/docs/refs/ThermaKin2Ds Code"}
out=${2:-"$repo/docs/refs/thermakin_build"}

[[ -f "$src/ThermaKin.cpp" ]] || { echo "ThermaKin source not found in $src" >&2; exit 1; }
rm -rf "$out"
mkdir -p "$out/src"
cp "$src"/*.cpp "$src"/*.h "$out/src/"

python3 - "$out/src/InOutput.cpp" <<'EOF'
import sys
path = sys.argv[1]
text = open(path).read()
old = "\tlocaltime_s(&tm_buf,&beg);\n\tasctime_s(tm_out,32,&tm_buf);\n"
new = ("#ifdef _MSC_VER\n" + old + "#else\n"
       "\tlocaltime_r(&beg,&tm_buf);\n\tasctime_r(&tm_buf,tm_out);\n#endif\n")
if text.count(old) != 1:
    sys.exit("unexpected InOutput.cpp; portability patch not applied")
open(path, "w").write(text.replace(old, new))
EOF

cxx=${CXX:-g++}
"$cxx" -std=c++11 -O2 -o "$out/thermakin" "$out"/src/*.cpp
{
    echo "compiler: $("$cxx" --version | head -1)"
    echo "flags: -std=c++11 -O2"
    echo "source_sha256: $(cd "$src" && cat $(ls *.cpp *.h | LC_ALL=C sort) | sha256sum | cut -d' ' -f1)"
    echo "executable_sha256: $(sha256sum "$out/thermakin" | cut -d' ' -f1)"
} | tee "$out/build_info.txt"
echo "built $out/thermakin"
