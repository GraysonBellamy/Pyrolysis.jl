#!/usr/bin/env bash
# Build a ThermaKin2Ds copy whose RAND absorption mode takes its seed from the
# THERMAKIN_SEED environment variable instead of srand(time(NULL)), so random-
# absorption runs are reproducible. Nothing from the source tree is committed:
# the patched copy and the executable go to an ignored directory under docs/refs/.
#
#   build_thermakin_seeded.sh [SOURCE_DIR] [BUILD_DIR]
#
# Patches, both numerically inert: the start-time banner (MSVC localtime_s/
# asctime_s -> POSIX, as test/verification/reference/matched/build_thermakin.sh)
# and the RAND seed in Boundary.cpp (unset THERMAKIN_SEED keeps time(NULL)).
set -euo pipefail

repo=$(cd "$(dirname "$0")/../../.." && pwd)
src=${1:-"$repo/docs/refs/ThermaKin2Ds Code"}
out=${2:-"$repo/docs/refs/thermakin_build_seeded"}

[[ -f "$src/ThermaKin.cpp" ]] || { echo "ThermaKin source not found in $src" >&2; exit 1; }
case "$out" in "$repo"/docs/refs/*) ;; *) echo "refusing to build outside docs/refs: $out" >&2; exit 1 ;; esac
rm -rf "$out"
mkdir -p "$out/src"
cp "$src"/*.cpp "$src"/*.h "$out/src/"

python3 - "$out/src" <<'PY'
import sys, pathlib
d = pathlib.Path(sys.argv[1])
p = d / "InOutput.cpp"
t = p.read_text()
old = "\tlocaltime_s(&tm_buf,&beg);\n\tasctime_s(tm_out,32,&tm_buf);\n"
new = ("#ifdef _MSC_VER\n" + old + "#else\n"
       "\tlocaltime_r(&beg,&tm_buf);\n\tasctime_r(&tm_buf,tm_out);\n#endif\n")
if t.count(old) != 1:
    sys.exit("unexpected InOutput.cpp; portability patch not applied")
p.write_text(t.replace(old, new))

p = d / "Boundary.cpp"
t = p.read_text()
old = "srand((unsigned)time(NULL));"
new = ("{ const char* tk_seed=getenv(\"THERMAKIN_SEED\"); "
       "srand(tk_seed ? (unsigned)strtoul(tk_seed,NULL,10) : (unsigned)time(NULL)); }")
if t.count(old) != 1:
    sys.exit("unexpected Boundary.cpp; seed patch not applied")
t = t.replace(old, new)
t = t.replace('#include "Boundary.h"', '#include "Boundary.h"\n#include <cstdlib>\n#include <ctime>', 1)
p.write_text(t)
PY

cxx=${CXX:-g++}
"$cxx" -std=c++11 -O2 -o "$out/thermakin" "$out"/src/*.cpp
{
    echo "compiler: $("$cxx" --version | head -1)"
    echo "flags: -std=c++11 -O2"
    echo "patches: InOutput.cpp banner (POSIX time); Boundary.cpp RAND seed from THERMAKIN_SEED"
    echo "source_sha256: $(cd "$src" && cat $(ls *.cpp *.h | LC_ALL=C sort) | sha256sum | cut -d' ' -f1)"
    echo "executable_sha256: $(sha256sum "$out/thermakin" | cut -d' ' -f1)"
} | tee "$out/build_info.txt"
echo "built $out/thermakin"
