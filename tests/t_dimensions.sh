#!/bin/bash
## tools/check_dimensions.sh reads src/registry.h by field layout, and broke
## silently once when the registry became per-file. The suite runs its parse
## offline, so the next layout change fails here and not for the next person
## who runs the script by hand. (The download half needs the network and is
## left to that person.)
set -uo pipefail
. "$(dirname "$0")/lib.sh"
root=$(cd "$(dirname "$0")/.." && pwd)

## Only in a tree that has the script and the registry: the conda test phase
## ships tests/ and docs/ only.
if [ ! -f "$root/tools/check_dimensions.sh" ] || [ ! -f "$root/src/registry.h" ]; then
  echo "  skip: no tools/check_dimensions.sh or src/registry.h here"
  done_testing
fi
n_coll=$(grep -c '^    { "[A-Za-z0-9]*", [0-9]*, [0-9]*, ' "$root/src/registry.h")
out=$(cd "$root" && YAME=${YAME:-yame} CHECK_DIMENSIONS_PARSE_ONLY=1 \
      sh tools/check_dimensions.sh 2>&1)
check_has "the registry parses" "parsed $n_coll collections" "$out"
done_testing
