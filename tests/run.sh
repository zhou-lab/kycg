#!/bin/bash
## Run kycg's whole test suite: the C unit binaries (pure functions over
## counts) and every tests/t_*.sh (the built binary, end to end).
##
##   make test                   # builds what is needed, then runs this
##   KYCG=/path/to/kycg YAME=/path/to/yame bash tests/run.sh
##
## Every shell test is self-contained: it packs its own fixtures from inline
## text with `yame pack`, makes its own temp directory, touches nothing
## shared, and needs no network and no data store. So this runs inside the
## conda test phase on every platform, against the installed package.
set -uo pipefail

here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)

## The binaries under test. Both are resolved from the environment first so
## the conda test phase can point at the INSTALLED package rather than a
## build tree -- that is the thing users get, and on macOS a freshly linked
## binary in the build directory cannot load libz yet.
export KYCG=${KYCG:-$root/kycg}
export YAME=${YAME:-$root/external/YAME/yame}
command -v "$KYCG" >/dev/null 2>&1 || [ -x "$KYCG" ] || {
  echo "no kycg at $KYCG (run make first)" >&2; exit 2; }
## yame is what packs the fixtures. kycg links libyame but does not ship the
## command, so a source tree needs `make yame-bin` and a conda test env needs
## the yame package -- both are how a USER gets it, which is the point.
command -v "$YAME" >/dev/null 2>&1 || [ -x "$YAME" ] || {
  echo "no yame at $YAME (run 'make yame-bin', or install the yame package)" >&2
  exit 2; }

pass=0; fail=0

## ---- the C unit tests: statistics and argv handling, no I/O ----
## Only in a build tree. In the conda test phase the binary under test is the
## INSTALLED kycg and tests/ arrives as source_files, so any test_* sitting
## there was compiled on somebody else's machine -- running it would report a
## failure about the packager's toolchain, not about this package.
for t in "$here"/test_*; do
  [ "$KYCG" = "$root/kycg" ] || continue
  case "$t" in *.c) continue;; esac
  [ -x "$t" ] || continue
  name=$(basename "$t")
  if out=$("$t" 2>&1); then
    pass=$((pass + 1)); echo "ok    $name"
  else
    fail=$((fail + 1)); echo "FAIL  $name"; echo "$out" | sed 's/^/      /'
  fi
done

## ---- the command-level tests ----
## Each makes its own temp directory and touches nothing shared, so they run
## JOBS at a time and the wall time is the slowest test rather than the sum.
## Output is captured per test and printed in name order, so the report reads
## the same as a serial run. JOBS=1 gives the serial run back.
ncpu=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
JOBS=${JOBS:-$(( ncpu < 8 ? ncpu : 8 ))}
logs=$(mktemp -d); trap 'rm -rf "$logs"' EXIT
running=0
## bash 4.3 brought `wait -n`. Asked by version rather than probed: `wait -n`
## with no children exits 127 even where it exists, so a probe reads as
## absent. macOS ships bash 3.2 and the conda test phase runs there.
have_wait_n=no
if [ "${BASH_VERSINFO[0]}" -gt 4 ] ||
   { [ "${BASH_VERSINFO[0]}" -eq 4 ] && [ "${BASH_VERSINFO[1]}" -ge 3 ]; }; then
  have_wait_n=yes
fi
## t_annotate builds a 27,722-row ordering and is by far the slowest; start it
## first or it begins after everything else and sets the wall time alone.
slow="t_annotate"
ordered=$(for n in $slow; do [ -f "$here/$n.sh" ] && echo "$here/$n.sh"; done
          for t in "$here"/t_*.sh; do
            b=$(basename "$t" .sh)
            [[ " $slow " == *" $b "* ]] || echo "$t"
          done)
for t in $ordered; do
  name=$(basename "$t" .sh)
  ( if bash "$t" > "$logs/$name.out" 2>&1; then : > "$logs/$name.ok"; fi ) &
  running=$((running + 1))
  if [ "$running" -ge "$JOBS" ]; then
    if [ "$have_wait_n" = yes ]; then wait -n; running=$((running - 1))
    else wait; running=0; fi
  fi
done
wait

for t in "$here"/t_*.sh; do
  name=$(basename "$t" .sh)
  if [ -e "$logs/$name.ok" ]; then
    pass=$((pass + 1)); echo "ok    $name"
  else
    fail=$((fail + 1)); echo "FAIL  $name"
    sed 's/^/      /' "$logs/$name.out"
  fi
done

echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ]
