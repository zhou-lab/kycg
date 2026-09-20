#!/bin/sh
# scripts/coverage.sh — line coverage of tests/run.sh over kycg's own sources.
#
#   scripts/coverage.sh            measure, print the table, rewrite docs/coverage.json
#   scripts/coverage.sh --print    measure and print only (leave the badge alone)
#   scripts/coverage.sh --check    measure and fail if docs/coverage.json is stale
#                                  by more than $TOLERANCE points (default 2.0)
#
# The number is gcov's own summary over src/*.c, counting EXECUTABLE lines. It
# is LINE coverage, not branch coverage -- branch coverage is lower, and for
# the argument parsers it is the more telling number, so measure it before
# quoting it.
#
# Only kycg's src/ is counted. YAME is a submodule with its own suite and its
# own badge; folding its ~20k lines in here would produce a number that says
# nothing about either.
#
# Everything happens in a scratch copy. An instrumented binary must never
# become the one in the repo: it is built -O0 and writes .gcda files beside
# itself on every run.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$here"
TOLERANCE=${TOLERANCE:-2.0}
mode=${1:---write}

command -v gcov >/dev/null || { echo "coverage: gcov not found (it ships with gcc)" >&2; exit 1; }
[ -f external/YAME/libyame.a ] || {
  echo "coverage: external/YAME/libyame.a is missing; run make first" >&2; exit 1; }
[ -x external/YAME/yame ] || {
  echo "coverage: external/YAME/yame is missing; run 'make yame-bin' first" >&2; exit 1; }

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
cp -r src tests tools Makefile "$work/"
# tools/ too: since the v1.50 bump kycg emits its own registry, and
# tests/t_docs.sh asserts the generator is present and that the committed
# header is regenerable from it. Without tools/ that check fails in the
# scratch copy and takes the whole suite with it.
#
# The docs tests read these and skip without them, which would quietly drop
# their lines out of the number.
cp -r docs "$work/"
cp README.md "$work/"
[ -d conda-recipe ] && cp -r conda-recipe "$work/"
# YAME is SYMLINKED, not copied: it is already built, it is not what we are
# measuring, and rebuilding htslib into the scratch tree costs minutes per
# run. Its objects carry no coverage flags, so nothing of YAME's lands in the
# gcov output even though every test links it.
ln -s "$here/external" "$work/external"
rm -f "$work"/src/*.o "$work"/src/*.gcda "$work"/src/*.gcno

# -O0 so line numbers map 1:1; optimisation merges and elides lines and the
# annotation stops meaning anything. The flags ride on CC because the Makefile
# ASSIGNS CFLAGS rather than appending to it -- a command-line CFLAGS would
# drop the -I flags the build needs. conda-recipe/build.sh does the same.
# `make` builds the CLI; `test-bins` adds the unit binaries, which link the
# statistics and digest objects directly. Without them those files show only
# what the CLI happens to reach and the number understates the suite badly.
( cd "$work" && make CC="${CC:-cc} -O0 -g --coverage" LDFLAGS="--coverage" \
             && make test-bins CC="${CC:-cc} -O0 -g --coverage" LDFLAGS="--coverage" ) \
    >"$work/build.log" 2>&1 || true
[ -x "$work/kycg" ] || { tail -20 "$work/build.log" >&2; echo "coverage: build failed" >&2; exit 1; }

# Counters accumulate across processes, so clear them and let ONLY the suite
# run: a stray `kycg --version` adds main.c's version path and moves the total.
rm -f "$work"/src/*.gcda
KYCG="$work/kycg" YAME="$here/external/YAME/yame" \
  CC="${CC:-cc} -O0 -g --coverage" bash "$work/tests/run.sh" >"$work/suite.log" 2>&1 || {
    tail -30 "$work/suite.log" >&2
    echo "coverage: the suite failed; coverage of a red suite is meaningless" >&2
    exit 1
}

( cd "$work" && gcov -n -o src src/*.c ) >"$work/gcov.txt" 2>/dev/null || true

pct=$(python3 - "$work/gcov.txt" <<'PY'
import re, sys
t = open(sys.argv[1]).read()
rows = re.findall(r"File '([^']+)'\nLines executed:([0-9.]+)% of (\d+)", t)
rows = [(f, float(p), int(n)) for f, p, n in rows if f.startswith("src/")]
if not rows:
    sys.exit("coverage: gcov produced no rows for src/")
tot = sum(n for _, _, n in rows)
cov = sum(p / 100 * n for _, p, n in rows)
for f, p, n in sorted(rows, key=lambda r: r[1]):
    print(f"  {f:<24}{n:>7}{p:>8.1f}%", file=sys.stderr)
print(f"  {'TOTAL':<24}{tot:>7}{100*cov/tot:>8.1f}%", file=sys.stderr)
print(f"{100*cov/tot:.1f}")
PY
)

echo "coverage: ${pct}% of executable lines (gcov, line coverage: tests/run.sh)"

badge="$here/docs/coverage.json"
case "$mode" in
  --print) exit 0 ;;
  --check)
      [ -f "$badge" ] || { echo "coverage: $badge is missing; run scripts/coverage.sh" >&2; exit 1; }
      old=$(sed -n 's/.*"message" *: *"\([0-9.]*\)%".*/\1/p' "$badge")
      awk -v a="$old" -v b="$pct" -v t="$TOLERANCE" 'BEGIN{d=a-b; if(d<0)d=-d; exit !(d>t)}' \
          && { echo "coverage: badge says ${old}% but the suite measures ${pct}% (> ${TOLERANCE} points); run scripts/coverage.sh and commit docs/coverage.json" >&2; exit 1; }
      echo "coverage: badge (${old}%) is within ${TOLERANCE} points"
      exit 0 ;;
esac

colour=$(awk -v p="$pct" 'BEGIN{
  print (p>=90)?"brightgreen":(p>=80)?"green":(p>=70)?"yellowgreen":(p>=60)?"yellow":(p>=50)?"orange":"red"}')
cat > "$badge" <<JSON
{
  "schemaVersion": 1,
  "label": "coverage",
  "message": "${pct}%",
  "color": "${colour}"
}
JSON
echo "coverage: wrote docs/coverage.json (${pct}%, ${colour})"
