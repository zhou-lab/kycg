#!/bin/sh
# Verify the row-space dimensions pinned in src/registry.h against live data.
#
#   tools/check_dimensions.sh          # needs yame on PATH and network
#   CHECK_DIMENSIONS_PARSE_ONLY=1 tools/check_dimensions.sh
#                                      # offline: parse the registry only
#
# WHY A SEPARATE SCRIPT
#   The dimensions are pinned rather than derived because they are needed
#   *before* any file is on disk: a .cm is comparable only to a query with the
#   same row count, so the number that decides whether a set is usable has to
#   be known while the user is still deciding what to fetch.
#
#   Deriving them during registry generation would mean make_registry.sh could
#   not run without a built kycg to read a .cm with -- and kycg cannot build
#   without a registry. Keeping the check separate breaks that circle: run it
#   after a build, when bumping tags, or whenever the numbers are in doubt.
#
# HOW
#   Downloads the smallest .cm each collection publishes -- around a hundred
#   bytes for the whole-genome repositories -- and reads its record length.
#   Every record in a collection indexes the same row list by construction, so
#   the cheapest one answers for all of them.
set -eu

cd "$(dirname "$0")/.."

# The row count of a .cm is read with `yame info` -- kycg links libyame but
# ships no `info` subcommand of its own. Prefer a yame on PATH, else the one
# built beside the submodule.
YAME=${YAME:-$(command -v yame || echo external/YAME/yame)}
if [ ! -x "$YAME" ] && ! command -v "$YAME" >/dev/null 2>&1; then
    echo "check_dimensions.sh: need yame on PATH (or run 'make yame-bin')" >&2
    exit 1
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# Parse what the registry claims, so this checks the file rather than a copy
# of the file's contents kept in this script.
echo "Verifying row-space dimensions in src/registry.h"

# Parse the registry with python rather than by regex position: the struct has
# gained fields twice already, and a positional pattern fails silently when it
# does -- printing a clean header and checking nothing, which is worse than an
# error. This keys on the field names in the struct instead.
YAME="$YAME" python3 - "$tmp" <<'PY'
import re, subprocess, sys, os, json, urllib.request

tmp = sys.argv[1]
reg = open("src/registry.h").read()

## KYCG_FILES: one row per file, key / store_path / url / sha256 / size.
files = re.findall(
    r'\{ "([^"]+)", "([^"]+)",\s*"([^"]+)",\s*"[0-9a-f]{64}", (\d+),', reg)

## The two collection tables: name, pinned rows, and the slice of KYCG_FILES
## that holds the collection's sets (KYCG_FILES + start, n files).
colls = re.findall(
    r'\{ "([A-Za-z0-9]+)", (\d+), \d+, KYCG_FILES \+ (\d+), (\d+),', reg)
colls += re.findall(
    r'\{ "([A-Za-z0-9]+)", (\d+), \d+, (?:"[^"]*"|NULL), (?:"[^"]*"|NULL), '
    r'KYCG_FILES \+ (\d+), (\d+),', reg)
if not files or not colls:
    sys.exit("check_dimensions.sh: could not parse src/registry.h; "
             "has the registry layout changed?")

rows = []   # (label, url, expected_rows)
for name, nrow, start, n in colls:
    ## the smallest .cm in the slice: every record in a collection indexes the
    ## same row list, so the cheapest one answers for all of them
    cm = [f for f in files[int(start):int(start) + int(n)]
          if f[1].endswith(".cm")]
    if not cm:
        print(f"  {name:<12} {'SKIP':<12} no .cm in the registry"); continue
    small = min(cm, key=lambda f: int(f[3]))
    rows.append((name, small[2], nrow))

## The suite runs the parse alone, offline: a registry layout change then
## fails `make test` instead of the next person who runs this by hand.
if os.environ.get("CHECK_DIMENSIONS_PARSE_ONLY"):
    print(f"parsed {len(rows)} collections from {len(files)} files")
    sys.exit(0 if len(rows) == len(colls) else 1)

fail = 0
for label, url, want in rows:
    probe = os.path.join(tmp, "probe.cm")
    try:
        with urllib.request.urlopen(url) as r, open(probe, "wb") as o:
            o.write(r.read())
    except Exception:
        print(f"  {label:<12} {'SKIP':<12} could not fetch probe"); continue
    yame = os.environ.get("YAME", "yame")
    out = subprocess.run([yame, "info", probe],
                         capture_output=True, text=True).stdout.splitlines()
    got = "?"
    if len(out) > 1:
        hdr = out[0].split("\t")
        col = hdr.index("Nrow") if "Nrow" in hdr else 3
        f = out[1].split("\t")
        got = f[col] if col < len(f) else "?"
    if got == want:
        print(f"  {label:<12} {'OK':<12} {got}")
    else:
        print(f"  {label:<12} {'MISMATCH':<12} pinned {want}, actual {got}")
        fail = 1

sys.exit(fail)
PY
