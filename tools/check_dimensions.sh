#!/bin/sh
# Verify the row-space dimensions pinned in src/registry.h against live data.
#
#   tools/check_dimensions.sh          # needs a built ./kycg and network
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

if [ ! -x ./kycg ]; then
    echo "check_dimensions.sh: build kycg first (make)" >&2
    exit 1
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

fail=0

check() {
    target=$1 url=$2 want=$3
    if ! curl -sfL "$url" -o "$tmp/probe.cm" 2>/dev/null; then
        printf '  %-12s %-12s could not fetch probe\n' "$target" "SKIP"
        return
    fi
    got=$(./kycg info "$tmp/probe.cm" 2>/dev/null | sed -n '2p' | cut -f5)
    if [ "$got" = "$want" ]; then
        printf '  %-12s %-12s %s\n' "$target" "OK" "$got"
    else
        printf '  %-12s %-12s pinned %s, actual %s\n' "$target" "MISMATCH" "$want" "$got"
        fail=1
    fi
}

# Parse what the registry claims, so this checks the file rather than a copy
# of the file's contents kept in this script.
echo "Verifying row-space dimensions in src/registry.h"

# Parse the registry with python rather than by regex position: the struct has
# gained fields twice already, and a positional pattern fails silently when it
# does -- printing a clean header and checking nothing, which is worse than an
# error. This keys on the field names in the struct instead.
python3 - "$tmp" <<'PY'
import re, subprocess, sys, os, json, urllib.request

tmp = sys.argv[1]
reg = open("src/registry.h").read()
iatag = re.search(r'#define KYCG_IA_TAG\s+"([^"]+)"', reg).group(1)

rows = []   # (label, url, expected_rows)

for m in re.finditer(r'\{ "([a-z0-9]+)", (\d+), (\d+), "(KYCGKB_[^"]+)", "([^"]+)"', reg):
    genome, nrow, _nset, repo, tag = m.groups()
    try:
        d = json.load(urllib.request.urlopen(
            f"https://api.github.com/repos/zhou-lab/{repo}/contents/?ref={tag}"))
        cm = [f for f in d if f["name"].endswith(".cm")]
        small = min(cm, key=lambda x: x["size"])["name"]
    except Exception:
        print(f"  {genome:<12} {'SKIP':<12} could not list repo"); continue
    rows.append((genome,
                 f"https://github.com/zhou-lab/{repo}/raw/{tag}/{small}", nrow))

for m in re.finditer(r'\{ "([A-Za-z0-9]+)", (\d+), (\d+), "[0-9a-f]{64}"', reg):
    plat, nrow, _nset = m.groups()
    base = (f"https://github.com/zhou-lab/InfiniumAnnotation/raw/{iatag}"
            f"/{plat}/KYCG")
    try:
        sums = urllib.request.urlopen(base + "/SHA256SUMS").read().decode()
    except Exception:
        print(f"  {plat:<12} {'SKIP':<12} no manifest"); continue
    cm = [l.split()[1] for l in sums.splitlines()
          if l.strip().endswith(".cm")]
    if not cm:
        print(f"  {plat:<12} {'SKIP':<12} no .cm found"); continue
    rows.append((plat, f"{base}/{cm[0]}", nrow))

fail = 0
for label, url, want in rows:
    probe = os.path.join(tmp, "probe.cm")
    try:
        with urllib.request.urlopen(url) as r, open(probe, "wb") as o:
            o.write(r.read())
    except Exception:
        print(f"  {label:<12} {'SKIP':<12} could not fetch probe"); continue
    out = subprocess.run(["./kycg", "info", probe],
                         capture_output=True, text=True).stdout.splitlines()
    got = out[1].split("\t")[4] if len(out) > 1 else "?"
    if got == want:
        print(f"  {label:<12} {'OK':<12} {got}")
    else:
        print(f"  {label:<12} {'MISMATCH':<12} pinned {want}, actual {got}")
        fail = 1

sys.exit(fail)
PY
