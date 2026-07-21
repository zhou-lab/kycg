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

# Whole genomes: { "hg38", 29401795, "KYCGKB_hg38", "v2", ... }
grep -oE '\{ "[a-z0-9]+", [0-9]+, "KYCGKB_[^"]+", "[^"]+"' src/registry.h |
while IFS= read -r line; do
    genome=$(echo "$line" | sed -E 's/.*\{ "([^"]+)".*/\1/')
    rows=$(echo   "$line" | sed -E 's/.*", ([0-9]+),.*/\1/')
    repo=$(echo   "$line" | sed -E 's/.*, "(KYCGKB_[^"]+)".*/\1/')
    tag=$(echo    "$line" | sed -E 's/.*"([^"]+)"$/\1/')
    # The smallest .cm in the repo is enough; ChromosomeXY is ~100 bytes.
    small=$(curl -sfL "https://api.github.com/repos/zhou-lab/$repo/contents/?ref=$tag" \
            | python3 -c '
import sys, json
d = json.load(sys.stdin)
cm = [f for f in d if f["name"].endswith(".cm")]
print(min(cm, key=lambda x: x["size"])["name"] if cm else "")' 2>/dev/null || true)
    [ -n "$small" ] || { printf '  %-12s %-12s no .cm found\n' "$genome" "SKIP"; continue; }
    check "$genome" "https://github.com/zhou-lab/$repo/raw/$tag/$small" "$rows"
done

# Arrays: { "MSA", 284309, "<sha>" } under the one InfiniumAnnotation tag.
iatag=$(sed -n 's/^#define KYCG_IA_TAG *"\([^"]*\)".*/\1/p' src/registry.h)
grep -oE '\{ "[A-Za-z0-9]+", [0-9]+, "[0-9a-f]{64}" \}' src/registry.h |
while IFS= read -r line; do
    plat=$(echo "$line" | sed -E 's/.*\{ "([^"]+)".*/\1/')
    rows=$(echo "$line" | sed -E 's/.*", ([0-9]+),.*/\1/')
    base="https://github.com/zhou-lab/InfiniumAnnotation/raw/$iatag/$plat/KYCG"
    small=$(curl -sfL "$base/SHA256SUMS" 2>/dev/null |
            awk '{print $2}' | grep '\.cm$' | head -1 || true)
    [ -n "$small" ] || { printf '  %-12s %-12s no .cm found\n' "$plat" "SKIP"; continue; }
    check "$plat" "$base/$small" "$rows"
done

exit $fail
