#!/bin/sh
# Write a SHA256SUMS manifest for a knowledgebase directory.
#
#   tools/make_kycgkb_sums.sh ~/repo/KYCGKB_hg38
#   tools/make_kycgkb_sums.sh ~/repo/KYCGKB_mm10
#
# WHY THIS FILE EXISTS
#   kycg fetches a collection by first fetching its manifest, checking that
#   manifest against a digest compiled into the binary, and then trusting every
#   file digest listed inside it. That is how the InfiniumAnnotation channel
#   already works: one pinned anchor per directory, everything else chaining
#   from it. The consequence worth having is that sets can be added upstream
#   without rebuilding kycg -- only the anchor moves, and only when the
#   maintainer decides it should.
#
#   The KYCGKB repositories carry the whole-genome knowledgebases but no
#   manifest, which is the one thing standing between them and that same
#   treatment. This script produces it.
#
# WHAT IT COVERS
#   Every data file: .cm knowledgebase sets, their .idx sidecars, and the
#   .cr whole-genome CpG reference. Not README.md, and not SHA256SUMS itself --
#   a manifest cannot list its own digest, and prose is not data.
#
#   Output format matches InfiniumAnnotation's exactly, so `shasum -a 256 -c
#   SHA256SUMS` verifies a directory with no kycg involved:
#
#     <64 hex>  <filename>
#
#   Sorted with LC_ALL=C, so the manifest is byte-reproducible on any machine.
#   That matters: the anchor is a digest *of this file*, so two people
#   regenerating it must agree to the byte.
set -eu

dir=${1:-.}

if [ ! -d "$dir" ]; then
    echo "make_kycgkb_sums.sh: no such directory: $dir" >&2
    exit 1
fi

cd "$dir"

# shasum is coreutils-independent and present on macOS; sha256sum on Linux.
if command -v shasum >/dev/null 2>&1; then
    hash_cmd="shasum -a 256"
elif command -v sha256sum >/dev/null 2>&1; then
    hash_cmd="sha256sum"
else
    echo "make_kycgkb_sums.sh: need shasum or sha256sum" >&2
    exit 1
fi

# Data files only, in a stable byte order.
files=$(LC_ALL=C ls | LC_ALL=C sort | grep -E '\.(cm|idx|cr)$' || true)

if [ -z "$files" ]; then
    echo "make_kycgkb_sums.sh: no .cm/.idx/.cr files in $dir" >&2
    exit 1
fi

tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT

for f in $files; do
    $hash_cmd "$f" >> "$tmp"
done

mv "$tmp" SHA256SUMS
trap - EXIT

n=$(wc -l < SHA256SUMS | tr -d ' ')
echo "wrote $dir/SHA256SUMS ($n files)" >&2
echo "anchor: $($hash_cmd SHA256SUMS | cut -d' ' -f1)" >&2
