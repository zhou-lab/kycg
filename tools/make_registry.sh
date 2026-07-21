#!/bin/sh
# Regenerate src/registry.h by pinning the digests each channel publishes.
#
#   tools/make_registry.sh [tag] > src/registry.h
#
# Two channels, two publishers, two digest schemes -- kycg verifies against
# whatever the publisher actually signs rather than inventing its own.
#
#   Arrays (InfiniumAnnotation, git raw at a tag)
#     Layout: <base>/<tag>/<PLATFORM>/KYCG/{SHA256SUMS, <Set>.<date>.cm[.idx]}
#     We pin sha256(<PLATFORM>/KYCG/SHA256SUMS). Every individual file digest
#     chains from that one anchor, so the pin is a single line per platform and
#     adding a set upstream does not require regenerating anything but the
#     anchor. A platform not published at the tag gets a NULL anchor and cannot
#     be fetched.
#
#   Whole genome (KYCGKB_<genome>, git raw at a tag)
#     Layout: <repo>/<tag>/{SHA256SUMS, cpg_nocontig.cr, <Set>.<date>.cm[.idx]}
#     Identical treatment to the arrays, since these repositories now publish a
#     SHA256SUMS of their own. Before that they were pinned file-by-file
#     against Zenodo's md5, which was the only digest that channel offered;
#     one anchor per genome replaces 82 pinned entries and lets sets be added
#     upstream without regenerating anything but the anchor.
#
#     Zenodo remains the archival record and keeps the DOI -- it is what the
#     paper cites -- so the record id and DOI stay in the registry as
#     provenance. They are no longer the fetch path.
#
#     File sizes are pinned separately, for display only. They let the browser
#     and the fetch plan say how large a download will be before it starts.
#     Nothing depends on them for correctness, and a file missing from the
#     table simply shows no size -- so adding a set upstream needs no rebuild.
#
# Regenerating requires network but downloads nothing large: only SHA256SUMS
# files and two JSON listings.
set -eu

# sha256 of stdin. shasum (perl Digest::SHA) is absent on Alpine and minimal
# RHEL images, where sha256sum is what exists; tools/make_kycgkb_sums.sh
# already probes for both and this script must not assume differently.
sha256_of() {
    if command -v shasum >/dev/null 2>&1; then
        shasum -a 256 | cut -d' ' -f1
    elif command -v sha256sum >/dev/null 2>&1; then
        sha256sum | cut -d' ' -f1
    else
        echo "make_registry.sh: need shasum or sha256sum" >&2
        exit 1
    fi
}

# python3 builds the per-genome size tables. Without it those tables come back
# empty, registry.h still compiles, and the browser silently shows no download
# sizes -- so check up front rather than degrading quietly.
command -v python3 >/dev/null 2>&1 || {
    echo "make_registry.sh: need python3" >&2
    exit 1
}
command -v curl >/dev/null 2>&1 || {
    echo "make_registry.sh: need curl" >&2
    exit 1
}

# Anything that could not be fetched. The script must not report success with
# NULL anchors in the output: that compiles, and the collection then simply
# cannot be verified at run time.
failures=0
note_failure() {
    failures=$((failures + 1))
    echo "make_registry.sh: FAILED to fetch $1" >&2
}

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT INT TERM

# Fetch a URL to a temp file and print its sha256. Prints nothing and returns
# 1 if the fetch fails.
#
# The digest is taken from the bytes on disk, not from $(curl ...): command
# substitution strips every trailing newline and printf adds exactly one back,
# so a manifest served with zero or two trailing newlines would be pinned to a
# digest kycg never reproduces, and every fetch from that collection would fail
# verification with no clue why.
fetch_to() {
    curl -sfL -o "$2" "$1" 2>/dev/null || return 1
    [ -s "$2" ] || return 1
    return 0
}

out=""
if [ "${1:-}" = "-o" ]; then
    out=${2:?-o needs a path}
    shift 2
fi

tag=${1:-v8}
base=${KYCG_IA_BASE_URL:-https://github.com/zhou-lab/InfiniumAnnotation/raw}

# With -o, write to a temp file and move it into place only on success. The
# documented `make_registry.sh > src/registry.h` truncates the header before
# the script runs, so any mid-run failure leaves a broken file and the previous
# contents recoverable only from git.
if [ -n "$out" ]; then
    exec > "$work/out.h"
fi

# Platforms carrying a KYCG/ directory in InfiniumAnnotation, with the size of
# the row space each one indexes -- the number of probes in its ordering.
#
# These are pinned rather than derived because they are the one fact a user
# needs *before* downloading anything: a .cm is only comparable to a query
# with the same row count, so showing the dimension next to each target is
# what makes it possible to pick the right one. Verified against the published
# .cm files.
#
# platform | rows
platforms="EPIC:866553 EPICv2:937690 HM27:27722 HM450:486427 \
Mammal40:38607 MM285:287692 MSA:284309"

# genome:repo:tag:zenodo_record:doi:rows   (rows = CpGs in cpg_nocontig.cr)
genomes="hg38:KYCGKB_hg38:v2:18175838:10.5281/zenodo.18175837:29401795 \
mm10:KYCGKB_mm10:v2:18175656:10.5281/zenodo.18175655:21867837"

cat <<EOF
/* registry.h -- GENERATED by tools/make_registry.sh. Do not edit.
 *
 * Regenerate with:  tools/make_registry.sh $tag > src/registry.h
 *
 * Arrays: sha256(<PLATFORM>/KYCG/SHA256SUMS) at the pinned tag is the trust
 * anchor; every file digest chains from it. NULL anchor = not published.
 *
 * Set counts are pinned too. A tag is immutable, so the number of sets it
 * publishes is a fixed fact about it -- there is nothing to discover at run
 * time, and the overview can show have/total without touching the network.
 *
 * Whole genome: the same, anchored on sha256(SHA256SUMS) in KYCGKB_<genome>
 * at its own tag. File sizes are pinned alongside for display only; nothing
 * depends on them, so a set added upstream needs no rebuild.
 *
 * Zenodo keeps the DOI and remains the citable archive. It is recorded here as
 * provenance and is no longer the fetch path.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026-present Wanding Zhou
 */
#ifndef _KYCG_REGISTRY_H
#define _KYCG_REGISTRY_H

#include <stdint.h>
#include <stddef.h>

#define KYCG_IA_BASE_URL   "$base"
#define KYCG_IA_TAG        "$tag"
#define KYCG_IA_SUMS_FILE  "SHA256SUMS"
#define KYCG_KB_BASE_URL   "https://github.com/zhou-lab"
#define KYCG_ZENODO_BASE   "https://zenodo.org/records"

/* One array platform's KYCG knowledgebase directory. */
typedef struct {
    const char *platform;
    uint64_t    rows;          /* probes in the platform's ordering */
    uint64_t    n_sets;        /* .cm sets published at this tag */
    const char *sums_sha256;   /* sha256 of <platform>/KYCG/SHA256SUMS */
    /* The probe ordering sits one level up, under the platform's own
     * manifest, and is fetched alongside any set: it is what gives a row
     * index a probe name, so a set without it is a column of anonymous bits. */
    const char *plat_sums_sha256; /* sha256 of <platform>/SHA256SUMS */
} kycg_array_reg_t;

static const kycg_array_reg_t KYCG_ARRAY_REGISTRY[] = {
EOF

for entry in $platforms; do
    plat=$(echo "$entry" | cut -d: -f1)
    nrow=$(echo "$entry" | cut -d: -f2)
    url="$base/$tag/$plat/KYCG/SHA256SUMS"
    # One fetch each, hashed from disk. Fetching twice and hashing the first
    # result while validating the second meant a transient failure on the first
    # request pinned sha256("") -- e3b0c442...b855 -- against a non-empty second
    # body, which compiles fine and permanently breaks verification.
    kf="$work/$plat.kycg.sums"
    pf="$work/$plat.plat.sums"
    if fetch_to "$url" "$kf"; then body=1; else body=""; note_failure "$url"; fi
    if fetch_to "$base/$tag/$plat/SHA256SUMS" "$pf"; then
        psha=$(sha256_of < "$pf")
    else
        psha=""
        note_failure "$base/$tag/$plat/SHA256SUMS"
    fi
    if [ -n "$body" ]; then
        sha=$(sha256_of < "$kf")
        nset=$(grep -c '\.cm$' < "$kf" || true)
        if [ -n "$psha" ]; then
            printf '    { "%s", %s, %s, "%s", "%s" },\n' \
                "$plat" "$nrow" "$nset" "$sha" "$psha"
        else
            printf '    { "%s", %s, %s, "%s", NULL },\n' \
                "$plat" "$nrow" "$nset" "$sha"
        fi
    else
        printf '    { "%s", %s, 0, NULL, NULL },  /* not published at %s */\n' \
            "$plat" "$nrow" "$tag"
    fi
done

cat <<'EOF'
    { NULL, 0, 0, NULL, NULL }
};

/* A published file size, for display only -- never for correctness. */
typedef struct {
    const char *name;
    uint64_t    size;
} kycg_fsize_t;

EOF

# One size table per genome, read from the GitHub tree at the pinned tag.
for g in $genomes; do
    genome=$(echo "$g" | cut -d: -f1)
    repo=$(echo "$g" | cut -d: -f2)
    gtag=$(echo "$g" | cut -d: -f3)
    printf 'static const kycg_fsize_t KYCG_SIZES_%s[] = {\n' "$genome"
    lf="$work/$genome.listing.json"
    if ! fetch_to "https://api.github.com/repos/zhou-lab/$repo/contents/?ref=$gtag" "$lf"; then
        note_failure "listing for $repo@$gtag (sizes will be absent)"
        : > "$lf"
    fi
    python3 -c '
import sys, json
try:
    d = json.load(sys.stdin)
except Exception:
    sys.exit(0)
if not isinstance(d, list):
    sys.exit(0)
for f in sorted(d, key=lambda x: x["name"]):
    if f.get("type") != "file":
        continue
    if f["name"] in ("README.md", "SHA256SUMS"):
        continue
    print("    { \"%s\", %d }," % (f["name"], f["size"]))
' < "$lf"
    printf '    { NULL, 0 }\n};\n\n'
done

cat <<'EOF'
/* One genome's whole-genome knowledgebase collection.
 *
 * `sums_sha256` is the trust anchor, exactly as for arrays. `record` and `doi`
 * point at the Zenodo deposit, which remains the citable archive but is no
 * longer where kycg fetches from. */
typedef struct {
    const char        *genome;
    uint64_t           rows;      /* CpGs in cpg_nocontig.cr */
    uint64_t           n_sets;    /* .cm sets published at this tag */
    const char        *repo;      /* GitHub repository name */
    const char        *tag;       /* pinned tag */
    const char        *sums_sha256;
    const char        *record;    /* Zenodo record id, provenance only */
    const char        *doi;
    const kycg_fsize_t *sizes;
} kycg_seq_reg_t;

static const kycg_seq_reg_t KYCG_SEQ_REGISTRY[] = {
EOF

for g in $genomes; do
    genome=$(echo "$g" | cut -d: -f1)
    repo=$(echo "$g" | cut -d: -f2)
    gtag=$(echo "$g" | cut -d: -f3)
    record=$(echo "$g" | cut -d: -f4)
    doi=$(echo "$g" | cut -d: -f5)
    nrow=$(echo "$g" | cut -d: -f6)
    url="https://github.com/zhou-lab/$repo/raw/$gtag/SHA256SUMS"
    gf="$work/$genome.sums"
    if fetch_to "$url" "$gf"; then body=1; else body=""; note_failure "$url"; fi
    if [ -n "$body" ]; then
        sha=$(sha256_of < "$gf")
        nset=$(grep -c '\.cm$' < "$gf" || true)
        printf '    { "%s", %s, %s, "%s", "%s", "%s", "%s", "%s", KYCG_SIZES_%s },\n' \
            "$genome" "$nrow" "$nset" "$repo" "$gtag" "$sha" "$record" "$doi" "$genome"
    else
        printf '    { "%s", %s, 0, "%s", "%s", NULL, "%s", "%s", KYCG_SIZES_%s },  /* no manifest at %s */\n' \
            "$genome" "$nrow" "$repo" "$gtag" "$record" "$doi" "$genome" "$gtag"
    fi
done

cat <<'EOF'
    { NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL }
};

#endif /* _KYCG_REGISTRY_H */
EOF

# A registry with NULL anchors compiles perfectly and simply cannot verify the
# collections it dropped. Refuse to look successful.
if [ "$failures" -gt 0 ]; then
    echo "make_registry.sh: $failures fetch(es) failed; registry NOT written" >&2
    exit 1
fi

if [ -n "$out" ]; then
    exec >&2
    mv "$work/out.h" "$out"
    echo "make_registry.sh: wrote $out"
fi
