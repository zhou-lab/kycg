#!/bin/sh
## Fail if docs/index.html has drifted from the version the build actually
## pins. The docs quote two versions in prose -- kycg's own and the coupled
## YAME -- and prose does not recompile, so it silently goes stale on a bump.
## This makes that a build failure instead: run it in CI (and before a release)
## so "forgot to update the docs" is caught, not shipped.
##
## Sources of truth (the same macros the binary prints via `kycg --version`):
##   src/kycg.h                    KYCG_VERSION  e.g. 0.1.0
##   external/YAME/src/yame_version.h  YAME_VERSION  e.g. v1.29
##
## Both must appear verbatim somewhere in docs/index.html. The check is a
## substring match, so it is agnostic to the surrounding HTML -- move or restyle
## the version however you like, it only has to be the right string.
set -eu

here=$(dirname "$0")
root=$(cd "$here/.." && pwd)
kycg_h="$root/src/kycg.h"
yame_h="$root/external/YAME/src/yame_version.h"
doc="$root/docs/index.html"

## Pull the quoted value out of a `#define NAME "value"` line.
define_of() {
    sed -n "s/.*#define $1 \"\\([^\"]*\\)\".*/\\1/p" "$2" | head -n1
}

[ -f "$yame_h" ] || {
    echo "check_docs_versions: $yame_h missing -- is the YAME submodule checked out?" >&2
    exit 2
}

kycg_ver=$(define_of KYCG_VERSION "$kycg_h")
yame_ver=$(define_of YAME_VERSION "$yame_h")

[ -n "$kycg_ver" ] && [ -n "$yame_ver" ] || {
    echo "check_docs_versions: could not read KYCG_VERSION / YAME_VERSION." >&2
    exit 2
}

fail=0
## -F: the versions carry '.' which is a regex metacharacter; match literally.
grep -qF "$kycg_ver" "$doc" || {
    echo "check_docs_versions: docs/index.html does not mention kycg $kycg_ver" >&2
    fail=1
}
grep -qF "$yame_ver" "$doc" || {
    echo "check_docs_versions: docs/index.html does not mention YAME $yame_ver" >&2
    fail=1
}

if [ "$fail" -ne 0 ]; then
    echo "check_docs_versions: docs/index.html is stale. Update it to kycg" \
         "$kycg_ver / YAME $yame_ver (the versions kycg --version prints)." >&2
    exit 1
fi

echo "check_docs_versions: docs/index.html in sync (kycg $kycg_ver, YAME $yame_ver)."
