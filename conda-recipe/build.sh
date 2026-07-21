#!/bin/bash
set -euo pipefail

# `source: path: ..` copies the working tree as it is on disk, including build
# products that are not tracked by git. external/YAME/libyame.a and its nested
# htslib/libhts.a are exactly that: built for the author's machine, and linking
# them into a package for another architecture would produce a binary that
# fails at load time rather than at build time. Drop them, and the stale .o
# they were made from, so they rebuild for THIS target.
#
# kycg's `make clean` deliberately does not reach into the submodule (that is
# `distclean`, which needs the submodule's own Makefile), so clean it by hand.
make clean >/dev/null 2>&1 || true
rm -f external/YAME/libyame.a external/YAME/htslib/libhts.a
find external/YAME -name '*.o' -delete 2>/dev/null || true

# htslib generates version.h itself and its own `clean` removes it, so there is
# nothing to stage here -- unlike config.h, which is committed upstream and
# means no autoconf run is needed at build time.

# conda supplies CC plus the sysroot / $PREFIX include+lib flags in CPPFLAGS /
# CFLAGS / LDFLAGS. The vendored YAME and htslib Makefiles set rigid CFLAGS and
# ignore CPPFLAGS, so bake conda's *compile* flags into CC -- every level of
# sub-make (kycg -> YAME -> htslib) uses CC verbatim, and a command-line CC
# propagates to all sub-makes. LDFLAGS is expanded explicitly by kycg's link
# rule -- make imports the environment into a variable, but an explicit recipe
# never expands it by itself, and conda-forge puts -L$PREFIX/lib only in
# LDFLAGS, so without that a linux-64 build fails to find -lz.
#
# CURL is left at its default of `auto`: curl-config comes from the libcurl
# host dependency, so detection succeeds and the package gets network support.
# If that ever silently failed the result would be a kycg that builds, tests
# clean, and cannot fetch -- which is why meta.yaml's test greps for it.
make -j"${CPU_COUNT:-1}" CC="${CC} ${CPPFLAGS} ${CFLAGS}"

make install PREFIX="${PREFIX}"
