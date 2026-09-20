#!/bin/bash
## The docs against the binary. Every claim checked here is one a reader acts
## on in their first five minutes, and each has been wrong at least once:
## a version that drifted from the recipe, an install line that produced a
## tree without the `yame` the very next section tells you to run.
##
## Skips cleanly where the sources are not present (a conda test environment
## copies only what meta.yaml's source_files names).
set -uo pipefail
. "$(dirname "$0")/lib.sh"
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)

## ---- one version, in four places ---------------------------------------- ##
v_bin=$("$KYCG" --version 2>&1 | head -1 | awk '{print $2}')
[ -n "$v_bin" ] || { echo "  FAIL --version printed no version"; exit 1; }

if [ -f "$root/src/kycg.h" ]; then
  v_hdr=$(sed -n 's/.*define KYCG_VERSION "\([^"]*\)".*/\1/p' "$root/src/kycg.h")
  check "src/kycg.h matches the binary" "$v_bin" "$v_hdr"
fi
if [ -f "$root/conda-recipe/meta.yaml" ]; then
  ## The recipe's version is what the package is published as; a drift here
  ## ships a 0.4 binary labelled 0.3 and there is no way to tell afterwards.
  v_rec=$(sed -n 's/^{% set version = "\([^"]*\)".*/\1/p' "$root/conda-recipe/meta.yaml")
  check "conda-recipe/meta.yaml matches the binary" "$v_bin" "$v_rec"
fi

## ---- every command the docs name is dispatchable ------------------------ ##
## The top-level usage is the list; each entry must actually be a subcommand.
for sub in $("$KYCG" 2>&1 | awk '/^Commands/{f=1; next} f && NF==0{exit} f{print $1}'); do
  out=$("$KYCG" "$sub" -h 2>&1)
  check_has "usage lists $sub, and it dispatches" "kycg $sub" "$out"
done

## kycg has no `info` subcommand -- the row-space check is `yame info`. If
## that ever changes, the docs below change with it.
out=$("$KYCG" info 2>&1 >/dev/null); rc=$?
check "there is no kycg info" 1 "$rc"

## ---- install instructions produce what the next section uses ------------ ##
## The docs tell readers to run `yame info` to check a row space. kycg links
## libyame but does not ship the yame COMMAND, so every install path must
## name it. This is the exact break that made the row-space section exit 127.
for f in "$root/README.md" "$root/docs/llms.txt" "$root/docs/index.html"; do
  [ -f "$f" ] || continue
  n=$(basename "$f")
  if grep -q "yame info" "$f"; then
    grep -q "conda install .*kycg yame" "$f" ||
      { echo "  FAIL $n tells the reader to run 'yame info' but never installs yame"
        fails=$((fails + 1)); }
  fi
done

## The source build has to be able to produce it too.
if [ -f "$root/Makefile" ]; then
  grep -q '^yame-bin:' "$root/Makefile" ||
    { echo "  FAIL the Makefile has no yame-bin target"; fails=$((fails + 1)); }
fi

## ---- the generator is kycg's own, and the docs say so ------------------ ##
## This guard is the inverse of the one it replaces. The emitter used to live
## in YAME (`--tool=kycg`), and the docs had to name the external/YAME path;
## since YAME v1.50 each tool emits its own registry, so kycg ships
## tools/make_registry.sh and a doc pointing into the submodule would send the
## reader to a script that no longer exists there.
if [ -d "$root/tools" ]; then
  [ -f "$root/tools/make_registry.sh" ] ||
    { echo "  FAIL kycg has no tools/make_registry.sh; it emits its own registry now"
      fails=$((fails + 1)); }
fi
for f in "$root/README.md" "$root/docs/llms.txt"; do
  [ -f "$f" ] || continue
  n=$(basename "$f")
  grep -q "external/YAME/tools/make_registry.sh" "$f" &&
    { echo "  FAIL $n still points make_registry.sh into the submodule"
      fails=$((fails + 1)); }
done

## and the generated headers really are regenerable from the checked-in table.
## Only in a tree that has the generators AND the catalog they read: a conda
## test environment gets tests/ and docs/ as source_files, not tools/ and not
## the YAME submodule, so there is nothing to regenerate from there.
if [ -x "$root/tools/make_registry.sh" ] &&
   [ -f "$root/external/YAME/tools/registry/lib.sh" ]; then
  ( cd "$root" && ./tools/make_registry.sh --check >/dev/null 2>&1 ) ||
    { echo "  FAIL src/registry.h is stale; run tools/make_registry.sh -o src/registry.h"
      fails=$((fails + 1)); }
  ( cd "$root" && ./tools/make_kbinfo.sh --check >/dev/null 2>&1 ) ||
    { echo "  FAIL src/kbinfo.h is stale; run tools/make_kbinfo.sh -o src/kbinfo.h"
      fails=$((fails + 1)); }
fi

## check_dimensions.sh reads the row count with `yame info`; kycg has none, so
## a `kycg info` (or `./kycg info`) call there prints MISMATCH for every row.
if [ -f "$root/tools/check_dimensions.sh" ]; then
  grep -Eq '(\./)?kycg["'"'"' ]+info|kycg", "info"' "$root/tools/check_dimensions.sh" &&
    { echo "  FAIL check_dimensions.sh calls a nonexistent 'kycg info'"
      fails=$((fails + 1)); }
fi

## ---- every code block on the page declares what it is ------------------ ##
## A block is either runnable -- data-example names the docs/examples script
## the docs gate runs (`make test-docs`) -- or it is an illustration and
## carries data-norun with the reason. A block that declares neither is a
## command nothing ever runs, which is how six documented lines on the
## published page came to fail as printed.
if [ -f "$root/docs/index.html" ] && [ -d "$root/docs/examples" ]; then
  n_pre=$(grep -o '<pre' "$root/docs/index.html" | wc -l | tr -d ' ')
  n_dec=$(grep -oE '<pre[^>]*data-(example|norun)=' "$root/docs/index.html" |
          wc -l | tr -d ' ')
  check "every <pre> declares example-or-norun" "$n_pre" "$n_dec"

  ## and every named script is really there
  for ex in $(grep -oE 'data-example="[^"]+"' "$root/docs/index.html" |
              sed 's/.*="//; s/"//'); do
    [ -f "$root/docs/examples/$ex.sh" ] ||
      { echo "  FAIL the page names docs/examples/$ex.sh, which does not exist"
        fails=$((fails + 1)); }
  done

  ## every example carries a title, and a skipped one says why
  for f in "$root"/docs/examples/*.sh; do
    b=$(basename "$f")
    grep -q '^## title:' "$f" ||
      { echo "  FAIL $b has no '## title:' line"; fails=$((fails + 1)); }
    bash -n "$f" ||
      { echo "  FAIL $b is not valid bash"; fails=$((fails + 1)); }
  done
fi

## ---- the collection count in the prose matches the binary -------------- ##
## The page and llms.txt state a number ("Ten collections", "Nine collections")
## that a reader counts on; a genome added upstream (mm39 was) makes it wrong.
ncoll=$(YAME_DATA_HOME="$(mktemp -d)" "$KYCG" fetch </dev/null 2>/dev/null |
        awk -F'\t' 'NR>1 && $2 ~ /genome|array/' | wc -l | tr -d ' ')
words="zero one two three four five six seven eight nine ten eleven twelve"
want=$(echo "$words" | cut -d' ' -f$((ncoll + 1)))
for f in "$root/docs/llms.txt" "$root/docs/index.html"; do
  [ -f "$f" ] || continue
  n=$(basename "$f")
  if grep -qiE "(nine|ten|eleven) collection" "$f"; then
    grep -qi "$want collection" "$f" ||
      { echo "  FAIL $n miscounts collections: binary lists $ncoll ($want)"
        fails=$((fails + 1)); }
  fi
done

## ---- the documented result columns are the ones emitted ----------------- ##
## llms.txt lists the output columns for agents reading it instead of the
## man page; a renamed column there is a silently broken downstream script.
if [ -f "$root/docs/llms.txt" ]; then
  d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
  awk 'BEGIN{for(i=0;i<40;i++) print i%3==0}'      > "$d/q.txt"
  awk 'BEGIN{for(i=0;i<40;i++) print (i<20)?"A":"B"}' > "$d/s.txt"
  pack_binary "$d/q.txt" "$d/q.cg"
  pack_state  "$d/s.txt" "$d/kb.cm"
  hdr=$("$KYCG" test -m "$d/kb.cm" "$d/q.cg" | head -1)
  for c in log10_p neglog10_fdr estimate overlap cf_jaccard beta depth; do
    case "$hdr" in *"$c"*) ;; *)
      echo "  FAIL the output has no $c column"; fails=$((fails + 1));; esac
    grep -q -- "$c" "$root/docs/llms.txt" ||
      { echo "  FAIL llms.txt never mentions the $c column"; fails=$((fails + 1)); }
  done
fi

done_testing
