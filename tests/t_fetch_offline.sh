#!/bin/bash
## `kycg fetch` without a network. Downloads happen in fetch and nowhere
## else, and every question is gated on an interactive terminal -- so off one,
## with an empty store, fetch must print its catalogue, touch nothing, and
## exit. That is the property a Nextflow job or a Docker build depends on.
##
## Nothing here reaches the network: the catalogue is compiled into the
## binary (the registry), and no target is named.
set -uo pipefail
. "$(dirname "$0")/lib.sh"
new_workdir

store=$PWD/store
mkdir -p "$store"
export YAME_DATA_HOME="$store"

## ---- the TSV catalogue -------------------------------------------------- ##
## Redirected stdout is not a terminal, so the browser must not open.
out=$("$KYCG" fetch < /dev/null 2>err.txt); rc=$?
check "fetch off a terminal exits 0" 0 "$rc"
check "fetch says nothing on stderr" "" "$(cat err.txt)"

## A comment line naming the store, then a header, then one row per target.
check_has "the store is stated" "$store" "$(printf '%s\n' "$out" | head -1)"
hdr=$(printf '%s\n' "$out" | sed -n '2p')
check "the catalogue header" "target	kind	rows	source	cached_sets" "$hdr"

## Every genome and array the build pins has to appear -- this is the list a
## user picks from, and it comes from the compiled registry.
for t in hg38 mm10 mm39 EPIC EPICv2 HM27 HM450 Mammal40 MM285 MSA; do
  check_has "the catalogue offers $t" "$t" "$out"
done

## Five columns on every row, so `kycg fetch | cut -f1` (the documented
## scripting use) cannot silently pick up a ragged line.
n_bad=$(printf '%s\n' "$out" | awk -F'\t' 'NR>2 && NF!=5' | wc -l | tr -d ' ')
check "every catalogue row has 5 fields" 0 "$n_bad"

## An empty store means nothing is cached, and the count must say so rather
## than report what upstream publishes.
n_nonzero=$(printf '%s\n' "$out" | awk -F'\t' 'NR>2 && $5 !~ /^0\//' | wc -l | tr -d ' ')
check "an empty store caches nothing" 0 "$n_nonzero"

## ---- nothing was downloaded -------------------------------------------- ##
check "the store is still empty" 0 "$(find "$store" -type f | wc -l | tr -d ' ')"

## ---- -d overrides the store -------------------------------------------- ##
other=$PWD/elsewhere
mkdir -p "$other"
out=$("$KYCG" fetch -d "$other" < /dev/null 2>/dev/null)
check_has "-d overrides YAME_DATA_HOME" "$other" "$(printf '%s\n' "$out" | head -1)"

## ---- -o narrows, and an unknown target is an error ---------------------- ##
out=$("$KYCG" fetch NoSuchTarget < /dev/null 2>&1 >/dev/null); rc=$?
awk -v rc="$rc" 'BEGIN{exit !(rc != 0)}' ||
  { echo "  FAIL an unknown fetch target was accepted (rc=$rc)"; fails=$((fails+1)); }
check_has "the unknown target is named" "NoSuchTarget" "$out"
check "still nothing downloaded" 0 "$(find "$store" -type f | wc -l | tr -d ' ')"

## ---- -o names sets, not files, and never an empty string --------------- ##
## The empty case is the unset-variable foot-gun: `-o "$SETS"` with SETS unset
## used to fetch the whole collection. Both are refused before any network.
out=$("$KYCG" fetch -f -o "" mm10 2>&1 >/dev/null); rc=$?
check "an empty -o exits non-zero" 1 "$rc"
check_has "and says -o needs set names" "set names" "$out"
out=$("$KYCG" fetch -f -o PMD.20220911.cm mm10 2>&1 >/dev/null); rc=$?
check "a filename in -o exits non-zero" 1 "$rc"
check_has "and names the confusion" "not a file name" "$out"
out=$("$KYCG" fetch -f -o "a/b" mm10 2>&1 >/dev/null); rc=$?
check "a path in -o exits non-zero" 1 "$rc"
check "the store is still empty after rejected -o" 0 \
      "$(find "$store" -type f | wc -l | tr -d ' ')"

done_testing
