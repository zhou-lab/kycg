#!/bin/bash
## The failure modes of `kycg test`. Everything is positional -- row i must
## mean the same CpG on both sides -- so the errors are the feature: each one
## has to name what went wrong and exit non-zero rather than produce a table
## of nonsense.
set -uo pipefail
. "$(dirname "$0")/lib.sh"
new_workdir

awk 'BEGIN{for(i=0;i<5;i++) print i%2}'  > q5.txt
awk 'BEGIN{for(i=0;i<7;i++) print (i%2)?"A":"B"}' > s7.txt
awk 'BEGIN{for(i=0;i<5;i++) print (i%2)?"A":"B"}' > s5.txt
pack_binary q5.txt query5.cg
pack_state  s7.txt  kb7.cm
pack_state  s5.txt  kb5.cm

## ---- the row-count mismatch -------------------------------------------- ##
## The one check kycg can make about row spaces, and the message users will
## actually meet. It must name BOTH files and BOTH counts: "row count
## mismatch" alone leaves them guessing which side to re-pack.
out=$("$KYCG" test -m kb7.cm query5.cg 2>&1 >/dev/null); rc=$?
check "mismatch exits non-zero" 1 "$rc"
check_has "mismatch names the query"  "query5.cg" "$out"
check_has "mismatch names the kb"     "kb7.cm"    "$out"
check_has "mismatch gives the query count" "5"    "$out"
check_has "mismatch gives the kb count"    "7"    "$out"
check_has "mismatch explains why" "different reference row lists" "$out"
## a user error, not an internal one: no [func:line] prefix leaks through.
check_lacks "mismatch has no [func] prefix" "[run_pair" "$out"
check_lacks "mismatch has no [ prefix at all" "[" "$out"

## Matching counts on the same fixture must then succeed -- otherwise the
## check above proves nothing about the mismatch specifically.
"$KYCG" test -m kb5.cm query5.cg > ok.tsv 2>/dev/null
check "matching row counts succeed" 0 "$?"
check "and produce a row per set" 2 "$(tail -n +2 ok.tsv | wc -l | tr -d ' ')"

## ---- inputs that are not there ----------------------------------------- ##
## A silent exit 1 is the worst answer: the user cannot tell a typo from an
## empty result. Every one of these must say something.
out=$("$KYCG" test -m no_such.cm query5.cg 2>&1 >/dev/null); rc=$?
check "missing knowledgebase exits non-zero" 1 "$rc"
check_has "missing knowledgebase is named"   "no_such.cm" "$out"
[ -n "$out" ] || { echo "  FAIL missing knowledgebase failed silently"
                   fails=$((fails + 1)); }

out=$("$KYCG" test -m kb5.cm no_such.cg 2>&1 >/dev/null); rc=$?
check "missing query exits non-zero" 1 "$rc"
check_has "missing query is named" "no_such.cg" "$out"

## A store name that matches nothing, off a terminal: an error with a way
## forward, never a wait for an answer nobody can give.
out=$(YAME_DATA_HOME="$PWD/empty_store" "$KYCG" test -m hg38:NoSuchSet query5.cg \
        < /dev/null 2>&1 >/dev/null); rc=$?
check "unknown set exits non-zero" 1 "$rc"
check_has "unknown set is named"   "NoSuchSet" "$out"
check_has "unknown set points at fetch" "kycg fetch" "$out"

## ---- input that is there but is not a store ---------------------------- ##
## A truncated or foreign file must be refused, not read as whatever the
## bytes happen to look like.
head -c 200 /dev/urandom > junk.cm
out=$("$KYCG" test -m junk.cm query5.cg 2>&1 >/dev/null); rc=$?
awk -v rc="$rc" 'BEGIN{exit !(rc != 0)}' ||
  { echo "  FAIL a junk knowledgebase was accepted (rc=$rc)"; fails=$((fails + 1)); }

: > empty.cm
out=$("$KYCG" test -m empty.cm query5.cg 2>&1 >/dev/null); rc=$?
awk -v rc="$rc" 'BEGIN{exit !(rc != 0)}' ||
  { echo "  FAIL an empty knowledgebase was accepted (rc=$rc)"; fails=$((fails + 1)); }

done_testing
