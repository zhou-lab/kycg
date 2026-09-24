#!/bin/bash
## The interactive browser, driven through a pty.
##
## `kycg fetch` with no redirect opens a full-screen tree, and `kycg test` /
## `kycg annotate` with no -m open a picker. None of it is reachable by a test
## that redirects -- which is every other test here -- so this one runs the
## binary under script(1) and types at it. Without this the whole of ui.c and
## the browser half of fetch.c are dark.
##
## Nothing is downloaded: no key sequence below ends in `f`.
set -uo pipefail
. "$(dirname "$0")/lib.sh"

have_pty || { echo "  skip: no usable pty here"; exit 0; }

new_workdir
store=$PWD/store
mkdir -p "$store"
export YAME_DATA_HOME="$store"

## pty_drive and the key names come from lib.sh.
drive() { pty_drive "$@"; }
ESC=$PTY_ESC
UP=$PTY_UP; DOWN=$PTY_DOWN; RIGHT=$PTY_RIGHT; LEFT=$PTY_LEFT
HOME_=$PTY_HOME; END_=$PTY_END; PGDN=$PTY_PGDN; PGUP=$PTY_PGUP

## ---- 1. the tree: open, move, unfold, check, quit ----------------------- ##
out=$(drive "$KYCG fetch hg38" "$DOWN$DOWN" "$UP" 'jjk' "$RIGHT" ' ' "$LEFT" 'q'); rc=$?
check "the browser exits 0 on q" 0 "$rc"
check_has "the tree names the store" "store:" "$out"
check_has "the tree lists the target" "hg38" "$out"
check_has "the tree shows a set" "ABCompartment" "$out"
## Opening with a target checked is the documented behaviour of `fetch <tgt>`.
check_has "the tree shows the set counts" "0/32" "$out"
check "nothing was downloaded" 0 "$(find "$store" -type f | wc -l | tr -d ' ')"

## ---- 2. every movement key, and the jumps ------------------------------- ##
## Home/End/PgUp/PgDn each have their own branch; a tree with one collection
## unfolded is long enough for the paging ones to do something.
out=$(drive "$KYCG fetch" "$RIGHT" "$END_" "$PGUP" "$PGDN" "$HOME_" "$DOWN" 'q'); rc=$?
check "movement keys leave the browser usable" 0 "$rc"
check_has "the whole catalogue is offered" "mm10" "$out"
check_has "arrays are offered too" "EPIC" "$out"

## ---- 3. the detail pane (i) and the recommend key (r) ------------------- ##
out=$(drive "$KYCG fetch hg38" "$RIGHT" "$DOWN" 'i' "$DOWN" 'i' 'q'); rc=$?
check "the detail pane opens and closes" 0 "$rc"
out=$(drive "$KYCG fetch hg38" 'r' "$ESC" 'q'); rc=$?
check "recommend then ESC exits cleanly" 0 "$rc"

## ---- 4. ESC leaves the browser without fetching ------------------------- ##
## A trailing q after the ESC: the key parser reads one byte and then waits
## to see whether more follow (an arrow key is ESC [ A), so a bare ESC as the
## LAST thing on the wire leaves it blocked on a read that never returns.
## Real terminals have the same ambiguity; here it would just hang the suite.
out=$(drive "$KYCG fetch hg38" "$RIGHT" ' ' "$ESC" 'q'); rc=$?
check "ESC leaves the browser" 0 "$rc"
check "ESC downloads nothing" 0 "$(find "$store" -type f | wc -l | tr -d ' ')"

## ---- 5. the picker: `kycg test` with no -m on a terminal ---------------- ##
## A different widget from the tree (a flat, filterable list), and the only
## way a user reaches it.
awk 'BEGIN{for(i=0;i<40;i++) print i%2}' > q.txt
pack_binary q.txt query.cg
out=$(drive "$KYCG test query.cg" 'j' 'q'); rc=$?
awk -v rc="$rc" 'BEGIN{exit !(rc != 0)}' ||
  { echo "  FAIL quitting the picker still ran a test (rc=$rc)"; fails=$((fails+1)); }
check_has "quitting the picker says nothing was selected" "No knowledgebase" "$out"

## `kycg annotate` offers the same picker, filtered to array platforms: a
## probe ID has no meaning without an ordering, so a genome cannot be offered.
printf 'Probe_ID\ncg00000001\n' > p.tsv
out=$(drive "$KYCG annotate p.tsv" 'q'); rc=$?
awk -v rc="$rc" 'BEGIN{exit !(rc != 0)}' ||
  { echo "  FAIL quitting the annotate picker still ran (rc=$rc)"; fails=$((fails+1)); }
check_lacks "the annotate picker offers no genome" "mm39" "$out"

## ---- 6. a terminal that cannot host the browser ------------------------- ##
## TERM=dumb is a terminal by isatty and not one by capability. The browser
## must decline in words rather than paint escape codes into a log.
out=$(TERM=dumb drive "$KYCG fetch hg38" 'q'); rc=$?
check_lacks "TERM=dumb draws no alternate screen" "$ESC[?1049h" "$out"

## ---- 7. NO_COLOR is not a knob ------------------------------------------ ##
## There was a block here asserting that NO_COLOR strips the colour and leaves
## the widget. It went when the UI became YAME's: the browser is coloured
## whenever it draws, and a reader who wants plain text has `-l`, `-y` and a
## pipe, all of which bypass the widget entirely. Nothing in kycg reads
## NO_COLOR any more, so there is nothing here to assert -- TERM=dumb and a
## redirected run, above and throughout, are the fallbacks that still hold.

done_testing
