#!/bin/bash
## The command-line surface: dispatch, help, version, and the store line.
## Nothing here reads data -- it is the contract a user meets in the first
## minute, and the part a refactor of main.c breaks silently.
set -uo pipefail
. "$(dirname "$0")/lib.sh"
new_workdir

## ---- dispatch --------------------------------------------------------- ##
out=$("$KYCG" 2>&1); rc=$?
check "bare kycg exits 1"        1 "$rc"
check_has "bare kycg lists fetch"    "fetch"    "$out"
check_has "bare kycg lists test"     "test"     "$out"
check_has "bare kycg lists annotate" "annotate" "$out"

## Usage goes to stderr, so `kycg > file` does not fill the file with help.
out=$("$KYCG" 2>/dev/null); check "usage is not on stdout" "" "$out"

out=$("$KYCG" bogus 2>&1); rc=$?
check "unknown command exits 1" 1 "$rc"
check_has "unknown command is named" "bogus" "$out"
## and the message carries no internal source location.
check_lacks "unknown command has no [func] prefix" "[main]" "$out"

## a bad OPTION (not a bad command) is an error too, and exits 1 -- the
## distinction from -h that makes help safe in a set -e script.
out=$("$KYCG" test -Z </dev/null 2>&1); rc=$?
check "an unknown option exits 1" 1 "$rc"

## -h and --help ask for help, which is success: exit 0, so `kycg -h` inside a
## `set -e` script does not kill it. (A usage ERROR -- no args, a bad option --
## still exits 1; those are checked above and below.)
h=$("$KYCG" -h 2>&1); rc=$?
check "-h exits 0" 0 "$rc"
check "--help exits 0" 0 "$("$KYCG" --help >/dev/null 2>&1; echo $?)"
check "--help matches -h" "$h" "$("$KYCG" --help 2>&1)"

for sub in test annotate fetch; do
  out=$("$KYCG" "$sub" -h 2>&1); rc=$?
  check "$sub -h exits 0" 0 "$rc"
  check_has "$sub -h names itself" "kycg $sub" "$out"
  check_has "$sub -h has an Options block" "Options" "$out"
done

## A subcommand with no arguments must say what is missing, not wait: the
## suite runs with stdin redirected, which is exactly the Nextflow/Docker
## case the browser is gated on.
out=$("$KYCG" test < /dev/null 2>&1); rc=$?
check "test with no args exits 1" 1 "$rc"
out=$("$KYCG" annotate < /dev/null 2>&1); rc=$?
check "annotate with no args exits 1" 1 "$rc"

## ---- --version -------------------------------------------------------- ##
## The one place that reports what this build is coupled to. The conda recipe
## greps it for the YAME version; check the shape it greps for.
v=$("$KYCG" --version 2>&1); rc=$?
check "--version exits 0" 0 "$rc"
check_has "--version names kycg" "kycg" "$v"
if ! printf '%s\n' "$v" | grep -qE "YAME[[:space:]]+v[0-9]+\.[0-9]+"; then
  echo "  FAIL --version reports a YAME version"
  echo "       actual: $v"; fails=$((fails + 1))
fi
## The store line must follow $YAME_DATA_HOME -- one shared store across the
## tools is the whole point of the libyame asset system.
v=$(YAME_DATA_HOME="$PWD/store_here" "$KYCG" --version 2>&1)
check_has "--version honours YAME_DATA_HOME" "$PWD/store_here" "$v"
## --version reports what it can verify and whether it can fetch (the README
## says it does): the registry tag and libcurl availability.
v=$("$KYCG" --version 2>&1)
check_has "--version reports the registry tag" "InfiniumAnnotation" "$v"
check_has "--version reports network support" "libcurl" "$v"
## -v is an alias for --version.
check "-v matches --version" "$v" "$("$KYCG" -v 2>&1)"

## ---- colour ----------------------------------------------------------- ##
## The help text is styled, but only on a terminal. Off one -- which is how
## every pipe, log and CI job sees it -- there must be no escape sequences.
esc=$(printf '\033')
out=$("$KYCG" 2>&1)
check_lacks "no ANSI escapes off a terminal" "$esc" "$out"

done_testing
