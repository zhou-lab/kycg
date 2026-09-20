## Shared helpers for the tests/t_*.sh command-level tests.
## Sourced, never run. Every test sets up its own temp directory and fixtures;
## this file only carries the assertions and the two or three pack recipes
## that would otherwise be copied into every script.
##
## The assertions print what was expected and what came back, because a test
## that only says "failed" costs a rerun to understand.

: "${KYCG:?export KYCG=/path/to/kycg}"
: "${YAME:?export YAME=/path/to/yame}"

fails=0

## A scratch directory for the caller, removed on exit. Tests cd into it, so
## every fixture path below is relative and nothing reaches outside.
new_workdir() {
  d=$(mktemp -d)
  trap 'rm -rf "$d"' EXIT
  cd "$d" || exit 1
}

## check <label> <expected> <actual>
check() {
  if [ "$2" = "$3" ]; then return 0; fi
  echo "  FAIL $1"
  echo "       expected: $2"
  echo "       actual:   $3"
  fails=$((fails + 1))
}

## check_close <label> <expected> <actual> <tolerance>  -- floating point
check_close() {
  if awk -v a="$2" -v b="$3" -v t="$4" \
      'BEGIN{d=a-b; if(d<0)d=-d; exit !(d<=t)}'; then return 0; fi
  echo "  FAIL $1"
  echo "       expected: $2 (+/- $4)"
  echo "       actual:   $3"
  fails=$((fails + 1))
}

## check_has <label> <needle> <haystack>  -- substring, for error messages
check_has() {
  case "$3" in
    *"$2"*) return 0;;
  esac
  echo "  FAIL $1"
  echo "       expected to contain: $2"
  echo "       actual:              $3"
  fails=$((fails + 1))
}

## check_lacks <label> <needle> <haystack>
check_lacks() {
  case "$3" in
    *"$2"*)
      echo "  FAIL $1"
      echo "       expected NOT to contain: $2"
      echo "       actual:                  $3"
      fails=$((fails + 1));;
  esac
}

## done_testing -- the exit status every test ends with
done_testing() {
  if [ "$fails" -eq 0 ]; then echo "  all checks passed"; exit 0; fi
  echo "  $fails check(s) failed"
  exit 1
}

## ---- fixture recipes -------------------------------------------------- ##
## kycg reads what yame writes, so the fixtures are packed rather than
## committed: a committed .cg would silently become a format the current YAME
## no longer writes, and the test would then be checking a museum piece.

## pack_binary <text-file> <out.cg>   one 0/1 per line: a CpG mask
pack_binary() { "$YAME" pack -f b "$1" "$2"; }

## pack_state <text-file> <out.cm>    one label per line: the knowledgebase.
## A state .cm carries every label as its own set, which is how one file gives
## `kycg test` many rows of output (chromhmm.cm works exactly this way).
pack_state() { "$YAME" pack -f s "$1" "$2"; }

## pack_mu <text-file> <out.cg>       two columns M and U: a sequencing query.
## This is what gives `kycg test` its beta and depth columns.
pack_mu() { "$YAME" pack -f m "$1" "$2"; }

## ---- driving the interactive widgets --------------------------------- ##
## The browser and the pickers only run when stdout is a terminal, so they are
## invisible to every test that redirects. script(1) gives us one.

## have_pty -- 0 if the terminal cases can be tested here
have_pty() {
  command -v script >/dev/null || return 1
  ## script(1) exists in some conda images but cannot allocate a pty there.
  script -qec true /dev/null >/dev/null 2>&1 || return 1
  return 0
}

## pty_drive <command> <keys...> -- run <command> on a pty and type at it.
## Each key argument is one group, sent a beat after the last, so a redraw
## happens in between (redraw paths are most of ui.c). The first key waits
## longer: anything sent before the program switches the terminal to raw mode
## is echoed and lost, and an opening screen may be reading a manifest.
##
## Never end with a bare ESC: the key parser reads one byte and waits to see
## whether more follow (an arrow is ESC [ A), so a trailing ESC blocks on a
## read that never returns. Send a q after it.
pty_drive() {
  _cmd=$1; shift
  ( sleep "${PTY_SETTLE:-1.5}"
    for _k in "$@"; do printf '%b' "$_k"; sleep "${PTY_BEAT:-0.4}"; done
    sleep "${PTY_TAIL:-1}" ) |
    timeout "${PTY_TIMEOUT:-60}" script -qec "$_cmd" /dev/null 2>&1
}

## The escape sequences the widgets parse.
PTY_ESC=$(printf '')
PTY_UP="$PTY_ESC[A";   PTY_DOWN="$PTY_ESC[B"
PTY_RIGHT="$PTY_ESC[C"; PTY_LEFT="$PTY_ESC[D"
PTY_HOME="$PTY_ESC[H"; PTY_END="$PTY_ESC[F"
PTY_PGUP="$PTY_ESC[5~"; PTY_PGDN="$PTY_ESC[6~"

## col <name> <tsv-file>  -- 1-based index of a named column, or empty
col() { awk -F'\t' -v n="$1" 'NR==1{for(i=1;i<=NF;i++) if($i==n){print i; exit}}' "$2"; }

## field <column-name> <db-name> <tsv-file>  -- one cell, by column and db
field() {
  awk -F'\t' -v n="$1" -v db="$2" '
    NR==1{for(i=1;i<=NF;i++){if($i==n)c=i; if($i=="db")d=i}; next}
    $d==db{print $c; exit}' "$3"
}
