#!/bin/bash
## `kycg test -p`: a query given as a list of probe IDs. The list is mapped
## through the platform ordering into the record `yame pack -f b` would write,
## so the one thing to prove is that both routes give the same table, plus the
## handling of the untidy lists people actually export.
set -uo pipefail
. "$(dirname "$0")/lib.sh"
new_workdir

## HM27's row count is pinned in the compiled registry, and -p refuses an
## ordering that disagrees, so the fixture has exactly this many rows.
N=27722
PLAT=HM27

store=$PWD/store
mkdir -p "$store/$PLAT"
{ echo "Probe_ID"
  awk -v n="$N" 'BEGIN{for(i=1;i<=n;i++) printf "cg%08d\n", i}'; } |
  gzip -c > "$store/$PLAT/$PLAT.ordering.tsv.gz"
export YAME_DATA_HOME="$store"

## Knowledgebase: rows 1-1000 CGI, 1001-2000 Shore, the rest OpenSea.
awk -v n="$N" 'BEGIN{for(i=1;i<=n;i++)
  print (i<=1000)?"CGI":((i<=2000)?"Shore":"OpenSea")}' > states.txt
pack_state states.txt cgi.cm

## Query: every 10th of the first 3000 rows, so 100 in CGI, 100 in Shore and
## 100 in OpenSea.
awk 'BEGIN{for(i=10;i<=3000;i+=10) printf "cg%08d\n", i}' > ids.txt
awk -v n="$N" 'BEGIN{for(i=1;i<=n;i++) print (i<=3000 && i%10==0)?1:0}' > bits.txt
pack_binary bits.txt q.cg

## ---- same table as the packed query ------------------------------------- ##
"$KYCG" test -m cgi.cm q.cg | cut -f2- > packed.tsv || exit 1
"$KYCG" test -p "$PLAT" -m cgi.cm ids.txt > list.tsv 2> err.txt || {
  echo "  FAIL test -p exited non-zero"; sed 's/^/       /' err.txt; exit 1; }
check "-p gives the packed query's table" "$(cat packed.tsv)" "$(cut -f2- list.tsv)"
check "query_file is the list" "ids.txt" "$(awk -F'\t' 'NR==2{print $1}' list.tsv)"
check "a clean list says nothing on stderr" "" "$(cat err.txt)"

## ---- an untidy list: header, comment, extra columns, duplicates, blanks,
## an ID the platform lacks, CRLF line ends ---------------------------------- ##
{ printf 'Probe_ID,logFC\r\n# exported from limma\r\n'
  sed 's/$/,0.5\r/' ids.txt
  head -20 ids.txt
  printf '\r\ncg99999999,0.1\r\n'; } > messy.csv
"$KYCG" test -p "$PLAT" -m cgi.cm messy.csv > messy.tsv 2> err.txt || exit 1
check "header, comments, extra columns and duplicates are ignored" \
      "$(cut -f2- list.tsv)" "$(cut -f2- messy.tsv)"
check_has "the missing ID is reported" "1 of 321 IDs" "$(cat err.txt)"
check_has "and named" "cg99999999" "$(cat err.txt)"

## ---- stdin ----------------------------------------------------------------- ##
"$KYCG" test -p "$PLAT" -m cgi.cm - < ids.txt > stdin.tsv || exit 1
check "- reads the list from stdin" "$(cut -f2- list.tsv)" "$(cut -f2- stdin.tsv)"

## ---- refusals -------------------------------------------------------------- ##
echo cg99999999 > none.txt
out=$("$KYCG" test -p "$PLAT" -m cgi.cm none.txt 2>&1 >/dev/null); rc=$?
check "a list with no probe on the platform exits non-zero" 1 "$rc"
check_has "and says to check -p" "-p" "$out"

out=$("$KYCG" test -p NOSUCH -m cgi.cm ids.txt 2>&1 >/dev/null); rc=$?
check "an unknown platform exits non-zero" 1 "$rc"
check_has "and names it" "NOSUCH" "$out"

out=$("$KYCG" test -p "$PLAT" -c -m cgi.cm ids.txt 2>&1 >/dev/null); rc=$?
check "-p with -c exits non-zero" 1 "$rc"

out=$(YAME_DATA_HOME="$PWD/nostore" "$KYCG" test -p "$PLAT" -m cgi.cm ids.txt \
        2>&1 >/dev/null); rc=$?
check "a missing ordering exits non-zero" 1 "$rc"
check_has "and names the file it wanted" "$PLAT.ordering.tsv.gz" "$out"

done_testing
