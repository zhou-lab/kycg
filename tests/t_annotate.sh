#!/bin/bash
## `kycg annotate`: a TSV of probe IDs in, the same TSV plus one column per
## knowledgebase out. The whole command rests on one fact -- a probe's row
## index IS its line number in the platform ordering -- so the fixtures build
## a real (synthetic) HM27 store and check membership against rows chosen by
## construction.
set -uo pipefail
. "$(dirname "$0")/lib.sh"
new_workdir

## HM27's row count is pinned in the compiled registry and annotate refuses an
## ordering that disagrees, so the fixture has to have exactly this many rows.
N=27722
PLAT=HM27

## ---- a synthetic store ------------------------------------------------- ##
## Zero-padded ids ascend in LC_ALL=C order, which is the ordering annotate's
## binary search assumes -- and verifies.
store=$PWD/store
mkdir -p "$store/$PLAT"
{ echo "Probe_ID"
  awk -v n="$N" 'BEGIN{for(i=1;i<=n;i++) printf "cg%08d\n", i}'; } |
  gzip -c > "$store/$PLAT/$PLAT.ordering.tsv.gz"

## A state knowledgebase over the same rows: the first 1000 probes are CGI,
## the next 1000 Shore, the rest OpenSea. So cg00000001 is CGI, cg00001500 a
## Shore, cg00027722 an OpenSea -- membership we know without reading a file.
awk -v n="$N" 'BEGIN{for(i=1;i<=n;i++)
  print (i<=1000)?"CGI":((i<=2000)?"Shore":"OpenSea")}' > states.txt
pack_state states.txt cgi.cm

cat > probes.tsv <<'T'
Probe_ID	effect	pval
cg00000001	0.51	1e-4
cg00001500	-0.22	3e-2
cg00027722	0.03	0.9
cg99999999	0.10	0.5
T

export YAME_DATA_HOME="$store"

## ---- the basic annotation ---------------------------------------------- ##
"$KYCG" annotate -m cgi.cm -p "$PLAT" probes.tsv > out.tsv 2> err.txt || {
  echo "  FAIL annotate exited non-zero"; sed 's/^/       /' err.txt; exit 1; }

check "every input row is kept" 4 "$(tail -n +2 out.tsv | wc -l | tr -d ' ')"
check "input columns are kept, in order" "Probe_ID	effect	pval" \
      "$(head -1 out.tsv | cut -f1-3)"
check "one column is added" 4 "$(head -1 out.tsv | awk -F'\t' '{print NF}')"

label() { awk -F'\t' -v p="$1" '$1==p{print $4}' out.tsv; }
check "first probe is CGI"        "CGI"     "$(label cg00000001)"
check "row 1500 is a Shore"       "Shore"   "$(label cg00001500)"
check "the last row is OpenSea"   "OpenSea" "$(label cg00027722)"
## A probe the ordering does not carry has no row, so it has no membership --
## and must be reported as unknown rather than dropped or guessed.
check "a probe absent from the ordering is NA" "NA" "$(label cg99999999)"
## The other columns must come back byte for byte.
check "the data columns are untouched" "$(cut -f1-3 probes.tsv)" "$(cut -f1-3 out.tsv)"

## ---- -i indicator columns ---------------------------------------------- ##
"$KYCG" annotate -i -m cgi.cm -p "$PLAT" probes.tsv > ind.tsv || exit 1
hdr=$(head -1 ind.tsv)
check_has "-i makes a column per set (CGI)"     "cgi:CGI"     "$hdr"
check_has "-i makes a column per set (Shore)"   "cgi:Shore"   "$hdr"
check_has "-i makes a column per set (OpenSea)" "cgi:OpenSea" "$hdr"
## The indicator columns are named <knowledgebase>:<set>, so the file a set
## came from stays visible when several are annotated at once.
cgi_col=$(col "cgi:CGI" ind.tsv)
[ -n "$cgi_col" ] || { echo "  FAIL no cgi:CGI column in: $hdr"; fails=$((fails+1)); }
check "-i marks the CGI probe 1" "1" \
      "$(awk -F'\t' -v c="$cgi_col" '$1=="cg00000001"{print $c}' ind.tsv)"
check "-i marks a non-member 0" "0" \
      "$(awk -F'\t' -v c="$cgi_col" '$1=="cg00001500"{print $c}' ind.tsv)"

## ---- -c: which column holds the probe ids ------------------------------- ##
## By name and by 1-based index, and on an input whose id column is not first.
awk -F'\t' 'BEGIN{OFS="\t"} {print $2, $1, $3}' probes.tsv > swapped.tsv
"$KYCG" annotate -m cgi.cm -p "$PLAT" -c Probe_ID swapped.tsv > byname.tsv || exit 1
check "-c by name finds the ids" "CGI" \
      "$(awk -F'\t' '$2=="cg00000001"{print $4}' byname.tsv)"
"$KYCG" annotate -m cgi.cm -p "$PLAT" -c 2 swapped.tsv > byindex.tsv || exit 1
if ! cmp -s byname.tsv byindex.tsv; then
  echo "  FAIL -c by name and by index disagree"; fails=$((fails + 1)); fi

## ---- -H: no header line ------------------------------------------------- ##
tail -n +2 probes.tsv > noheader.tsv
"$KYCG" annotate -H -m cgi.cm -p "$PLAT" noheader.tsv > outH.tsv || exit 1
check "-H keeps every row (no header consumed)" 4 "$(wc -l < outH.tsv | tr -d ' ')"
check "-H annotates the first row too" "CGI" \
      "$(awk -F'\t' '$1=="cg00000001"{print $4}' outH.tsv)"

## ---- -s: the separator for a probe in several sets ---------------------- ##
## A state knowledgebase gives each row exactly one label, so membership in
## two sets at once needs a multi-record .cm: two binary masks concatenated
## and named with `yame index -s`, which is how the real ones are built.
awk -v n="$N" 'BEGIN{for(i=1;i<=n;i++) print (i<=100)?1:0}' > mA.txt
awk -v n="$N" 'BEGIN{for(i=1;i<=n;i++) print (i<=50||i>27000)?1:0}' > mB.txt
pack_binary mA.txt a.cm
pack_binary mB.txt b.cm
cat a.cm b.cm > two.cm
printf 'setA\nsetB\n' > names.txt
"$YAME" index -s names.txt two.cm
"$KYCG" annotate -m two.cm -p "$PLAT" probes.tsv > multi.tsv || exit 1
check "a probe in both sets joins with the default comma" "setA,setB" \
      "$(awk -F'\t' '$1=="cg00000001"{print $4}' multi.tsv)"
"$KYCG" annotate -s '|' -m two.cm -p "$PLAT" probes.tsv > multi_sep.tsv || exit 1
check "-s changes the separator" "setA|setB" \
      "$(awk -F'\t' '$1=="cg00000001"{print $4}' multi_sep.tsv)"

## ---- a set named the way fetch names it --------------------------------- ##
## `-m HM27:CGI` resolves against the store instead of naming a path, and the
## platform comes from the spec itself -- no -p, because the name already says
## which ordering the row indices belong to.
mkdir -p "$store/$PLAT/KYCG"
cp cgi.cm "$store/$PLAT/KYCG/CGI.20220904.cm"
"$KYCG" annotate -m "$PLAT:CGI" probes.tsv > byspec.tsv 2>spec_err.txt || {
  echo "  FAIL annotate by spec exited non-zero"; sed 's/^/       /' spec_err.txt
  fails=$((fails + 1)); }
if [ -s byspec.tsv ]; then
  check "a store set resolves by name" "CGI" \
        "$(awk -F'\t' '$1=="cg00000001"{print $4}' byspec.tsv)"
  check "and the platform came from the spec" 4 \
        "$(head -1 byspec.tsv | awk -F'\t' '{print NF}')"
fi

## ---- the refusals ------------------------------------------------------- ##
## A platform the build does not pin: annotate cannot invent an ordering.
out=$("$KYCG" annotate -m cgi.cm -p NotAPlatform probes.tsv 2>&1 >/dev/null); rc=$?
check "an unknown platform exits non-zero" 1 "$rc"
check_has "the unknown platform is named" "NotAPlatform" "$out"

## No platform, and a plain path cannot imply one.
out=$("$KYCG" annotate -m cgi.cm probes.tsv < /dev/null 2>&1 >/dev/null); rc=$?
check "a path without -p exits non-zero" 1 "$rc"
check_has "and says to name the platform" "-p" "$out"

## An ordering that is not C-sorted: the binary search would return confident
## wrong answers, so it must be refused rather than used.
bad=$PWD/badstore
mkdir -p "$bad/$PLAT"
{ echo "Probe_ID"
  awk -v n="$N" 'BEGIN{for(i=n;i>=1;i--) printf "cg%08d\n", i}'; } |
  gzip -c > "$bad/$PLAT/$PLAT.ordering.tsv.gz"
out=$(YAME_DATA_HOME="$bad" "$KYCG" annotate -m cgi.cm -p "$PLAT" probes.tsv \
        2>&1 >/dev/null); rc=$?
check "an unsorted ordering exits non-zero" 1 "$rc"
check_has "and says it is not sorted" "not sorted" "$out"

## An ordering of the wrong length is a different row space, not a smaller one.
short=$PWD/shortstore
mkdir -p "$short/$PLAT"
{ echo "Probe_ID"
  awk 'BEGIN{for(i=1;i<=100;i++) printf "cg%08d\n", i}'; } |
  gzip -c > "$short/$PLAT/$PLAT.ordering.tsv.gz"
out=$(YAME_DATA_HOME="$short" "$KYCG" annotate -m cgi.cm -p "$PLAT" probes.tsv \
        2>&1 >/dev/null); rc=$?
check "a short ordering exits non-zero" 1 "$rc"
check_has "and gives both counts" "$N" "$out"

## No ordering at all: the message has to say where it looked and how to get it.
none=$PWD/nostore
mkdir -p "$none"
out=$(YAME_DATA_HOME="$none" "$KYCG" annotate -m cgi.cm -p "$PLAT" probes.tsv \
        2>&1 >/dev/null); rc=$?
check "a missing ordering exits non-zero" 1 "$rc"
check_has "and names the path it wanted" "$PLAT.ordering.tsv.gz" "$out"

done_testing
