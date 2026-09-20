#!/bin/bash
## `kycg test`: the options around the statistics -- output plumbing, FDR
## grouping, naming. Each one is a flag a user reaches for once and never
## checks again, so a silent regression in any of them is invisible.
set -uo pipefail
. "$(dirname "$0")/lib.sh"
new_workdir

## Two knowledgebases, so FDR grouping has something to group BY: the default
## corrects within each knowledgebase file, -G corrects across all of them.
awk 'BEGIN{srand(7); for(i=0;i<300;i++) print (rand()<0.3)?1:0}' > q.txt
awk 'BEGIN{srand(11); s[0]="Enh"; s[1]="Prom"; s[2]="Quies"; s[3]="Tx";
           for(i=0;i<300;i++) print s[int(rand()*4)]}' > s1.txt
awk 'BEGIN{srand(13); s[0]="CGI"; s[1]="Shore"; s[2]="OpenSea";
           for(i=0;i<300;i++) print s[int(rand()*3)]}' > s2.txt
pack_binary q.txt query.cg
pack_state  s1.txt chrom.cm
pack_state  s2.txt cgi.cm

## ---- -m is repeatable; the group column names the knowledgebase --------- ##
"$KYCG" test -m chrom.cm -m cgi.cm query.cg > both.tsv || exit 1
check "rows from both knowledgebases" 7 "$(tail -n +2 both.tsv | wc -l | tr -d ' ')"
check "group is the knowledgebase file" "chrom.cm" "$(field group Enh both.tsv)"
check "group for the second file"       "cgi.cm"   "$(field group CGI both.tsv)"

## ---- -H drops the header ----------------------------------------------- ##
n_with=$(wc -l < both.tsv)
"$KYCG" test -H -m chrom.cm -m cgi.cm query.cg > noh.tsv || exit 1
check "-H drops exactly one line" "$((n_with - 1))" "$(wc -l < noh.tsv | tr -d ' ')"
check_lacks "-H leaves no header" "query_file" "$(head -1 noh.tsv)"

## ---- -o writes to a file and nothing to stdout -------------------------- ##
out=$("$KYCG" test -o via_o.tsv -m chrom.cm query.cg 2>/dev/null)
check "-o leaves stdout empty" "" "$out"
"$KYCG" test -m chrom.cm query.cg > via_redirect.tsv || exit 1
if ! cmp -s via_o.tsv via_redirect.tsv; then
  echo "  FAIL -o and a redirect differ"; diff via_o.tsv via_redirect.tsv | head -5
  fails=$((fails + 1))
fi

## ---- -F reports the path as given, not the basename -------------------- ##
"$KYCG" test    -m "$PWD/chrom.cm" query.cg > plain.tsv || exit 1
"$KYCG" test -F -m "$PWD/chrom.cm" query.cg > fullp.tsv || exit 1
check "default reports the basename" "chrom.cm" "$(field db_file Enh plain.tsv)"
check "-F reports the full path" "$PWD/chrom.cm" "$(field db_file Enh fullp.tsv)"

## ---- -M loads the knowledgebase into memory: same numbers -------------- ##
"$KYCG" test -M -m chrom.cm query.cg > mem.tsv || exit 1
if ! cmp -s mem.tsv via_redirect.tsv; then
  echo "  FAIL -M changes the output"; diff mem.tsv via_redirect.tsv | head -5
  fails=$((fails + 1))
fi

## ---- -G: global FDR vs within-knowledgebase ---------------------------- ##
## -G moves the correction from "within each knowledgebase" to "across all of
## them", so the FDR column must change while nothing upstream of it does.
## (The direction is NOT checked: BH is a step-up procedure, and its
## cumulative minimum runs over every test in the correction set, so pooling
## two groups can move a given row either way. The values themselves are
## tests/test_enrich.c's job, against R's p.adjust.)
"$KYCG" test -G -m chrom.cm -m cgi.cm query.cg > global.tsv || exit 1
n_diff=0
for db in Enh Prom Quies Tx CGI Shore OpenSea; do
  check "$db log10_p is untouched by -G" \
        "$(field log10_p "$db" both.tsv)" "$(field log10_p "$db" global.tsv)"
  a=$(field neglog10_fdr "$db" both.tsv)
  b=$(field neglog10_fdr "$db" global.tsv)
  awk -v a="$a" -v b="$b" 'BEGIN{exit (a==b)}' && n_diff=$((n_diff + 1))
done
awk -v n="$n_diff" 'BEGIN{exit !(n > 0)}' ||
  { echo "  FAIL -G changed no FDR at all"; fails=$((fails + 1)); }

## The group column names the FDR stratum, so it must track what -G did: the
## default is the knowledgebase file; under -G the stratum is global, and the
## column should say so rather than keep naming a file the correction pooled.
for db in Enh Prom Quies Tx; do
  check "$db group is the file by default" "chrom.cm" "$(field group "$db" both.tsv)"
  check "$db group is (global) under -G"   "(global)" "$(field group "$db" global.tsv)"
done

## ---- -s names the query samples ---------------------------------------- ##
printf 'my_sample\n' > names.txt
"$KYCG" test -s names.txt -m chrom.cm query.cg > named.tsv || exit 1
check "-s renames the query" "my_sample" "$(field query Enh named.tsv)"

## ---- the group column is the knowledgebase file ------------------------- ##
## It is what FDR is stratified by, so it has to identify the file exactly --
## a prettified name would make two different files look like one stratum.
cp chrom.cm KYCG.MM285.ChromHMM.20220414.cm
"$KYCG" test -m KYCG.MM285.ChromHMM.20220414.cm query.cg > named_kb.tsv || exit 1
check "the group is the file name, verbatim" "KYCG.MM285.ChromHMM.20220414.cm" \
      "$(field group Enh named_kb.tsv)"

## ---- an option that is not one ------------------------------------------ ##
out=$("$KYCG" test -Z -m chrom.cm query.cg 2>&1 >/dev/null); rc=$?
awk -v rc="$rc" 'BEGIN{exit !(rc != 0)}' ||
  { echo "  FAIL an unknown option was accepted"; fails=$((fails + 1)); }

## ---- several queries in one call --------------------------------------- ##
awk 'BEGIN{srand(23); for(i=0;i<300;i++) print (rand()<0.5)?1:0}' > q2.txt
pack_binary q2.txt query2.cg
"$KYCG" test -m chrom.cm query.cg query2.cg > multi.tsv || exit 1
check "two queries give two blocks" 8 "$(tail -n +2 multi.tsv | wc -l | tr -d ' ')"
check "both query files are named" 2 \
  "$(tail -n +2 multi.tsv | cut -f1 | sort -u | wc -l | tr -d ' ')"

done_testing
