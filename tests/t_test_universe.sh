#!/bin/bash
## `kycg test -u`: restrict the universe. Dropping a row from the universe is
## the same as never having had it, so the reference is `yame rowsub` cutting
## the query and the knowledgebase down to the kept rows first: both routes
## must give the same table. The rest is the spec forms and the refusals.
set -uo pipefail
. "$(dirname "$0")/lib.sh"
new_workdir

## A multiple of 8 rows, so the last restriction byte is exactly full: the
## case where ANDing two restrictions once read past the buffer.
N=4000

## Knowledgebase: three states in blocks, so every state overlaps every
## restriction differently.
awk -v n="$N" 'BEGIN{for(i=1;i<=n;i++)
  print (i%7<2)?"A":((i%7<4)?"B":"C")}' > states.txt
pack_state states.txt kb.cm

## Query: every 3rd row of the first 3000.
awk -v n="$N" 'BEGIN{for(i=1;i<=n;i++) print (i<=3000 && i%3==0)?1:0}' > q.txt
pack_binary q.txt q.cg

## Two restrictions: drop the last 1200 rows, and drop every 5th row.
awk -v n="$N" 'BEGIN{for(i=1;i<=n;i++) print (i<=n-1200)?1:0}' > r1.txt
awk -v n="$N" 'BEGIN{for(i=1;i<=n;i++) print (i%5)?1:0}' > r2.txt
pack_binary r1.txt r1.cg
pack_binary r2.txt r2.cg

## A categorical restriction: keep rows whose label is "keep".
awk -v n="$N" 'BEGIN{for(i=1;i<=n;i++) print (i%4)?"keep":"drop"}' > r3.txt
pack_state r3.txt r3.cm

## The comparable columns: everything but the file names (group is the
## knowledgebase file too).
cols() { awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i!="query_file" && $i!="db_file" && $i!="group") k[i]=1}
                     {s=""; for(i=1;i<=NF;i++) if(k[i]) s=s (s?"\t":"") $i; print s}' "$1" | sort; }

## reference <keep-bits.txt> <label> -- rowsub both files to the kept rows
reference() {
  awk '$1==1{print NR}' "$1" > idx.txt
  "$YAME" rowsub -l idx.txt q.cg > q_sub.cg || exit 1
  "$YAME" rowsub -l idx.txt kb.cm > kb_sub.cm || exit 1
  "$KYCG" test -m kb_sub.cm q_sub.cg > "ref_$2.tsv" || exit 1
}

## ---- one restriction ------------------------------------------------------- ##
"$KYCG" test -u r1.cg -m kb.cm q.cg > u1.tsv || exit 1
reference r1.txt r1
check "-u matches rowsub on every column" "$(cols ref_r1.tsv)" "$(cols u1.tsv)"
check "the universe is the kept rows" 2800 "$(field nU A u1.tsv)"

## ---- two restrictions, ANDed --------------------------------------------- ##
"$KYCG" test -u r1.cg -u r2.cg -m kb.cm q.cg > u12.tsv || exit 1
paste -d' ' r1.txt r2.txt | awk '{print ($1 && $2)?1:0}' > r12.txt
reference r12.txt r12
check "two -u are ANDed, as rowsub to the intersection" "$(cols ref_r12.tsv)" "$(cols u12.tsv)"

## ---- a state of a categorical record ------------------------------------- ##
"$KYCG" test -u r3.cm:keep -m kb.cm q.cg > u3.tsv || exit 1
awk '{print ($1=="keep")?1:0}' r3.txt > r3bits.txt
reference r3bits.txt r3
check "FILE:STATE keeps that state's rows" "$(cols ref_r3.tsv)" "$(cols u3.tsv)"

## ---- a set knowledgebase: the popcount path, not YAME's fmt6 summary ---- ##
awk -v n="$N" 'BEGIN{for(i=1;i<=n;i++) print (i%11<4)?1:0}' > m.txt
pack_binary m.txt m.cm
"$KYCG" test -u r1.cg -u r2.cg -m m.cm q.cg > ub.tsv || exit 1
awk '$1==1{print NR}' r12.txt > idx.txt
"$YAME" rowsub -l idx.txt q.cg > q_sub.cg || exit 1
"$YAME" rowsub -l idx.txt m.cm > m_sub.cm || exit 1
"$KYCG" test -m m_sub.cm q_sub.cg > ref_b.tsv || exit 1
check "a set knowledgebase matches rowsub too" "$(cols ref_b.tsv)" "$(cols ub.tsv)"

## ---- refusals -------------------------------------------------------------- ##
out=$("$KYCG" test -u r3.cm -m kb.cm q.cg 2>&1 >/dev/null); rc=$?
check "a categorical file without a state exits non-zero" 1 "$rc"
check_has "and asks for one" "name the state" "$out"
check_has "in the form it was given" "r3.cm:<state>" "$out"

out=$("$KYCG" test -u r3.cm:nosuch -m kb.cm q.cg 2>&1 >/dev/null); rc=$?
check "an unknown state exits non-zero" 1 "$rc"
check_has "and names it" "nosuch" "$out"

out=$("$KYCG" test -u nosuch.cg -m kb.cm q.cg 2>&1 >/dev/null); rc=$?
check "a missing file exits non-zero" 1 "$rc"
check_has "and names it" "nosuch.cg" "$out"

awk 'BEGIN{for(i=1;i<=16;i++) print 1}' > short.txt
pack_binary short.txt short.cg
out=$("$KYCG" test -u short.cg -m kb.cm q.cg 2>&1 >/dev/null); rc=$?
check "a restriction of another row count exits non-zero" 1 "$rc"
check_has "and says so" "row count mismatch" "$out"

## A valued query: YAME's summary takes its universe from every row, so
## blanking rows in the query alone would give wrong statistics.
awk -v n="$N" 'BEGIN{for(i=1;i<=n;i++) print (i%3)?"3\t1":"0\t0"}' > mu.txt
pack_mu mu.txt mu.cg
out=$("$KYCG" test -u r1.cg -m kb.cm mu.cg 2>&1 >/dev/null); rc=$?
check "a format 3 query with -u exits non-zero" 1 "$rc"
check_has "and says to binarize" "binarize" "$out"

out=$("$KYCG" test -c -u r1.cg -m kb.cm mu.cg 2>&1 >/dev/null); rc=$?
check "-u with -c exits non-zero" 1 "$rc"
check_has "and says they do not combine" "does not combine with -c" "$out"

done_testing
