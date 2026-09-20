#!/bin/bash
## `kycg test`: the numbers. A 40-row fixture small enough that the 2x2
## counts can be recounted in awk and the hypergeometric tail computed
## EXACTLY in integer arithmetic (math.comb), so every expected value here
## comes from an independent implementation rather than a pasted constant.
set -uo pipefail
. "$(dirname "$0")/lib.sh"
new_workdir

## ---- fixture: 40 CpGs, a query mask and a 3-state knowledgebase --------- ##
## Column 1 is the query (1 = in the set), column 2 the chromatin-like state.
cat > fixture.txt <<'T'
1	Enh
1	Enh
0	Enh
1	Enh
0	Enh
1	Prom
0	Prom
1	Prom
0	Prom
0	Prom
1	Quies
0	Quies
0	Quies
0	Quies
0	Quies
0	Quies
1	Quies
0	Quies
0	Quies
0	Quies
1	Enh
0	Enh
1	Enh
0	Prom
1	Prom
0	Quies
0	Quies
1	Quies
0	Quies
0	Quies
1	Enh
0	Enh
0	Prom
1	Prom
0	Quies
0	Quies
0	Quies
1	Quies
0	Quies
0	Quies
T
cut -f1 fixture.txt > q.txt
cut -f2 fixture.txt > s.txt
pack_binary q.txt query.cg
pack_state  s.txt  kb.cm

"$KYCG" test -m kb.cm query.cg > res.tsv 2> err.txt || {
  echo "  FAIL kycg test exited non-zero"; sed 's/^/       /' err.txt; exit 1; }

## ---- the header is the documented one ---------------------------------- ##
want_hdr="query_file	query	db_file	db	group	nU	nQ	nD	overlap	estimate	log10_p	p_value	fdr	neglog10_fdr	cf_jaccard	cf_mcc	cf_overlap	cf_npmi	cf_dice	beta	depth"
check "header columns" "$want_hdr" "$(head -1 res.tsv)"

## One row per (query sample, knowledgebase record): 1 sample x 3 states.
check "row count" 3 "$(tail -n +2 res.tsv | wc -l | tr -d ' ')"
check "file column is a basename" "kb.cm" "$(field db_file Enh res.tsv)"

## ---- the 2x2 counts, recounted from the same text ---------------------- ##
for db in Enh Prom Quies; do
  eval "$(awk -F'\t' -v db="$db" '
    {nU++; if($1==1) nQ++; if($2==db){nD++; if($1==1) nDQ++}}
    END{printf "enU=%d; enQ=%d; enD=%d; enDQ=%d\n", nU, nQ, nD, nDQ}' fixture.txt)"
  check "$db nU"      "$enU"  "$(field nU      "$db" res.tsv)"
  check "$db nQ"      "$enQ"  "$(field nQ      "$db" res.tsv)"
  check "$db nD"      "$enD"  "$(field nD      "$db" res.tsv)"
  check "$db overlap" "$enDQ" "$(field overlap "$db" res.tsv)"

  ## ---- effect sizes: the three that are pure arithmetic on the counts --- ##
  ## (MCC and NPMI are checked by tests/test_enrich.c against their own
  ## reference values; repeating their formulas here would only re-test the
  ## same algebra twice.)
  ej=$(awk -v a="$enDQ" -v d="$enD" -v q="$enQ" 'BEGIN{printf "%.9g", a/(d+q-a)}')
  eo=$(awk -v a="$enDQ" -v d="$enD" -v q="$enQ" 'BEGIN{m=(d<q?d:q); printf "%.9g", a/m}')
  ed=$(awk -v a="$enDQ" -v d="$enD" -v q="$enQ" 'BEGIN{printf "%.9g", 2*a/(d+q)}')
  check_close "$db cf_jaccard" "$ej" "$(field cf_jaccard "$db" res.tsv)" 1e-9
  check_close "$db cf_overlap" "$eo" "$(field cf_overlap "$db" res.tsv)" 1e-9
  check_close "$db cf_dice"    "$ed" "$(field cf_dice    "$db" res.tsv)" 1e-9

  ## ---- the odds ratio ---------------------------------------------------- ##
  ee=$(awk -v a="$enDQ" -v d="$enD" -v q="$enQ" -v u="$enU" \
      'BEGIN{printf "%.9g", log((a*(u-q-d+a))/((q-a)*(d-a)))/log(2)}')
  check_close "$db estimate" "$ee" "$(field estimate "$db" res.tsv)" 1e-8
done

## ---- the hypergeometric tail, exactly ---------------------------------- ##
## 40 rows means the sum over the tail is a ratio of binomial coefficients
## with no floating point in it at all. This is the check that the log-space
## saddle-point implementation is not merely self-consistent.
if command -v python3 >/dev/null 2>&1; then
  for db in Enh Prom Quies; do
    eval "$(awk -F'\t' -v db="$db" '
      {nU++; if($1==1) nQ++; if($2==db){nD++; if($1==1) nDQ++}}
      END{printf "enU=%d; enQ=%d; enD=%d; enDQ=%d\n", nU, nQ, nD, nDQ}' fixture.txt)"
    ep=$(python3 -c '
import sys
from math import comb, log10
from fractions import Fraction
nU, nQ, nD, nDQ = (int(a) for a in sys.argv[1:])
tot = comb(nU, nD)
p = sum(Fraction(comb(nQ, x) * comb(nU - nQ, nD - x), tot)
        for x in range(nDQ, min(nQ, nD) + 1))
print("%.12g" % log10(p))' "$enU" "$enQ" "$enD" "$enDQ")
    check_close "$db log10_p (exact)" "$ep" "$(field log10_p "$db" res.tsv)" 1e-9
  done
else
  echo "  skip log10_p exact check (no python3)"
fi

## ---- alternatives ------------------------------------------------------ ##
## greater and less are the two tails of one distribution: at the observed
## count they must together cover more than the whole mass (both include the
## point itself), and two.sided is twice the smaller, capped at 1.
"$KYCG" test -a less      -m kb.cm query.cg > less.tsv || exit 1
"$KYCG" test -a two.sided -m kb.cm query.cg > two.tsv  || exit 1
for db in Enh Prom Quies; do
  g=$(field log10_p "$db" res.tsv)
  l=$(field log10_p "$db" less.tsv)
  t=$(field log10_p "$db" two.tsv)
  awk -v g="$g" -v l="$l" 'BEGIN{exit !(10^g + 10^l > 1)}' ||
    { echo "  FAIL $db greater+less should exceed 1 (got $g, $l)"; fails=$((fails+1)); }
  et=$(awk -v g="$g" -v l="$l" \
      'BEGIN{m=(g<l?g:l)+log(2)/log(10); printf "%.12g", (m>0?0:m)}')
  check_close "$db two.sided" "$et" "$t" 1e-9
done

## An alternative that is not one of the three must be refused, not guessed.
out=$("$KYCG" test -a sideways -m kb.cm query.cg 2>&1); rc=$?
check "bad -a exits non-zero" 1 "$rc"

## ---- beta and depth come from a sequencing query ----------------------- ##
## A binary mask has no methylation level to report; an M/U query does, and
## those two columns are what the metagene example reads.
awk '{m=($1==1)?8:1; u=($1==1)?2:9; print m"\t"u}' q.txt > mu.txt
pack_mu mu.txt seq.cg
"$KYCG" test -m kb.cm seq.cg > seq.tsv || exit 1
b=$(field beta Enh seq.tsv)
awk -v b="$b" 'BEGIN{exit !(b>=0 && b<=1)}' ||
  { echo "  FAIL beta out of [0,1]: $b"; fails=$((fails+1)); }
d=$(field depth Enh seq.tsv)
check "depth of a 10-read fixture" "10.000" "$d"

done_testing
