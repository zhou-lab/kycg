#!/bin/bash
## `kycg test -c`: the set enrichment score of a continuous query. Checks the
## sign, the symmetry zhou-lab/knowYourCG#6 found broken in R, that every mask
## format scores the same set the same way, and the refusals.
set -uo pipefail
. "$(dirname "$0")/lib.sh"
new_workdir

## 200 rows of values, every tenth one missing. Set "top" is the 20 largest
## values, "bot" the 20 smallest, "rnd" every seventh row.
awk 'BEGIN{srand(5); for(i=0;i<200;i++)
       if(i%10==9) print "NA"; else printf "%.6f\n", rand()}' > v.txt
"$YAME" pack -f n v.txt q.cg
awk '{print ($1=="NA") ? "NA" : sprintf("%.6f", 1-$1)}' v.txt > flip.txt
"$YAME" pack -f n flip.txt flip.cg

awk '$1!="NA"{print NR-1"\t"$1}' v.txt | sort -k2,2g > ranked.txt
awk 'NR<=20{print $1}' ranked.txt > bot.rows
tail -n 20 ranked.txt | cut -f1 > top.rows
mark() { awk -v f="$1" 'BEGIN{while((getline l < f)>0) s[l]=1}
                        {print (NR-1 in s) ? 1 : 0}' v.txt; }
mark top.rows > top.txt; mark bot.rows > bot.txt
awk '{print (NR % 7 == 0) ? 1 : 0}' v.txt > rnd.txt
awk '{print (NR % 10 == 0) ? 1 : 0}' v.txt > gone.txt   # only missing rows
for s in top bot rnd gone; do pack_binary $s.txt $s.cm; done
cat top.cm bot.cm rnd.cm gone.cm > sets.cm
printf 'top\nbot\nrnd\ngone\n' > sets.names
"$YAME" index -s sets.names sets.cm

## ---- the sign: GSEA's ---------------------------------------------------- ##
"$KYCG" test -c -m sets.cm q.cg > res.tsv || exit 1
check "same columns as a count test" \
  "$("$KYCG" test -m sets.cm q.cg 2>/dev/null | head -1)" "$(head -1 res.tsv)"
check "four sets, four rows" 4 "$(tail -n +2 res.tsv | wc -l | tr -d ' ')"
check_close "largest values: estimate +1" 1 "$(field estimate top res.tsv)" 1e-12
check_close "smallest values: estimate -1" -1 "$(field estimate bot res.tsv)" 1e-12
check "universe is the rows with a value" 180 "$(field nU top res.tsv)"
check "nD counts the set" 20 "$(field nD top res.tsv)"
awk -v p="$(field log10_p top res.tsv)" 'BEGIN{exit !(p < -10)}' ||
  { echo "  FAIL a perfect set is not significant"; fails=$((fails + 1)); }

## ---- a set with no member that has a value: estimate 0, p 1 ------------ ##
check "no valued members: overlap 0" 0 "$(field overlap gone res.tsv)"
check "no valued members: estimate 0" 0 "$(field estimate gone res.tsv)"
check "no valued members: log10_p 0" 0 "$(field log10_p gone res.tsv)"
check "no valued members: no beta" NA "$(field beta gone res.tsv)"

## ---- reversing the query negates the estimate and keeps the p-value ----- ##
## The property knowYourCG#6 found broken: two tails with separate nulls
## gave the same set a different p-value when the query was negated.
"$KYCG" test -c -m sets.cm flip.cg > flip.tsv || exit 1
for s in top bot rnd; do
  check_close "$s estimate negates" \
    "$(field estimate $s res.tsv)" \
    "$(awk -v e="$(field estimate $s flip.tsv)" 'BEGIN{printf "%.12f", -e}')" 1e-9
  check "$s log10_p is unchanged" \
    "$(field log10_p $s res.tsv)" "$(field log10_p $s flip.tsv)"
done

## ---- deterministic: the null is seeded from the set size ---------------- ##
"$KYCG" test -c -m sets.cm q.cg > again.tsv || exit 1
cmp -s res.tsv again.tsv ||
  { echo "  FAIL two runs differ"; fails=$((fails + 1)); }

## ---- every mask format scores the same set the same way ----------------- ##
awk '{print ($1) ? "top" : "rest"}' top.txt > state.txt
pack_state state.txt state.cm
"$KYCG" test -c -m state.cm q.cg > st.tsv || exit 1
check "a state mask gives a row per state" 2 "$(tail -n +2 st.tsv | wc -l | tr -d ' ')"
check "format 2 matches format 0" \
  "$(field estimate top res.tsv)" "$(field estimate top st.tsv)"
check "format 2 rest is the mirror" \
  "$(field estimate top res.tsv)" \
  "$(awk -v e="$(field estimate rest st.tsv)" 'BEGIN{print -e}')"
awk '{print $1 "\t1"}' top.txt > d.txt
"$YAME" pack -f d d.txt top6.cm
"$KYCG" test -c -m top6.cm q.cg > f6.tsv || exit 1
check "format 6 matches format 0" \
  "$(field estimate top res.tsv)" "$(tail -n +2 f6.tsv | cut -f"$(col estimate f6.tsv)")"

## ---- a format 3 query is scored as beta --------------------------------- ##
awk '{if($1=="NA") print "0\t0";
      else {m=int($1*1000+0.5); print m "\t" 1000-m}}' v.txt > mu.txt
pack_mu mu.txt mu.cg
"$KYCG" test -c -m sets.cm mu.cg > mu.tsv || exit 1
check_close "format 3: largest betas, estimate +1" 1 "$(field estimate top mu.tsv)" 1e-12
check "format 3: uncovered rows leave the universe" 180 "$(field nU top mu.tsv)"
check_close "beta is the members' mean value" \
  "$(awk -v f=top.rows 'BEGIN{while((getline l<f)>0) s[l+1]=1}
       NR in s{t+=$1; n++} END{printf "%.4f", t/n}' v.txt)" \
  "$(field beta top mu.tsv)" 1e-3

## ---- the refusals ------------------------------------------------------- ##
err=$("$KYCG" test -c -a less -m sets.cm q.cg 2>&1 >/dev/null)
check_has "-a with -c is refused" "-a does not apply with -c" "$err"
err=$("$KYCG" test -c -m sets.cm top.cm 2>&1 >/dev/null)
check_has "a format 0 query is refused" "must be format 3" "$err"
err=$("$KYCG" test -c -m mu.cg q.cg 2>&1 >/dev/null)
check_has "a format 3 knowledgebase is refused" "must be format 0, 1, 2 or 6" "$err"
check_has "usage names -c" "continuous query" "$("$KYCG" test -h 2>&1)"

done_testing
