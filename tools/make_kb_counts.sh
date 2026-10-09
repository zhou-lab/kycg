#!/usr/bin/env bash
## Count the rows in every set of every knowledgebase the registry offers,
## for the Knowledgebases tab of docs/index.html.
##
##   tools/make_kb_counts.sh -o docs/knowledgebase_counts.tsv
##
## The sets live in the .cm files, not in the registry, so this reads the
## store ($YAME_DATA_HOME). Each .cm is summarized against an all-ones
## universe of its collection's row count, and N_mask is the set size.
## A genome .cm is tens of millions of rows: run it under sbatch. The output
## is committed, so `make test` and a docs rebuild never need the store.
set -euo pipefail
export LC_ALL=C

cd "$(dirname "$0")/.."
out=${2:-docs/knowledgebase_counts.tsv}
[ "${1:-}" = "-o" ] || { echo "usage: $0 -o <out.tsv>" >&2; exit 2; }
store=${YAME_DATA_HOME:?set YAME_DATA_HOME to the store}
yame=${YAME:-external/YAME/yame}
tsv=external/YAME/tools/registry/files.tsv
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

## collection -> row count, from the files that define a row space
grep -v '^#' "$tsv" |
  awk -F'\t' '$5 != "" {split($2, p, "/"); print p[1] "\t" $5}' \
  > "$work/rows.tsv"

## A knowledgebase with more than this many sets is a tiling (genomic bins,
## k-mers) rather than a list, so it gets one summary row -- set "*", n = the
## number of sets, then the min / median / max set size and the first and last
## set names -- instead of a row per set.
max_sets=1500
printf 'collection\tfile\tset\tn\tmin\tmedian\tmax\tfirst\tlast\n' \
  > "$work/out.tsv"
while IFS=$'\t' read -r coll rows; do
  ## one all-ones universe per row space
  ## (awk, not `yes | head`: under pipefail yes's SIGPIPE fails the script)
  awk -v n="$rows" 'BEGIN {for (i = 0; i < n; i++) print 1}' |
    "$yame" pack -f b - "$work/$coll.all.cg"
  grep -v '^#' "$tsv" | cut -f2 |
    grep "^$coll/KYCG/.*\.cm$" |
    while read -r path; do
      f=$(basename "$path")
      ## Mask column is the record name (format 0) or the state (format 2)
      "$yame" summary -m "$store/$path" "$work/$coll.all.cg" |
        awk -F'\t' 'NR > 1 {print $4 "\t" $7}' > "$work/sets.tsv"
      k=$(wc -l < "$work/sets.tsv")
      if [ "$k" -le "$max_sets" ]; then
        awk -F'\t' -v c="$coll" -v f="$f" -v OFS='\t' \
          '{print c, f, $1, $2, "", "", "", "", ""}' "$work/sets.tsv"
      else
        first=$(head -n 1 "$work/sets.tsv" | cut -f1)
        last=$(tail -n 1 "$work/sets.tsv" | cut -f1)
        cut -f2 "$work/sets.tsv" | sort -n |
          awk -v c="$coll" -v f="$f" -v a="$first" -v z="$last" -v OFS='\t' \
            '{v[NR] = $1}
             END {print c, f, "*", NR, v[1], v[int((NR + 1) / 2)], v[NR], a, z}'
      fi >> "$work/out.tsv"
    done
  rm -f "$work/$coll.all.cg"
done < "$work/rows.tsv"

mv "$work/out.tsv" "$out"
echo "wrote $out: $(($(wc -l < "$out") - 1)) sets"
