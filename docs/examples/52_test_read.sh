## title: read the result
## Sorted by significance. p_value and fdr underflow to 0 for strong results
## by design -- threshold on log10_p / neglog10_fdr, which stay in log space.
cut -f3,4,10,13 res_all.tsv | head -6

## the underflow the page warns about is real: a top hit's fdr reads 0
## while its neglog10_fdr is large.
awk -F'\t' 'NR==2 && $13=="0" && $14>100 {found=1} END{exit !found}' res_all.tsv
