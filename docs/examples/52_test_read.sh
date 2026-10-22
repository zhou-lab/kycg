## title: read the result
## Sorted by significance. p_value and fdr underflow to 0 for strong results
## by design -- threshold on log10_p / neglog10_fdr, which stay in log space.
cd kycg-examples/data
kycg test -m mm10:ChromHMM,CGI,PMD onecell.cg > res3.tsv
cut -f3,4,10,13 res3.tsv | head -6

## the rows the page prints, in order
test "$(cut -f4 res3.tsv | sed -n 2,6p | tr '\n' ' ')" = \
     "commonHMD Quies QuiesG Tx Shore "
## the underflow the page warns about is real: a top hit's fdr reads 0
## while its neglog10_fdr is large.
awk -F'\t' 'NR==2 && $13=="0" && $14>100 {found=1} END{exit !found}' res3.tsv
