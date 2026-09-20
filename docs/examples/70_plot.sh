## title: plot the result with cinderplot
## kycg emits TSV; cinderplot renders it as its own process. Drop the
## zero-overlap sentinel first: estimate is log2(OR), and log2(0) clamps near
## -1022, which would squash the whole x-axis.
awk -F'\t' 'NR==1||($10>-1000&&$10<1000)' res_all.tsv > v.tsv
cinderplot 'v.tsv
  + aes(estimate, neglog10_fdr)
  + geom_point() + theme_minimal()' -o v.pdf

test -s v.pdf
