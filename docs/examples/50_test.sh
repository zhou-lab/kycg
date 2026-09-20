## title: the core test
## One row per (query sample, knowledgebase record): the 2x2 counts, a
## hypergeometric tail in log10, six effect sizes, and a grouped FDR.
kycg fetch -f mm10:ChromHMM,PMD

kycg test -m "${YAME_DATA_HOME:?}/mm10/KYCG/ChromHMM.20220414.cm" \
  kycg-examples/data/onecell.cg > res.tsv

## the documented columns are the ones emitted
head -1 res.tsv | grep -q 'log10_p'
head -1 res.tsv | grep -q 'neglog10_fdr'
test "$(wc -l < res.tsv)" -gt 1
