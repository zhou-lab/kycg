## title: the core test
## One row per (query sample, knowledgebase record): the 2x2 counts, a
## hypergeometric tail in log10, six effect sizes, and a grouped FDR.
## Runs in kycg-examples/data, as the page says after the install block.
cd kycg-examples/data
kycg fetch -f mm10:ChromHMM,CGI,PMD
kycg test -m mm10:ChromHMM onecell.cg > res.tsv

## the documented columns are the ones emitted
head -1 res.tsv | grep -q 'log10_p'
head -1 res.tsv | grep -q 'neglog10_fdr'
test "$(wc -l < res.tsv)" -gt 1
