## title: test against everything in a collection
## -m mm10 names every set fetched for mm10; pooling them into ONE invocation
## decompresses the query once instead of once per set.
cd kycg-examples/data
kycg test -m mm10 onecell.cg > res_all.tsv

## more knowledgebases than the single-set run above
test "$(wc -l < res_all.tsv)" -gt "$(wc -l < res.tsv)"
