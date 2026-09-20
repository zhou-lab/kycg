## title: test against everything in a collection
## -m is repeatable, so a whole collection is a shell loop -- and pooling them
## into ONE invocation decompresses the query once instead of once per set.
kycg test $(for f in "${YAME_DATA_HOME:?}"/mm10/KYCG/*.cm; do
        printf -- '-m %s ' $f; done) \
  kycg-examples/data/onecell.cg > res_all.tsv

## more knowledgebases than the single-set run above
test "$(wc -l < res_all.tsv)" -gt "$(wc -l < res.tsv)"
