## title: annotate a table of probe IDs
## The lookup half of the same data: where does THIS probe fall? Arrays only,
## and the platform is named, never inferred.
kycg fetch -f EPIC:CGI,ProbeType

printf 'Probe_ID\teffect\ncg00000029\t0.51\ncg00000109\t-0.22\n' > hits.tsv
kycg annotate -m EPIC:CGI,ProbeType hits.tsv > annotated.tsv
cat annotated.tsv

## every input column and row kept, in order, plus one column per set
test "$(wc -l < annotated.tsv)" = "$(wc -l < hits.tsv)"
test "$(head -1 annotated.tsv | awk -F'\t' '{print NF}')" -eq 4

## -i gives a 0/1 indicator column per set instead of joined names. The
## columns are named <knowledgebase>:<set>, so one CGI knowledgebase becomes
## CGI:Island, CGI:OpenSea, CGI:Shelf, CGI:Shore.
kycg annotate -i -m EPIC:CGI hits.tsv > indicator.tsv
head -1 indicator.tsv
grep -q 'CGI:Island' indicator.tsv
grep -q 'CGI:Shore'  indicator.tsv
## and the values are 0/1, not names
awk -F'\t' 'NR==2{for(i=3;i<=NF;i++) if($i!="0" && $i!="1") exit 1}' indicator.tsv
