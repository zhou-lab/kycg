## title: fetch a knowledgebase into the store
## -f skips the browser: plan, then fetch. The row list (cpg_nocontig.cr) comes
## too -- it is what gives a row its identity -- so this is 2 files, 22.8 MB.
kycg fetch -f mm10:CGI

## it landed where the page says it does, and a re-fetch is a no-op
test -f "${YAME_DATA_HOME:?}/mm10/KYCG/CGI.20220904.cm"
test -f "${YAME_DATA_HOME:?}/mm10/cpg_nocontig.cr"
kycg fetch -f mm10:CGI
