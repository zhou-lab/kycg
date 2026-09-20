## title: the catalogue as TSV
## Redirected stdout means "give me data": no browser, no questions. This is
## the documented scripting form, and it must keep working for `cut`.
kycg fetch | cut -f1

## every collection the build pins shows up
kycg fetch | tail -n +2 | grep -q '^mm10'
kycg fetch | tail -n +2 | grep -q '^EPIC'
