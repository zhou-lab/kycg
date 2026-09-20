## title: check the row space before testing
## Everything is positional: row i must mean the same CpG on both sides.
## `yame info` is how you confirm it -- kycg has no `info` subcommand.
cd kycg-examples/data
yame info onecell.cg chromhmm.cm

## both index the same 21,867,837-row mm10 reference
test "$(yame info onecell.cg  | tail -1 | cut -f4)" \
   = "$(yame info chromhmm.cm | tail -1 | cut -f4)"
