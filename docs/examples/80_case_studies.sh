## title: the case studies, end to end
## The two scripts behind the figure gallery. They fetch what they need, take
## their tools from PATH, and regenerate results/ and figures/ in the clone.
bash kycg-examples/scripts/tissue_tfbs_enrichment.sh
bash kycg-examples/scripts/metagene_profile.sh

## the biology the page claims: fetal brain returns neuronal regulators
grep -qi 'DLX\|OTX2\|LHX2\|EOMES\|PAX6' kycg-examples/results/fetal_brain_2s.tsv
## and every advertised figure was written
ls kycg-examples/figures/tissue_tfbs.* kycg-examples/figures/brain_volcano.* \
   kycg-examples/figures/hep_bar.*
