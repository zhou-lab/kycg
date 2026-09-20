## title: the sample queries and knowledgebases
## Neither the conda package nor the main repo ships query data; the worked
## examples live in kycg-examples, and every block after this one uses them.
## Re-running is a no-op, so the gate does not re-clone.
[ -d kycg-examples ] || git clone https://github.com/zhou-lab/kycg-examples
ls kycg-examples/data/onecell.cg kycg-examples/data/chromhmm.cm
