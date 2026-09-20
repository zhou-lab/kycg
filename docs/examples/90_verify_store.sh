## title: re-verify the store with no kycg code
## kycg writes a SHA256SUMS derived from its own registry rows, so the store is
## re-checkable by hand with no kycg code. --ignore-missing because you rarely
## hold every set the directory could hold -- the manifest describes the
## directory, not the subset you fetched.
## GNU coreutils calls it sha256sum; macOS ships shasum.
SUM=$(command -v sha256sum || command -v shasum)
case $SUM in *shasum) SUM="$SUM -a 256";; esac

cd "${YAME_DATA_HOME:?}/mm10/KYCG"
$SUM -c --ignore-missing SHA256SUMS
