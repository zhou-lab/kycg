## title: install kycg and yame
## norun: installs packages into the reader's environment
## kycg ships one binary. `yame` is a separate package -- kycg links libyame
## statically, but the row-space check below is the yame COMMAND.
conda install -c zhou-lab -c conda-forge kycg yame
