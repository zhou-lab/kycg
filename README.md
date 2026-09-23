<div align="center">

# KnowYourCG 2

[![build](https://github.com/zhou-lab/kycg/actions/workflows/conda-build.yml/badge.svg)](https://github.com/zhou-lab/kycg/actions/workflows/conda-build.yml)
[![conda](https://img.shields.io/conda/vn/zhou-lab/kycg?label=conda)](https://anaconda.org/zhou-lab/kycg)
[![license](https://img.shields.io/badge/license-BSD--2--Clause%20(academic)%20%2F%20commercial-blue.svg)](LICENSE)
[![coverage](https://img.shields.io/endpoint?url=https%3A%2F%2Fzhou-lab.github.io%2Fkycg%2Fcoverage.json)](tests/run.sh)
[![docs](https://img.shields.io/badge/docs-online-blueviolet)](https://zhou-lab.github.io/kycg/)
[![language](https://img.shields.io/badge/C-C99-00599C)](src/kycg.h)
[![row spaces](https://img.shields.io/badge/row%20spaces-arrays%20%7C%2029M%20CpG%20genomes-brightgreen)](https://zhou-lab.github.io/kycg/)
[![vs R](https://img.shields.io/badge/stats%20vs%20R-machine%20precision-success)](tests/validate_vs_R.R)

**CpG-set enrichment as a single C binary** — fetch verified knowledgebases,
test a query set, annotate a probe table. Arrays and whole genomes, no R.

📖 **[Documentation](https://zhou-lab.github.io/kycg/)** ·
🧬 [KnowYourCG v1 (R/Bioconductor)](https://bioconductor.org/packages/knowYourCG/) ·
🔬 [Validated against it](tests/validate_vs_R.R)

</div>

KnowYourCG 2 is the second implementation of
[KnowYourCG](https://bioconductor.org/packages/knowYourCG/), built directly on
[YAME](https://github.com/zhou-lab/YAME)'s bit-packed CpG formats and run from
the command line. The binary is `kycg`.

The Bioconductor R package is KnowYourCG v1 and continues on the 1.x series;
this is v2. They are parallel implementations rather than an upgrade path —
one is an R package, the other a command-line tool — and their statistics
agree to floating-point tolerance.

**Cite:** Goldberg *et al.* KnowYourCG. *Sci Adv* 2025;11(43):eadw3027.
[doi:10.1126/sciadv.adw3027](https://doi.org/10.1126/sciadv.adw3027)

## Install

```bash
conda install -c zhou-lab -c conda-forge kycg yame
```

Install `yame` alongside kycg. kycg links libyame statically, so the backend
*library* is inside the binary — but the `yame` *command* is a separate
package, and it is what inspects a file (`yame info`). kycg has no `info`
subcommand of its own.

Knowledgebases are not in the package. `kycg fetch` pulls them into your own
store on demand, verified against digests compiled into the binary.

## Use

```bash
kycg fetch                        # browse and download knowledgebases
kycg test -m mm10:CGI query.cg    # set enrichment -> TSV on stdout
kycg annotate -m EPIC:CGI in.tsv  # label a table of probe IDs
kycg --version                    # build, coupled YAME, registry tags
```

Each subcommand takes `-h`. Everything else — what the store holds, how to
read a result, how the digests work, the case studies — is on the
[documentation page](https://zhou-lab.github.io/kycg/), one self-contained
page whose every runnable block is tested against the binary on each release.

## Build from source

```bash
git clone --recurse-submodules https://github.com/zhou-lab/kycg
cd kycg
make              # kycg, against YAME's libyame.a
make yame-bin     # the yame command itself, for `yame info`
make test         # 6 C unit tests + 9 command-level tests
```

Dependencies are YAME's: vendored htslib, zlib, libm, pthreads, librt.
libcurl is optional and used only by `kycg fetch` (`make CURL=0` forces it
off). `make install PREFIX=/usr/local` installs the binary. After a submodule
bump, `make registry` regenerates the compiled catalogue.

## License

Academic and non-profit research use is under the 2-Clause BSD License;
commercial use or transfer goes through Dr. Wanding Zhou (zhouw3@chop.edu).
Copyright (C) 2026-present The Children's Hospital of Philadelphia. See
[`LICENSE`](LICENSE) for the terms.

This is not a standard SPDX identifier, so GitHub reports the repository as
"Other" and shows no license in the sidebar. That is expected.

kycg links libyame statically from the pinned YAME submodule, which carries
the same terms as of YAME v1.53, so the distributed binary and this
repository's sources are governed alike.
