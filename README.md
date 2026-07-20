# kycg

Functional analysis of DNA methylation at CpG resolution — a C command-line
reimplementation of the [knowYourCG](https://bioconductor.org/packages/knowYourCG/)
R/Bioconductor package, with [YAME](https://github.com/zhou-lab/YAME) as its
computational backend.

See [DESIGN.md](DESIGN.md) for the architecture, the statistical formulas, and
the phasing. This README covers what is implemented and how to run it.

## Status

Phase 0 (foundation) and Phase 1 (`kycg test`) are implemented and validated.
`kycg test` supersedes the knowYourCG sequencing workflow.

| Phase | Scope | State |
|---|---|---|
| 0 | Submodule, build, dispatch, `kycg info`, row-count assertion | done |
| 1 | `kycg test` — hypergeometric, group-stratified BH, effect sizes | done |
| 2 | `kycg fetch` + registry | not started |
| 3 | Plot recipes (with cinderplot) | not started |
| 4 | `proximity`, `sea`, `anno`, `bed2cg` | not started |

## Building

YAME is a submodule; the build drives its `make lib` target automatically.

```bash
git clone --recurse-submodules <this repo>
cd kycg
make
```

If you already cloned without `--recurse-submodules`:

```bash
git submodule update --init --recursive
```

Dependencies are YAME's: vendored htslib, zlib, libm, pthreads. Nothing else.

## Usage

### `kycg test` — set enrichment

```bash
kycg test -m knowledgebase.cm query.cg > results.tsv
```

For every (query sample, knowledgebase record) pair this emits the four
contingency counts, a one-sided hypergeometric tail probability in log10, six
effect-size coefficients, and a false discovery rate corrected within
knowledgebase.

```
    -m FILE   knowledgebase (.cm) to test against [required]
    -a STR    alternative: greater|less|two.sided [greater]
    -G        correct FDR globally instead of within knowledgebase
    -s FILE   sample names for the query
    -M        load the knowledgebase into memory
    -F        report full file paths instead of basenames
    -H        suppress the header line
    -o FILE   write to FILE instead of stdout
```

Output columns:

```
query_file  query  db_file  db  group  nU  nQ  nD  overlap
estimate  log10_p  p_value  fdr  neglog10_fdr
cf_jaccard  cf_mcc  cf_overlap  cf_npmi  cf_dice
beta  depth
```

The first eight are `yame summary`'s own columns, renamed; everything from
`estimate` onward is what kycg adds. `Log2OddsRatio` becomes `estimate` and
gains full precision — `yame` prints it at `%1.2f`, which is fine for eyeballing
and too coarse to sort or plot by.

**FDR is corrected within one query sample and one knowledgebase.** That is the
family knowYourCG corrects over: `testEnrichment()` takes a single query and
calls `p.adjust()` on that query's result frame, within `group`. Passing more
samples or more knowledgebases in one invocation therefore does not change any
FDR already reported. `-G` widens the stratum to all knowledgebases for a
sample, matching `mtc_by_group=FALSE`.

`neglog10_fdr` is precomputed because cinderplot has no `-log10()` transform in
aesthetics and it is the y-axis of three separate plots.

> **`p_value` and `fdr` are convenience columns and will read `0` for any
> strong result.** They are recovered as `10^log10_p`, which underflows in
> double precision somewhere around `log10_p = -308`; real enrichments run to
> `-1500` and beyond. Sort, threshold, and plot on `log10_p` and
> `neglog10_fdr`, which are carried in log space end to end.

### `kycg info` — describe a file

```bash
$ kycg info onecell.cg chromhmm.cm
file          record  name  format         n_rows    n_set
onecell.cg    1       NA    3:mu           21867837  NA
chromhmm.cm   1       NA    2:categorical  21867837  NA
```

Use this to confirm two files index the same row space before testing them.

### Plotting

kycg emits TSV; [cinderplot](https://github.com/zhou-lab/cinderplot) renders it
as a separate process. There is no linkage between them.

```bash
kycg test -m ChromHMM.cm query.cg > res.tsv
cinderplot 'res.tsv + aes(estimate, neglog10_fdr) + geom_point()' -o volcano.pdf
```

## Row spaces

Everything here is positional: a CpG has no identity beyond its index into a
reference row list. A `.cg` and a `.cm` are comparable only if row *i* means
the same CpG in both.

kycg asserts that the query and knowledgebase record lengths agree and fails
loudly otherwise:

```
[run_pair:148] Row count mismatch: query 'small.cg' record '1' has 5 rows but
knowledgebase 'chromhmm.cm' record '1' has 21867837 rows.
These files index different reference row lists and cannot be compared.
```

It does no more than that, by design (DESIGN.md §2). kycg does not infer
platforms and cannot detect two files that share a row count but come from
different row spaces — sequencing uses a whole-genome `.cr`, arrays use a
per-platform `ordering.tsv.gz`, and a `.cm` from one is meaningless in the
other. Keeping them matched is the user's responsibility.

## Correctness

The statistics are validated at two levels.

**Unit** — `make test` checks the hypergeometric tail and the BH correction
against values produced by R's own `phyper()` and `p.adjust()`, across small,
array-scale (~486K rows), whole-genome-scale (29M rows), and deep-tail
(`log10 p ≈ -584141`) regimes. Worst observed relative error: **2.6e-16**.

**Unit** — `tests/test_enrich.c` covers FDR stratification and the effect-size
clamping edge cases (empty overlap, perfect containment, empty knowledgebase).

**End to end** — `tests/validate_vs_R.R` re-derives every column of a real
`kycg test` run using the same R calls knowYourCG makes, and compares:

```bash
$ kycg test -m chromhmm.cm onecell.cg > res.tsv
$ Rscript tests/validate_vs_R.R res.tsv
Comparing 19 rows against R
  log10_p        OK     worst abs diff 4.09e-12
  estimate       OK     worst abs diff 4.75e-10
  log10 fdr      OK     worst abs diff 2.96e-12
  ...
PASS: kycg agrees with R on every column.
```

DESIGN.md sets ~1e-10 agreement on `log10.p.value` as the Phase 1 acceptance
bar; the residual above is dominated by output text precision, not by the
computation. The counts themselves are byte-identical to `yame summary`.

Regenerate the unit-test reference values with `Rscript tests/ref.R`.

### Two implementation notes

**Everything is log space.** knowYourCG's `p.value` column underflows to 0 for
exactly the results users care about, and `log10.p.value` is its primary sort
key throughout. kycg computes the tail in log space and applies BH in log space
too — a linear-space BH would assign FDR 0 to every underflowed row.

**The point mass uses Loader's saddle-point expansion**, not a difference of
`lgamma` calls, which is also what R's `dhyper()` does internally. Over a 29M-row
universe the naive form loses about seven digits to cancellation — coarser than
the acceptance bar. htslib's `kt_fisher_exact()` (vendored in YAME) is not usable
here for a separate reason: it takes `int` counts and accumulates in linear
space, looping "until underflow" by construction.

## Layout

```
src/hypergeo.{c,h}   log-space hypergeometric tail + BH   (pure, no YAME)
src/enrich.{c,h}     effect sizes, FDR strata, ordering, TSV emission
src/test.c           kycg test
src/info.c           kycg info
src/main.c           subcommand dispatch (the only main())
tests/               unit tests + R cross-validation
external/YAME/       submodule
```

`hypergeo.c` and `enrich.c` link without YAME, which is what lets the test
binary build against them directly.

## License

AGPL-3.0-or-later, matching YAME.
