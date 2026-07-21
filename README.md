# kycg

Functional analysis of DNA methylation at CpG resolution — a C command-line
reimplementation of the [knowYourCG](https://bioconductor.org/packages/knowYourCG/)
R/Bioconductor package, with [YAME](https://github.com/zhou-lab/YAME) as its
computational backend.

**Docs site:** [`docs/index.html`](docs/index.html) — a single self-contained
page covering the whole workflow. See [DESIGN.md](DESIGN.md) for the
architecture, the statistical formulas, and the phasing. This README covers
what is implemented and how to run it.

## Status

Phase 0 (foundation) and Phase 1 (`kycg test`) are implemented and validated.
`kycg test` supersedes the knowYourCG sequencing workflow.

| Phase | Scope | State |
|---|---|---|
| 0 | Submodule, build, dispatch, `kycg info`, row-count assertion | done |
| 1 | `kycg test` — hypergeometric, group-stratified BH, effect sizes | done |
| 2 | `kycg fetch` + registry, `kycg list` | done |
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

### `kycg fetch` — build the knowledgebase store

```bash
kycg list                     # browse, check what you want, fetch it
kycg fetch mm10               # every set for a target
kycg fetch hg38:CGI,ChromHMM  # just those sets
kycg list hg38                # the individual sets (works offline)
```

`kycg list` is the interactive surface: browsing the catalogue, choosing from
it, and fetching are one activity on one screen. `→` unfolds a target, `space`
checks a set, `f` fetches everything checked, `r` refreshes the catalogues,
`d` points the browser at a different store. Sets already present show a green ✓ and cannot be checked —
there is nothing to ask for. `kycg fetch` with no target simply opens it.

Fetching does not leave the browser. The plan, the confirmation and the
progress bar render in a panel across the bottom rows while the catalogue
stays visible above; when it finishes, the counts and check marks refresh in
place with your folds and cursor where you left them. Only `q` exits.

```
    target    kind          source                 cached_sets
❯ ▾ hg38      whole genome  zenodo:18175838        4/33
   ├ [x] ABCompartment   ABCompartment.20220911.cm   9.5 KB  -
   ├  ✓  Blacklist       Blacklist.20220304.cm       3.1 KB  cached
   ├ [ ] CTCFbind        CTCFbind.20220911.cm        159 KB  -
  row 4 of 42  •  1 selected  •  → open  ← close  space select  f fetch  q quit
```

Either way, nothing downloads until a plan has been shown — which files, how
large, where they land, and what is already present — followed by a `Proceed?`
confirmation.

**Every prompt is gated on an interactive terminal.** DESIGN.md's original rule
was "never prompt", because a prompt hangs a Nextflow job or a Docker build with
no indication of why. That guarantee is preserved: off a TTY an explicit target
proceeds without asking, and a missing target is an error rather than a wait.
`-y` forces the same on a TTY. kycg still downloads in `fetch` and nowhere else.

```
    -d DIR    store directory [$KYCG_DATA_DIR, else ~/.cache/kycg]
    -o SETS   comma-separated subset, by set name (CGI,ChromHMM,TFBS)
    -n        dry run: list what would be fetched, download nothing
    -f        re-download even if present and verified
    -t TAG    InfiniumAnnotation tag, arrays only [v8]
```

Nine collections: `hg38` (32 sets) and `mm10` (29) from `KYCGKB_<genome>`,
plus `MSA` `EPICv2` `EPIC` `HM450` `HM27` `MM285` `Mammal40` from
InfiniumAnnotation. Fetched sets are ordinary files — pass one to
`kycg test -m`.

**One trust model for both channels.** Each publishes a `SHA256SUMS` at a
pinned tag; kycg pins `sha256(SHA256SUMS)` in the binary, verifies that file
after downloading it, then trusts every digest listed inside. Downloads land on
a `.part` sibling and are renamed only after their digest matches, so an
interrupted fetch cannot leave a file that later reads as valid. Afterwards the
store is re-verifiable with `shasum -a 256 -c SHA256SUMS` and no kycg code at
all. Because the anchor is the manifest rather than the file list, sets can be
added upstream without rebuilding kycg.

The Zenodo deposits remain the citable archive and keep the DOIs
([hg38](https://doi.org/10.5281/zenodo.18175837),
[mm10](https://doi.org/10.5281/zenodo.18175655)); they are recorded in the
registry as provenance and are no longer the fetch path.

**libcurl is optional.** `test`, `info`, and `list` build and run without it;
only `fetch` needs it, and it says so plainly if the build lacks it. `CURL=0`
forces it off.

### A build is coupled to a generation of the data

kycg pins a specific tag per collection and can verify only those. `--version`
says which:

```
$ kycg --version
kycg 0.1.0
knowledgebases pinned by this build:
  hg38       KYCGKB_hg38 v2
  mm10       KYCGKB_mm10 v1
  arrays     InfiniumAnnotation v8
```

Nothing updates on its own, deliberately: the digest a download is checked
against is compiled in, so a tag this build does not pin is one it cannot
verify, and `-t` on such a tag is refused rather than silently fetched.
Following an upstream tag is two commands and one generated file:

```bash
tools/make_registry.sh v9 > src/registry.h && make
tools/check_dimensions.sh          # confirm the pinned row counts still hold
```

The trade is reproducibility for immediacy — a given kycg release means an
exact, known set of knowledgebases, and data updates ride software releases.

### `kycg test` — set enrichment

```bash
kycg test -m knowledgebase.cm query.cg > results.tsv
```

For every (query sample, knowledgebase record) pair this emits the four
contingency counts, a one-sided hypergeometric tail probability in log10, six
effect-size coefficients, and a false discovery rate corrected within
knowledgebase.

`-m` is repeatable, which is what makes a store worth having — testing one
query against 30 knowledgebases in 30 processes decompresses the query 30
times, and pooling them decompresses it once:

```bash
kycg test $(for f in ~/.cache/kycg/mm10/*.cm; do printf -- '-m %s ' $f; done) \
  query.cg > res.tsv
```

Omitting `-m` on a terminal offers the store instead, **filtered to the sets
whose row count matches the query**. A `.cm` from another row space is not a
worse choice but a meaningless one — `kycg test` would refuse it on the
row-count check anyway — so the picker turns that late error into a list you
cannot pick wrong from:

```
Knowledgebases matching 21,867,837 rows (29 of 35 in the store)

  [ ] mm10/Blacklist.20220304.cm
  [x] mm10/CGI.20220904.cm
❯ [x] mm10/ChromHMM.20220414.cm
  [ ] mm10/ChromHMMfullStack.20231222.cm
  ...
  29/29 shown · 2 selected · arrows move  space toggles  a/n all/none
  / filter  enter accept  esc cancel
```

The picker and `kycg list` both render in place: arrows or `j`/`k` to move,
`/` to search, `q` to quit.

`kycg list` is a tree — `→` unfolds a target to show the sets it holds, `←`
folds it back, so the `cached_sets` count and the sets it counts are one
keystroke apart:

```
    target    kind          source                 cached_sets
  ▸ hg38      whole genome  zenodo:18175838        3/33
❯ ▾ mm10      whole genome  zenodo:18175656        29/29
    ├ CGI                CGI.20220904.cm            120 KB  cached
    ├ ChromHMM           ChromHMM.20220414.cm       857 KB  cached
    ├ PMD                PMD.20220911.cm           16.3 KB  cached
```

Whole-genome targets unfold offline, since their file list is compiled in.
Array platforms keep their file list in `SHA256SUMS` — that is what lets
upstream add a set without a kycg rebuild — so unfolding an unfetched platform
pulls that manifest (a couple of KB, verified against the compiled anchor) to
show the catalogue. That is the only implicit request kycg makes, it happens
only on a terminal, and never when stdout is redirected.

Every target carries the size of the row space it indexes, which is the fact
you need before fetching anything: a `.cm` is comparable only to a query with
the same row count, so `kycg test` will refuse a mismatch.

`cached_sets` reads *have/total* once the catalogue is known. For whole
genomes that is always, since the file list is compiled in. For an array
platform the total lives in its manifest, so it appears as soon as that is to
hand — locally, after unfolding the platform, or after `r`, which pulls them
all. Until then only the count is shown, because drawing the overview happens
on every keystroke and must never reach the network.

The picker, the browser and the tree are full-screen: they take the alternate
screen buffer, scroll a fixed-height viewport, and hand the terminal back
exactly as they found it. On a terminal that cannot support that — `NO_COLOR`,
`TERM=dumb` — they are skipped entirely in favour of a numbered prompt and
plain text carrying the same information. `kycg list` writes plain TSV whenever
stdout is redirected, so piping into `cut` or `awk` is unaffected.

```
    -m FILE   knowledgebase (.cm) [required]; repeatable
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
