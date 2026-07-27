# kycg

Functional analysis of DNA methylation at CpG resolution — a C command-line
reimplementation of the [knowYourCG](https://bioconductor.org/packages/knowYourCG/)
R/Bioconductor package, with [YAME](https://github.com/zhou-lab/YAME) as its
computational backend.

**Cite:** Goldberg *et al.* KnowYourCG. *Sci Adv* 2025;11(43):eadw3027.
[doi:10.1126/sciadv.adw3027](https://doi.org/10.1126/sciadv.adw3027)

**Docs site:** [`docs/index.html`](docs/index.html) — a single self-contained
page covering the whole workflow. This README covers what is implemented and
how to run it.

## Status

Phases 0-2 are implemented and validated. `kycg test` supersedes the
knowYourCG sequencing workflow.

| Phase | Scope | State |
|---|---|---|
| 0 | Submodule, build, dispatch, row-count assertion | done |
| 1 | `kycg test` — hypergeometric, group-stratified BH, effect sizes | done |
| 2 | `kycg fetch` — registry, catalogue browser, verified store | done |
| 3 | Plot recipes (with cinderplot) | not started |
| 4 | `annotate` (done), `proximity`, `sea`, `bed2cg` | partial |

## Installing

```bash
conda install -c zhou-lab -c conda-forge kycg
```

Knowledgebases are not in the package — `kycg fetch` pulls them into your own
store on demand. The compendium is hundreds of megabytes and versioned
separately from the software, so bundling it would turn every data update into
a software release. What the binary does carry is the registry of tags and
digests it can verify; `kycg --version` prints it.

See [`conda-recipe/`](conda-recipe/) to build the package yourself.

## Building from source

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
libcurl is optional and used only by `kycg fetch`; `make CURL=0` forces it off,
and `kycg --version` reports whether a build has it.

`make install PREFIX=/usr/local` installs the binary.

## Usage

### `kycg fetch` — browse and build the knowledgebase store

```bash
kycg fetch                    # browse everything
kycg fetch hg38               # open with hg38 checked; press f to start
kycg fetch hg38:CGI,ChromHMM  # open with just those checked
kycg fetch -f hg38            # download it, no browser, no questions
kycg fetch | grep -v '^#'     # TSV catalogue, for scripts
```

There is one interactive path and one scripted one. On a terminal, naming a
target opens the catalogue with that target already checked — you see exactly
what will be downloaded and how large it is, can narrow it in the same screen,
and press `f` to start. `-f` skips all of that, and so does a redirected
stdout, which is what keeps scripts working unchanged. The redirected form
leads with a `# store: <path>` comment naming the store it described, so
strip comment lines before parsing. `→` unfolds a target, `space`
checks a set, `f` fetches everything checked, `d` points the browser at a different store. Sets already present show a green ✓ and cannot be checked —
there is nothing to ask for. `kycg fetch` with no target simply opens it.

Two keys answer the questions the columns cannot. `r` checks that collection's
recommended selection in one keypress — the dozen or so sets worth having
before you know what you are looking for, rather than all forty. `i` describes
the set under the cursor: what the annotation means for methylation, the
upstream database, the publication, and what was done to it on the way in.
The pane is open by default and follows the cursor, so arrowing down walks the
catalogue with each set explained as you reach it. It collapses on collection
rows, where the columns already say what there is to say. `i` hides it when you
want the full screen for scanning.
That last one is the answer to "is `TFBSrm` the same as `TFBS`, differently
filtered?" — it is not; `rm` is ReMap. The text lives in
[`data/knowledgebases.tsv`](data/knowledgebases.tsv) and is compiled into the
binary, so both keys work with no network and no data files to find. Fields
nobody has been able to establish read `not recorded` rather than being hidden.

Fetching does not leave the browser. The plan, the confirmation and the
progress bar render in a panel across the bottom rows while the catalogue
stays visible above; when it finishes, the counts and check marks refresh in
place with your folds and cursor where you left them. Only `q` exits.

```
    target    kind          rows        source           cached_sets
❯ ▾ hg38      whole genome  29,401,795  KYCGKB_hg38 v2    4/32
   ├ [x] ABCompartment   ABCompartment.20220911.cm   9.5 KB  -
   ├  ✓  Blacklist       Blacklist.20220304.cm       3.1 KB  cached
   ├ [ ] CTCFbind        CTCFbind.20220911.cm        159 KB  -
  row 4 of 42  •  1 selected  •  → open  ← close  space select  r recommended  i hide  f fetch  d store   q quit
```

**Nothing is asked when nobody can answer.** The original rule was
"never prompt", because a prompt hangs a Nextflow job or a Docker build with no
indication of why. That guarantee holds: off a TTY a named target downloads
without asking, and a missing target prints the catalogue rather than waiting.
kycg still downloads in `fetch` and nowhere else.

```
    -d DIR    store directory [$YAME_DATA_HOME, else the shared store]
    -o SETS   comma-separated subset, by set name (CGI,ChromHMM,TFBS)
    -f        download now: no browser, no questions
    -r        re-download even what is present and verified
    -t TAG    InfiniumAnnotation tag, arrays only [v8]
```

Nine collections: `hg38` (32 sets) and `mm10` (28) from `KYCGKB_<genome>`,
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

**libcurl is optional.** `test` and `info` build and run without it;
only `fetch` needs it, and it says so plainly if the build lacks it. `CURL=0`
forces it off.

### A build is coupled to a generation of the data

kycg pins a specific tag per collection, compiled into `src/registry.h`, and
can verify only those. `--version` reports the build and the coupled YAME:

```
$ kycg --version
kycg 0.4
    built against  YAME v1.33
    store          ~/.local/share/yame   ($YAME_DATA_HOME unset; -d overrides)
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

### `kycg annotate` — label a table of probe IDs

```bash
kycg annotate -m EPIC:CGI,ProbeType hits.tsv > annotated.tsv
```

The lookup half of the same data `test` aggregates over: testing asks whether a
group of probes is enriched somewhere, annotating asks where one probe actually
is. It reads a TSV, keeps every column and row in order, and appends one column
per knowledgebase.

```
Probe_ID    beta  CGI      ProbeType
cg00000029  0.5   Shore    cg
cg00000103  0.5   OpenSea  cg
cg99999999  0.1   NA       NA
```

`-i` gives an indicator matrix instead — one 0/1 column per set, which is what
you want feeding a model rather than reading.

Omitting `-m` on a terminal opens the same browser `fetch` and `test` use,
listing the array platforms: check what you want, `f` to fetch anything
missing, `a` to annotate. Off a terminal it stays an error rather than a wait.

Naming a set that exists but is not downloaded opens that browser too, with the
set already checked — `f` to fetch it, `q` to carry on. A name that matches
nothing is a typo and stays an error. Both `annotate` and `test` behave this
way, and both refuse rather than prompt off a terminal.

```
    -m SPEC   path or platform[:sets]; repeatable [required]
    -p PLAT   platform, when -m is a plain path
    -c COL    probe ID column: a name or 1-based index [Probe_ID, else 1]
    -i        indicator columns (0/1) per set
    -s SEP    separator when a probe is in several sets [,]
    -H        the input has no header line
    -o FILE   write to FILE instead of stdout
```

**Arrays only, and the platform is named rather than inferred.** A probe ID
means nothing without the ordering that gives it a row, so `annotate` reads the
platform's `ordering.tsv.gz` — fetched with any set for that platform — and
binary-searches it. knowYourCG guesses the platform from probe ID patterns;
kycg does not, for the same reason it never guesses row spaces: a wrong guess
produces a full table of confident, wrong answers.

A probe absent from the ordering and a probe in no set both read `NA`, which
are different facts, so the count of the first goes to stderr where it cannot
corrupt the table.

### `kycg test` — set enrichment

```bash
kycg test -m mm10:CGI query.cg > results.tsv
```

`-m` takes a path, or a set named the way `kycg fetch` names it — so a
knowledgebase is written the same way whether you are downloading it or
testing against it:

```bash
kycg test -m mm10:CGI           query.cg    # one set from the store
kycg test -m mm10:CGI,ChromHMM  query.cg    # two of them
kycg test -m mm10               query.cg    # everything cached for mm10
kycg test -m path/to/some.cm    query.cg    # a path still works
```

A spec that is an existing file is taken literally, so plain paths are
unaffected and nothing is ambiguous.

For every (query sample, knowledgebase record) pair this emits the four
contingency counts, a one-sided hypergeometric tail probability in log10, six
effect-size coefficients, and a false discovery rate corrected within
knowledgebase.

`-m` is repeatable, and `-m mm10` expands to everything cached for that
target — which is what makes a store worth having. Testing one query against
30 knowledgebases in 30 processes decompresses the query 30 times; pooling
them decompresses it once.

Omitting `-m` on a terminal opens the same browser `kycg fetch` uses, showing
only the collections whose row count matches the query — a `.cm` from another
row space is not a worse choice but a meaningless one, and `kycg test` would
refuse it anyway, so the list is one you cannot pick wrong from:

```
Knowledgebases for 21,867,837 rows -- space to choose, t to test
    target  kind          rows        cached_sets
❯ ▾ mm10    whole genome  21,867,837  28
   ├ [x] CGI          CGI.20220904.cm            cached
   ├ [ ] EvoCons      EvoCons.20220314.cm        -
  row 3 of 30 · 1 selected · → open  ← close  space select  r recommended  i hide  f fetch  t test  q quit
```

It lists everything a collection publishes, not just what you have — so if the
set you want is missing, check it, press `f` to fetch it, then `t` to test
against it, without leaving. Targets are filtered by the row counts pinned in
the registry, so this is a table lookup rather than a scan of the store.

The picker and the catalogue both render in place: arrows or `j`/`k` to move,
`/` to search, `q` to quit.

The catalogue is a tree — `→` unfolds a target to show the sets it holds, `←`
folds it back, so the `cached_sets` count and the sets it counts are one
keystroke apart:

```
    target    kind          rows        source           cached_sets
  ▸ hg38      whole genome  29,401,795  KYCGKB_hg38 v2    3/32
❯ ▾ mm10      whole genome  21,867,837  KYCGKB_mm10 v2   28/28
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
hand — locally, or after unfolding the platform, which pulls it. Until then
only the count is shown, because drawing the overview happens on every
keystroke and must never reach the network.

The picker, the browser and the tree are full-screen: they take the alternate
screen buffer, scroll a fixed-height viewport, and hand the terminal back
exactly as they found it. On a terminal that cannot support that — `NO_COLOR`,
`TERM=dumb` — they are skipped entirely in favour of a numbered prompt and
plain text carrying the same information. `kycg fetch` writes plain TSV whenever
stdout is redirected, so piping into `cut` or `awk` is unaffected.

```
    -m SPEC   path or target[:sets]; repeatable [required]
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

### Inspecting a file

kycg builds on YAME's `.cx` formats, so file inspection is `yame info` — it
describes a `.cg` / `.cm`'s format, row count, and (for a categorical set) its
state keys:

```bash
$ yame info onecell.cg chromhmm.cm
```

That is how you confirm two files index the same row space before testing.
`kycg test` also fails loudly on any row-count mismatch, naming both files, so
the check is a convenience rather than a guardrail.

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

It does no more than that, by design. kycg does not infer platforms and
cannot detect two files that share a row count but come from different row
spaces — sequencing uses a whole-genome `.cr`, arrays use a
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

The acceptance bar is ~1e-10 agreement on `log10.p.value`; the residual above
is dominated by output text precision, not by the computation. The counts
themselves are byte-identical to `yame summary`.

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
src/args.{c,h}       argv permutation so options may follow operands
src/store.{c,h}      store paths, enumeration, member-name safety
src/digest.{c,h}     sha256
src/registry.h       generated: pinned tags, anchors, row counts, sizes
src/kbinfo.h         generated: per-set provenance from data/knowledgebases.tsv
src/ui.{c,h}         terminal layer: tree browser, panels, progress, TTY gating
src/fetch.c          kycg fetch — catalogue, plan, verified download
src/test.c           kycg test
src/main.c           subcommand dispatch (the only main())
data/                knowledgebase metadata (source of truth for kbinfo.h)
tools/               generators for registry.h and kbinfo.h, dimension checks
conda-recipe/        conda package
tests/               unit tests + R cross-validation
external/YAME/       submodule
```

`hypergeo.c` and `enrich.c` link without YAME, which is what lets the test
binary build against them directly.

## License

AGPL-3.0-or-later, matching YAME.
