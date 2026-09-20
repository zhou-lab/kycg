# The test suite

`make test` builds everything and runs `tests/run.sh`, which is the whole
suite: the C unit binaries plus every `t_*.sh`. The conda recipe runs the same
`tests/run.sh` in its `test:` phase, against the **installed** package on every
platform it builds, so what is verified is what a user gets.

```sh
make test                                   # everything
KYCG=./kycg YAME=external/YAME/yame bash tests/t_annotate.sh   # one test
scripts/coverage.sh                         # line coverage + the badge
```

## The two layers

**C unit tests** (`test_*.c`) — the pure functions, linked against
`hypergeo.o enrich.o args.o store.o digest.o`. `test_hypergeo.c` checks the
log-space tail and BH against values `tests/ref.R` generated from R's own
`phyper` / `p.adjust`, from small counts to a 29M-row deep tail;
`test_digest.c` checks SHA-256 and MD5 against the published FIPS/RFC vectors
rather than against themselves; `test_store.c` states the traversal attack the
manifest-name guard exists to stop, and `test_storewalk.c` the enumeration
rules (a symlinked file is kept, a symlinked directory is not descended --
a duplicate `.cm` would inflate BH's m and shift every FDR in its stratum).
They run only in a build tree: in a conda test environment the binary under
test is the installed one, and a `test_*` copied along as a source file was
compiled on somebody else's machine.

**Command-level tests** (`t_*.sh`) — the built binary, end to end. Each is
self-contained: it packs its own fixtures from inline text with `yame pack`,
makes its own temp directory, touches nothing shared, and needs no network and
no data store. Fixtures are packed rather than committed because kycg reads
what YAME writes — a committed `.cg` would quietly become a format the current
YAME no longer produces, and the test would then be checking a museum piece.

| test | what it pins |
|---|---|
| `t_usage` | dispatch, `-h`/`--help`, `--version`, `$YAME_DATA_HOME`, no ANSI off a terminal |
| `t_test_stats` | the 2x2 counts, effect sizes and odds ratio recounted in awk; the hypergeometric tail computed EXACTLY in integer arithmetic on a 40-row fixture; the three alternatives; `beta`/`depth` from an M/U query |
| `t_test_options` | `-H -o -F -M -G -s`, repeated `-m`, several queries, the group column |
| `t_rowspace` | the row-count mismatch message, missing inputs, junk and empty stores |
| `t_annotate` | a synthetic 27,722-row HM27 store: membership, `-i -c -s -H`, a set named `HM27:CGI`, and the four refusals (unknown platform, no `-p`, unsorted ordering, wrong-length ordering) |
| `t_fetch_offline` | the TSV catalogue off a terminal, `-d`, and that nothing is downloaded |
| `t_fetch_http` | a loopback mirror: download and verify, skip, `-r`, 404, corrupted bytes, a tampered manifest, an unpinned tag, the companion, the "published but not here" offer, and the same through the browser |
| `t_browser` | the tree and the pickers through a pty: every movement key, the detail pane, ESC, `TERM=dumb`, quitting a picker |
| `t_docs` | one version in binary/header/recipe, every advertised subcommand dispatches, every install path that says `yame info` also installs `yame`, documented columns exist |

`lib.sh` carries the assertions and the pack recipes. Nothing else is shared;
a little duplication between tests is deliberate, so each reads top to bottom.

## The documented workflows — `make test-docs`

A third layer, outside `make test`. Every runnable block on `docs/index.html`
is a script under `docs/examples/`, and `tests/docs_gate.py` runs them all in
order in one sandbox, against this checkout's binaries — the proof that what a
reader copies off the page still works.

```sh
make test-docs                                  # ~25s cold, <10s warm
KYCG_DOCS_BIN=~/.conda/envs/kycgtest/bin make test-docs   # test an installed copy
```

It is deliberately **not** part of `make test` and never runs in CI or a conda
build: it needs the network, a `git clone` and cinderplot. The sandbox
persists (`KYCG_DOCS_SANDBOX`, default `~/tmp/kycg/docs_gate`) so a second run
downloads nothing; everything the examples *built* is cleared first, so a
stale artifact can never stand in for a command that stopped working.

Each script carries `## title:`, and one the reader should not run carries
`## norun: <reason>` — the gate skips it and says why. The page declares the
same thing per block: `data-example="NN_name"` names the script that runs it,
`data-norun="<reason>"` marks an illustration. `t_docs.sh` fails if a block
declares neither or names a script that is not there, so a new command cannot
reach the page ungated without an explicit reason.

Two limits, both deliberate. The gate checks **exit codes**, not output — every
pasted result on the page is prose (release SOP 4c). And a `norun` block's
contents are checked by nothing, which is why they rot; the SOP says to read
them by hand.

## Coverage

`scripts/coverage.sh` measures gcov line coverage of the suite over `src/`
in a scratch copy (YAME is symlinked, not rebuilt) and writes
`docs/coverage.json`, the shields endpoint behind the README badge. CI runs
`--check`, which fails when the badge and the measurement differ by more than
two points in either direction.

**84.6% as of 2026-09-18**, from 11.2% of lines even being reachable before
the `t_*.sh` layer existed. Where the rest is:

| file | lines | cov | why |
|---|---:|---:|---|
| `args.c`, `main.c`, `store.c` | 122 | 100% | |
| `hypergeo.c` | 162 | 98% | the deep-tail fallback |
| `enrich.c` | 113 | 97% | |
| `digest.c` | 77 | 96% | |
| `test.c` | 256 | 93% | |
| `annotate.c` | 325 | 85% | allocation-failure paths |
| `fetch.c` | 1309 | 84% | redraw-after-refresh inside the browser |
| `ui.c` | 719 | 74% | terminal-capability fallbacks, resize and signal paths |

What is still dark is mostly unreachable by design: malloc failure, terminals
that cannot host a widget, and redraws that need a resize. The pty harness
(`pty_drive` in `lib.sh`) and the loopback mirror are the tools for anything
that is not.

Two rules keep the number honest. **A release may not lower it** (see the
release SOP), and code no command can reach is deleted rather than tested:
writing this suite found 377 lines of unreachable widgets in `ui.c`
(`kycg_ui_choose`, `kycg_ui_multiselect`, `kycg_ui_confirm`, `kycg_ui_ask`)
and a dead group-naming helper in `enrich.c`. Removing them is why the
denominator fell from 3,335 lines to 3,105.
