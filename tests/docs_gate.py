#!/usr/bin/env python3
"""Run every runnable docs/examples/*.sh, in order, as a reader would.

This is the documented-workflow gate: the proof that the page's workflows run
on the binary that ships. It is NOT part of `make test` -- it needs the
network (about 60 MB of knowledgebases the first time), a git clone, and a few
minutes -- and it must never run in CI or a conda build. `make test-docs` runs
it; the release SOP runs that under sbatch.

Each script runs with `bash -eo pipefail` from one working directory, so a
block can use what an earlier block built, exactly as on the page. Scripts
with a `## norun:` line are skipped and reported.

The sandbox PERSISTS between runs (KYCG_DOCS_SANDBOX, default
~/tmp/kycg/docs_gate): the store and the kycg-examples clone stay, so a second
run downloads nothing and only re-verifies digests. Everything the examples
BUILT is deleted before a run, so "output exists" can never be the reason a
block passes.

PATH puts the working tree first -- ./kycg and external/YAME/yame -- so the
gate tests the code in this checkout, never an installed copy. Point
KYCG_DOCS_BIN at a DIRECTORY (a conda env's bin/) to test what a reader
installs instead.
"""
import glob, os, shutil, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES = os.path.join(ROOT, "docs", "examples")
SANDBOX = os.environ.get("KYCG_DOCS_SANDBOX") or os.path.expanduser("~/tmp/kycg/docs_gate")
WORK = os.path.join(SANDBOX, "work")     # where every example runs
STORE = os.path.join(SANDBOX, "store")   # $YAME_DATA_HOME, as a reader's would be
LOGS = os.path.join(SANDBOX, "logs")
BIN = os.environ.get("KYCG_DOCS_BIN") or (ROOT + ":" + os.path.join(ROOT, "external", "YAME"))

# Inputs, not outputs: a fetched store and a cloned examples repo are what the
# reader already has by the time the later blocks run. Everything else in WORK
# is something an example built, and is cleared so a stale artifact can never
# stand in for a command that no longer works.
KEEP = {"kycg-examples"}


def header(path):
    title, norun = "", None
    for line in open(path, encoding="utf-8"):
        if line.startswith("## title:"):
            title = line[9:].strip()
        elif line.startswith("## norun:"):
            norun = line[9:].strip()
        elif not line.startswith("#"):
            break
    return title, norun


def clear_outputs():
    for name in os.listdir(WORK):
        if name in KEEP:
            continue
        p = os.path.join(WORK, name)
        shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)


def main():
    for d in (WORK, STORE, LOGS):
        os.makedirs(d, exist_ok=True)
    env = dict(os.environ,
               HOME=os.path.join(SANDBOX, "home"),
               YAME_DATA_HOME=STORE,
               PATH=BIN + ":" + os.environ["PATH"])
    os.makedirs(env["HOME"], exist_ok=True)
    clear_outputs()

    scripts = sorted(glob.glob(os.path.join(EXAMPLES, "*.sh")))
    if not scripts:
        sys.exit("docs_gate: no scripts in %s" % EXAMPLES)

    results, t0 = [], time.time()
    for path in scripts:
        name = os.path.basename(path)[:-3]
        title, norun = header(path)
        if norun:
            results.append((name, "SKIP", norun))
            print("%-24s SKIP   %s" % (name, norun), flush=True)
            continue
        with open(os.path.join(LOGS, name + ".log"), "w") as log:
            t1 = time.time()
            rc = subprocess.run(["bash", "-eo", "pipefail", path], cwd=WORK,
                                env=env, stdout=log, stderr=subprocess.STDOUT,
                                timeout=7200).returncode
        state = "ok" if rc == 0 else "FAIL rc=%d" % rc
        results.append((name, state, title))
        print("%-24s %-10s %s  (%ds)" % (name, state, title, time.time() - t1),
              flush=True)

    failed = [r for r in results if r[1].startswith("FAIL")]
    skipped = [r for r in results if r[1] == "SKIP"]
    print("\n%d blocks, %d failed, %d skipped  (%ds; sandbox %s)"
          % (len(results), len(failed), len(skipped), time.time() - t0, SANDBOX))
    for name, state, _ in failed:
        print("  %s: %s -- see %s/%s.log" % (name, state, LOGS, name))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
