#!/bin/bash
## `kycg fetch` against a loopback mirror: the download half of fetch.c.
##
## Nothing here touches the network. YAME_ASSETS_MIRROR replaces the scheme
## and host of every URL and keeps the path, so a local server standing at
## mirror/zhou-lab/InfiniumAnnotation/<tag>/... is served exactly what
## github.com would be asked for. Verification is untouched: the manifests in
## tests/fixtures are the REAL published ones, whose sha256 is the anchor
## compiled into this build, and the one .cm fixture is the real file whose
## digest that manifest lists. A mirror can therefore fail to serve the right
## bytes but cannot pass wrong ones -- which is what makes the failure cases
## below meaningful rather than staged.
##
## The server reads its behaviour from a file on every request, so one server
## covers success, 404, corruption and a tampered manifest.
set -uo pipefail
. "$(dirname "$0")/lib.sh"
here=$(cd "$(dirname "$0")" && pwd)

command -v python3 >/dev/null || { echo "  skip: no python3 for the mirror"; exit 0; }
[ -d "$here/fixtures/mirror" ] || { echo "  skip: no mirror fixtures"; exit 0; }

new_workdir

## ---- the mirror tree ---------------------------------------------------- ##
## The tag follows the submoduled YAME, never a literal here. src/registry.h
## is projected from that submodule's tools/registry/files.tsv, so the tag the
## binary requests moves with every catalog bump -- a pinned v8.1 turns the
## NEXT bump into seven mirror 404s that read as a fetch bug. --version prints
## what THIS binary compiled, which is exactly the tag it will ask for.
TAG=$("$KYCG" --version | grep -oE 'InfiniumAnnotation@v[0-9.]+' | cut -d@ -f2)
[ -n "$TAG" ] || { echo "  skip: --version reports no InfiniumAnnotation tag"; exit 0; }
## Since YAME v1.50 the URL rule is raw.githubusercontent.com/<org>/<repo>/
## <tag>/<path>, and YAME_ASSETS_MIRROR keeps everything after the host.
tree=mirror/zhou-lab/InfiniumAnnotation/$TAG
mkdir -p "$tree"
cp -r "$here/fixtures/mirror/EPIC" "$tree/"
SET=Blacklist.20220304.cm
SETSHA=3482751b12d943a0ea42ceba244c7c9abdafeccea2492a00026f680bce14fce5

cat > srv.py <<'PY'
import os, sys
from http.server import BaseHTTPRequestHandler, HTTPServer
port, root = int(sys.argv[1]), sys.argv[2]
class H(BaseHTTPRequestHandler):
    def do_GET(s):
        mode = open("mode").read().strip()
        f = root + s.path
        ## The probe ordering is a planned item on every EPIC fetch (it is
        ## what makes a set interpretable) and the mirror does not carry it,
        ## so every run here exercises the partial-failure path too.
        if s.path.endswith(".ordering.tsv.gz"):
            s.send_response(404); s.end_headers(); return
        if s.path.endswith("Blacklist.20220304.cm"):
            if mode == "404":
                s.send_response(404); s.end_headers(); return
            if mode == "corrupt":
                b = b"not the bytes the registry pinned\n"
                s.send_response(200); s.send_header("Content-Length", str(len(b)))
                s.end_headers(); s.wfile.write(b); return
        if not os.path.isfile(f):
            s.send_response(404); s.end_headers(); return
        b = open(f, "rb").read()
        s.send_response(200); s.send_header("Content-Length", str(len(b)))
        s.end_headers(); s.wfile.write(b)
    def log_message(s, *a): pass
HTTPServer(("127.0.0.1", port), H).serve_forever()
PY

port=$(python3 -c 'import socket; s=socket.socket(); s.bind(("127.0.0.1",0)); print(s.getsockname()[1])')
echo ok > mode
python3 srv.py "$port" "$PWD/mirror" >server.log 2>&1 &
srv=$!
## The workdir trap from lib.sh removes the directory; add the server to it.
trap 'kill $srv 2>/dev/null; rm -rf "$d"' EXIT
for i in $(seq 1 50); do
  python3 -c "import socket; socket.create_connection(('127.0.0.1', $port), timeout=1)" 2>/dev/null && break
  sleep 0.2
done

export YAME_ASSETS_MIRROR="http://127.0.0.1:$port"
export YAME_DATA_HOME="$PWD/store"

## A build without libcurl cannot do any of this, and says so.
probe=$("$KYCG" fetch -f EPIC:Blacklist </dev/null 2>&1)
case $probe in
  *"no network support"*) echo "  skip: built without libcurl"; exit 0;;
esac

sha_of() { python3 -c "import hashlib,sys; print(hashlib.sha256(open(sys.argv[1],'rb').read()).hexdigest())" "$1"; }
dest="$YAME_DATA_HOME/EPIC/KYCG/$SET"

## ---- 1. the happy path (the probe above already ran it) ----------------- ##
check_has "the fetch names the set" "Blacklist" "$probe"
[ -f "$dest" ] || { echo "  FAIL nothing landed at $dest"; sed 's/^/       /' <<<"$probe"
                    fails=$((fails + 1)); }
if [ -f "$dest" ]; then
  check "the fetched bytes verify" "$SETSHA" "$(sha_of "$dest")"
fi
## The manifest is written into the store, which is what makes the next run
## able to answer "what do I have" without the network.
check "the manifest is kept in the store" 1 \
      "$([ -f "$YAME_DATA_HOME/EPIC/KYCG/SHA256SUMS" ] && echo 1 || echo 0)"
## A download that failed half way leaves a .part; a finished one must not.
check "no .part files are left" 0 \
      "$(find "$YAME_DATA_HOME" -name '*.part' | wc -l | tr -d ' ')"

## ---- 2. a second run skips what is present and verified ----------------- ##
## Not asserted on the exit code: the probe ordering rides along with every
## EPIC fetch and this mirror does not carry it (case 8), so every run here
## ends non-zero by design. What matters is that the SET was not re-fetched --
## the inode is the evidence.
before=$(stat -c %i "$dest" 2>/dev/null || stat -f %i "$dest")
out=$("$KYCG" fetch -f EPIC:Blacklist </dev/null 2>&1)
check_has "a redundant fetch reports it as current" "already current" "$out"
after=$(stat -c %i "$dest" 2>/dev/null || stat -f %i "$dest")
check "the file was not rewritten" "$before" "$after"

## ---- 3. -r re-downloads it anyway --------------------------------------- ##
out=$("$KYCG" fetch -f -r EPIC:Blacklist </dev/null 2>&1)
check "-r still verifies" "$SETSHA" "$(sha_of "$dest")"
## -r means fetch it again even though it is current, so the run must NOT
## report it as already current the way case 2 does.
check_lacks "-r re-fetches rather than skipping" "already current" "$out"

## ---- 3b. stale .part files are swept, live ones are not ----------------- ##
## A download lands on a per-pid ".part" and is renamed once it verifies, so a
## killed fetch leaves one behind. The next fetch into that directory clears
## the ones older than a day and leaves anything newer alone -- another
## process may be part way through writing it, and deleting that would turn a
## concurrent fetch into a corrupt file. mtime is the only thing separating
## the two cases, so the test sets it rather than waiting a day.
partdir="$YAME_DATA_HOME/EPIC/KYCG"
: > "$partdir/stale.part"
: > "$partdir/fresh.part"
touch -d '3 days ago' "$partdir/stale.part" 2>/dev/null ||
  touch -t "$(date -v-3d +%Y%m%d%H%M 2>/dev/null)" "$partdir/stale.part"
"$KYCG" fetch -f -r EPIC:Blacklist </dev/null >/dev/null 2>&1
check "a day-old .part is swept" 0 \
      "$([ -e "$partdir/stale.part" ] && echo 1 || echo 0)"
check "a .part being written now is left alone" 1 \
      "$([ -e "$partdir/fresh.part" ] && echo 1 || echo 0)"
/bin/rm -f "$partdir/fresh.part"

## ---- 4. the catalogue counts what is cached, offline -------------------- ##
## No mirror needed: this is read from the store and the compiled set counts.
cat=$(YAME_ASSETS_MIRROR="http://127.0.0.1:1" "$KYCG" fetch </dev/null 2>/dev/null)
line=$(printf '%s\n' "$cat" | grep '^EPIC	')
check_has "the catalogue counts the cached set" "1/17" "$line"

## ---- 5. a 404 is an error, and leaves nothing behind -------------------- ##
/bin/rm -rf "$YAME_DATA_HOME"
echo 404 > mode
out=$("$KYCG" fetch -f EPIC:Blacklist </dev/null 2>&1); rc=$?
awk -v rc="$rc" 'BEGIN{exit !(rc != 0)}' ||
  { echo "  FAIL a 404 fetch reported success"; fails=$((fails + 1)); }
check_has "the 404 is reported" "404" "$out"
check "nothing landed after a 404" 0 \
      "$(find "$YAME_DATA_HOME" -name "$SET" 2>/dev/null | wc -l | tr -d ' ')"
check "no .part after a 404" 0 \
      "$(find "$YAME_DATA_HOME" -name '*.part' 2>/dev/null | wc -l | tr -d ' ')"

## ---- 6. bytes that do not match the digest are refused ------------------ ##
## The case the whole anchor chain exists for: the server is reachable, the
## manifest is genuine, and the file is wrong.
/bin/rm -rf "$YAME_DATA_HOME"
echo corrupt > mode
out=$("$KYCG" fetch -f EPIC:Blacklist </dev/null 2>&1); rc=$?
awk -v rc="$rc" 'BEGIN{exit !(rc != 0)}' ||
  { echo "  FAIL corrupted bytes were accepted"; fails=$((fails + 1)); }
printf '%s\n' "$out" | grep -qiE 'sha|digest|checksum|verif' ||
  { echo "  FAIL the corruption was not described:"; sed 's/^/       /' <<<"$out"
    fails=$((fails + 1)); }
check "corrupted bytes are not kept" 0 \
      "$(find "$YAME_DATA_HOME" -name "$SET" 2>/dev/null | wc -l | tr -d ' ')"
check "no .part after corruption" 0 \
      "$(find "$YAME_DATA_HOME" -name '*.part' 2>/dev/null | wc -l | tr -d ' ')"

## ---- 7. the manifest is DERIVED, and is not a source of trust ---------- ##
## This case replaces "a tampered manifest fails the compiled anchor". There
## is no anchor and no fetched manifest any more: kycg verifies each file
## against the digest compiled into its registry, and then WRITES a manifest
## derived from those same rows. So the mirror cannot influence verification
## at all -- the stronger guarantee, tested in case 6 -- and what lands in the
## store is kycg's own description of the directory.
echo ok > mode
/bin/rm -rf "$YAME_DATA_HOME"
"$KYCG" fetch -f EPIC:Blacklist </dev/null >/dev/null 2>&1
man="$YAME_DATA_HOME/EPIC/KYCG/SHA256SUMS"
check "a derived manifest is written" 1 "$([ -f "$man" ] && echo 1 || echo 0)"
if [ -f "$man" ]; then
  check_has "it lists the set" "$SET" "$(cat "$man")"
  check_has "with the pinned digest" "$SETSHA" "$(cat "$man")"
  ## the companion lives in another directory, so it is not this one's business
  check_lacks "it does not list the ordering" "ordering" "$(cat "$man")"
  ## and it is a real sha256sum file: the store re-verifies with no kycg code.
  ## GNU coreutils calls the tool sha256sum; macOS ships shasum -a 256.
  SUM=$(command -v sha256sum || command -v shasum)
  case $SUM in *shasum) SUM="$SUM -a 256";; esac
  ( cd "$YAME_DATA_HOME/EPIC/KYCG" && $SUM -c --ignore-missing SHA256SUMS ) \
      >/dev/null 2>&1 ||
    { echo "  FAIL the derived manifest does not verify with $SUM -c"
      fails=$((fails + 1)); }
fi

## ---- 8. the companion is planned, and its failure is reported ----------- ##
## The probe ordering rides along with any set -- asking for one means asking
## for the thing that makes it interpretable. The mirror does not carry it, so
## the run must report a failure even though the set itself landed.
/bin/rm -rf "$YAME_DATA_HOME"
out=$("$KYCG" fetch -f EPIC:Blacklist </dev/null 2>&1); rc=$?
check_has "the ordering is part of the plan" "ordering" "$out"
awk -v rc="$rc" 'BEGIN{exit !(rc != 0)}' ||
  { echo "  FAIL a failed companion reported overall success"; fails=$((fails + 1)); }
check "the set itself still landed" 1 \
      "$([ -f "$dest" ] && echo 1 || echo 0)"

## ---- 9. the registry is the authority for what exists ------------------ ##
## This replaces "a tag this build holds no digest for" (-t is gone: tags are
## per file now, so there is no global tag to override). The property that
## took its place: a name with no row carries no digest, so it cannot be
## verified, so it is not offered -- whatever the mirror is willing to serve.
out=$("$KYCG" fetch -f EPIC:NoSuchSet </dev/null 2>&1); rc=$?
awk -v rc="$rc" 'BEGIN{exit !(rc != 0)}' ||
  { echo "  FAIL a set with no registry row was fetched"; fails=$((fails + 1)); }
check_has "the unknown set is named" "NoSuchSet" "$out"

## ---- 10. what was fetched is usable ------------------------------------- ##
## The point of all of it: a set from the store resolves by name in `kycg
## test`, against a query on the same row space.
echo ok > mode
/bin/rm -rf "$YAME_DATA_HOME"
"$KYCG" fetch -f EPIC:Blacklist </dev/null >/dev/null 2>&1
## `yame info` is TAB separated; Nrow is the fourth column.
rows=$("$YAME" info "$dest" 2>/dev/null | tail -1 | cut -f4)
if [ -n "$rows" ] && [ "$rows" -gt 0 ] 2>/dev/null; then
  awk -v n="$rows" 'BEGIN{for(i=0;i<n;i++) print (i%7==0)}' > q.txt
  pack_binary q.txt query.cg
  out=$("$KYCG" test -m EPIC:Blacklist query.cg 2>&1); rc=$?
  check "a fetched set resolves by name in kycg test" 0 "$rc"
  check_has "and reports the set" "Blacklist" "$out"
else
  echo "  skip: could not read the row count from the fetched set"
fi

## ---- 10b. a set that exists upstream but is not here yet ---------------- ##
## Off a terminal there is nobody to answer an offer, so this has to be an
## error -- but the useful kind: it distinguishes "no such set" from "not
## downloaded", and prints the command that fixes it.
## The store must already hold the collection: "published but not here" is a
## statement about a collection kycg can see, and an empty store cannot say
## more than "nothing matches".
/bin/rm -rf "$YAME_DATA_HOME"
"$KYCG" fetch -f EPIC:Blacklist </dev/null >/dev/null 2>&1
out=$("$KYCG" test -m EPIC:CGI query.cg </dev/null 2>&1 >/dev/null); rc=$?
awk -v rc="$rc" 'BEGIN{exit !(rc != 0)}' ||
  { echo "  FAIL an undownloaded set was treated as usable"; fails=$((fails+1)); }
check_has "the set is named" "CGI" "$out"
check_has "and it says where it is" "not in the store" "$out"
check_has "and how to get it" "kycg fetch -f EPIC:CGI" "$out"
## A name nothing publishes is the other mistake, and reads differently.
out=$("$KYCG" test -m EPIC:NoSuchSet query.cg </dev/null 2>&1 >/dev/null)
check_has "a name nothing publishes is a typo" "no set named" "$out"

## ---- 11. the browser, fetching for real --------------------------------- ##
## The interactive half of the same machinery: the tree opens with EPIC's sets
## already checked, `f` builds the plan, the panel asks before spending the
## bandwidth, and the answer drives the download. Nothing here is reachable
## from a redirected run, and it is most of what ui.c does.
if have_pty; then
  echo ok > mode
  /bin/rm -rf "$YAME_DATA_HOME"
  ## Only Blacklist is on the mirror, so this is also the mixed outcome: one
  ## file fetched, the rest 404. The tally line has to say both.
  out=$(PTY_SETTLE=4 PTY_BEAT=2 PTY_TAIL=6 \
        pty_drive "$KYCG fetch EPIC" 'f' 'y' ' ' 'q'); rc=$?
  check "the browser exits cleanly after a fetch" 0 "$rc"
  check "the browser fetched the one set the mirror has" 1 \
        "$([ -f "$dest" ] && echo 1 || echo 0)"
  [ -f "$dest" ] && check "and it verifies" "$SETSHA" "$(sha_of "$dest")"
  ## The panel reports the outcome; a failure count of zero here would mean
  ## the 404s were swallowed.
  printf '%s\n' "$out" | grep -qi "fetched" ||
    { echo "  FAIL the panel never reported a tally"; fails=$((fails + 1)); }
  printf '%s\n' "$out" | grep -qi "failed" ||
    { echo "  FAIL the 16 missing files were not reported as failures"
      fails=$((fails + 1)); }

  ## Answering no must cost nothing: the plan is built and then dropped.
  /bin/rm -rf "$YAME_DATA_HOME"
  out=$(PTY_SETTLE=4 PTY_BEAT=2 PTY_TAIL=4 \
        pty_drive "$KYCG fetch EPIC" 'f' 'n' 'q'); rc=$?
  check "declining the prompt exits cleanly" 0 "$rc"
  check "declining the prompt downloads nothing" 0 \
        "$(find "$YAME_DATA_HOME" -name '*.cm' 2>/dev/null | wc -l | tr -d ' ')"

  ## ---- 11b. the picker `kycg test` opens when -m is omitted ------------ ##
  ## Only sets whose row count matches the query are offered -- the one check
  ## kycg can make about row spaces, applied before the user can get it wrong.
  echo ok > mode
  /bin/rm -rf "$YAME_DATA_HOME"
  "$KYCG" fetch -f EPIC:Blacklist </dev/null >/dev/null 2>&1
  nrow=$("$YAME" info "$dest" 2>/dev/null | tail -1 | cut -f4)
  if [ -n "$nrow" ] && [ "$nrow" -gt 0 ] 2>/dev/null; then
    awk -v n="$nrow" 'BEGIN{for(i=0;i<n;i++) print (i%9==0)}' > big.txt
    pack_binary big.txt big.cg
    out=$(PTY_SETTLE=3 PTY_BEAT=1.2 PTY_TAIL=4 \
          pty_drive "$KYCG test big.cg" \
                    "$PTY_RIGHT" "$PTY_DOWN$PTY_DOWN" ' ' 't'); rc=$?
    check "the picker runs the test it was given" 0 "$rc"
    check_has "and the result names the set" "Blacklist" "$out"
    ## A query on a different row space has nothing to offer.
    awk 'BEGIN{for(i=0;i<40;i++) print i%2}' > tiny.txt
    pack_binary tiny.txt tiny.cg
    out=$(PTY_SETTLE=3 PTY_BEAT=1 PTY_TAIL=3 pty_drive "$KYCG test tiny.cg" 'q')
    check_lacks "a mismatched query is offered nothing" "Blacklist" "$out"
  else
    echo "  skip: could not read the row count for the picker case"
  fi

  ## ---- 11c. -m names a set that is published but not here -------------- ##
  ## The interactive half of case 10b. Off a terminal that is an error naming
  ## the fetch command; on one, kycg opens the catalogue with exactly that set
  ## checked, and the analysis carries on afterwards rather than making the
  ## user retype it. The proof it carried on is the test RESULT at the end --
  ## the spec is re-resolved after the browser closes, so a report naming the
  ## set can only mean the fetch landed and the resolve then succeeded.
  ## Priming matters: "published but not here" is a claim about a collection
  ## kycg can see, so the collection is fetched first and only the set removed.
  echo ok > mode
  /bin/rm -rf "$YAME_DATA_HOME"
  "$KYCG" fetch -f EPIC:Blacklist </dev/null >/dev/null 2>&1
  nrow=$("$YAME" info "$dest" 2>/dev/null | tail -1 | cut -f4)
  if [ -n "$nrow" ] && [ "$nrow" -gt 0 ] 2>/dev/null; then
    awk -v n="$nrow" 'BEGIN{for(i=0;i<n;i++) print (i%11==0)}' > pre.txt
    pack_binary pre.txt pre.cg
    /bin/rm -f "$dest"                 # the set goes, the collection stays
    out=$(PTY_SETTLE=4 PTY_BEAT=2 PTY_TAIL=6 \
          pty_drive "$KYCG test -m EPIC:Blacklist pre.cg" 'f' 'y' ' ' 'q'); rc=$?
    check_has "kycg says it is opening the catalogue" "Opening the catalogue" "$out"
    check "the offered set was fetched" 1 \
          "$([ -f "$dest" ] && echo 1 || echo 0)"
    check "and the test then ran" 0 "$rc"
    check_has "and the result names the set" "Blacklist" "$out"
  else
    echo "  skip: could not read the row count for the prefill case"
  fi

  ## ---- 12. regression: a checked row that carries no name -------------- ##
  ## The companion (the probe ordering) is shown as a row but is not a choice,
  ## so it has no key. Space on a COLLECTION checks every row under it, that
  ## one included, and the action then handed a NULL name to strdup: `kycg
  ## fetch EPIC`, down to another collection, space, f -- two keystrokes from
  ## a fresh browser -- died with SIGSEGV. It must now simply work.
  ## Two declines, not one: this selects sets in TWO collections (EPIC, which
  ## the target preselects, and the one the cursor lands on), and fetch_picked
  ## confirms once per collection. Answering only the first leaves the second
  ## prompt waiting, which reads as a hang.
  /bin/rm -rf "$YAME_DATA_HOME"
  out=$(PTY_SETTLE=4 PTY_BEAT=1.5 PTY_TAIL=4 \
        pty_drive "$KYCG fetch EPIC" \
                  "$PTY_RIGHT" "$PTY_HOME" "$PTY_DOWN$PTY_DOWN" ' ' 'f' \
                  'n' 'n' 'q'); rc=$?
  check "checking a keyless row and acting on it does not crash" 0 "$rc"
  awk -v rc="$rc" 'BEGIN{exit !(rc != 139)}' ||
    { echo "  FAIL SIGSEGV is back"; fails=$((fails + 1)); }
else
  echo "  skip: no usable pty for the browser cases"
fi

done_testing
