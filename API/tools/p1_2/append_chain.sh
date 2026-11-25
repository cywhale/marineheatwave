# P1-2: two consecutive monthly appends (538 -> 539 -> 540) on a WRITABLE COPY,
# with fixtures/stubs (no DB, NOAA, network, .env).  Everything under W.
# usage: chain.sh W PY TPY API SRC PORT
#   PY  python of the production-lock venv (updater + API + gunicorn)
#   TPY python of the tools venv (candidate runner)
#   API checkout dir with tools/, mhw_app.py, src/, update_mhw.py
#   SRC read-only source store (never written)
#
# Exit status (Codex r22): 0 and "CHAIN OK" in logs/steps.txt + out/chain_status.txt
# only when every step is proven: run / verify JSON ok=true (exit 0 and
# check-result), server and recover exit 0, the negative control failed exactly
# as expected (verify exit 3, ok=false, the wrong-month check among the failed
# ones), and the source snapshots before / after are identical.  Anything else:
# "FAIL: <reason>", diagnostics (recover --mode list for every API dir; --mode
# stop only when list reports a candidate ALIVE, and stop signals only processes
# whose identity is proven), "CHAIN FAILED", exit 1.  Nothing later runs.
set -u
W=$1; PY=$2; TPY=$3; API=$4; SRC=$5; PORT=$6
cd "$W" || exit 1
mkdir -p logs out home tmp mpl xdg/cache xdg/config xdg/data xdg/state
E="env -i PATH=/usr/bin:/bin LANG=C.UTF-8 HOME=$W/home TMPDIR=$W/tmp XDG_CACHE_HOME=$W/xdg/cache XDG_CONFIG_HOME=$W/xdg/config XDG_DATA_HOME=$W/xdg/data XDG_STATE_HOME=$W/xdg/state MPLCONFIGDIR=$W/mpl MPLBACKEND=Agg PYTHONNOUSERSITE=1 POLARS_FORCE_PKG=compat"
T="$E $PY -I -B $API/tools/p1_2/updater_check.py"
ts() { date -u +%FT%TZ; }
step() { echo "== $(ts) $*" | tee -a logs/steps.txt; }
APIDIRS=""
diagnose() {  # lists candidates; never signals a process whose identity is not proven
  for d in $APIDIRS; do
    (cd "$API" && $E $TPY -B -m tools.bench.candidate recover --dir "$W/$d" --mode list --port $PORT) \
      > "logs/fail_recover_${d}_list.json" 2>&1
    lrc=$?
    echo "diagnose $d: recover --mode list rc=$lrc" | tee -a logs/steps.txt
    if [ $lrc -eq 1 ]; then  # EXIT_FAILED: candidate processes alive (identity proven)
      (cd "$API" && $E $TPY -B -m tools.bench.candidate recover --dir "$W/$d" --mode stop --port $PORT) \
        > "logs/fail_recover_${d}_stop.json" 2>&1
      echo "diagnose $d: recover --mode stop rc=$?" | tee -a logs/steps.txt
    fi
  done
}
fail() {
  echo "FAIL: $*" | tee -a logs/steps.txt
  diagnose
  echo "FAILED: $*" > out/chain_status.txt
  echo "== $(ts) CHAIN FAILED" | tee -a logs/steps.txt
  exit 1
}
# a result JSON must say ok=true (the command's exit status alone is not enough)
expect_ok() { $T check-result --file "$1" --expect ok >> logs/check_results.jsonl 2>&1 || fail "$2: $1 is not ok=true"; }
cmpsnap() { $E $PY -I -B -c "import json,sys; a=json.load(open(sys.argv[1])); b=json.load(open(sys.argv[2])); print('identical' if a==b else 'DIFFERENT', len(a), sum(v[1] for v in a.values()))" "$1" "$2"; }
same_snap() {  # $1 $2 label: fails unless the two snapshots are identical
  R=$(cmpsnap "$1" "$2") || fail "$3: comparison crashed"
  echo "$3: $R" | tee -a logs/steps.txt
  case "$R" in identical*) ;; *) fail "$3: $R" ;; esac
}
rm -f out/chain_status.txt

step "source snapshot (read-only)"
[ -w "$SRC/.zmetadata" ] && fail "source is writable"
$T snapshot --store "$SRC" --out out/src_snap_before.json > logs/src_snap_before.log 2>&1 || fail "snapshot-src rc=$?"
step "copy"
SECONDS=0
if [ "$(uname)" = Darwin ]; then cp -c -R "$SRC" store.rw; else cp -a "$SRC" store.rw; fi || fail copy
chmod -R u+w store.rw || fail chmod
echo "copy seconds: $SECONDS" | tee -a logs/steps.txt
$T snapshot --store store.rw --out out/snap0.json > logs/snap0.log 2>&1 || fail "snap0 rc=$?"
same_snap out/src_snap_before.json out/snap0.json "copy vs source"

for N in 1 2; do
  P=$((N-1))
  step "append $N: metadata copy, fixtures"
  $T meta-copy --store store.rw --out out/meta$P > logs/meta$P.log 2>&1 || fail "meta$P rc=$?"
  $T fixtures --store store.rw --out out/fix$N > logs/fix$N.log 2>&1 || fail "fix$N rc=$?"
  step "append $N: run (update_mhw.append_to_zarr)"
  $T run --store store.rw --fixtures out/fix$N --api "$API" --out out/run$N.json > logs/run$N.log 2>&1 \
    || fail "run$N rc=$? (3: completed / audit check failed, see failed_checks)"
  expect_ok out/run$N.json run$N
  $T snapshot --store store.rw --out out/snap$N.json > logs/snap$N.log 2>&1 || fail "snap$N rc=$?"
  $T diff --before out/snap$P.json --after out/snap$N.json --out out/diff$N.json > logs/diff$N.log 2>&1 \
    || fail "diff$N rc=$?"
  EARLIER=""; [ $N -eq 2 ] && EARLIER="--earlier-fixtures out/fix1"
  step "append $N: verify"
  $T verify --store store.rw --fixtures out/fix$N --before out/snap$P.json --before-store out/meta$P $EARLIER \
     --api "$API" --no-api --out out/verify$N.json > logs/verify$N.log 2>&1
  rc=$?; echo "verify$N rc=$rc" | tee -a logs/steps.txt
  [ $rc -eq 0 ] || fail "verify$N rc=$rc (3: ok=false, see failed_checks)"
  expect_ok out/verify$N.json verify$N
  FX="--fixtures out/fix1"; [ $N -eq 2 ] && FX="--fixtures out/fix1 --fixtures out/fix2"
  $T expect-points $FX --out out/points$N.json > logs/points$N.log 2>&1 || fail "points$N rc=$?"
  step "append $N: loopback API (gunicorn + uvicorn_worker) on the copy"
  rm -rf "api$N" && mkdir -p "api$N/app/API/data" && cp -a "$API/mhw_app.py" "$API/src" "api$N/app/API/" \
    && ln -s "$W/store.rw" "api$N/app/API/data/mhw.zarr" || fail "api dir $N"
  APIDIRS="$APIDIRS api$N"
  (cd "$API" && $E $TPY -B -m tools.p1_2.append_server --dir "$W/api$N" --port $PORT --python $PY \
     --gunicorn "$(dirname $PY)/gunicorn" --points "$W/out/points$N.json" --step after_append$N) > logs/server$N.log 2>&1
  rc=$?; echo "server$N rc=$rc" | tee -a logs/steps.txt
  [ $rc -eq 0 ] || fail "server$N rc=$rc"
  (cd "$API" && $E $TPY -B -m tools.bench.candidate recover --dir "$W/api$N" --mode list --port $PORT) > logs/recover$N.json 2>&1
  rc=$?; echo "recover$N rc=$rc" | tee -a logs/steps.txt
  [ $rc -eq 0 ] || fail "recover$N rc=$rc (candidate alive, port bound or identity undetermined)"
done

step "negative control: the twice-appended store checked against the FIRST append's expectations (must fail)"
$T verify --store store.rw --fixtures out/fix1 --before out/snap1.json --before-store out/meta1 --api "$API" --no-api \
   --out out/verify_negative_wrong_n.json > logs/verify_negative.log 2>&1
rc=$?; echo "negative verify rc=$rc (expected 3: a complete result with ok=false)" | tee -a logs/steps.txt
[ $rc -eq 3 ] || fail "negative control: verify rc=$rc, expected 3 (0 = it passed, 1 = it crashed)"
$T check-result --file out/verify_negative_wrong_n.json --expect fail --expect-failed dates_last_is_fixture_month \
  >> logs/check_results.jsonl 2>&1 || fail "negative control: not the expected failure"
echo "negative control: failed as expected" | tee -a logs/steps.txt
step "source snapshot again"
$T snapshot --store "$SRC" --out out/src_snap_after.json > logs/src_snap_after.log 2>&1 || fail "snapshot-src-after rc=$?"
same_snap out/src_snap_before.json out/src_snap_after.json "source before vs after"
echo "OK" > out/chain_status.txt
step "CHAIN OK"
