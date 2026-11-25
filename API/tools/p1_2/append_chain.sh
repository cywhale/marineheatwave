# P1-2: two consecutive monthly appends (538 -> 539 -> 540) on a WRITABLE COPY,
# with fixtures/stubs (no DB, NOAA, network, .env).  Everything under W.
# usage: chain.sh W PY TPY API SRC PORT
#   PY  python of the production-lock venv (updater + API + gunicorn)
#   TPY python of the tools venv (candidate runner)
#   API checkout dir with tools/, mhw_app.py, src/, update_mhw.py
#   SRC read-only source store (never written)
set -u
W=$1; PY=$2; TPY=$3; API=$4; SRC=$5; PORT=$6
cd "$W" || exit 1
mkdir -p logs out home tmp mpl xdg/cache xdg/config xdg/data xdg/state
E="env -i PATH=/usr/bin:/bin LANG=C.UTF-8 HOME=$W/home TMPDIR=$W/tmp XDG_CACHE_HOME=$W/xdg/cache XDG_CONFIG_HOME=$W/xdg/config XDG_DATA_HOME=$W/xdg/data XDG_STATE_HOME=$W/xdg/state MPLCONFIGDIR=$W/mpl MPLBACKEND=Agg PYTHONNOUSERSITE=1 POLARS_FORCE_PKG=compat"
T="$E $PY -I -B $API/tools/p1_2/updater_check.py"
ts() { date -u +%FT%TZ; }
step() { echo "== $(ts) $*" | tee -a logs/steps.txt; }
fail() { echo "FAIL: $*" | tee -a logs/steps.txt; exit 1; }
cmpsnap() { $E $PY -I -B -c "import json,sys; a=json.load(open(sys.argv[1])); b=json.load(open(sys.argv[2])); print('identical' if a==b else 'DIFFERENT', len(a), sum(v[1] for v in a.values()))" "$1" "$2"; }

step "source snapshot (read-only)"
[ -w "$SRC/.zmetadata" ] && fail "source is writable"
$T snapshot --store "$SRC" --out out/src_snap_before.json > logs/src_snap_before.log 2>&1 || fail snapshot-src
step "copy"
SECONDS=0
if [ "$(uname)" = Darwin ]; then cp -c -R "$SRC" store.rw; else cp -a "$SRC" store.rw; fi || fail copy
chmod -R u+w store.rw || fail chmod
echo "copy seconds: $SECONDS" | tee -a logs/steps.txt
$T snapshot --store store.rw --out out/snap0.json > logs/snap0.log 2>&1 || fail snap0
echo "copy vs source: $(cmpsnap out/src_snap_before.json out/snap0.json)" | tee -a logs/steps.txt
cmpsnap out/src_snap_before.json out/snap0.json | grep -q '^identical' || fail "copy differs from source"

for N in 1 2; do
  P=$((N-1))
  step "append $N: metadata copy, fixtures"
  $T meta-copy --store store.rw --out out/meta$P > logs/meta$P.log 2>&1 || fail meta$P
  $T fixtures --store store.rw --out out/fix$N > logs/fix$N.log 2>&1 || fail fix$N
  step "append $N: run (update_mhw.append_to_zarr)"
  $T run --store store.rw --fixtures out/fix$N --api "$API" --out out/run$N.json > logs/run$N.log 2>&1 || fail run$N
  $T snapshot --store store.rw --out out/snap$N.json > logs/snap$N.log 2>&1 || fail snap$N
  $T diff --before out/snap$P.json --after out/snap$N.json --out out/diff$N.json > logs/diff$N.log 2>&1 || fail diff$N
  EARLIER=""; [ $N -eq 2 ] && EARLIER="--earlier-fixtures out/fix1"
  step "append $N: verify"
  $T verify --store store.rw --fixtures out/fix$N --before out/snap$P.json --before-store out/meta$P $EARLIER \
     --api "$API" --no-api --out out/verify$N.json > logs/verify$N.log 2>&1; echo "verify$N rc=$?" | tee -a logs/steps.txt
  FX="--fixtures out/fix1"; [ $N -eq 2 ] && FX="--fixtures out/fix1 --fixtures out/fix2"
  $T expect-points $FX --out out/points$N.json > logs/points$N.log 2>&1 || fail points$N
  step "append $N: loopback API (gunicorn + uvicorn_worker) on the copy"
  rm -rf "api$N" && mkdir -p "api$N/app/API/data" && cp -a "$API/mhw_app.py" "$API/src" "api$N/app/API/" \
    && ln -s "$W/store.rw" "api$N/app/API/data/mhw.zarr" || fail "api dir $N"
  (cd "$API" && $E $TPY -B -m tools.p1_2.append_server --dir "$W/api$N" --port $PORT --python $PY \
     --gunicorn "$(dirname $PY)/gunicorn" --points "$W/out/points$N.json" --step after_append$N) > logs/server$N.log 2>&1
  echo "server$N rc=$?" | tee -a logs/steps.txt
  (cd "$API" && $E $TPY -B -m tools.bench.candidate recover --dir "$W/api$N" --mode list --port $PORT) > logs/recover$N.json 2>&1
  echo "recover$N rc=$?" | tee -a logs/steps.txt
done

step "negative control: the twice-appended store checked against the FIRST append's expectations (must fail)"
$T verify --store store.rw --fixtures out/fix1 --before out/snap1.json --before-store out/meta1 --api "$API" --no-api \
   --out out/verify_negative_wrong_n.json > logs/verify_negative.log 2>&1; echo "negative rc=$?" | tee -a logs/steps.txt
step "source snapshot again"
$T snapshot --store "$SRC" --out out/src_snap_after.json > logs/src_snap_after.log 2>&1 || fail snapshot-src-after
echo "source before vs after: $(cmpsnap out/src_snap_before.json out/src_snap_after.json)" | tee -a logs/steps.txt
step "CHAIN DONE"
