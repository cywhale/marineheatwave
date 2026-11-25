#!/usr/bin/env bash
# Shell fallback for the isolated candidate runner (doc 24 r4 §3.5).  Linux only.
# Use only when `python -B -m tools.bench.candidate recover` cannot run.
#
#   bash candidate_recover.sh list|stop RUN_DIR PORT
#
# RUN_DIR is one step directory ($D/run/<step>) holding candidate.id and
# candidate.pids.  Same identity rules as tools/bench/procident.py:
#   V  recorded processes alive with the recorded start (/proc/PID/stat field 22)
#   T  V + members of sessions/groups of V + children of T (recursively)
#   F  provably foreign: sid/pgid = a recorded number whose leader is alive with another start
#   U  undetermined: sid/pgid = a recorded number, not in T, not provably foreign
# Any U stops everything before a signal is sent.  No `kill -- -PGID`: every
# kill re-reads that PID's start first (bash has no pidfd, so a tiny window
# between the check and the kill remains).
#
# A /proc entry of our uid that exists but cannot be read or parsed might be
# a candidate process, so it makes the state undetermined too (no signal,
# never reported clean).  Unreadable entries of other uids are not ours.
#
# Exit: 0 clean | 1 candidates alive (list) | 2 usage/format | 5 left after stop
#       or port bound or port check failed | 6 undetermined | 7 recorded number
#       now used by another process (nothing signalled for it).
#
# Real use, with a clean environment:
#   env -i PATH=/usr/bin:/bin bash candidate_recover.sh list|stop RUN_DIR PORT
#
# Test-only overrides, honoured only with MHW_RECOVER_TEST_MODE=1 (any of them
# without it is refused): MHW_RECOVER_PROC (proc root), MHW_RECOVER_KILL (kill
# command, called as: CMD -s SIG PID), MHW_RECOVER_PORTCHECK (command called
# with PORT; prints listeners), MHW_RECOVER_TERM_WAIT, MHW_RECOVER_KILL_WAIT,
# MHW_RECOVER_POLL (seconds).
set -u -o pipefail

MODE=${1:-}
RUN=${2:-}
PORT=${3:-}
PROC=${MHW_RECOVER_PROC:-/proc}
KILL=${MHW_RECOVER_KILL:-kill}
TERM_WAIT=${MHW_RECOVER_TERM_WAIT:-30}
KILL_WAIT=${MHW_RECOVER_KILL_WAIT:-5}
POLL=${MHW_RECOVER_POLL:-1}
NUM='^(0|[1-9][0-9]*)$'

fail() { echo "$2"; exit "$1"; }

if [ -n "${MHW_RECOVER_PROC+x}${MHW_RECOVER_KILL+x}${MHW_RECOVER_PORTCHECK+x}${MHW_RECOVER_TERM_WAIT+x}${MHW_RECOVER_KILL_WAIT+x}${MHW_RECOVER_POLL+x}" ]; then
  [ "${MHW_RECOVER_TEST_MODE:-}" = 1 ] \
    || fail 2 "REFUSED: MHW_RECOVER_* overrides are test-only; run with env -i (see header) or set MHW_RECOVER_TEST_MODE=1 in tests"
  echo "TEST MODE: proc=$PROC kill=$KILL"
fi

case "$MODE" in list|stop) ;; *) fail 2 "usage: $0 list|stop RUN_DIR PORT" ;; esac
[[ "$PORT" =~ $NUM ]] || fail 2 "FORMAT: port must be a number"
[ -d "$RUN" ] || fail 2 "FORMAT: run dir $RUN does not exist"
for v in "$TERM_WAIT" "$KILL_WAIT" "$POLL"; do [[ "$v" =~ $NUM ]] || fail 2 "FORMAT: wait values must be integers"; done
[ "$POLL" -ge 1 ] || fail 2 "FORMAT: poll interval must be >= 1"

# start value of one PID ("" when gone, zombie or unparsable)
start_of() {
  local s rest
  s=$(cat "$PROC/$1/stat" 2>/dev/null) || return 0
  [[ "$s" == "$1 ("* ]] || return 0
  rest=${s##*) }
  set -- $rest
  [ $# -ge 20 ] || return 0
  case "$1" in Z|X|x) return 0 ;; esac
  [[ "${20}" =~ $NUM ]] && printf '%s' "${20}"
  return 0
}

# an entry that exists but cannot be used: "B pid" when it is ours (blocks), else a warning
bad() {
  [ -d "$PROC/$1" ] || return 0  # vanished while being read: gone
  if [ -O "$PROC/$1" ]; then
    echo "B $1"; echo "warning: $2 for $1 (owned by us: undetermined)" >&2
  else
    echo "warning: $2 for $1 (another uid: not a candidate)" >&2
  fi
}

# one line per live process: "S pid ppid pgid sid start"; "B pid" for unusable entries of ours
snapshot() {
  local f p s rest
  for f in "$PROC"/[0-9]*/stat; do
    p=${f#"$PROC"/}; p=${p%/stat}
    [[ "$p" =~ $NUM ]] || continue
    s=$(cat "$f" 2>/dev/null) || { bad "$p" "unreadable stat"; continue; }
    [[ "$s" == "$p ("* ]] || { bad "$p" "unparsable stat"; continue; }
    rest=${s##*) }
    set -- $rest
    if [ $# -lt 20 ]; then bad "$p" "short stat"; continue; fi
    case "$1" in Z|X|x) continue ;; esac
    if ! [[ "$2" =~ $NUM && "$3" =~ $NUM && "$4" =~ $NUM && "${20}" =~ $NUM ]]; then
      bad "$p" "non-numeric stat"; continue
    fi
    echo "S $p $2 $3 $4 ${20}"
  done
}

port_bound() {  # 0 bound, 1 free, 2 check failed
  local out
  if [ -n "${MHW_RECOVER_PORTCHECK:-}" ]; then
    out=$($MHW_RECOVER_PORTCHECK "$PORT") || return 2
  else
    out=$(ss -ltnH "sport = :$PORT") || return 2
  fi
  [ -n "$out" ] && return 0
  return 1
}

# ---------------------------------------------------------------- state files
if [ ! -e "$RUN/candidate.id" ]; then
  # gunicorn never ran without candidate.id (the child writes it before exec);
  # list processes whose cwd is the candidate app dir, never signal them.
  D=${RUN%/run/*}
  hits=""
  for f in "$PROC"/[0-9]*; do
    [ "$(readlink "$f/cwd" 2>/dev/null)" = "$D/app/API" ] && hits="$hits ${f#"$PROC"/}"
  done
  if [ -n "$hits" ]; then fail 6 "UNDETERMINED: no candidate.id; processes with the candidate cwd:$hits (not signalled)"; fi
  bads=$(snapshot 2>/dev/null | awk '$1=="B" {printf " %s", $2}') || fail 6 "UNDETERMINED: process table could not be read"
  if [ -n "$bads" ]; then fail 6 "UNDETERMINED: no candidate.id; unreadable entries of our uid:$bads (not signalled)"; fi
  fail 0 "no candidate.id and no process with the candidate cwd: no candidate process"
fi
idtext=$(cat "$RUN/candidate.id") || fail 2 "FORMAT: cannot read candidate.id"
[[ "$idtext" != *$'\n'* ]] || fail 2 "FORMAT: candidate.id must be one line"
re='^mhw-candidate-id v1 linux-ticks (0|[1-9][0-9]*) (0|[1-9][0-9]*) (0|[1-9][0-9]*) (0|[1-9][0-9]*)$'
[[ "$idtext" =~ $re ]] || fail 2 "FORMAT: candidate.id: unexpected content or unit (shell accepts linux-ticks only)"
PGID=${BASH_REMATCH[1]}; SID=${BASH_REMATCH[2]}; MPID=${BASH_REMATCH[3]}; MSTART=${BASH_REMATCH[4]}
[ "$PGID" = "$SID" ] && [ "$SID" = "$MPID" ] && [ "$MPID" -gt 1 ] || fail 2 "FORMAT: candidate.id pgid/sid/pid must be equal and > 1"

[ -e "$RUN/candidate.pids" ] || fail 2 "FORMAT: candidate.pids is missing"
RECS=""
n=0
while IFS= read -r line || [ -n "$line" ]; do
  n=$((n + 1))
  if [ "$n" -eq 1 ]; then
    [ "$line" = "# mhw-candidate-pids v1 linux-ticks" ] || fail 2 "FORMAT: candidate.pids header"
    continue
  fi
  [[ "$line" =~ ^(0|[1-9][0-9]*)\ (0|[1-9][0-9]*)\ (0|[1-9][0-9]*)\ (0|[1-9][0-9]*)\ (0|[1-9][0-9]*)$ ]] \
    || fail 2 "FORMAT: candidate.pids line $n"
  RECS+="R $line"$'\n'
done < "$RUN/candidate.pids" || fail 2 "FORMAT: cannot read candidate.pids"
[ "$n" -ge 1 ] || fail 2 "FORMAT: candidate.pids is empty"

# ------------------------------------------------------------------ resolution
AWK='
$1=="M" { G=$2; S=$3; MP=$4; MS=$5; num[G]=1; num[S]=1; ls[S]=ls[S] " " MS " "; nr++; rp[nr]=MP; rs[nr]=MS; next }
$1=="R" { nr++; rp[nr]=$2; rs[nr]=$3; num[$5]=1; num[$6]=1
          if ($2==$6 || $2==$5) ls[$2]=ls[$2] " " $3 " "; next }
$1=="B" { print "X", $2; next }
$1=="S" { n++; P[n]=$2; PP[$2]=$3; PG[$2]=$4; SS[$2]=$5; ST[$2]=$6; alive[$2]=1; next }
END {
  for (i=1; i<=nr; i++) { p=rp[i]; if ((p in alive) && ST[p]==rs[i]) { T[p]="recorded"; as[SS[p]]=1; ag[PG[p]]=1 } }
  ms = !(MP in alive) ? "gone" : (ST[MP]==MS ? "alive" : "reused")
  for (i=1; i<=n; i++) { p=P[i]; if (p in T) continue
    if (SS[p] in as) T[p]="session"; else if (PG[p] in ag) T[p]="group" }
  do { ch=0; for (i=1; i<=n; i++) { p=P[i]; if (!(p in T) && (PP[p] in T)) { T[p]="child"; ch=1 } } } while (ch)
  if (MP in T) print "T", MP, ST[MP], T[MP], PP[MP], PG[MP], SS[MP]
  for (i=1; i<=n; i++) { p=P[i]; if ((p in T) && p!=MP) print "T", p, ST[p], T[p], PP[p], PG[p], SS[p] }
  for (i=1; i<=n; i++) { p=P[i]; if (p in T) continue
    hit=0; foreign=0
    for (k=1; k<=2; k++) { x=(k==1) ? SS[p] : PG[p]
      if (!(x in num)) continue
      hit=1
      if ((x in alive) && (x in ls) && index(ls[x], " " ST[x] " ")==0) foreign=1 }
    if (!hit) continue
    print (foreign ? "F" : "U"), p, ST[p] }
  print "MASTER", ms
}'

RES=""
resolve_now() {
  local snap
  snap=$(snapshot) || return 1
  RES=$( { echo "M $PGID $SID $MPID $MSTART"; printf '%s' "$RECS"; printf '%s\n' "$snap"; } | awk "$AWK" ) || return 1
  # carry every proven target into later rounds as recorded (pid + start), so
  # stopping the master first cannot turn unrecorded workers into undetermined ones
  RECS+=$(printf '%s\n' "$RES" | awk '$1=="T" {print "R", $2, $3, $5, $6, $7}')$'\n'
  return 0
}
count() { printf '%s\n' "$RES" | awk -v t="$1" '$1==t {c++} END {print c+0}'; }
state() {  # none | alive | undetermined
  if [ "$(count U)" -gt 0 ] || [ "$(count X)" -gt 0 ]; then echo undetermined
  elif [ "$(count T)" -gt 0 ]; then echo alive
  else echo none; fi
}
show() {
  echo "identity: pgid=$PGID sid=$SID master=$MPID start=$MSTART"
  printf '%s\n' "$RES" | awk '
    $1=="T" { print "candidate " $2 " start=" $3 " via=" $4 }
    $1=="F" { print "foreign   " $2 " start=" $3 " (recorded number reused; not signalled)" }
    $1=="U" { print "UNDETERMINED " $2 " start=" $3 " (not signalled)" }
    $1=="X" { print "UNREADABLE " $2 " (our uid; cannot be checked; nothing signalled)" }
    $1=="MASTER" { print "master: " $2 }'
}
mismatch() { [ "$(count F)" -gt 0 ] || printf '%s\n' "$RES" | grep -qx 'MASTER reused'; }
send_all() {  # $1 = TERM | KILL; re-check each start right before its kill
  local tag p st how rest cur
  while read -r tag p st how rest; do
    [ "$tag" = T ] || continue
    cur=$(start_of "$p")
    if [ "$cur" != "$st" ]; then echo "skip $1 $p: start changed or process gone"; continue; fi
    if $KILL -s "$1" "$p" 2>/dev/null; then echo "sent $1 $p"; else echo "kill $1 $p failed (gone?)"; fi
  done <<< "$RES"
}
finish() {  # final exit for a resolved, non-undetermined state (caller has shown it)
  local st=$1 pb
  if [ "$st" = alive ]; then
    [ "$MODE" = list ] && fail 1 "CANDIDATE PROCESSES ALIVE"
    fail 5 "CANDIDATE STILL ALIVE"
  fi
  port_bound; pb=$?
  [ "$pb" -eq 0 ] && fail 5 "PORT $PORT STILL BOUND"
  [ "$pb" -eq 2 ] && fail 5 "PORT CHECK FAILED: cannot confirm port $PORT is free"
  if mismatch; then fail 7 "no candidate process; a recorded number now belongs to another process (not signalled)"; fi
  fail 0 "no candidate process; port $PORT free"
}

resolve_now || fail 6 "UNDETERMINED: process table could not be read"
st=$(state)
if [ "$st" = undetermined ]; then show; fail 6 "UNDETERMINED: nothing signalled"; fi
[ "$MODE" = list ] && { show; finish "$st"; }

for sig in TERM KILL; do
  [ "$st" = alive ] || break
  show
  send_all "$sig"
  limit=$TERM_WAIT; [ "$sig" = KILL ] && limit=$KILL_WAIT
  waited=0
  while :; do
    resolve_now || fail 6 "UNDETERMINED: process table could not be read"
    st=$(state)
    if [ "$st" = undetermined ]; then show; fail 6 "UNDETERMINED: stopped before further signals"; fi
    [ "$st" = alive ] || break
    [ "$waited" -ge "$limit" ] && break
    sleep "$POLL"; waited=$((waited + POLL))
  done
done
show
finish "$st"
