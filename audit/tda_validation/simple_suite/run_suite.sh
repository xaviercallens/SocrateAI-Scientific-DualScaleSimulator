#!/usr/bin/env bash
# Runs simple_suite.py cases, one process per case, under a 6 GiB address-space
# cap and a per-command timeout (shared machine; see expectations.json).
# Usage: bash run_suite.sh <timeout_sec> <simple_suite.py args...>
# Example: bash run_suite.sh 590 --case P7 --n 500
set -u
PY=/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python
HERE="$(cd "$(dirname "$0")" && pwd)"
T="$1"; shift
export OMP_NUM_THREADS=1
export SUITE_WRAPPER="OMP_NUM_THREADS=1 timeout $T prlimit --as=6442450944 --"
cd "$HERE"
timeout "$T" prlimit --as=6442450944 -- "$PY" simple_suite.py "$@"
rc=$?
echo "exit=$rc args=$*"
if [ $rc -ne 0 ]; then
  # killed (timeout=124, or a signal) before the case could write its own record: write one
  tag=$(echo "$*" | tr -c 'A-Za-z0-9' '_')
  mkdir -p cases
  printf '{"status": "killed_or_crashed", "exit_code": %d, "args": "%s", "command": "%s", "pass": false}\n' \
    "$rc" "$*" "$SUITE_WRAPPER $PY simple_suite.py $*" > "cases/killed_${tag}.json"
fi
exit $rc
