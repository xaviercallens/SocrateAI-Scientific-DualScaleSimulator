#!/usr/bin/env bash
# Read-only re-orientation after a crash, spend-limit stop or long pause. Changes nothing.
# Usage: scripts/session_status.sh
set -u
cd "$(dirname "$0")/.."

section() { printf '\n=== %s ===\n' "$1"; }

section "main vs origin"
git fetch -q origin 2>/dev/null || echo "(fetch failed: offline or no credentials)"
git status -sb | head -1
echo "unpushed commits: $(git rev-list --count origin/main..main 2>/dev/null)"
echo "commits behind:   $(git rev-list --count main..origin/main 2>/dev/null)"
git log --oneline -5

section "working tree (tracked changes, then untracked outside .claude/worktrees)"
git status --porcelain | grep -v '^??' || echo "(no tracked changes)"
git status --porcelain | grep '^??' | grep -v '\.claude/worktrees' || echo "(no untracked files)"

section "stashes (other sessions' too -- never pop blindly)"
git stash list --format='%h %gs' | head -5

section "loop/* branches not merged into main"
for b in $(git branch --format='%(refname:short)' --list 'loop/*' 'papers/*' 'haiku/*'); do
  n=$(git rev-list --count "main..$b")
  [ "$n" -gt 0 ] && echo "$b: $n commit(s) not in main"
done
echo "(done)"

section "last CI runs"
if command -v gh >/dev/null; then gh run list -L 3 2>&1 | cut -f1-5,7; else echo "gh not installed"; fi

section "Zenodo (public API, read-only)"
curl -sL --max-time 15 "https://zenodo.org/api/records/22683564/versions/latest" \
  | python3 -c 'import json,sys; r=json.load(sys.stdin); print(r.get("id"), r.get("doi"), r["metadata"].get("version"), r["metadata"].get("publication_date"))' \
  2>/dev/null || echo "(unreachable)"

section "Lean axiom audit (cached report; regenerate with python3 scripts/lean_axiom_audit.py)"
if [ -f audit/lean_axiom_report.json ]; then
  python3 - <<'PY'
import json
r = json.load(open("audit/lean_axiom_report.json"))
bad = [(lib, t) for lib, v in r.items() if isinstance(v, dict)
       for t, d in v.get("theorems", {}).items() if d.get("classification") != "STANDARD"]
print(f"{len(r)} libraries, {len(bad)} non-standard theorem(s)")
for lib, t in bad: print(" ", lib, t)
PY
else
  echo "(no cached report)"
fi

section "next"
echo "Read the memory index first: ~/.claude/projects/-home-callensxavier-gmail-com-SocrateAI-Scientific-DualScaleSimulator/memory/MEMORY.md"
echo "Push, reset --hard, rm and any Zenodo deposit still need explicit user confirmation."
