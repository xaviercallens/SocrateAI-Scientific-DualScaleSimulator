# Stream 2 → all streams — state of 2026-10-10 (FYI; nothing required)

**From:** Stream 2, after T0's instructions of 2026-10-08/10 (Zenodo deposit; consider scoring; start T3; validate K3×T²; re-run Home
numerics; improve the Simulator). **To:** Stream 1, LeanMaster, DualScaleSimulator, Stream 3.

## What changed

| Item | Result | Where |
|---|---|---|
| Zenodo, Stream 1 | Version 5 deposited: **10.5281/zenodo.23248582** (concept 10.5281/zenodo.22853238). The Lean development is identical to v4; the paper has a dated addendum §12 | LeanProposal v0.27, PRs #4–#8 |
| T3 | **Hard gate** (AM-9, D30′): a `T3_DISAGREE` certificate triggers F1, an absent one does not; never scored; no candidate removed | K3 `v0.3.20-am9-t3-hard-gate` |
| Scoring | **Considered and not adopted** (D31′): no ranking before the v1.0 freeze, C4/C5 have no checkers, s7 and s10 tie on every hard gate, ledger items 8/9 apply. `autoevolve/run_ranking.py` does not exist | K3 TODO / ledger item 13 |
| WP-TW2 step 2b-ii | The section at z = 1/27 is built and read exactly: P̄·Ō = 5, contr = 0, h = 14, agreeing with steps 0 and 2a | K3 `v0.3.21-tw2-section-z127` |
| K3_CRITERIA mirrors | Re-pinned (Stream 1, PR #8) and re-synced (Home, PR #16) to the AM-9 text | — |
| Home numerics | Full suite **678 passed** under `~/venv` (the earlier "iminuit missing" note was wrong). WP-E6 synthetic validation re-run: **bit-identical**. **No real-data run**: the real-data mode is blocked pending T0's rulings | Home TODO |
| Simulator | Branch `stream2-results-leanflow-2026-10-08` (commit aba1ef6, **unpushed**): Stream 2's results recorded as mathematics, a Reading S cross-check test, a stale level test fixed, false LeanFlow badges withdrawn. Suite 190 passed, 5 skipped | Simulator worktree |

## Standing answers, for every stream

- **The physical K3×T² dual-scale theory is not validated and cannot be today.** It is Tier C: no EFT matching, no derived
  observable, F5b blocks the flux route. What is validated is the mathematics (Reading S lattice, T3, the TW2 section). See
  `STREAM2_K3xT2_DUAL_SCALE_SCOPING_2026_10_10.md`.
- Home has **no** "K3 numeric selection" pipeline. The K3 selection is Stream 2's and feeds no observable.
- Nothing from these results is a simulation input or a prior (ledger item 4; D2 of the Simulator is unchanged).

## Owned by T0 (open)

1. The Simulator branch: push it, and note that its local `main` is 11 commits ahead of `origin/main` (earlier sessions' work, never
   pushed). Its ROADMAP requires explicit confirmation for any push.
2. Home real-data WP-E6: the C1–C5 theory-corrections ruling, the eligibility rule (E-A/B/C), and setting `REAL_DATA_RULING_PIN`.
3. WP-TW3: twist data over a threefold base (no T0 text specifies it).
4. Home CI: the `workflow` token scope for the prepared install fix.
5. Proposal decisions 3–4: a `K3_CRITERIA.md` edit for the score formula and soft weights (§7 freeze sign-offs).

## Owned by Stream 1's sessions

The kernel check of the Reading S lattice and of the TW2 lattice arithmetic (`STREAM2_TO_STREAM1_TW2_LATTICE_ATTESTATION_REQUEST_2026_10_07.md`
and its 2026-10-10 addendum). The statement must carry set equality of the complement (`LL.md` §1).

*Generated-by: Claude (Sonnet 5.5), Stream 2 | Verified-by: the releases, PRs and certificates named | Reviewed-by: N*
