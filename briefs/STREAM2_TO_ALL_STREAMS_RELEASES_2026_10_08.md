# Stream 2 → all streams — the three papers are released (2026-10-08)

**From:** Stream 2, on T0's instruction *"complete the work and communicate to others streams and publish the stream 1, stream 2
and stream 3 results as tex and pdf. commit, push, merge and release"*. **To:** Stream 1, LeanMaster, DualScaleSimulator,
Stream 3. Nothing is required of anyone. Nothing was deposited on Zenodo or submitted anywhere; GitHub releases only.

| Stream | Release | PDF (tex in the repo) | Gate the tag waited for |
|---|---|---|---|
| 1 — Theory (LeanProposal) | `v0.27-paper-addendum-2026-10-08` @ 8dbf438 | `dual-scale-stream1-paper.pdf` (`paper/main.tex`, 49 pp.) | `release_gates.sh` **ALL GATES OK** on the tagged commit |
| 2 — Selection & Geometry (K3-DarkMatter) | `v0.3.19-k3xt2-reading-s` @ 5fe0f98 | `stream2_k3_identification_twisted_route_2026_10_08.pdf` (`papers/`) | Agora CI Gate post-merge run 37832068928, green |
| 3 — Experimentation (Home) | `v5.13.0-stream3-note-2026-10-08` @ e1502cc | `stream3_status_note_2026_10_08.pdf` (`papers/`) | Home CI is **red at its install step** (pre-existing; fix waits on a `workflow`-scoped token); local tests pass |

## What each says, in one line

- **Stream 1:** the version-4 text unchanged, plus a dated addendum §12 on companion results. It adds no kernel-checked statement, and
  the rescaling gap of its conjecture stays open.
- **Stream 2:** the selected surfaces are X₃ (s7) and X₄ (s10). "K3×T²" means the Shioda–Inose pairing (D29′), with T(E×E′) ≅ U⊕⟨2n⟩
  computed for n = 7 and n = 10. At the selected point z = ∞ there is no twisted-route section.
- **Stream 3:** no dark-matter result. It contains calibrations and the Phase 0 recount (221 → 182 / 152 open cells).

## Standing limits, for every stream

No physical reading of any of this (Tier C stays blocked). No ranking of s7 against s10. No Kodaira labels. A compactification's
T² is not E or E′. No number from these papers is a simulation input or a prior.

## Open, owned elsewhere

- **T0:** Zenodo version-5 deposit for Stream 1; proposal decisions 2–4 (T3 hard gate, K3_CRITERIA edit, scoring); the `workflow`
  token scope for Home CI.
- **Stream 1's own sessions:** the TW2 lattice-attestation request of 2026-10-07.

*Generated-by: Claude (Opus 5.5), Stream 2 | Verified-by: the three release pages and their gate records | Reviewed-by: N*
