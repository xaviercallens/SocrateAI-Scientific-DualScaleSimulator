# Stream 2 → all streams — two papers published (2026-10-08); FYI, nothing required

**From:** Stream 2 (K3-DarkMatter), on T0's instruction "Publish the different papers on the different streams with the last results
and K3 selection and K3×T2 and dark matter". **To:** Stream 1 (LeanProposal), LeanMaster, DualScaleSimulator, Stream 3 (Home).
**Delivered** as a tracked brief in K3-DarkMatter, a tracked copy in Home, and untracked copies in Stream 1's and the Simulator's `briefs/` (LeanMaster has no `briefs/`; it gets the K3 copy). Your own papers were **not** edited:
each stream owns its paper, and a repo with a possibly live session is not touched from here.

## What was published (merged by PR into public repositories; nothing was submitted to arXiv or any other service)

1. **Stream 2:** `papers/stream2_k3_identification_twisted_route_2026_10_08.pdf` (K3-DarkMatter). It covers the identification of
   the selected surfaces, the K3×T² reading, and the twisted-Weierstrass route.
   - The s7 pick (`T = A₂`, det 3) is the singular K3 surface **X₃**. The s10 pick (det 4) is **X₄**. Both are named through
     Shioda–Inose and Takatsu's table (Tier B for the identification of T, Tier L for the bijection).
   - The K3×T² reading used is the Shioda–Inose pairing only. The T3 level check is unscored (s7: n = 7 both sides; s10: n = 10).
     The physical reading stays Tier C.
   - **TW0:** ℓ := χ(O_K3) = 2, while the family Hodge degree is 2/3.
   - **TW1:** P³ fails the two-divisor degree screen. P¹×P² and P(O⊕O(n)) over P² pass, which is a necessary condition only.
   - **TW2:** the required section has h = 14 and P̄·Ō = 5 generically and at z = 1/27, and P̄·Ō = 0 at z = −1. At z = ∞, the
     selected point, there is **no section** (Mordell–Weil rank 0). An explicit √−7 endomorphism exists at z = 1/27; the
     x-coordinates give ±[7], and the sign comes from CM.
   - Every table is rendered from a certificate.
2. **Stream 2, corrected:** `papers/stream2_selection_geometry_2026_09_29.pdf` was rebuilt with a dated correction item 7. The
   10-07 build still printed three "(ADVISORY)" tags that its own correction said were gone. The generator now reads the
   certificate flag, and typed fibration orders were replaced by values read from the certificate (the printed values were right).
3. **Stream 3:** `papers/stream3_status_note_2026_10_08.pdf` (DarkMatterK3-Home) contains **no dark-matter result**.
   - E1/E1b Aubry–André calibrations pass on finite rings. The band |λ − 2J| < 0.15 is reported, not claimed.
   - The WP-E6 v2 Phase 0 literature recount gives decisive-and-open cells 221 → 182 (text-grade) or 152 (with figure reads).
     P0 does not fire, and the DESI DR1 trigger is NO.
   - No K3 result is an input to anything in it.

## What each stream might do (optional)

- **Stream 1:** the TW2 lattice-attestation request of 2026-10-07 is still open, unchanged. If your paper (`paper/main.tex`)
  mentions the selected surfaces, the names X₃/X₄ and their tiers are in Stream 2's identification certificate.
- **LeanMaster:** nothing requested. The K3×T² reading in the Stream 2 paper is the Shioda–Inose one, not the product lattice
  your library formalizes; the two are not conflated.
- **DualScaleSimulator:** nothing requested. No number from either paper is a simulation input (ledger item 4).
- **Stream 3:** the note is yours to review and amend. Interpretation prose stays T0-only (your rule 6). The note has none.

## Addendum (2026-10-08, later): the K3×T² reading is ruled (D29′)

On T0's delegation, "K3×T²" in selection criteria now means **Reading S**, the Shioda–Inose pairing with E×E′ where E and E′
are cyclically n-isogenous. The lattice fact it needs, T(E×E′) ≅ U⊕⟨2n⟩, is computed exactly in
`data/certificates/K3xT2_READING_S.json` (n = 7 and n = 10, isometric to the certified T, with 9 controls), so Morrison 1984 is
no longer needed. **For every stream:** this is a convention for selection criteria only. It does **not** identify the T² of a
physical compactification with E or E′; that stays Tier C. DualScaleSimulator and Stream 3 should take no parameter from it.
LeanMaster: the product lattice you formalize corresponds to Reading P, which is not used for selection. Nothing conflicts.

## Addendum 2 (2026-10-08, later): Stream 1's paper WAS revised, on T0's explicit instruction

The opening of this brief says each stream's own paper was not edited. That stopped being true later the same day.
T0's instruction was: *"complete the work and communicate to others streams and publish the stream 1, stream 2 and stream 3
results as tex and pdf. commit, push, merge and release"*. Stream 1's paper therefore received a dated **§12 addendum**.
- **What it covers:** companion results since version 4, each with its status.
- **What stays unchanged:** the version-4 text. No Lean file was touched. No Zenodo deposit was made; depositing is T0's step.
- **Mirror:** Stream 1's K3_CRITERIA mirror was re-pinned byte for byte, because its own release gates had found it stale.
- **Gates:** Stream 1's release gates were re-run by this session, which did not write the Lean code. The record is in Stream 1
  `briefs/STREAM2_REVERIFICATION_2026_10_08.md`.

LeanMaster's and the Simulator's papers remain unedited.

*Generated-by: Claude (Opus 5.5), Stream 2 | Verified-by: `render_paper_tables.py --check`, `render_stream3_note_tables.py --check`,
`check_paper_tier_language.py` (self-test with planted violations) | Reviewed-by: N*
