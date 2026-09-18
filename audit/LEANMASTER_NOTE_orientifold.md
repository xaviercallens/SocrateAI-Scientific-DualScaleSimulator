# Note for the LeanMaster session: orientifold and tadpole files

**From:** the DualScaleSimulator session, 2026-09-18.
**LeanMaster checked at:** `v3.19.0`, read-only; nothing in LeanMaster was edited.
**Tier:** literature (L). No Lean was built for this note.

## Sources (pinned)
Pinned on DualScaleSimulator branch `loop/k3t2-rigidity`, `audit/k3t2_rigidity_v2/sources/`, with sha256 in `SHA256SUMS`:
- **Tripathy & Trivedi (TT)**, *Compactification with flux on K3 and tori*, JHEP 03 (2003) 028, arXiv:hep-th/0301139.
- **Berkooz, Leigh, Polchinski, Schwarz, Seiberg & Witten (BLPSSW)**, *Anomalies, dualities, and topology of D=6 N=1 superstring vacua*, Nucl. Phys. B475 (1996) 115–148, arXiv:hep-th/9605184.

What they say (line numbers refer to the pdftotext `.txt` files):
- **TT §2.2, lines 160–163:** "The Z2 orientifold symmetry has 4 fixed points on the T², an O7-plane is located at each of these fixed points. To cancel the resulting 7−brane charge 16 D7-branes need to be added."
- **TT lines 164–171:** O7s and D7s wrap K3. This induces 2 units of D3 charge per O7 and 1 per D7, a total of 24. Their eq. (2.3) is **½ N_flux + N_D3 = 24**, where N_flux = (2π)⁻⁴(α′)⁻² ∫_{K3×T²} H₃∧F₃ is integrated over the covering space. The footnote says that in F-theory there are 24 (p,q) 7-branes, each with one unit.
- **TT lines 312–314:** in their symmetric configuration "each O7 plane has 4 D7 branes on it". So an O7-plane has charge −4 in D7 units. That value is arithmetic (16/4), not a quote.
- **BLPSSW, text lines 227–229:** at the Gimon–Polchinski orbifold point, "there must be one instanton hidden at each fixed point".
- **BLPSSW p. 10:** after blow-up, "eight 5-branes on a smooth K3". The total 16 + 8 = 24 is our inference.

Consequence: in the IIB orientifold on K3 × T²/ℤ₂, the orientifold acts on T², not on K3.
- The 4 O7-planes sit at the fixed points of T²/ℤ₂ and wrap K3.
- The 16 fixed points of T⁴/ℤ₂ belong to the orbifold limit of the K3 factor, and **no O7-plane sits there**.
- In the T-dual Type I / Gimon–Polchinski description, those 16 points carry O5-planes.

## Files that disagree with these sources

1. **`StringTheoryFoundation/StringTheory/TadpoleCancellation.lean`**
   - Its docstring and definitions put **16 O7⁻ planes at the T⁴/ℤ₂ fixed points** (`numFixedPointsT4Z2 = 16`, `chargeO7Minus = −4`, `total_O7_charge_is_minus_64`) and **32 D7-branes of charge +2** (`total_D7_charge_is_64`). The arithmetic is correct, but it describes no construction.
   - Its D3 condition (docstring near `d3TadpoleTarget`, line 63ff) reads "N_D3 + ½∫H₃∧F₃ = χ(K3)/24 = 24/24 = 1", and `d3_tadpole_target_is_one` proves 24/24 = 1. TT's eq. (2.3) has **24** on the right-hand side, not 1.
2. **`DualScaleM24Formalization/Moonshine/KummerTadpole.lean`**, lines 25–26 and 100–127
   - It uses 16 D7 stacks of charge **+4** and 4 O7-planes of charge **−16**, cites Gimon–Polchinski eq. (3.12) and Sen, and proves `rr_tadpole_cancellation` (64 − 64 = 0).
   - The plane count (4 O7, 16 D7) matches TT. The per-object charges are 4 times TT's counting (D7 = +1, O7 = −4). The total 16·4 − 4·16 = 0 still holds.
   - Please check which normalisation Gimon–Polchinski (3.12) actually uses, and state it in the docstring.
3. **Consistent:** `DualScaleStream2/Flux/Tadpole.lean:116` `tadpole_budget` (flux + n = 24). Its docstring sources this to Dasgupta–Rajesh–Sethi, M-theory on K3 × K3 (½∫G∧G + n = 24).
   - This has the same form as TT's eq. (2.3) with flux := ½ N_flux.
   - M-theory on K3 × K3 is dual to the IIB orientifold on K3 × T²/ℤ₂, so two independent literature routes give the same 24.
   - Suggest citing TT eq. (2.3) next to DRS.

## Suggested Lean-side changes (for the LeanMaster owner to decide)
- **In `TadpoleCancellation.lean`:**
  - Replace the 16-O7 model with `o7Count := 4` (the fixed points of T²/ℤ₂), `d7Count := 16`, D7 charge +1 and O7 charge −4, so that `16 + 4 * (-4) = 0`.
  - Add the induced D3 count `4 * 2 + 16 * 1 = 24`.
  - Delete or re-document `d3_tadpole_target_is_one`, since the target is 24.
- **Keep** `num_fixed_points_is_16` (2⁴ = 16 for T⁴/ℤ₂). Re-document it as the orbifold limit of K3 (the Kummer construction), not as orientifold planes.
- **In `KummerTadpole.lean`:** state the charge normalisation explicitly, or switch to TT's.

## In this repository
- `audit/PAPER_FACTS.md` F1 now records the resolution and the withdrawal.
- `T_duality_Alone.tex` and `leanflow_engine.tex` are corrected.
- The local `proofs/` files (`KummerTDAAnomalyCertification.lean`, `TadpoleCancellation.lean`) still encode the withdrawn assignment. `proofs/` is review-only here, because LeanMaster is the replacement.
