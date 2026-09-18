export const meta = {
  name: 'k3t2-rigidity-loop-v2',
  description: 'K3 x T2 blind re-derivation v2: sealed targets outside the worktree, structural (non-circular) rigidity chained across tracks, all LeanMaster targets compared, polar bug fixed, new flux-vacua counting track (Tripathy-Trivedi tadpole 24), two independent skeptics',
  whenToUse: 'Second independent route to LeanMaster K3 x T2 results; rerun whenever LeanMaster advances.',
  phases: [
    { title: 'Targets', detail: 'seal every LeanMaster K3xT2 statement outside the worktree' },
    { title: 'Blind compute', detail: '5 tracks, blind to LeanMaster and to v1 comparison files' },
    { title: 'Chain', detail: 'structural rigidity chain across tracks (no literal answers)' },
    { title: 'Compare', detail: 'every sealed target, per-row' },
    { title: 'Reverse', detail: 'Lean statements beyond their checked range' },
    { title: 'Skeptic', detail: 'two lenses: blindness/circularity and mathematical correctness' },
    { title: 'Synthesize', detail: 'ledger + report text' },
  ],
}

const WT = '/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2'
const OUT = `${WT}/audit/k3t2_rigidity_v2`
const SEALED = '/mnt/disks/disk-socrateai-local-1/k3t2-sealed-v2'
const PY = '/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python'
const LM = '/home/callensxavier_gmail_com/SocrateAI-Scientific-Agora-LeanMaster'
const SRC = `${OUT}/sources`

const RULES = `
GROUND RULES (strict):
- Work ONLY in git worktree ${WT} (branch loop/k3t2-rigidity), outputs under ${OUT}/<your-dir>/. Never touch proofs/, main, other worktrees; never push. Commit your own files with a message ending: Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
- Python: ${PY} (gudhi 3.13, sympy, mpmath, numpy, scipy, networkx, pytest). EXACT arithmetic (int, fractions.Fraction, sympy.Rational) for every integer/rational claim.
- Evidence-bound: every number comes from a script you wrote and ran; save script + JSON; state the exact command INCLUDING all arguments (v1 lost reproducibility because 'theta_forms.py 22' was undocumented). Scripts must be runnable from a clean checkout with the stated command.
- NO LITERAL ANSWERS IN A COMPUTATION PATH (the main v1 failure). A known value may appear only in an "expected" field filled after computing, with its source. A rigidity scan whose selecting condition is built from the true value (e.g. "== 24", a target series built with exponent 24) is NOT a rigidity result: label it "normalisation" and say which input fixes it. A scan counts as RIGID only if its selecting condition is structural (consistency, integrality, modularity, exactness, agreement of two independent computations) and a neighbouring value demonstrably fails.
- Tiers: B = exact arithmetic + negative control; L = literature identification; C = conjecture. Never write "proved".
`
const BLIND = `
BLINDNESS: do NOT open, grep, list or cat anything under ${LM}, ${SEALED}, proofs/, lean_foundation/, audit/PAPER_FACTS.md, or these v1 files that contain Lean values: ${WT}/audit/k3t2_rigidity/{targets_sealed.json,comparison.json,agreement_ledger.json,build_ledger.py,REPORT.md,reverse/,skeptic/}. You MAY read and reuse v1 blind track code in ${WT}/audit/k3t2_rigidity/{genus-moonshine,dyons,lattices-duality,tda-gudhi}/ (it was written blind) but must fix its known defects listed in your task. A later auditor greps your scripts and shell history.`

const TRACK = {
  type: 'object',
  properties: {
    track: { type: 'string' },
    results: { type: 'array', items: { type: 'object', properties: {
      id: { type: 'string' }, quantity: { type: 'string' }, computed: { type: 'string' }, method: { type: 'string' },
      script: { type: 'string' }, command: { type: 'string' } }, required: ['id', 'quantity', 'computed', 'script', 'command'] } },
    rigidity: { type: 'array', items: { type: 'object', properties: {
      parameter_inserted: { type: 'string' }, selecting_condition: { type: 'string' }, condition_uses_true_value: { type: 'boolean' },
      solution_set: { type: 'string' }, negative_control: { type: 'string' }, classification: { type: 'string', description: 'RIGID | NORMALISATION | CONDITIONAL_ON_INPUT' },
      input_it_depends_on: { type: 'string' } },
      required: ['parameter_inserted', 'selecting_condition', 'condition_uses_true_value', 'solution_set', 'classification'] } },
    exports: { type: 'object', description: 'computed quantities other tracks may consume, with file path', additionalProperties: true },
    could_not_do: { type: 'array', items: { type: 'string' } }, commit: { type: 'string' },
  },
  required: ['track', 'results', 'rigidity', 'could_not_do'],
}
const TARGETS = {
  type: 'object',
  properties: { tag: { type: 'string' }, commit: { type: 'string' }, file: { type: 'string' }, count: { type: 'number' },
    by_track: { type: 'object', additionalProperties: { type: 'number' } } },
  required: ['tag', 'file', 'count'],
}
const COMPARE = {
  type: 'object',
  properties: { rows: { type: 'array', items: { type: 'object', properties: {
    target_id: { type: 'string' }, lean_file_line: { type: 'string' }, lean_value: { type: 'string' }, blind_value: { type: 'string' }, blind_source: { type: 'string' },
    status: { type: 'string', description: 'AGREE | DISAGREE | NOT_COMPUTED | NOT_COMPARABLE | OUT_OF_SCOPE' }, independent_routes: { type: 'number' }, note: { type: 'string' } },
    required: ['target_id', 'status'] } },
    counts: { type: 'object', additionalProperties: { type: 'number' } }, blind_only: { type: 'array', items: { type: 'string' } } },
  required: ['rows', 'counts'],
}
const REVERSE = {
  type: 'object',
  properties: { predictions: { type: 'array', items: { type: 'object', properties: {
    from_theorem: { type: 'string' }, prediction: { type: 'string' }, computed: { type: 'string' }, holds: { type: 'boolean' },
    proposed_lean_statement: { type: 'string' }, lean_reachable_now: { type: 'boolean' }, script: { type: 'string' } },
    required: ['from_theorem', 'prediction', 'computed', 'holds'] } }, commit: { type: 'string' } },
  required: ['predictions'],
}
const AUDIT = {
  type: 'object',
  properties: {
    lens: { type: 'string' },
    violations: { type: 'array', items: { type: 'string' } },
    reruns: { type: 'array', items: { type: 'object', properties: { command: { type: 'string' }, reproduced: { type: 'boolean' }, note: { type: 'string' } }, required: ['command', 'reproduced'] } },
    disputed_rows: { type: 'array', items: { type: 'object', properties: { target_id: { type: 'string' }, reason: { type: 'string' } }, required: ['target_id', 'reason'] } },
    rigidity_verdicts: { type: 'array', items: { type: 'object', properties: { parameter: { type: 'string' }, genuinely_rigid: { type: 'boolean' }, reason: { type: 'string' } }, required: ['parameter', 'genuinely_rigid', 'reason'] } },
    math_errors: { type: 'array', items: { type: 'string' } }, commit: { type: 'string' },
  },
  required: ['lens', 'violations', 'reruns', 'disputed_rows', 'rigidity_verdicts'],
}

const TRACKS = [
  { key: 'A-genus', prompt: `TRACK A - K3 elliptic genus, Mathieu moonshine, twining.
Reuse/repair v1 genus-moonshine code. Compute exactly (through q^10 at least): Z_K3 = 2 phi_{0,1} from thetas (phi_{0,1} = 4 sum_{i=2,3,4} (theta_i(tau,z)/theta_i(tau,0))^2); Z(tau,0); c(4n-l^2) and its D-dependence; H(tau) = 2q^{-1/8}(-1 + sum A_n q^n) via the Appell-Lerch decomposition, A_1..A_10.
Twining: for classes 2A, 3A, 5A, 7A use H_g = (chi(g)/24) H - F_g / eta^3 with Lambda_N = N q d/dq log(eta(N tau)/eta(tau)); the literature forms (Cheng-Duncan-Harvey Table 3, FROM MEMORY, treat as tier-L inputs and say so): F_2A = 16 Lambda_2, F_3A = 6 Lambda_3, F_5A = 2 Lambda_5, F_7A = Lambda_7, with chi(2A)=8, chi(3A)=6, chi(5A)=4, chi(7A)=3 (chi = number of fixed points on 24 letters). STRUCTURAL CHECK that does not use these as answers: the twined coefficients must be integers and, together with the untwined A_n, must be consistent with class functions of M24 - at minimum check integrality, the sign pattern, and that A_n^(g) = A_n mod the order of g (congruence) for n <= 8; report failures. Test whether A_2*60 = 4*A_1*77 (the '27720 lock', with the convention you actually compute; state it) survives each twining.
Rigidity: (a) mu-term coefficient N in 20..28 selected by z-independence of H (two y-slices agree) - structural; (b) replace the overall factor k of Z = k phi_{0,1}: select k by a STRUCTURAL condition only (e.g. the q^0 y^{+-1} coefficient must equal 2 = number of ground states of a sigma model with h^{0,0}=h^{2,0}=1?) - if you cannot find a condition that avoids the known answer, classify NORMALISATION and say so. Export: chi_from_genus = Z(tau,0) to ${OUT}/A-genus/exports.json.` },
  { key: 'B-dyons', prompt: `TRACK B - quarter-BPS dyons on K3 x T2 (DMVV product; Dabholkar-Murthy-Zagier).
Reuse v1 dyons code BUT fix: (1) part3_polar.py:141 - the s<=-1 branch of A_{2,m} must keep the q^{m s^2 + s}... prefactor (v1 used n = -2s - sk; correct for m=1 is n = s^2 - s - s*k; derive the general-m form yourself from A_{2,m} = sum_s q^{m s^2+s} y^{2ms+1}/(1-q^s y)^2 expanded in |q|<|y|<1); (2) document every argument (theta_forms.py QMAX etc.); (3) scans must be JOINT, never coordinate-wise at the true value; (4) the Goettsche comparison must NOT be a table built from exponent 24: instead obtain the exponent from Track-A-style computation of c(0)+2c(-1) inside your own code, or compare G_k(z=0) with Euler numbers of Hilb^k(K3) computed from Goettsche's formula using chi(K3) taken from YOUR computed c(0)+2c(-1).
Compute: G_k(z=0), k=0..8; Delta*psi_m checks for m=-1..2 (DMZ (5.16) forms as tier-L ansatz: Delta psi_{-1}=1/A, Delta psi_0 = 2B/A, 4 Delta psi_1 = 9B^2/A + 3 E4 A); Hurwitz class numbers H(D) for D<=40 by reduced-form counting; immortal m=1: remainder of Delta psi_1 - N A_{2,1} vs 3E4A - M Hhat, JOINT (N,M) grid wide enough to show uniqueness on the FULL window; also show that with M*Hhat replaced by 0 the identity fails (i.e. class numbers are needed).
Rigidity: joint (N,M); c(3), c(4) shifted by -2..+2 each (not only +1); exponent multiplier kappa selected ONLY by integrality/structure (if only the Euler-number target selects it, classify NORMALISATION).` },
  { key: 'C-lattices', prompt: `TRACK C - lattices, T-duality, dual-scale bound, with a CHAINED structural signature.
Reuse v1 lattices-duality code. Add: symbolic sympy proofs-by-computation for d=1..4 (generic symbols, not samples) of: factorized T-duality squares to 1; basis change [[A,0],[0,A^{-T}]] preserves eta; eta H(G,0) eta = H(G^{-1},0); dualScale(G) := tr G + tr G^{-1} and its minimiser: for SPD G, tr G + tr G^{-1} - 2d = sum_i (x_i - 1)^2/x_i >= 0 with equality iff G = 1 (show symbolically in eigenvalues + a d=2 fully symbolic SPD check). E8 Cartan: det, positive definiteness (Sylvester), 240 roots. Mukai (4,20), Gamma^{6,22}. Tadpole budget: for naturals flux, n with flux + n = 24, enumerate and show n <= 24 and flux = 0 -> n = 24 (trivial but state it).
CHAINED RIGIDITY (the key improvement): do NOT type (3,19). Derive the K3 signature from topology: read chi = 24 and b2 = 22 from ${OUT}/D-tda/exports.json if it exists (poll the file for up to 20 minutes; else compute chi from Track A's ${OUT}/A-genus/exports.json; say which), then Noether's formula chi(O) = (c1^2 + c2)/12 with c1 = 0 (Calabi-Yau: state as definitional input), c2 = chi_top, giving chi(O) = 2 = 1 - h^{0,1} + h^{0,2}; with h^{0,1} = b1/2 = 0 get h^{2,0} = 1; Hodge index b2^+ = 2 h^{2,0} + 1; so (b2^+, b2^-) follows. Then the lattice enumeration over (m,n) with 2m+8n = 22 is selected by this DERIVED signature. Classify: RIGID-CONDITIONAL on the inputs you name (c1=0, Kahler) - honest.` },
  { key: 'D-tda', prompt: `TRACK D - topology of the Kummer K3 by GUDHI (INRIA), repaired.
Reuse v1 tda-gudhi code. Fix: (1) only use grid sizes where the singular closed stars are pairwise DISJOINT (v1: N=4 invalid, N=6 valid) - check the premise programmatically and add N=8 as a second valid data point if it runs in < 15 min; (2) read the number of fixed points and all Euler characteristics used in Mayer-Vietoris from computed results, no literals (v1 typed 16 and 32); (3) translation control: must pass the same 2-to-1 halving check as the negation quotient, or be dropped; (4) compute chi from the f-vector AND from Betti numbers and compare; (5) compute Betti numbers over Z/2, Z/3, Z/5 for the link (RP^3 torsion control).
Report b0..b4 of the resolved K3 (hybrid MV, labelled), chi, and K3 x T2 by Kunneth, and the rigidity scan over the number k of resolved points with the SELECTING CONDITION structural (e.g. the result must be a closed orientable 4-manifold whose b1 = 0 and whose middle homology has the parity needed for an even intersection form?) - if only chi == 24 selects it, classify NORMALISATION. Export chi and b2 to ${OUT}/D-tda/exports.json as soon as they are computed (Track C waits for it).` },
  { key: 'E-flux', prompt: `TRACK E - flux vacua on the Type IIB orientifold K3 x T2/Z2 (Tripathy-Trivedi, pinned at ${SRC}/hep-th_0301139_TripathyTrivedi.txt, sha256 in ${SRC}/SHA256SUMS; BLPSSW hep-th/9605184 also pinned there). This track is literature-driven (tier L inputs), not blind to LeanMaster by construction but still must not open ${LM} or ${SEALED}.
1. Quote (with line numbers of the .txt) TT's orientifold content (O7-planes, D7-branes, their location) and the D3 tadpole condition (their eq. 2.3) and the D3-charge bookkeeping that gives 24. Check the arithmetic 4*2 + 16*1 = 24 and the F-theory count (24 7-branes x 1).
2. Quote BLPSSW's statement about instantons hidden at the Gimon-Polchinski fixed points and the five-brane count after blow-up; check 16 + 8 = 24 and state clearly whether BLPSSW literally states the sum or it is an inference.
3. COUNT flux vacua in the simplest tractable setting TT treat: read their parametrisation of G = F - tau H (flux vectors in H^2(K3) x H^1(T2)-type lattice, their sections 3-4, e.g. the alpha/beta vectors and the tadpole contribution N_flux = alpha^2 x ... as in their eq. (4.14)) and write it down explicitly with line references. Then, within an explicitly stated finite truncation (e.g. flux vectors restricted to a small rank-2 or rank-4 sublattice such as U + U or U + (-E8)-rank-2 slice, entries bounded), ENUMERATE integer flux choices satisfying their supersymmetry conditions and N_flux + N_D3 = 24 with N_D3 >= 0, modulo the obvious symmetries you can justify (state them). Report counts per N_flux value. This is a finite discrete landscape: the point is to show the continuous moduli are traded for a finite list of integers, and to measure how large that list is in the truncation. Negative control: drop the tadpole bound and show the count grows without limit as the entry bound grows.
4. State which moduli remain unfixed per TT (their abstract: 'some of the Kahler moduli').
Tier: counts are tier B within the stated truncation; the physics identification tier L. Put the explicit truncation in every reported number.` },
]

phase('Targets')
const targetsP = agent(`${RULES}
You are the ONLY agent allowed to read LeanMaster (${LM}, READ-ONLY). Write the sealed target list to ${SEALED}/targets_sealed_v2.json (this directory is OUTSIDE the worktree on purpose; do NOT write it anywhere under ${WT}; do not commit it). Run git -C ${LM} describe --tags and rev-parse HEAD.
Include EVERY theorem (verbatim statement, file:line, the concrete numbers it asserts, its checked order/range, committed yes/no) from: DualScaleMoonshine/*, DualScaleDyons/*, DualScaleStream2/{Lattice,TDuality,DFT,DualScale,Flux}/*, DoubleFieldTheory/K3Topology.lean, and any Kummer/Betti/Euler statements elsewhere (grep). Skip statements that are only about the simulator's Mapper constants (list them separately as out_of_scope). Assign each to a track: A-genus, B-dyons, C-lattices, D-tda, E-flux (E: tadpole statements). Also record, as a separate list 'orientifold_statements', any statement about O7/O5/D7/D5 charges or D3 tadpoles, verbatim, for the orientifold audit.
Use grep -n "^theorem" then read each statement fully (multi-line). Do not paraphrase. Return only the summary (tag, file, count, by_track).`,
  { label: 'targets:seal-all', phase: 'Targets', model: 'sonnet', schema: TARGETS })

phase('Blind compute')
const blindP = parallel(TRACKS.map(t => () => agent(`${RULES}${BLIND}
${t.prompt}
Write to ${OUT}/${t.key}/ (scripts, results.json, exports.json when asked). Keep each script under ~15 min (lower truncation and say so rather than hang). Commit. Put everything you could not do in could_not_do.`,
  { label: `blind:${t.key}`, phase: 'Blind compute', model: 'sonnet', schema: TRACK })))

const [targets, blindRaw] = await Promise.all([targetsP, blindP])
const blind = (blindRaw || []).filter(Boolean)
log(`sealed targets: ${targets ? targets.count : 0}; blind tracks: ${blind.length}/5`)
if (!targets || blind.length === 0) return { stopped: 'no targets or blind results', targets, blind }

phase('Chain')
const chain = await agent(`${RULES}${BLIND}
CHAIN CHECK. The v2 design is that rigidity should come from chaining tracks, not from typed answers:
 D (GUDHI) -> chi, b2  ->  C (Noether + Hodge index) -> signature -> lattice (3U + 2(-E8))  ; A (elliptic genus) -> Z(tau,0) must equal D's chi ; B (DMVV) -> Euler numbers of Hilb^k must follow from A's c(0)+2c(-1) and D's chi.
Track results: ${JSON.stringify(blind.map(b => ({ track: b.track, exports: b.exports, rigidity: b.rigidity })))}
Open the exports.json files under ${OUT}/*/ and verify each link numerically with a small script ${OUT}/chain/chain_check.py (read values from the files, no literals). Report for each link: consistent or not, and which inputs are still definitional (e.g. c1 = 0). Commit. Return a short plain-text summary.`,
  { label: 'chain:cross-track', phase: 'Chain', model: 'sonnet' })

phase('Compare')
const compare = await agent(`${RULES}
COMPARATOR (you may read LeanMaster read-only and ${SEALED}). Sealed file: ${targets.file} (${targets.count} targets). Blind results: ${JSON.stringify(blind.map(b => ({ track: b.track, results: b.results.map(r => ({ id: r.id, quantity: r.quantity, computed: r.computed, script: r.script })) })))}
CHAIN: ${chain}
Make a row for EVERY sealed target (v1 rowed only 12 of 32 - that is the defect to fix). Open the blind results.json files, not only the summaries. Statuses: AGREE / DISAGREE / NOT_COMPUTED / NOT_COMPARABLE (convention differs; give the conversion and whether it matches under it) / OUT_OF_SCOPE (simulator constants). independent_routes = number of genuinely independent blind computations behind a row (tracks computing the same theta quantity count as one). Investigate any DISAGREE with a third minimal computation. Write ${OUT}/comparison.json and commit it (it contains Lean values; that is fine now).`,
  { label: 'compare:all-targets', phase: 'Compare', schema: COMPARE })

phase('Reverse')
const reverse = await agent(`${RULES}
REVERSE PASS (you may read LeanMaster and all of ${OUT}). Take 5-8 Lean theorems limited to a finite range and test them beyond it with exact arithmetic, EXCLUDING the five v1 already did (p24 k=6..8, DMZ m=4, polar m=4,5, trace bound d=7..10). Candidates: immortal dyons m=4 structure (-p24(5) H|V_4 or whatever the Lean pattern generalises to - for non-prime m say what the Hecke-like operator should be and test), twined Goettsche characters for classes beyond those Lean covers, elliptic genus index-1 property to q^20, umbral exactness at more levels, tadpole statements vs Tripathy-Trivedi (4 O7 at T2/Z2 fixed points, 16 D7, N_flux + N_D3 = 24). For each: prediction, computation, holds or the minimal counterexample, and a decide-checkable Lean-style statement. Save under ${OUT}/reverse/ and commit.
COMPARISON: ${JSON.stringify(compare && compare.counts)}`,
  { label: 'reverse:beyond-range', phase: 'Reverse', model: 'sonnet', schema: REVERSE })

phase('Skeptic')
const ctx = `BLIND: ${JSON.stringify(blind)}
CHAIN: ${chain}
COMPARISON: ${JSON.stringify(compare)}
REVERSE: ${JSON.stringify(reverse)}`
const skeptics = await parallel([
  () => agent(`${RULES}
SKEPTIC, lens = BLINDNESS AND CIRCULARITY. You produced none of this; default to distrust.
${ctx}
1. grep all scripts under ${OUT}/{A-genus,B-dyons,C-lattices,D-tda,E-flux,chain}/ and 'git log -p' of the v2 commits for LeanMaster paths, ${SEALED}, the v1 forbidden files, and for literal answer values in computation paths (24, 22, 16, 32, (3,19), 324, 648, 90, 462, 45, 231, 1, 24, 324, 3200, 25650, 176256, -128, 216, 240, ...). A literal only in an 'expected' field filled after computing is fine. List file:line.
2. For every rigidity entry, decide genuinely rigid or not, reading the code of the selecting condition. The chained signature is rigid only conditionally: say on what.
3. Rerun every stated command from a clean state (git stash nothing; run in a copy under /tmp if needed) and report reproduced or not.
Write ${OUT}/skeptic_blindness/verdict.json and commit.`, { label: 'skeptic:blindness', phase: 'Skeptic', schema: AUDIT }),
  () => agent(`${RULES}
SKEPTIC, lens = MATHEMATICAL CORRECTNESS. You produced none of this; default to distrust.
${ctx}
Check the mathematics, not the bookkeeping: (a) conventions (A_n vs 2A_n, H(0) = -1/12, strips for A_{2,m}, signature sign conventions, twining forms F_g and chi(g) values used in Track A, which were given FROM MEMORY - verify them structurally or flag); (b) Noether/Hodge-index chain in Track C (is c1 = 0 and b1 = 0 used correctly? is h^{0,1} = b1/2 legitimate here?); (c) the Mayer-Vietoris step in Track D (premises checked at each N?); (d) Track E: is the flux parametrisation faithful to Tripathy-Trivedi (compare to the pinned text line by line), are the supersymmetry conditions and tadpole contribution correct, is the truncation stated, is the symmetry quotient justified; does the BLPSSW 16+8 inference hold as stated; (e) every DISAGREE or NOT_COMPARABLE row: who is right. Re-derive at least one number per track by an independent minimal script.
Write ${OUT}/skeptic_math/verdict.json and commit.`, { label: 'skeptic:math', phase: 'Skeptic', schema: AUDIT }),
])

phase('Synthesize')
const report = await agent(`${RULES}
Produce the v2 REPORT as your final text (the harness blocks subagents from writing .md files; the orchestrator will save it). Also write ${OUT}/agreement_ledger_v2.json with a small builder script and commit those two.
Inputs: TARGETS ${JSON.stringify(targets)}; CHAIN ${chain}; COMPARISON ${JSON.stringify(compare)}; REVERSE ${JSON.stringify(reverse)}; SKEPTICS ${JSON.stringify(skeptics)}; BLIND rigidity/could_not_do ${JSON.stringify(blind.map(b => ({ track: b.track, rigidity: b.rigidity, could_not_do: b.could_not_do })))}
Sections: 1 Bottom line (counts after removing skeptic-disputed rows; compare with v1: 5 clean / 3 qualified / 0 disagree of 12 rowed out of 32; rigidity v1 3/11). 2 Per-track table. 3 Rigidity: RIGID vs NORMALISATION vs CONDITIONAL, with the inputs each depends on; one sentence on what 'zero free parameters' means here. 4 The chain D->C->lattice and A/B/D consistency. 5 Flux vacua (Track E): the finite list, its size in the truncation, what stays unfixed, the orientifold facts with line references. 6 Reverse predictions (held/failed, Lean statements ready to hand over). 7 Everything not done, disputed, from memory. Header: "Generated by workflow k3t2-rigidity-loop-v2; computations tier B at best; physical identifications tier L/C."`,
  { label: 'synthesize:report-v2', phase: 'Synthesize' })

return { targets, chain, counts: compare && compare.counts, reverse, skeptics, report }
