export const meta = {
  name: 'k3t2-rigidity-loop-v3',
  description: 'K3 x T2 blind re-derivation v3: reseal at LeanMaster HEAD (v3.25+, Stream 8 included), v2 skeptic corrections applied (declared normalisations, route independence, realizable controls, TT unimodular lattice, script-generated results, relative paths), new blind Track F (which K3: T2 self-dual point, T(A), Kummer code vs Golay), reverse pass reviewed by skeptics',
  whenToUse: 'Second independent route to LeanMaster K3 x T2 results; rerun whenever LeanMaster advances. Supersedes k3t2-rigidity-loop-v2 (kept as the record of the 2026-09-19 run).',
  phases: [
    { title: 'Targets', detail: 'seal every LeanMaster K3xT2 statement at HEAD, outside the worktree; diff against the v2 seal' },
    { title: 'Blind compute', detail: '6 tracks (A-F), blind to LeanMaster and to all v1/v2 comparison files' },
    { title: 'Chain', detail: 'cross-track links, each labelled CROSS_METHOD or POINTER, runnable from a clean clone' },
    { title: 'Compare', detail: 'every sealed target, per row, independence counted over disjoint typed inputs' },
    { title: 'Reverse', detail: 'Lean statements beyond their checked range (excluding the 8 already done)' },
    { title: 'Skeptic', detail: 'blindness/circularity/reproducibility and mathematical correctness, both also review the reverse pass' },
    { title: 'Synthesize', detail: 'ledger builder + report text (orchestrator saves REPORT.md)' },
  ],
}

// args (optional): { only?: string[] track keys, e.g. ['F-whichk3'] }
const ONLY = (args && Array.isArray(args.only)) ? args.only : null

const WT = '/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2'
const OUT = `${WT}/audit/k3t2_rigidity_v3`
const V2 = `${WT}/audit/k3t2_rigidity_v2`
const SEALED = '/mnt/disks/disk-socrateai-local-1/k3t2-sealed-v3'
const SEALED_V2 = '/mnt/disks/disk-socrateai-local-1/k3t2-sealed-v2'
const PY = '/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python'
const LM = '/home/callensxavier_gmail_com/SocrateAI-Scientific-Agora-LeanMaster'
const SRC = `${V2}/sources`

const RULES = `
GROUND RULES (strict; v3 adds items 6-12, each one a defect the v2 skeptics found):
1. Work ONLY in git worktree ${WT} (branch loop/k3t2-rigidity), outputs under ${OUT}/<your-dir>/. Never touch proofs/, main, other worktrees; never push.
2. Python: ${PY} (gudhi 3.13, sympy, mpmath, numpy, scipy, networkx, pytest). EXACT arithmetic (int, fractions.Fraction, sympy.Rational / algebraic numbers) for every integer/rational claim.
3. Evidence-bound: every number comes from a script you wrote and ran; save script + JSON.
4. NO LITERAL ANSWERS IN A COMPUTATION PATH. A known value may appear only in an "expected" field filled after computing, with its source.
5. Tiers: B = exact arithmetic + negative control; L = literature identification; C = conjecture.
6. WORDING: never write "proved", "proof", "proof-by-computation" or "theorem" for your own work (v2: E, C and B did). Write "verified by exact computation (tier B)" or "holds for every case enumerated".
7. DECLARED INPUTS: every typed normalisation or structural input (e.g. an overall factor k in Z = k*phi_{0,1}, a typed Gram matrix, a fixed-point count, a quoted literature number) goes in ${OUT}/<your-dir>/inputs.json with {name, value, tier, why}. Every result lists shared_inputs: the names of the declared inputs it depends on. A value obtained as (typed input) x (computed number) is NOT an independent derivation of that value (v2: 24 = typed 2 x computed 12 was reported as independent). Downstream code must use the solution a scan returns, never re-type it (v2: A typed 24 at run_track_a_v2.py:322 instead of the scan's solution).
8. RIGIDITY LABELS, exactly one of: RIGID (structural selecting condition that does not contain the target, and a neighbouring value demonstrably fails) | RIGID_GIVEN_DEFINITION (RIGID once named defining inputs are granted) | CONDITIONAL_ON_INPUT | NORMALISATION (the selecting condition contains the target) | VERIFIED_IDENTITY (a universal identity; nothing is selected) | INVARIANCE (the answer does not change across a nuisance choice) | NON_DISCRIMINATING (every value passes). The negative control must perturb THE SAME parameter the selection is about (v2: E's rank-2 emptiness control tested the signature, not the factor 8). Scan grids must not be centred on the answer and must step through every integer (v2: B's grid was centred on (324,648) and sampled even M only); prefer an exact solve over a grid.
9. CONTROLS MUST BE REALIZABLE: a counterfactual input used as a control must correspond to an object that exists, or be labelled HYPOTHETICAL and not counted as discrimination (v2: C's chi in {36,48} controls used surfaces that do not exist).
10. REPRODUCIBILITY: results.json and exports.json must be WRITTEN BY A SCRIPT, never hand-authored. Paths must be relative to the repo root, found from the script location with pathlib (Path(__file__).resolve().parents[k]). No absolute /mnt/... paths inside scripts. State each command as: cd <dir relative to repo root> && ${PY} <script> <all arguments>. Docstrings must state the same arguments as the committed run (v2: docstrings said 16, the run used 40).
11. COMMITS: stage ONLY your own files by explicit path (git add <paths>), run git diff --cached --name-only and check it before committing; never git add -A or git add . (v2: a shared-index sweep put Track C's files in a Track D commit). Message ends with: Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
12. Tier-L inputs quoted from memory must say FROM MEMORY. Where a structural re-derivation is possible, do it (for example, the twining data from frame shapes, as below).
`

const BLIND = `
BLINDNESS: do NOT open, grep, list or cat anything under ${LM}, ${SEALED}, ${SEALED_V2}, proofs/, lean_foundation/, audit/PAPER_FACTS.md, or files that contain Lean values:
  v1: ${WT}/audit/k3t2_rigidity/{targets_sealed.json,comparison.json,agreement_ledger.json,build_ledger.py,REPORT.md,reverse/,skeptic/};
  v2: ${V2}/{comparison.json,comparison/,agreement_ledger_v2.json,build_ledger_v2.py,REPORT.md,reverse/,skeptic_blindness/,skeptic_math/}.
You MAY read and reuse the v2 blind track code in ${V2}/{A-genus,B-dyons,C-lattices,D-tda,E-flux,chain}/ and the v1 blind code in ${WT}/audit/k3t2_rigidity/{genus-moonshine,dyons,lattices-duality,tda-gudhi}/ (written blind). You must fix the defects listed in your task. A later auditor greps your scripts and shell history.`

const TRACK = {
  type: 'object',
  properties: {
    track: { type: 'string' },
    inputs: { type: 'array', items: { type: 'object', properties: { name: { type: 'string' }, value: { type: 'string' }, tier: { type: 'string' }, why: { type: 'string' } }, required: ['name', 'value', 'tier'] } },
    results: { type: 'array', items: { type: 'object', properties: {
      id: { type: 'string' }, quantity: { type: 'string' }, computed: { type: 'string' }, method: { type: 'string' },
      shared_inputs: { type: 'array', items: { type: 'string' } },
      script: { type: 'string' }, command: { type: 'string' } }, required: ['id', 'quantity', 'computed', 'shared_inputs', 'script', 'command'] } },
    rigidity: { type: 'array', items: { type: 'object', properties: {
      parameter_inserted: { type: 'string' }, selecting_condition: { type: 'string' }, condition_uses_true_value: { type: 'boolean' },
      solution_set: { type: 'string' }, negative_control: { type: 'string' }, control_perturbs_same_parameter: { type: 'boolean' },
      classification: { type: 'string', description: 'RIGID | RIGID_GIVEN_DEFINITION | CONDITIONAL_ON_INPUT | NORMALISATION | VERIFIED_IDENTITY | INVARIANCE | NON_DISCRIMINATING' },
      input_it_depends_on: { type: 'string' } },
      required: ['parameter_inserted', 'selecting_condition', 'condition_uses_true_value', 'solution_set', 'classification'] } },
    literature_checks: { type: 'array', description: 'quotation arithmetic on literature numbers (tier L), kept apart from computed results', items: { type: 'string' } },
    exports: { type: 'object', description: 'computed quantities other tracks may consume, with file path', additionalProperties: true },
    could_not_do: { type: 'array', items: { type: 'string' } }, commit: { type: 'string' },
  },
  required: ['track', 'inputs', 'results', 'rigidity', 'could_not_do'],
}
const TARGETS = {
  type: 'object',
  properties: { tag: { type: 'string' }, commit: { type: 'string' }, file: { type: 'string' }, count: { type: 'number' },
    by_track: { type: 'object', additionalProperties: { type: 'number' } },
    changed_since_v2_seal: { type: 'array', items: { type: 'string' }, description: 'statement ids added/removed/changed relative to the v2 seal (ids only, no values)' } },
  required: ['tag', 'file', 'count'],
}
const COMPARE = {
  type: 'object',
  properties: { rows: { type: 'array', items: { type: 'object', properties: {
    target_id: { type: 'string' }, lean_file_line: { type: 'string' }, lean_value: { type: 'string' }, blind_value: { type: 'string' }, blind_source: { type: 'string' },
    status: { type: 'string', description: 'AGREE | DISAGREE | NOT_COMPUTED | NOT_COMPARABLE | OUT_OF_SCOPE' },
    independent_routes: { type: 'number' }, shared_inputs: { type: 'array', items: { type: 'string' } }, post_seal_changed: { type: 'boolean' }, note: { type: 'string' } },
    required: ['target_id', 'status'] } },
    counts: { type: 'object', additionalProperties: { type: 'number' } }, blind_only: { type: 'array', items: { type: 'string' } } },
  required: ['rows', 'counts'],
}
const REVERSE = {
  type: 'object',
  properties: { predictions: { type: 'array', items: { type: 'object', properties: {
    from_theorem: { type: 'string' }, prediction: { type: 'string' }, computed: { type: 'string' }, holds: { type: 'boolean' },
    independence: { type: 'string', description: 'what the check shares with the Lean statement or with literature quotes; do not overstate' },
    proposed_lean_statement: { type: 'string' }, lean_reachable_now: { type: 'boolean' }, script: { type: 'string' } },
    required: ['from_theorem', 'prediction', 'computed', 'holds', 'independence'] } }, commit: { type: 'string' } },
  required: ['predictions'],
}
const AUDIT = {
  type: 'object',
  properties: {
    lens: { type: 'string' },
    violations: { type: 'array', items: { type: 'string' } },
    reruns: { type: 'array', items: { type: 'object', properties: { command: { type: 'string' }, reproduced: { type: 'boolean' }, clean_clone: { type: 'boolean' }, note: { type: 'string' } }, required: ['command', 'reproduced'] } },
    disputed_rows: { type: 'array', items: { type: 'object', properties: { target_id: { type: 'string' }, reason: { type: 'string' } }, required: ['target_id', 'reason'] } },
    rigidity_verdicts: { type: 'array', items: { type: 'object', properties: { parameter: { type: 'string' }, genuinely_rigid: { type: 'boolean' }, final_label: { type: 'string' }, reason: { type: 'string' } }, required: ['parameter', 'genuinely_rigid', 'reason'] } },
    reverse_verdicts: { type: 'array', items: { type: 'object', properties: { item: { type: 'string' }, stands: { type: 'boolean' }, reason: { type: 'string' } }, required: ['item', 'stands', 'reason'] } },
    math_errors: { type: 'array', items: { type: 'string' } }, commit: { type: 'string' },
  },
  required: ['lens', 'violations', 'reruns', 'disputed_rows', 'rigidity_verdicts', 'reverse_verdicts'],
}

const TRACKS = [
  { key: 'A-genus', prompt: `TRACK A - K3 elliptic genus, Mathieu moonshine, twining.
Reuse ${V2}/A-genus. Compute exactly (through q^12 at least): phi_{0,1} = 4 sum_{i=2,3,4} (theta_i(tau,z)/theta_i(tau,0))^2; Z_K3 = k*phi_{0,1}; Z(tau,0); c(4n-l^2) and its D-dependence; H(tau) = 2q^{-1/8}(-1 + sum A_n q^n) via the Appell-Lerch decomposition, A_1..A_12.
V2 DEFECTS TO FIX: (a) k is a DECLARED INPUT (rule 7). Look for a structural selector for k that does not contain Z(tau,0) or the ground-state count; if none exists, keep NORMALISATION and propagate k into shared_inputs of every downstream value (chi_from_genus, c(D), A_n, twinings). (b) The mu-term scan: report its solution as a function of k (run k = 1, 2, 3), and label it honestly; downstream H must use the scan's returned solution, not a typed number. (c) Twining data: DO NOT take chi(g) and F_g from memory. Derive them structurally: from the cycle shape (frame shape) of g as a permutation of 24 points, chi(g) = number of 1-cycles; the eta-product eta_g(tau) = prod eta(a tau)^{m_a}; use dim M_2(Gamma_0(N)) and the requirement that the twined genus has the right ground-state term to fix F_g inside the weight-2 space. Obtain the cycle shapes from M24 as built in Track F if its exports exist (poll ${OUT}/F-whichk3/exports.json up to 20 min), else generate the permutation group yourself from the extended Golay code (quadratic residues mod 23; see Track F) and read cycle shapes of elements of orders 2, 3, 5, 7. (d) Do the sign-pattern check that v2 left open. (e) Twined coefficients for 2A, 3A, 5A, 7A (and 4B, 11A if the structural derivation reaches them): integrality and A_n^(g) = A_n mod ord(g) for n <= 12. (f) Test the A_2*60 = 4*A_1*77 lock under each twining, with the convention you compute.
Export chi_from_genus WITH its shared_inputs (k) to ${OUT}/A-genus/exports.json.` },
  { key: 'B-dyons', prompt: `TRACK B - quarter-BPS dyons on K3 x T2 (DMVV product; Dabholkar-Murthy-Zagier).
Reuse ${V2}/B-dyons. V2 DEFECTS TO FIX: (1) results.json must be written by a script (v2 was hand-authored). (2) The typed factor 2 in the DMVV input (theta_forms.py:135, twoB) is a DECLARED INPUT; propagate it. (3) Euler numbers of Hilb^k(K3): EVALUATE THE DMVV PRODUCT AT z = 0 yourself (1/prod (1-p^m q^n)^{c(4mn)} truncated), and compare with Goettsche's formula prod (1-p^m)^{-chi} using chi from Track D's exports. v2 compared the Goettsche formula with itself (part1_euler_numbers.py:88-90), which is circular. (4) Immortal m=1 identity Delta*psi_1 - N*A_{2,1} = 3E4A - M*Hhat: solve for (N, E4A coefficient, M) EXACTLY as an over-determined linear system (report rank and residual), then show a neighbouring integer fails. No centred grid. (5) Extend: Delta*psi_m for m = 3 (DMZ-style line) if reachable; Hurwitz class numbers H(D) by reduced-form counting for D <= 400, with an independent Dirichlet L-function cross-check for fundamental discriminants (h(D) = -(1/|D|) sum chi_D(a) a).
Commands must use ${PY}, with every argument stated and the docstrings in agreement.` },
  { key: 'C-lattices', prompt: `TRACK C - lattices, T-duality, dual-scale bound, chained signature.
Reuse ${V2}/C-lattices. V2 DEFECTS TO FIX: (1) T-duality and generalized-metric identities are VERIFIED_IDENTITY, not RIGID. (2) The typed Gram matrix 3U + 2(-E8) is a DECLARED INPUT; its signature and rank are not independent routes. (3) chi controls must be REALIZABLE (rule 9). A compact complex surface with c1 = 0 and b1 = 0 has K = O, so chi(O) = 2 and chi_top = 24 follows by Noether. Present the chain as: K = O (the definition of K3), b1 = 0, Serre duality and Noether (tier L), giving chi_top; then Hodge index gives (b2+, b2-); the only discriminating controls are the signatures of candidate lattices mU + n(-E8) that do exist. Label the result RIGID_GIVEN_DEFINITION. (4) The symbolic tie-back for the dual-scale bound tr G + tr G^-1 >= 2d at d = 1..4, entry-level, not only d = 2. (5) Redo the lambda-in-R T-duality check instead of carrying it from v1.
Export chi_top and signature, each with shared_inputs, to ${OUT}/C-lattices/exports.json.` },
  { key: 'D-tda', prompt: `TRACK D - topology of the Kummer K3 by GUDHI (INRIA).
Reuse ${V2}/D-tda. V2 DEFECTS TO FIX: (1) Use the OPEN-star disjointness criterion for the singular points. Mayer-Vietoris needs disjoint open neighbourhoods; v2's closed-star criterion was stricter than needed and wrongly rejected N = 4. Check the premise at N = 4, 6, 8 and verify that the link of each singular point has RP^3 homology over Z/2, Z/3 and Z/5. (2) COMPUTE b3 and b4 from the MV sequence, including the rank of the H_3 map, instead of stating them. (3) The number of singular points and every Euler characteristic come from computed results; no literals. (4) Try to evaluate the parity of the intersection form: build a chain-level cup/intersection pairing on H_2 over Z/2 if feasible within 30 min, else record it in could_not_do. (5) Relative paths only.
Report b0..b4, chi (from the f-vector and from Betti numbers), K3 x T2 by Kunneth, and the k-resolved-points scan with honest labels. Export chi, b2 and n_singular to ${OUT}/D-tda/exports.json AS SOON AS they are computed (Tracks B, C, E and F read them).` },
  { key: 'E-flux', prompt: `TRACK E - flux vacua on the Type IIB orientifold K3 x T2/Z2 (Tripathy-Trivedi, pinned at ${SRC}/hep-th_0301139_TripathyTrivedi.txt, sha256 in ${SRC}/SHA256SUMS; BLPSSW hep-th/9605184 is pinned there too). This track is literature-driven (tier L inputs). It must still not open ${LM} or any sealed file.
Reuse ${V2}/E-flux. V2 DEFECTS TO FIX: (1) Work PRIMARILY in TT's unimodular presentation of the flux lattice ((A.2)-(A.3), U^3 type). Compute |det| of the span of TT's (A.6) forms relative to it. If the index is > 1, the (A.6) counts are counts of a SUBFAMILY: say so, and give counts in the unimodular lattice as the headline. Any divisibility statement about alpha^2 must be verified in the unimodular lattice, with a control that perturbs the lattice's evenness (e.g. an odd lattice). (2) TT line ~941 ("N_flux = 24 and no D3-branes") against eq. (2.3) with the 1/2: decide which is right by computation from TT's own worked examples (4.13)-(4.32), with line numbers. (3) Quotation arithmetic (4*2 + 16*1, 24*1, the BLPSSW 24 - 8) goes in literature_checks (tier L), NOT in results. Read the number of T^4/Z2 singular points from Track D's exports (poll up to 20 min). (4) Symmetry quotient: name the group and justify that it is a subgroup of O(Gamma) preserving the supersymmetry conditions and the orientation of the positive block. Call orbit counts "counts modulo G", never "upper bounds on physically distinct vacua". (5) Extend to TT's section 4.2 and 5 families if the parametrisation is explicit; state the truncation in every number. (6) The negative control (no tadpole bound implies unbounded growth) stays.
State which moduli remain unfixed, quoting TT with line numbers.` },
  { key: 'F-whichk3', prompt: `TRACK F (NEW) - "which K3": the self-dual T2, abelian surfaces and the Kummer code. All exact arithmetic. Every identification with a named object is tier L; every computed number is tier B.
1. T2 MODULI. Moduli (tau, rho) in H x H. Duality group: SL(2,Z)_tau x SL(2,Z)_rho, extended by the exchange tau <-> rho and by (tau, rho) -> (-tau-bar, -rho-bar). (a) Find the points of the fundamental domain whose stabiliser in SL(2,Z) is larger than {+-1}, by exact search over reduced matrices with bounded entries, and give each stabiliser order. (b) Among the points (tau, rho), find those whose stabiliser in the full duality group is maximal, and give its order. (c) At each such point, build the Narain lattice Gamma^{2,2}, enumerate the vectors with p_R = 0 and p_L^2 = 2, and identify the enhanced gauge root system from the Cartan matrix of a simple system (exact). Nothing typed: the point and the algebra are outputs.
2. ABELIAN SURFACE. For A = E_tau x E_tau' at the tau values found in 1(a), compute NS(A) and T(A) inside H^2(A,Z) = U^3 by exact linear algebra over Q(sqrt(-3)) or Q(i). The period is the wedge of the holomorphic 1-forms; NS = integral classes orthogonal to it; T = the orthogonal complement of NS. Report the Gram matrix of T(A), its discriminant, and whether it is (a scaling of) a root lattice. Then T(Km A) = T(A)(2) (Nikulin, tier L) gives the Kummer surface's transcendental form and discriminant. Enumerate reduced positive-definite even binary forms of discriminant <= 50 and locate T(A) and T(Km A) among them.
3. T4 LATTICES. Among rank-4 root lattices (A4, D4, A3+A1, A2+A2, A2+A1+A1, A1^4), compute |Aut(L)| exactly by enumerating automorphisms that preserve the minimal-vector set, and report which is maximal. For the maximal one, find complex structures on R^4/L compatible with an order-3 or order-4 automorphism, and compute T(A) for the resulting abelian surface as in 2.
4. KUMMER CODE AND GOLAY CODE. (a) The 16 nodes = F_2^4. The Kummer code K is spanned by the indicator vectors of the affine hyperplanes of F_2^4. Compute dim K and its weight enumerator. (b) Build the extended binary Golay code G24 from the quadratic residues mod 23. Compute its weight enumerator and the number of weight-8 words. (c) Build M24 as the permutation group preserving G24. Generate it by PSL(2,23) acting on the projective line P^1(F_23) plus one more permutation that preserves G24, found by search and checked. Give its order via Schreier-Sims (sympy.combinatorics) and the cycle shapes of elements of orders 2, 3, 4, 5, 6, 7, 8, 11, 12, 14, 15, 21, 23. Export the cycle shapes for Track A. (d) Test whether 16 of the 24 coordinates can be chosen so that the weight-8 words of G24 supported inside those 16 restrict to a code equivalent to K (exact search over octads and complements). Report the answer and the number of such choices.
5. SYMPLECTIC ACTIONS. On H^*(T^4, Z), compute the Lefschetz numbers of the translations by 2-torsion points and of -1, and hence the fixed-point counts on the Kummer surface (the Lefschetz fixed-point formula is tier L). Compare with the orders and cycle shapes from 4(c) where an identification is natural, and label that identification tier C.
Keep each script under 20 min; lower the bounds and say so rather than hang. Export to ${OUT}/F-whichk3/exports.json as soon as each part finishes (Track A reads the cycle shapes).` },
]
const ACTIVE = ONLY ? TRACKS.filter(t => ONLY.includes(t.key)) : TRACKS
if (ONLY) log(`running only tracks: ${ACTIVE.map(t => t.key).join(', ')} (others skipped)`)

phase('Targets')
const targetsP = agent(`${RULES}
You are the ONLY agent allowed to read LeanMaster (${LM}, READ-ONLY) and the v2 seal. Write the sealed target list to ${SEALED}/targets_sealed_v3.json. That directory is OUTSIDE the worktree on purpose: do NOT write it anywhere under ${WT}, and do not commit it. Record git -C ${LM} describe --tags and rev-parse HEAD in the file. Expect v3.25.0 or later.
Include EVERY theorem (verbatim statement, file:line, the concrete numbers it asserts, its checked order/range, committed yes/no) from: DualScaleMoonshine/*, DualScaleDyons/* (this includes Stream 8: WhichK3, SelfDualT2, KummerE3, ForgerE4, KummerD4, and any newer files), DualScaleStream2/{Lattice,TDuality,DFT,DualScale,Flux}/*, DualScaleM24Formalization/Moonshine/KummerTadpole.lean, StringTheoryFoundation/StringTheory/{TadpoleCancellation,TadpoleConstraint}.lean, DoubleFieldTheory/K3Topology.lean, and any Kummer/Betti/Euler/Golay/octad statements elsewhere (grep). List statements that are only about the simulator's Mapper constants separately as out_of_scope.
Assign each to a track: A-genus, B-dyons, C-lattices, D-tda, E-flux (tadpoles, fluxes), F-whichk3 (T2 self-dual point, transcendental lattices, binary forms, Kummer code, Golay/octads, M24 group structure, symplectic automorphisms, fixed points). Record a separate list 'orientifold_statements' (O7/O5/D7/D5 charges, D3 tadpoles), verbatim.
Diff against the v2 seal ${SEALED_V2}/targets_sealed_v2.json: for each statement id, mark it added, removed, changed or unchanged. Return the ids that are not unchanged in changed_since_v2_seal. Return ids only, never values: the orchestrator's prompts to blind agents must stay value-free.
Use grep -n "^theorem" and read each statement in full (they can span several lines). Do not paraphrase. Return only the summary.`,
  { label: 'targets:seal-all', phase: 'Targets', model: 'sonnet', schema: TARGETS })

phase('Blind compute')
const blindP = parallel(ACTIVE.map(t => () => agent(`${RULES}${BLIND}
${t.prompt}
Write to ${OUT}/${t.key}/: scripts, inputs.json, and results.json and exports.json generated by the scripts. Keep each script under ~20 min; lower the truncation and say so rather than hang. Commit only your own paths (rule 11). Put everything you could not do in could_not_do.`,
  { label: `blind:${t.key}`, phase: 'Blind compute', model: 'sonnet', schema: TRACK })))

const [targets, blindRaw] = await Promise.all([targetsP, blindP])
const blind = (blindRaw || []).filter(Boolean)
log(`sealed targets: ${targets ? targets.count : 0}; changed since v2 seal: ${targets && targets.changed_since_v2_seal ? targets.changed_since_v2_seal.length : '?'}; blind tracks: ${blind.length}/${ACTIVE.length}`)
if (!targets || blind.length === 0) return { stopped: 'no targets or blind results', targets, blind }

phase('Chain')
const chain = await agent(`${RULES}${BLIND}
CHAIN CHECK. Links:
- D (GUDHI) gives chi, b2 and n_singular.
- C (definition of K3 + Noether + Hodge index) gives chi_top, then the signature, then the lattice.
- A gives Z(tau,0), which must equal D's chi.
- B: the DMVV product at z=0 must give Goettsche with D's chi.
- E reads n_singular from D.
- F supplies M24 cycle shapes to A. Its T(A) and T(Km A) discriminants are consistent with its binary-form enumeration.
Track results: ${JSON.stringify(blind.map(b => ({ track: b.track, inputs: b.inputs, exports: b.exports, rigidity: b.rigidity })))}
Write ${OUT}/chain/chain_check.py. It reads every value from the exports files through paths relative to the repo root, with no literals. It must run from a CLEAN CLONE: test it with git clone ${WT} /tmp/k3t2v3_clone && cd /tmp/k3t2v3_clone && ${PY} audit/k3t2_rigidity_v3/chain/chain_check.py, then delete the clone. v2's chain failed from a clean checkout.
For each link report: CONSISTENT or not; CROSS_METHOD (two different methods) or POINTER (reads a value produced by the same method, or a declared input); and the declared inputs shared by both ends. A link whose ends share a declared normalisation is not independent corroboration. Commit only your files. Return a short plain-text summary.`,
  { label: 'chain:cross-track', phase: 'Chain', model: 'sonnet' })

phase('Compare')
const compare = await agent(`${RULES}
COMPARATOR (you may read LeanMaster read-only and ${SEALED}). Sealed file: ${targets.file} (${targets.count} targets; ids changed since the v2 seal: ${JSON.stringify(targets.changed_since_v2_seal || [])}).
Blind results: ${JSON.stringify(blind.map(b => ({ track: b.track, results: b.results.map(r => ({ id: r.id, quantity: r.quantity, computed: r.computed, shared_inputs: r.shared_inputs, script: r.script })) })))}
CHAIN: ${chain}
Make a row for EVERY sealed target. Open the blind results.json and inputs.json files, not only the summaries. Statuses: AGREE / DISAGREE / NOT_COMPUTED / NOT_COMPARABLE (the conventions differ; give the conversion and say whether the values match under it) / OUT_OF_SCOPE.
independent_routes counts only blind routes whose declared inputs are DISJOINT, with the shared inputs listed in shared_inputs. v2 overcounted 13 rows: GUDHI plus a normalisation-conditional genus is 1 route, and a typed Gram matrix is 0 routes.
Tautological routes (the same function on both sides) count as NOT_COMPUTED.
Before finishing, run git -C ${LM} log <seal commit>..HEAD -- <file> for each file behind a DISAGREE row. If the statement changed after the seal, set post_seal_changed and compare against both versions.
Investigate every DISAGREE with a third minimal computation, and say who is right. Write ${OUT}/comparison.json and commit it (it contains Lean values; that is fine now).`,
  { label: 'compare:all-targets', phase: 'Compare', schema: COMPARE })

phase('Reverse')
const reverse = await agent(`${RULES}
REVERSE PASS (you may read LeanMaster and all of ${OUT}). Take 5-8 Lean theorems that are limited to a finite range and test them beyond it with exact arithmetic.
EXCLUDE the 8 already done:
- v1: p24 k=6..8, DMZ m=4, polar m=4,5, trace bound d=7..10;
- v2: c_depends_only_on_D to q^20, moonshine twined_div24/twinedAll_div/modules_decompose to n=20, immortal m=4 Hecke form.
Candidates:
- umbral exactness at more levels;
- DMZ m=3 lines;
- Moore's binary-form count identity beyond Lean's range;
- twined Goettsche characters for more classes;
- Stream 8 statements (WhichK3, SelfDualT2, KummerE3, ForgerE4, KummerD4) beyond their checked cases, e.g. other attractive K3s or other CM points;
- tadpole statements against TT.
For each: the prediction, the computation, whether it holds or the minimal counterexample, and a decide-checkable Lean-style statement.
In the independence field, state exactly what the check shares with the Lean statement or with literature quotes. v2's item 3 overstated its independence: it matched two DMZ quotes, not a DMVV computation.
Save under ${OUT}/reverse/ and commit.
COMPARISON counts: ${JSON.stringify(compare && compare.counts)}`,
  { label: 'reverse:beyond-range', phase: 'Reverse', model: 'sonnet', schema: REVERSE })

phase('Skeptic')
const ctx = `BLIND: ${JSON.stringify(blind)}
CHAIN: ${chain}
COMPARISON: ${JSON.stringify(compare)}
REVERSE: ${JSON.stringify(reverse)}`
const skeptics = await parallel([
  () => agent(`${RULES}
SKEPTIC, lens = BLINDNESS, CIRCULARITY AND REPRODUCIBILITY. You produced none of this; default to distrust.
${ctx}
1. grep all scripts under ${OUT}/{A-genus,B-dyons,C-lattices,D-tda,E-flux,F-whichk3,chain}/, and git log -p of the v3 commits, for: LeanMaster paths, ${SEALED}, ${SEALED_V2}, the forbidden v1/v2 files, and literal answer values in computation paths. A literal is allowed only in an 'expected' field or in inputs.json as a declared input. List file:line.
2. Check rule 7: every typed normalisation is declared, and its shared_inputs propagate to every downstream value and to the comparison's independent_routes.
3. For every rigidity entry, decide its final label (rule 8), reading the code of the selecting condition and of the negative control.
4. Clone the worktree to /tmp (git clone ${WT} /tmp/k3t2v3_skb) and rerun every stated command there. Report reproduced/not and clean_clone = true, then delete the clone. Flag hand-authored JSON and absolute paths.
5. Review the REVERSE items for independence claims (reverse_verdicts).
Write ${OUT}/skeptic_blindness/verdict.json and commit.`, { label: 'skeptic:blindness', phase: 'Skeptic', schema: AUDIT }),
  () => agent(`${RULES}
SKEPTIC, lens = MATHEMATICAL CORRECTNESS. You produced none of this; default to distrust.
${ctx}
Check the mathematics, not the bookkeeping:
(a) Track A: the structural derivation of chi(g) and F_g from cycle shapes; the conventions (A_n vs 2A_n, H(0) = -1/12, the strips for A_{2,m}); the sign pattern.
(b) Track B: the DMVV product evaluated at z=0 (genuinely, not via Goettsche); the exact solve for (N, M); the Hurwitz table against the Dirichlet cross-check.
(c) Track C: the chain K = O, b1 = 0 => chi(O) = 2 => chi_top => (b2+, b2-); the realizability of every control.
(d) Track D: the open-star premise at each N; the MV computation of b3 and b4; the parity of the intersection form, if attempted.
(e) Track E: the unimodular lattice against TT (A.2)-(A.3) line by line; the index of the (A.6) span; the line-941 resolution; the symmetry group justification.
(f) Track F: the stabiliser orders and enhanced root system; the NS/T computation; the |Aut| enumeration; the Golay/M24 construction (order check); the Kummer-code embedding search; the Lefschetz numbers.
(g) Every DISAGREE or NOT_COMPARABLE row: who is right.
(h) Every REVERSE item.
Re-derive at least one number per track with an independent minimal script.
Write ${OUT}/skeptic_math/verdict.json and commit.`, { label: 'skeptic:math', phase: 'Skeptic', schema: AUDIT }),
])

phase('Synthesize')
const report = await agent(`${RULES}
Produce the v3 REPORT as your final text (the harness blocks subagents from writing .md files; the orchestrator will save it). Also write ${OUT}/agreement_ledger_v3.json with a builder script ${OUT}/build_ledger_v3.py. The builder reads only committed files inside the worktree, and it applies both skeptics' disputes by explicit, documented rules. Commit those two.
Inputs:
- TARGETS ${JSON.stringify(targets)}
- CHAIN ${chain}
- COMPARISON ${JSON.stringify(compare)}
- REVERSE ${JSON.stringify(reverse)}
- SKEPTICS ${JSON.stringify(skeptics)}
- BLIND inputs/rigidity/could_not_do ${JSON.stringify(blind.map(b => ({ track: b.track, inputs: b.inputs, rigidity: b.rigidity, could_not_do: b.could_not_do })))}
Sections:
1. Bottom line: counts after the skeptics. Compare with v2 (312 rows at LeanMaster v3.20.0: 77 clean, 15 qualified, 0 DISAGREE against the blind tracks, 3 sealed-Lean errors since fixed in v3.21.0, 207 NOT_COMPUTED, 0 rows with two independent routes; rigidity 4/14, of which both skeptics agree on 2) and with v1 (5/3/0 of 12 rowed out of 32; rigidity 3/11). Say what changed in LeanMaster since the v2 seal.
2. Per-track table, including Track F.
3. Declared inputs: the list, and which results depend on each.
4. Rigidity table with final labels. One sentence on what "zero free parameters" means here.
5. The chain, with each link marked CROSS_METHOD or POINTER.
6. Flux vacua (Track E) in the unimodular lattice, and what stays unfixed.
7. Which K3 (Track F), with every identification tiered.
8. Reverse predictions, with the skeptics' verdicts.
9. Everything not done, disputed, or from memory.
Header: "Generated by workflow k3t2-rigidity-loop-v3; computations tier B at best; physical identifications tier L/C."`,
  { label: 'synthesize:report-v3', phase: 'Synthesize' })

return { targets, chain, counts: compare && compare.counts, reverse, skeptics, report }
