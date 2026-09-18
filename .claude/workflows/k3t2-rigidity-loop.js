export const meta = {
  name: 'k3t2-rigidity-loop',
  description: 'Second, independent route to the zero-parameter mathematics of K3 x T2: blind Python/GUDHI computations, deformation (rigidity) tests, then comparison with LeanMaster kernel-checked theorems, and a reverse pass from Lean statements to new computed predictions',
  whenToUse: 'When LeanMaster advances: re-derive its K3 x T2 numbers independently (no access to the Lean code), show no free parameter can be inserted, and report agreement or discrepancy.',
  phases: [
    { title: 'Targets', detail: 'sealed list of LeanMaster statements (comparator only)' },
    { title: 'Blind compute', detail: '4 independent tracks, forbidden to open LeanMaster' },
    { title: 'Compare', detail: 'match blind numbers against Lean statements' },
    { title: 'Reverse', detail: 'Lean statements -> new predictions -> compute beyond Lean range' },
    { title: 'Skeptic', detail: 'blindness audit, reruns, hard-coding hunt' },
    { title: 'Synthesize', detail: 'agreement ledger + report' },
  ],
}

const WT = '/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2'
const PY = '/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python'
const LM = '/home/callensxavier_gmail_com/SocrateAI-Scientific-Agora-LeanMaster'
const OUT = `${WT}/audit/k3t2_rigidity`

const RULES = `
GROUND RULES (strict):
- Work ONLY in the git worktree ${WT} (branch loop/k3t2-rigidity), outputs under ${OUT}/ . Never touch proofs/, main, or any other worktree; never push. 
- Python: ${PY} (gudhi 3.13, sympy, mpmath, numpy, scipy, networkx). Use EXACT arithmetic (python int / fractions.Fraction / sympy Rational) for every integer or rational claim; floats only for eigenvalue signs and the trace bound.
- Evidence-bound: every reported number comes from a script you wrote and ran; save script + JSON output; give the command. Never type a target number into a script as a constant to compare with itself: numbers must be COMPUTED from the defining formula. Known literature values may appear only in a separate "expected" field that you fill AFTER the computation, stating your source (formula or paper eq.), and mark "from memory, unverified" when that is the case.
- Tiers: B = exact arithmetic with a negative control; L = literature identification (that the number means dyons, moonshine, etc.); C = conjecture. Your computations are tier B at best. Never write "proved".
- Commit on loop/k3t2-rigidity; message ends with: Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
`
const BLIND = `
BLINDNESS (this is the point of the exercise): you must NOT open, grep, or list anything under ${LM} or any path containing "LeanMaster", nor lean_foundation/, proofs/, audit/PAPER_FACTS.md, audit/*.md of this repo, nor workshopcosmo.py / axioms.json / simulation_results.json (they contain the answers). Derive everything from the mathematical definitions given below. A later agent audits your shell history and scripts for violations.
RIGIDITY TEST (the "zero free parameter" experiment): for each structure, insert a would-be free parameter (examples are given per track), scan it over a range including the true value, and record for which values ALL consistency conditions hold. Report the solution set. "Zero free parameters" is supported only if the solution set is a single point and at least one nearby value demonstrably fails (negative control).`

const TRACK = {
  type: 'object',
  properties: {
    track: { type: 'string' },
    results: { type: 'array', items: { type: 'object', properties: {
      id: { type: 'string' }, quantity: { type: 'string' }, computed: { type: 'string' }, method: { type: 'string' },
      script: { type: 'string' }, command: { type: 'string' }, exact: { type: 'boolean' } },
      required: ['id', 'quantity', 'computed', 'method', 'script'] } },
    rigidity: { type: 'array', items: { type: 'object', properties: {
      parameter_inserted: { type: 'string' }, scanned: { type: 'string' }, conditions: { type: 'string' },
      solution_set: { type: 'string' }, negative_control: { type: 'string' } },
      required: ['parameter_inserted', 'scanned', 'solution_set', 'negative_control'] } },
    could_not_do: { type: 'array', items: { type: 'string' } },
    commit: { type: 'string' },
  },
  required: ['track', 'results', 'rigidity', 'could_not_do'],
}
const TARGETS = {
  type: 'object',
  properties: { tag: { type: 'string' }, targets: { type: 'array', items: { type: 'object', properties: {
    id: { type: 'string' }, track: { type: 'string' }, theorem: { type: 'string' }, file_line: { type: 'string' },
    statement_verbatim: { type: 'string' }, numbers: { type: 'string' }, order_or_range: { type: 'string' }, committed: { type: 'boolean' } },
    required: ['id', 'track', 'theorem', 'file_line', 'statement_verbatim', 'numbers'] } } },
  required: ['tag', 'targets'],
}
const COMPARE = {
  type: 'object',
  properties: { rows: { type: 'array', items: { type: 'object', properties: {
    target_id: { type: 'string' }, lean_value: { type: 'string' }, blind_value: { type: 'string' },
    status: { type: 'string', description: 'AGREE | DISAGREE | NOT_COMPUTED | NOT_COMPARABLE' }, note: { type: 'string' } },
    required: ['target_id', 'status'] } },
    agree: { type: 'number' }, disagree: { type: 'number' }, not_computed: { type: 'number' },
    blind_only_findings: { type: 'array', items: { type: 'string' }, description: 'things the blind route computed that Lean does not state (candidates for new theorems)' },
    lean_only: { type: 'array', items: { type: 'string' } } },
  required: ['rows', 'agree', 'disagree', 'not_computed'],
}
const REVERSE = {
  type: 'object',
  properties: { predictions: { type: 'array', items: { type: 'object', properties: {
    from_theorem: { type: 'string' }, prediction: { type: 'string' }, computed: { type: 'string' }, holds: { type: 'boolean' },
    proposed_lean_statement: { type: 'string' }, script: { type: 'string' } },
    required: ['from_theorem', 'prediction', 'computed', 'holds'] } }, commit: { type: 'string' } },
  required: ['predictions'],
}
const AUDIT = {
  type: 'object',
  properties: {
    blindness_violations: { type: 'array', items: { type: 'string' } },
    hardcoded_answers: { type: 'array', items: { type: 'string' } },
    reruns: { type: 'array', items: { type: 'object', properties: { script: { type: 'string' }, reproduced: { type: 'boolean' }, note: { type: 'string' } }, required: ['script', 'reproduced'] } },
    disputed_rows: { type: 'array', items: { type: 'string' } },
    rigidity_verdicts: { type: 'array', items: { type: 'object', properties: { parameter: { type: 'string' }, genuinely_rigid: { type: 'boolean' }, reason: { type: 'string' } }, required: ['parameter', 'genuinely_rigid', 'reason'] } },
    tda_valid: { type: 'boolean' }, tda_note: { type: 'string' } },
  required: ['blindness_violations', 'hardcoded_answers', 'reruns', 'rigidity_verdicts', 'tda_valid'],
}

const TRACKS = [
  { key: 'genus-moonshine', prompt: `TRACK A - K3 elliptic genus and Mathieu moonshine.
Definitions: Jacobi thetas theta_1..4(tau,z); phi_{0,1} = 4 * sum_{i=2,3,4} (theta_i(tau,z)/theta_i(tau,0))^2 ; Z_K3 = 2*phi_{0,1}; phi_{-2,1} = theta_1(tau,z)^2 / eta(tau)^6. Work with truncated q-series whose coefficients are Laurent polynomials in y, exact integers, through at least q^8.
Compute: (1) Z_K3(tau, z=0) (must be constant = Euler number); the coefficients c(n,l) and check they depend only on 4n - l^2; list c(D) for D = -1,0,3,4,7,8,... (2) decompose Z_K3 into N=4 c=6 characters (massless h=1/4,l=0 ; massive h=1/4+n, l=1/2), or equivalently extract H(tau) = 2 q^{-1/8}(-1 + sum A_n q^n) via Z_K3 = 24*mu-term + H*theta_1^2/eta^3 (Appell-Lerch mu) ; report A_1..A_8 computed, and try to decompose A_1..A_5 into dimensions of M24 irreducibles (dims 1,23,45,45,231,231,252,253,483,770,770,990,990,1035,1035,1035',1265,1771,2024,2277,3312,3520,5313,5544,5796,10395 - this list is a given). (3) the twined series for class 2A: H_2A = (2 eta-quotient formula: H_g = (chi(g)/24) H - F_g/eta^3 with F_2A = 16 Lambda_2, Lambda_N = N q d/dq log(eta(N tau)/eta(tau))) ; report its first 6 coefficients, and test whether the identity "A_2 * 60 = 4 * A_1 * 77" that holds untwined (check it) survives when A_n is replaced by the 2A-twined coefficients.
Rigidity parameters to insert: (a) replace the 24 multiplying the mu-term by N in 20..28 and ask for which N the resulting H has integer coefficients AND polar term -2; (b) replace the overall 2 in Z_K3 = 2 phi_{0,1} by a rational k and require Z(tau,0) = chi of a compact hyperkahler 4-fold (24) with integrality.` },
  { key: 'dyons', prompt: `TRACK B - quarter-BPS dyons on K3 x T2 (Dijkgraaf-Moore-Verlinde-Verlinde product; Dabholkar-Murthy-Zagier).
Definitions: let c(4n-l^2) be the Fourier coefficients of 2*phi_{0,1}(tau,z) (compute them yourself from thetas: phi_{0,1} = 4 sum_{i=2,3,4} (theta_i(tau,z)/theta_i(tau,0))^2). DMVV: sum_k G_k p^k = prod_{r>=1, s>=0, t in Z} (1 - p^r q^s y^t)^{-c(4rs - t^2)} is the generating function of elliptic genera of Sym^k(K3). 
Compute exactly: (1) G_k at z=0 (y=1), q^0, for k = 0..5 : these must equal the Euler numbers of Hilb^k(K3); compare AFTER computing with the coefficients of prod (1-p^n)^{-24}. (2) With A = phi_{-2,1} = theta_1^2/eta^6, B = phi_{0,1}, E4, E6 from divisor sums, Delta = eta^24: compute psi_m via Delta*psi_m = G_{m+1}/A for m = -1,0,1 and check whether Delta*psi_0 = 2B/A and 4*Delta*psi_1 = 9 B^2/A + 3 E4 A hold to the order you can reach. (3) m=1 finite part: remove the polar part p24(2) * A_{2,1}, A_{2,m} = sum_s q^{m s^2 + s} y^{2ms+1}/(1 - q^s y)^2 expanded in |q|<|y|<1, and test whether the remainder equals 3 E4 A - 648 * Hhat where Hhat = sum H(4n - l^2) q^n y^l with H the Hurwitz class numbers, which you must compute independently by counting reduced binary quadratic forms (H(0) = -1/12). Report H(D) for D = 3,4,7,8,11,12,15.
Rigidity parameters: (a) exponent multiplier: use -kappa*c(...) in the product and find which kappa gives integer Euler numbers matching prod(1-p^n)^{-24}; (b) polar multiplicity: remove N * A_{2,1} with N in 320..328 and record for which N the remainder has no pole at z=0 (y=1); (c) replace 648 by M in a window and record for which M the identity holds.` },
  { key: 'lattices-duality', prompt: `TRACK C - lattices and T-duality.
Compute from explicit integer Gram matrices (build them yourself): (1) E8 Cartan matrix: determinant, even, positive definite (exact: leading principal minors via sympy), number of roots (norm-2 vectors, by enumeration) ; (2) U = [[0,1],[1,0]]; K3 lattice = 3U + 2(-E8): rank, determinant, signature COMPUTED from the 22x22 Gram matrix (exact inertia via sympy LDL or Sturm on the characteristic polynomial, not floats), evenness, unimodularity; Mukai lattice = K3 + U (4,20); K3 x T2 charge lattice = Mukai + U + U... report signature (p,q) you obtain for Gamma^{6,22} = 6U + 2(-E8). (3) O(d,d;Z): with eta = [[0,1],[1,0]] (blocks), verify for d=2,3 on 200 seeded random integer samples each that theta-shifts [[1,Theta],[0,1]] with antisymmetric integer Theta, basis changes [[A,0],[0,A^{-T}]] with A in GL(d,Z), and factorized dualities preserve eta; that factorized duality squares to 1; and that eta*H(G,0)*eta = H(G^{-1},0) for the generalized metric H(G,B). Negative control: symmetric Theta must FAIL. (4) dual-scale bound: for random SPD G (d = 1..6, 2000 seeded samples) min of tr G + tr G^{-1} - 2d; show equality iff G = 1 (perturbation scan), exact proof sketch via eigenvalues x + 1/x >= 2. (5) flux tadpole arithmetic: chi(K3)=24 from Betti numbers 1,0,22,0,1; chi(K3 x K3)/24; chi(K3 x T2) via Kunneth.
Rigidity parameters: (a) number n of E8 copies and m of U copies: for which (m,n) with rank 22 is the lattice even unimodular of signature (3,19)? (enumerate; recall even unimodular requires p - q = 0 mod 8 - CHECK it computationally, do not assume); (b) T-duality radius map R -> lambda * alpha'/R: for which lambda is it an involution preserving the spectrum m^2 = (n/R)^2 + (w R/alpha')^2 ? (scan lambda, exact).` },
  { key: 'tda-gudhi', prompt: `TRACK D - topology of K3 x T2 by ACTUAL topological data analysis with INRIA GUDHI (import gudhi; SimplexTree / PeriodicCubicalComplex / RipsComplex; persistent homology with coefficient field Z/3 or Z/5 to avoid 2-torsion; say which).
Goal: obtain the Betti numbers of the Kummer construction K3 = resolution of T^4/Z_2 and of K3 x T^2 from computation, not from memory.
(1) Known-answer controls: T^2 -> (1,2,1); T^4 -> (1,4,6,4,1) via gudhi PeriodicCubicalComplex or a simplicial torus.
(2) T^4/Z_2 (x -> -x): build a Z_2-invariant triangulation of T^4 (e.g. Freudenthal/Kuhn triangulation of the periodic cubical grid Z_N^4 with N = 4 or 6 so that the 16 fixed points are vertices and the action is simplicial and regular - subdivide if needed so the quotient is a simplicial complex), form the quotient complex by identifying simplices with their images, insert into a gudhi SimplexTree and compute Betti numbers over Z/3. Expected rational answer is for you to discover; report it together with the Euler characteristic computed from simplex counts (independent check) and the count of fixed vertices.
(3) Resolution: each of the 16 singular points (link RP^3) is replaced by a (-2)-curve neighbourhood (T^*S^2). Either (i) do it combinatorially (remove open star of each fixed vertex, glue a triangulated disk bundle over S^2 with boundary RP^3) and recompute with gudhi, or, if gluing is out of reach, (ii) compute with gudhi the homology of the complement and of the link (must be RP^3: Betti (1,0,0,1) over Z/3, and over Z/2 (1,1,1,1) as torsion control) and finish with an explicit Mayer-Vietoris calculation, clearly labelled as hybrid. Report b_0..b_4, chi, and hence the signature ingredients b_2 = b_2^+ + b_2^-.
(4) K3 x T^2 by Kunneth from your computed numbers: Betti numbers and chi.
(5) A sampled-point-cloud persistent-homology run as a SECOND TDA route on something feasible: e.g. 1500 seeded points on a flat T^2 embedded in R^4 (Clifford torus) with RipsComplex / alpha complex -> two long H1 bars, one H2 bar; and a null (uniform ball) control.
Rigidity parameter: the order of the orbifold group. Repeat (2) for T^4/Z_3 is hard; instead vary the NUMBER of resolved points k = 0..16 in the Mayer-Vietoris/Euler count and record for which k the result is a smooth closed manifold with chi = 24 and vanishing first Chern class indicator you can actually compute (Noether: chi = 24 <=> 12*(chi(O)) with c1 = 0; state what is computed vs quoted).` },
]

phase('Targets')
const targetsP = agent(`${RULES}
TASK (you are the only agent allowed to read LeanMaster, READ-ONLY): build the sealed target list ${OUT}/targets_sealed.json (do not commit it until told; just write it) from ${LM}. Run git -C ${LM} describe --tags and git -C ${LM} status --short (mark targets living in uncommitted files as committed=false).
For each of these areas quote 4-10 theorems VERBATIM with file:line and extract the concrete numbers they assert and the order/range they cover:
 track genus-moonshine: DualScaleMoonshine (elliptic genus Z(tau,0)=24, index-1 property, EOT coefficients / A_n, twining at 2A, ForgerTest lock_at_identity and ratio_fails_*, Shadow polar multiplicity 24 vs 23/25);
 track dyons: DualScaleDyons (goettsche / p24, c_first, DMZ (5.16) identities, polar_part_removes_pole, hurwitz_values, immortal_m1, immortal_m1_needs_H, immortal_m2/m3, twined_goettsche_*);
 track lattices-duality: DualScaleStream2 + StringTheoryFormalization (cartanE8_posDef, sigK3_eq, sigMukai_eq, sigK3T2_eq, thetaShift_isODD, basisChange_isODD, factorized_mul_self, tduality_inverts_metric, dual-scale trace bound, dualScale_eq_iff, tadpole_conservation, euler numbers);
 track tda-gudhi: any statement about Betti numbers / Euler characteristic / Kummer (b2 = 22, chi = 24, 16 fixed points) wherever it lives.
Use #check-free reading: grep -n "^theorem" then sed the statement. Do not build Lean. Do not paraphrase statements.`,
  { label: 'targets:sealed-from-LeanMaster', phase: 'Targets', model: 'haiku', schema: TARGETS })

phase('Blind compute')
const blindP = parallel(TRACKS.map(t => () => agent(`${RULES}${BLIND}
${t.prompt}
Write scripts to ${OUT}/${t.key}/ , run them (keep each run under ~10 minutes; lower the truncation order rather than hang, and say so), save results.json, commit. List honestly in could_not_do anything you did not manage; a smaller exact result beats a larger uncertain one.`,
  { label: `blind:${t.key}`, phase: 'Blind compute', model: 'sonnet', schema: TRACK })))

const [targets, blindRaw] = await Promise.all([targetsP, blindP])
const blind = (blindRaw || []).filter(Boolean)
log(`targets: ${targets ? targets.targets.length : 0}; blind tracks returned: ${blind.length}/4`)
if (!targets || blind.length === 0) return { stopped: 'no targets or no blind results', targets, blind }

phase('Compare')
const compare = await agent(`${RULES}
You are the COMPARATOR. You computed nothing. Inputs:
SEALED LEAN TARGETS: ${JSON.stringify(targets)}
BLIND RESULTS: ${JSON.stringify(blind)}
For every target decide AGREE / DISAGREE / NOT_COMPUTED / NOT_COMPARABLE by opening the blind results.json files under ${OUT}/ (not just the summaries) and the Lean statement in ${LM} (read-only). Be strict about conventions (A_n vs 2A_n, H(0) = -1/12, which strip for the polar part, signature sign conventions): a convention difference is NOT_COMPARABLE until you state the exact conversion and both sides match under it. A DISAGREE is the most valuable outcome: investigate which side is wrong by a third, minimal computation of your own and say so. Also list what the blind route found that Lean does not state (e.g. signatures computed from Gram matrices where Lean's are bookkeeping) and what Lean has that nobody recomputed. Write ${OUT}/comparison.json and commit it together with targets_sealed.json.`,
  { label: 'compare:lean-vs-blind', phase: 'Compare', schema: COMPARE })

phase('Reverse')
const reverse = await agent(`${RULES}
REVERSE PASS: Lean statements -> predictions -> experiment. You MAY read LeanMaster (read-only) and the blind scripts under ${OUT}/.
LEAN TARGETS: ${JSON.stringify(targets && targets.targets.map(t => ({ id: t.id, theorem: t.theorem, numbers: t.numbers, range: t.order_or_range })))}
COMPARISON: ${JSON.stringify(compare)}
Pick 5-8 Lean theorems whose statement is limited to a finite order or a finite list (e.g. "through q^9", m = 1,2,3, k <= 4, 26 classes at levels 1-7) and treat each as a HYPOTHESIS making predictions beyond its checked range: next q-orders, m = 4 and m = 5 immortal dyons (H|V_m structure for prime m), Goettsche numbers k = 5..8 as twined characters at class 2A (non-negative integer multiplicities need the character table - only if you can build it; otherwise restrict to integrality and the 1A value), signatures of the lattices for which Lean only does bookkeeping, the trace bound in d = 7..10. Extend or reuse the blind scripts (exact arithmetic) to TEST each prediction. For each that holds, write the precise Lean-style statement that LeanMaster could add (same style as the quoted theorems, decide-checkable). For any that FAILS, report it prominently with the minimal counterexample. Save under ${OUT}/reverse/ and commit.`,
  { label: 'reverse:lean-to-predictions', phase: 'Reverse', model: 'sonnet', schema: REVERSE })

phase('Skeptic')
const audit = await agent(`${RULES}
You are the SKEPTIC / auditor; you produced none of this. Default to distrust.
BLIND RESULTS: ${JSON.stringify(blind)}
COMPARISON: ${JSON.stringify(compare)}
REVERSE: ${JSON.stringify(reverse)}
1. BLINDNESS AUDIT: grep every script under ${OUT}/ (git log -p on the branch too) for "LeanMaster", "lean_foundation", "PAPER_FACTS", "axioms.json"; then hunt for hard-coded answers: target numbers (45, 231, 770, 2277, 324, 3200, 25650, 176256, -128, 216, 648, 27720, 22, (3,19) ...) appearing as literals that feed a comparison rather than being computed. A literal inside an "expected" field filled after computation is fine; a literal inside the computation path is a violation. List file:line.
2. RERUN at least one script per track plus the reverse scripts; report reproduced or not (byte-compare results where deterministic).
3. RIGIDITY: for each inserted parameter judge whether the scan genuinely shows a one-point solution set (conditions independent of the parameter's definition, a nearby value really fails) or is circular (the condition was built from the true value). 
4. TDA: confirm gudhi is imported and used for the Betti numbers claimed, the T^2/T^4 controls pass, the quotient complex is a valid simplicial complex (check: every simplex's faces present, vertex count vs orbit count, Euler characteristic from f-vector equals alternating Betti sum), and that any Mayer-Vietoris step is labelled hybrid.
5. Dispute any comparison row you think is wrongly marked AGREE.`,
  { label: 'skeptic:audit', phase: 'Skeptic', schema: AUDIT })

phase('Synthesize')
const report = await agent(`${RULES}
Write ${OUT}/REPORT.md (plain language, short sentences, for a physicist in a hurry) and ${OUT}/agreement_ledger.json, then commit. Inputs:
TARGETS tag: ${targets.tag}; COMPARISON: ${JSON.stringify(compare)}
BLIND: ${JSON.stringify(blind.map(b => ({ track: b.track, rigidity: b.rigidity, could_not_do: b.could_not_do, n_results: b.results.length })))}
REVERSE: ${JSON.stringify(reverse)}
AUDIT: ${JSON.stringify(audit)}
Sections: 1. Bottom line: how many Lean statements were independently reproduced (AGREE / DISAGREE / NOT_COMPUTED), after removing rows the skeptic disputed or that came from scripts with blindness/hard-coding violations. 2. Table per track: quantity, Lean theorem, Lean value, blind value, status. 3. Rigidity: for each inserted parameter, the solution set and the skeptic's verdict; state in one sentence what "zero free parameters" does and does not mean here (it is a statement about the K3 x T2 mathematics - integers fixed by topology and modularity - NOT about the cosmological sectors, which LeanMaster's own Stream 3/5 notes keep at tier C). 4. What TDA (GUDHI) contributed, with controls. 5. Reverse pass: new predictions that held (with proposed Lean statements, ready to hand to the LeanMaster session) and any that failed. 6. Everything not done, disputed, or from memory, verbatim. Header: "Generated by workflow k3t2-rigidity-loop; computations are tier B at best; physical identifications tier L/C." Return the report markdown as your final text.`,
  { label: 'synthesize:report', phase: 'Synthesize' })

return { tag: targets.tag, agree: compare && compare.agree, disagree: compare && compare.disagree, not_computed: compare && compare.not_computed, audit, reverse, report }
