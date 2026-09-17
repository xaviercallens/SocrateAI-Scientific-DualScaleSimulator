export const meta = {
  name: 'physics-claims-audit',
  description: 'Haiku fleet: extract every Lean/physics/simulation/manuscript claim, assign Mathesis tier (A/B/L/C/X), verify with 3 evidence-bound skeptics, map local Lean files to LeanMaster replacements, write audit ledgers',
  whenToUse: 'After editing proofs/, rust_simulator/, workshopcosmo.py or the manuscript; before submission; when LeanMaster (the string-theory formalization) advances. Caller writes returned `ledger` → audit/physics_claims_ledger.md and `replacementMap` → audit/lean_replacement_map.md. Args: {judgeModel?: "haiku"|"sonnet", only?: string[] group ids, skipReplacementMap?: boolean}',
  phases: [
    { title: 'Extract', detail: 'one Haiku reader per Lean group / simulator / manuscript slice', model: 'haiku' },
    { title: 'Verify', detail: '3 Haiku skeptics per finding: quote fidelity, domain (decisive-experiment), steelman', model: 'haiku' },
    { title: 'Replacement', detail: 'per local Lean file: read candidate LeanMaster module (read-only), judge statement adequacy', model: 'haiku' },
  ],
}

const REPO = '/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator'
const LM = '/home/callensxavier_gmail_com/SocrateAI-Scientific-Agora-LeanMaster'
const PAPER = `${REPO}/papers/T-dulaity alone/T_duality_Alone.tex`
const A = args || {}
const JUDGE = A.judgeModel || 'haiku'

const GROUNDING = `
GROUNDING RULES (mandatory — violating findings are discarded automatically):
- Read files with \`cat -n\` / \`sed -n 'A,Bp'\` BEFORE reporting. Never report on a file you did not open in this session. A previous Haiku audit of this repo invented findings with zero tool calls — do not repeat that.
- Every finding needs the 1-indexed line and a VERBATIM one-line quote copied from your output (no line-number prefix).
- No numbers or citations from memory: compute (python3 -c) or quote a file.
- Read the Lean STATEMENT, not the docstring that summarises it. Record SOUND claims too.
- LLMs are not arbiters of mathematical truth: for Lean, only the kernel/#print axioms output counts as proof status.
Repo root: ${REPO}. The sibling repo ${LM} is being edited by another session: READ-ONLY, never write there.`

const TIERS = `
Mathesis epistemic tiers (ordering X < C < L < B < A; a claim may not sit above what it rests on):
- A: kernel-verified, 0 sorry, \`#print axioms\` shows only propext / Classical.choice / Quot.sound, AND the statement adequately expresses the claim (not vacuous).
- B: exact Q/Z arithmetic or exact computation with a negative control that fails.
- L: literature result, quoted with a source.
- C: conjecture / physical interpretation.
- X: exploratory (floats, sampling, unvalidated simulation, LLM output).
Forbidden phrasing for anything below A: "formally verified", "certified", "proves", "kernel-certified", "zero axioms", "100% rigor", "first fully certified", "formalizes string theory". Tier C sentences may not use predicts/establishes/shows/implies/governs/determines/demonstrates unless prefixed "we conjecture" / "if".
Benchmarks/speedups without an execution log in the repo must read "pending hardware verification". A fitted parameter is not a prediction.`

const CATEGORIES = `
Categories (exactly one per finding):
- vacuous_by_definition: true only because a def was chosen to force it (function hard-coded to 0, identical if/else branches, field := true checked by rfl, entropy := 0.0, c(N) chosen monotone).
- axiom_assumes_conclusion: content supplied by an \`axiom\` or structure/class field stating the conclusion.
- inconsistent_axiom: axiom quantified over all types (or trivial typeclass) that is false for a concrete instance, so it proves False (e.g. 1/(1/x)=x at Nat, x=2).
- mislabelled_theorem: name/docstring claims more than the statement proves.
- physics_error: wrong vs standard references (O-plane dimension for the orbifold, dilaton shift factor 1/2, charge units, summing different-degree cohomology into one integer, …).
- numerical_bug: dead branch, clamp disabling a check, hard-coded root disagreeing with the potential, floors/slack factors hiding failure, integrator mislabelled.
- unsupported_physical_claim: physical consequence asserted with no derivation/code.
- unreproducible_number: quantitative claim with no generating code/log in the repo.
- inconsistent_numbers: same quantity differs across files (Rule 8 cross-consistency).
- overclaimed_verification: tier-A wording for something below tier A.
- sound: correct and honestly stated.`

const GROUPS = [
  { id: 'lean-duality', files: ['proofs/BuscherRules.lean', 'proofs/LeanscratchDB/DoubleScaleT2.lean', 'proofs/GysinBEMSequence.lean', 'proofs/LeanscratchDB/TopologicalTDuality.lean', 'proofs/FourierMukai.lean'],
    focus: 'Buscher conventions (sign of g_{mu y}; dilaton Phi - (1/2) log g_yy) and consistency between the two Buscher files; global axioms and whether they are inconsistent; Gysin "exactness" fields (im = ker vs composition-zero); theorem names vs statements (euler_class_pullback_zero).' },
  { id: 'lean-k3-tadpole', files: ['proofs/TadpoleCancellation.lean', 'proofs/KummerTDAAnomalyCertification.lean', 'proofs/KummerOrbifoldResolution.lean', 'proofs/MukaiLatticeK3.lean'],
    focus: 'O-planes from reflecting 4 torus coordinates (Op, p = 9 - #reflected); Gimon-Polchinski vs Sen T2/Z2 limit (4 O7- + 16 D7) charge units; D3 tadpole chi/24 context and sign; Green-Schwarz 6D condition (n_H - n_V + 29 n_T = 273) vs struct hard-coded to 0; classTadpoleAnomaly returning 0 for every class; Kummer b2/Euler/signature arithmetic; Mukai pairing/vector conventions (likely sound).' },
  { id: 'lean-cosmo', files: ['proofs/SwamplandDistanceConjecture.lean', 'proofs/FluxVacuumDecayCTheorem.lean', 'proofs/TachyonCondensationKTheory.lean', 'proofs/DbraneBoundaryStates.lean', 'proofs/FTheoryCosmology.lean', 'lakefile.lean'],
    focus: 'discreteMassBound defined as 0 for d>0; dualTowerMass identical branches; c(N)=100(N+1) and bounceAction=42(n+1) hand-picked; is_minimum_action := true; rrCharge = rank + c1 + c2; Callan-Harvey terms defined as -x and +x; FTheoryCosmology.lean LaTeX tokens (\\forall) and absence from lakefile targets.' },
  { id: 'lean-moonshine-scratch', files: ['proofs/MathieuVertexOperators.lean', 'proofs/LeanscratchDB/HoloAlg.lean', 'proofs/LeanscratchDB/QuantumEntanglement.lean', 'proofs/LeanscratchDB/EmergentCosmology.lean', 'proofs/ExportAxioms.lean', 'axioms.json'],
    focus: 'EOT multiplicities 90=45+45bar, 462=231+231bar, 1540=770+770bar (compare with the table in ' + LM + '/DualScaleStream2/Moonshine/EOT.lean, read-only); 3-point exponent algebra; is R_NL = 462/(4*90) = 77/60 derived as a bispectrum/f_NL quantity anywhere, or numerology (tier B arithmetic + tier C interpretation)?; numSupercharges = numFixedPointsT2Z2 as a "congruence" theorem; HoloAlg universe-polymorphic axioms (atiyah_singer_k3 over any R with empty FieldIndex class); placeholders := 0.0; axioms.json named "verified axioms" though it holds constants; "sorry" only inside comments.' },
  { id: 'sim-rust', files: ['rust_simulator/src/vacuum_decay_cdl.rs', 'rust_simulator/src/swampland_geodesic.rs', 'rust_simulator/src/tachyon_condensation.rs', 'rust_simulator/src/kummer_langevin.rs', 'rust_simulator/src/main.rs'],
    focus: 'CDL: gravity included or flat-space Coleman bounce?; compute roots of V\'(phi)=m2 phi - 3 kappa phi^2 + 4 lambda phi^3 with the default config and compare with hard-coded phi_true=2.05; can the overshoot branch (final_phi < phi_false) ever fire given the clamp?; .max(1.0) action floor. Swampland: is the update really Velocity-Verlet/symplectic; mass formula completeness (T-modulus winding); metric normalisation vs hard-coded alpha=1/sqrt2; 1.5x slack in bound_satisfied. Tachyon/Langevin: what k_theory_charge_conserved is computed from; RNG seeding. Count #[test].' },
  { id: 'sim-python', files: ['workshopcosmo.py', 'scripts/frontier_loops_simulation.py', 'scripts/tda_mapper.py', 'scripts/exact_math_middleware.py', 'tests/test_workshopcosmo.py'],
    focus: 'How zero_sorry / theorems_verified_count / lean_certification_verified are computed (compile? grep? hard-coded?) — a checker that never fails on purpose is untested; absolute paths /home/xavkal; rusty-SUNDIALS dependency on an absent sibling repo yet "sundials_integrated" reported; whether tests assert physics or only status strings; Mapper beta_1 (graph cycle rank vs persistent homology); symmetron/Cassini numbers; w0/wa, Delta chi2, Bayes-factor provenance (distinguishable from LambdaCDM?).' },
  { id: 'paper-front', files: ['papers/T-dulaity alone/T_duality_Alone.tex'], range: [1, 310],
    focus: 'Abstract + Buscher/Mukai/Kummer sections: map each "formally verify/certify/prove" sentence to the Lean theorem it relies on and assign the justified tier; Mapper numbers (187/557/376) vs simulation_results.json and tda_mapper_skeleton.json; Symmetron/Cassini numbers vs simulation_results.json.' },
  { id: 'paper-dynamics', files: ['papers/T-dulaity alone/T_duality_Alone.tex'], range: [311, 873],
    focus: 'LeanFlow architecture, swampland/tachyon/tunneling, flops/index theorem, topological T-duality, moduli cosmology: compare every number with *_summary.json and frontier_loops_summary.json (bounce action / bubble radius differ between vacuum_decay_cdl_summary.json and frontier_loops_summary.json); flag claims resting on vacuous/axiom-backed theorems.' },
  { id: 'paper-results', files: ['papers/T-dulaity alone/T_duality_Alone.tex', 'papers/T-dulaity alone/submission/COVER_LETTER_CPC.tex', '.zenodo.json', 'deployment/huggingface/README.md', 'deployment/huggingface/space/README.md'], range: [874, 1278],
    focus: 'Benchmarks (1,520x vs SciPy, 7.1x vs C++ CVODE, $ costs, 99.1%): grep -rn the repo for any harness or timing log; GPU/TPU claims vs code; observational claims (NANOGrav, w0/wa, Delta chi2, Bayes factor); falsifiability (does any prediction state OBSERVABLE / PRECISION / FALSIFY IF?); conclusion wording vs forbidden phrases. ALSO read fully and audit the PUBLIC metadata: .zenodo.json (description field), deployment/huggingface/README.md and deployment/huggingface/space/README.md (report their true file lines).' },
]

const FINDINGS = {
  type: 'object',
  required: ['findings', 'files_read'],
  properties: {
    files_read: { type: 'array', items: { type: 'string' } },
    findings: {
      type: 'array',
      items: {
        type: 'object',
        required: ['file', 'line', 'quote', 'claim', 'category', 'severity', 'claimed_status', 'justified_tier', 'explanation', 'suggested_fix'],
        properties: {
          file: { type: 'string', description: 'repo-relative path' },
          line: { type: 'integer' },
          quote: { type: 'string' },
          claim: { type: 'string' },
          category: { type: 'string', enum: ['vacuous_by_definition', 'axiom_assumes_conclusion', 'inconsistent_axiom', 'mislabelled_theorem', 'physics_error', 'numerical_bug', 'unsupported_physical_claim', 'unreproducible_number', 'inconsistent_numbers', 'overclaimed_verification', 'sound'] },
          severity: { type: 'string', enum: ['critical', 'high', 'medium', 'low'] },
          claimed_status: { type: 'string', description: 'how the repo/paper presents it, e.g. "formally verified", "certified", "PASS", "prediction"' },
          justified_tier: { type: 'string', enum: ['A', 'B', 'L', 'C', 'X'] },
          explanation: { type: 'string' },
          suggested_fix: { type: 'string', description: 'minimal honest change: corrected statement, code edit, or tier-appropriate rewording' },
        },
      },
    },
  },
}

const VERDICT = {
  type: 'object',
  required: ['refuted', 'confidence', 'evidence'],
  properties: {
    refuted: { type: 'boolean' },
    confidence: { type: 'string', enum: ['low', 'medium', 'high'] },
    corrected_tier: { type: 'string', enum: ['A', 'B', 'L', 'C', 'X', 'unchanged'] },
    evidence: { type: 'string', description: 'commands run and output, or the reference/convention relied on' },
  },
}

function extractPrompt(g) {
  const read = g.range
    ? `Read the manuscript slice: sed -n '${g.range[0]},${g.range[1]}p' "${PAPER}" | cat -n  — displayed number + ${g.range[0] - 1} = true file line; report TRUE lines. Read any JSON/summary files needed to cross-check numbers.`
    : `Read each file fully: ${g.files.map(f => `${REPO}/${f}`).join(', ')}.`
  return `Strict referee audit of a string-cosmology research repo (Lean 4 + Rust/Python simulations + manuscript).
${read}
Focus checklist: ${g.focus}
${TIERS}
${CATEGORIES}
${GROUNDING}
Kernel facts you may re-check (the build is green): cd ${REPO} && write a scratch file under /tmp importing the module and run \`timeout 600 lake env lean -DmaxHeartbeats=1000000 /tmp/x.lean\` with \`#print axioms <theorem>\`. Do not start \`lake build\` (another agent may be compiling; this VM tolerates 2-3 compile jobs).
Return every distinct finding, defects AND sound claims.`
}

const resolveFile = f => `${REPO}/${f.file.replace(/^\/?/, '').replace(REPO + '/', '')}`
const LENSES = [
  { key: 'fidelity', model: 'haiku', effort: 'low', prompt: f => `QUOTE FIDELITY check only. Run: sed -n '${Math.max(1, f.line - 3)},${f.line + 3}p' "${resolveFile(f)}"
Does this text appear (whitespace-insensitive) within ±3 lines: ${JSON.stringify(f.quote)} ?
Does the source actually say what is claimed: ${JSON.stringify(f.claim)} ?
refuted=true if the quote is absent or the finding misreads the source.` },
  { key: 'domain', model: JUDGE, prompt: f => `Decisive-experiment referee: try to REFUTE this finding — check it as hard as you would a positive claim.
[${f.category}, ${f.severity}, justified tier ${f.justified_tier}] ${f.file}:${f.line}
Quote: ${JSON.stringify(f.quote)}
Claim: ${f.claim}
Explanation: ${f.explanation}
Open the source yourself (${REPO}). Lean: does truth follow only from unfolding a def / an axiom / a structure field? Use #print axioms if needed (scratch file in /tmp, \`lake env lean\`). Physics: standard results (Polchinski II; Gimon-Polchinski hep-th/9601038; Sen hep-th/9605150; Bouwknegt-Evslin-Mathai hep-th/0306062; Coleman 1977; Coleman-De Luccia 1980; Eguchi-Ooguri-Tachikawa 1004.0956; Ooguri-Vafa hep-th/0605264) — do not quote equations from memory as fact; if you rely on one, say so and lower confidence. Numerics: compute with python3.
If category is "sound", hunt for a defect instead. Set corrected_tier if the tier is wrong.
refuted=true if wrong, overstated, or mis-categorised; if genuinely unsure, refuted=true with confidence low.` },
  { key: 'steelman', model: 'haiku', prompt: f => `Argue the AUTHORS' side. Is there a standard convention, reading, or explicitly stated scope under which the flagged item is acceptable as written?
[${f.category}] ${f.file}:${f.line} — ${f.claim}
Quote: ${JSON.stringify(f.quote)}
Critic: ${f.explanation}
Open the source (${REPO}). refuted=true ONLY if you name a concrete convention/reference/scope in the repo text that makes it correct as stated. Otherwise refuted=false.` },
]

const groups = A.only ? GROUPS.filter(g => A.only.includes(g.id)) : GROUPS
if (A.only) log(`Restricted to groups: ${groups.map(g => g.id).join(', ')}`)

const perGroup = await pipeline(
  groups,
  g => agent(extractPrompt(g), { label: `extract:${g.id}`, phase: 'Extract', schema: FINDINGS, model: 'haiku' }),
  (res, g) => {
    if (!res) return []
    const valid = res.findings.filter(f => f.quote && f.quote.trim().length > 3 && f.line >= 1)
    if (valid.length < res.findings.length) log(`${g.id}: dropped ${res.findings.length - valid.length} finding(s) lacking line/quote`)
    return parallel(valid.map((f, i) => async () => {
      const [fid, dom, steel] = await parallel(LENSES.map(l => () =>
        agent(l.prompt(f), { label: `${l.key}:${g.id}#${i}`, phase: 'Verify', schema: VERDICT, model: l.model, effort: l.effort })))
      let status
      if (!fid || fid.refuted) status = 'rejected_quote'
      else if (dom && dom.refuted && steel && steel.refuted) status = 'rejected'
      else if ((dom && dom.refuted) || (steel && steel.refuted)) status = 'disputed'
      else status = 'confirmed'
      const tier = dom && dom.corrected_tier && dom.corrected_tier !== 'unchanged' ? dom.corrected_tier : f.justified_tier
      return { group: g.id, ...f, justified_tier: tier, status, votes: { fidelity: fid, domain: dom, steelman: steel } }
    }))
  },
)

const all = perGroup.filter(Boolean).flat().filter(Boolean)

// Replacement map: local Lean file -> LeanMaster candidate (content from the 2026-09-17 methodology brief; adequacy re-judged here)
const REPLACEMENTS = [
  { local: ['proofs/BuscherRules.lean', 'proofs/LeanscratchDB/DoubleScaleT2.lean'], candidates: ['DualScaleStream2/DFT/GeneralizedMetric.lean', 'DualScaleStream2/TDuality/ODD.lean', 'DoubleFieldTheory/TDualityBuscher.lean'] },
  { local: ['proofs/MukaiLatticeK3.lean', 'proofs/FourierMukai.lean'], candidates: ['DualScaleStream2/Lattice/Mukai.lean', 'DualScaleStream2/Lattice/Reflection.lean'] },
  { local: ['proofs/TadpoleCancellation.lean', 'proofs/KummerTDAAnomalyCertification.lean'], candidates: ['DualScaleStream2/Flux/Tadpole.lean', 'StringTheoryFormalization/StringDynamics/KummerBlowup.lean'] },
  { local: ['proofs/KummerOrbifoldResolution.lean'], candidates: ['StringTheoryFormalization/StringDynamics/KummerBlowup.lean', 'DualScaleStream2/Lattice/K3T2Signature.lean'] },
  { local: ['proofs/MathieuVertexOperators.lean'], candidates: ['DualScaleStream2/Moonshine/EOT.lean'] },
  { local: ['proofs/GysinBEMSequence.lean', 'proofs/LeanscratchDB/TopologicalTDuality.lean'], candidates: ['StringTheoryFormalization/StringDynamics/TDualityGysin.lean'] },
  { local: ['proofs/SwamplandDistanceConjecture.lean'], candidates: ['StringTheoryFormalization/StringDynamics/SwamplandSafe.lean', 'StringTheoryFormalization/Frontier/ModuliGeodesics.lean'] },
  { local: ['proofs/FluxVacuumDecayCTheorem.lean', 'proofs/TachyonCondensationKTheory.lean', 'proofs/DbraneBoundaryStates.lean', 'proofs/FTheoryCosmology.lean', 'proofs/LeanscratchDB/HoloAlg.lean', 'proofs/LeanscratchDB/QuantumEntanglement.lean', 'proofs/LeanscratchDB/EmergentCosmology.lean'], candidates: [] },
]

const REPL = {
  type: 'object',
  required: ['entries'],
  properties: {
    entries: {
      type: 'array',
      items: {
        type: 'object',
        required: ['local_file', 'local_verdict', 'recommended_action', 'candidate', 'candidate_adequacy', 'evidence', 'toolchain_note'],
        properties: {
          local_file: { type: 'string' },
          local_verdict: { type: 'string', description: 'tier-honest one-liner: e.g. "vacuous: mass bound defined as 0"' },
          recommended_action: { type: 'string', enum: ['replace_with_candidate', 'keep_as_tier_B_arithmetic', 'downgrade_to_tier_C_scaffold', 'delete', 'needs_new_formalization'] },
          candidate: { type: 'string', description: 'LeanMaster module path, or "none"' },
          candidate_adequacy: { type: 'string', enum: ['adequate', 'partial', 'inadequate', 'none'] },
          theorem_mapping: { type: 'array', items: { type: 'string' }, description: 'local theorem -> LeanMaster theorem (verbatim names from files)' },
          evidence: { type: 'string', description: 'quoted Lean statements (not docstrings) from both sides' },
          toolchain_note: { type: 'string' },
        },
      },
    },
  },
}

let replacement = []
if (!A.skipReplacementMap) {
  phase('Replacement')
  replacement = (await parallel(REPLACEMENTS.map((r, i) => () => agent(
`Statement-adequacy review for replacing local Lean files with the maintained string-theory formalization.
Local files (in ${REPO}): ${r.local.join(', ')} — read them fully.
Candidate modules (in ${LM}, READ-ONLY, another session is working there): ${r.candidates.length ? r.candidates.join(', ') : 'none listed — search with: grep -rln "<key concept>" ' + LM + '/DualScaleStream2 ' + LM + '/StringTheoryFormalization'}.
For each local file: read the local theorem STATEMENTS; read the candidate STATEMENTS (not docstrings); decide whether the candidate states the same mathematical content non-vacuously. Check the candidate's axiom status in ${LM}/LEDGER.md or by grep for "axiom " / "sorry" / "native_decide" in the candidate file.
WARNING: ${LM}/DualScaleM24Formalization/FrontierTriad/* and Moonshine/KummerTadpole.lean copy this repo's vacuous definitions (c := 100*(n+1), rrCharge = rank+c1+c2, hard-coded V/E/beta1) — never propose those as replacements.
Downstream import: the tested pattern is ${LM}/examples/consumer_demo (toolchain v4.33.1, path require, packagesDir reuse). Toolchains: this repo ${'`'}cat ${REPO}/lean-toolchain${'`'} (mathlib unpinned) vs ${'`'}cat ${LM}/lean-toolchain${'`'} and the mathlib rev in ${LM}/lakefile.lean — state compatibility for a lake ${'`'}require${'`'}.
FIRST load the skills \`string-theory-foundation\` and \`leanmaster-theorem-search\`. The authoritative list of verified LeanMaster statements is ${LM}/docs/VERIFIED_FOUNDATION.md §3 (verbatim #check + #print axioms, release v2.1.0) — prefer candidates that appear there, quote them from that file, and mark candidates absent from it as partial at best. Note: an axiom audit of this repo's proofs/ (83 theorems, 5 failing on 2026-09-17) shows most local theorems use only standard axioms yet are vacuous — judge adequacy on statements, not on the audit. with the Skill tool (if available) and use their entry points / search index to find counterparts; the candidate list above is only a starting hint.
Rules: quote statements verbatim; no memory-sourced math; if a local theorem has no counterpart say needs_new_formalization.`,
    { label: `replace#${i}:${r.local[0].split('/').pop()}`, phase: 'Replacement', schema: REPL, model: JUDGE === 'haiku' ? 'haiku' : JUDGE }))))
    .filter(Boolean).flatMap(x => x.entries)
}

const order = { critical: 0, high: 1, medium: 2, low: 3 }
const pick = s => all.filter(f => f.status === s).sort((a, b) => order[a.severity] - order[b.severity])
const confirmed = pick('confirmed'), disputed = pick('disputed')
const rejected = all.filter(f => f.status.startsWith('rejected'))
const esc = s => String(s == null ? '' : s).replace(/\|/g, '\\|').replace(/\n/g, ' ')
const row = f => `| ${f.severity} | ${f.justified_tier} | ${f.category} | \`${esc(f.file)}:${f.line}\` | ${esc(f.claimed_status)} | ${esc(f.claim)} | ${esc(f.explanation)} | ${esc(f.suggested_fix)} |`
const table = rows => rows.length ? ['| sev | tier | category | location | presented as | claim | why | fix |', '|---|---|---|---|---|---|---|---|', ...rows.map(row)].join('\n') : '_none_'
const counts = {}
for (const f of confirmed) counts[f.category] = (counts[f.category] || 0) + 1
const footer = `\n---\nGenerated-by: physics-claims-audit workflow (Haiku; domain judge ${JUDGE}) | Verified-by: 3-lens Haiku skeptics (not a kernel proof) | Reviewed-by: T0 N\n`

const ledger = `# Physics & Formal-Claims Ledger

Negative results first. Tiers follow Mathesis (A kernel-verified & adequate, B exact arithmetic w/ negative control, L literature, C conjecture, X exploratory).
${all.length} findings — ${confirmed.length} confirmed, ${disputed.length} disputed, ${rejected.length} rejected.

## Confirmed defects
${table(confirmed.filter(f => f.category !== 'sound'))}

## Confirmed by category
${Object.entries(counts).map(([k, v]) => `- ${k}: ${v}`).join('\n') || '_none_'}

## Disputed — needs a T0/physicist ruling
${table(disputed)}
${disputed.map(f => `- \`${esc(f.file)}:${f.line}\` domain: ${esc(f.votes.domain && f.votes.domain.evidence)} / steelman: ${esc(f.votes.steelman && f.votes.steelman.evidence)}`).join('\n')}

## Confirmed sound claims
${table(confirmed.filter(f => f.category === 'sound'))}

## Rejected (transparency)
${rejected.map(f => `- [${f.status}] \`${esc(f.file)}:${f.line}\` ${esc(f.claim)}`).join('\n') || '_none_'}
${footer}`

const replMd = `# Local Lean → LeanMaster Replacement Map

Local \`proofs/\` is review-only; the string-theory formalization is maintained in SocrateAI-Scientific-Agora-LeanMaster. Adequacy judged on statements, not docstrings.

| local file | local verdict | action | candidate | adequacy | toolchain |
|---|---|---|---|---|---|
${replacement.map(e => `| \`${esc(e.local_file)}\` | ${esc(e.local_verdict)} | ${e.recommended_action} | \`${esc(e.candidate)}\` | ${e.candidate_adequacy} | ${esc(e.toolchain_note)} |`).join('\n')}

## Theorem mappings and evidence
${replacement.map(e => `### ${esc(e.local_file)}\n${(e.theorem_mapping || []).map(m => `- ${esc(m)}`).join('\n')}\n\n${esc(e.evidence)}`).join('\n\n')}
${footer}`

// Files are written by the caller from the returned markdown (avoids lossy copy through an agent).
return { total: all.length, confirmed: confirmed.length, disputed: disputed.length, rejected: rejected.length, byCategory: counts, replacement, ledger, replacementMap: replMd, findings: all }
