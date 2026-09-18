export const meta = {
  name: 'zero-param-loop',
  description: 'Forward loop (theory -> experiment -> GUDHI TDA -> reduce) then reverse loop (reduced -> TDA -> experiment -> theory), incrementally cutting free parameters 6 -> 4 -> 2 -> 0 with independent skeptics',
  whenToUse: 'Iterate the dual-scale model toward zero free parameters using real public data, INRIA GUDHI, and LeanMaster theorems as fixed constants.',
  phases: [
    { title: 'Ground', detail: 'Lean constants inventory, real data fetch, parametrized simulator harness' },
    { title: 'Check harness', detail: 'independent check of harness + data before any round' },
    { title: 'Forward', detail: 'sweep vs real data, GUDHI persistent homology, propose one reduction step' },
    { title: 'Reverse', detail: 'reduced model -> predicted topology -> refit real data -> compare' },
    { title: 'Skeptic', detail: 'independent refutation of the round' },
    { title: 'Adapt theory', detail: 'three angles on how the hypothesis must change' },
    { title: 'Synthesize', detail: 'ledger + report' },
  ],
}

const WT = '/mnt/disks/disk-socrateai-local-1/dualscale-wt-loop'
const PY = '/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python'
const LM = '/home/callensxavier_gmail_com/SocrateAI-Scientific-Agora-LeanMaster'
const MAX_ROUNDS = (args && args.maxRounds) || 3

const RULES = `
GROUND RULES (strict):
- Work ONLY inside the git worktree ${WT} (branch loop/zero-param). Never touch proofs/, never touch main, never push, never deposit to Zenodo. LeanMaster (${LM}) is READ-ONLY.
- Python: always use ${PY} (has gudhi 3.13, numpy, scipy, pandas, scikit-learn, networkx, matplotlib). Seed every random call (seed 42).
- Evidence-bound: every number you report must come from a command you ran; give the command and the file it wrote. Never invent experimental bounds, DOIs or citations. If you cannot fetch or verify something, write "ABSENT" and continue.
- Tiers (Mathesis): A kernel-checked Lean statement exactly as written; B exact arithmetic with negative control; L literature; C conjecture; X exploratory numerics. A Lean theorem fixes a NUMBER; the claim that this number IS the physical parameter is tier C unless derived.
- A parameter counts as removed only if (a) a tier-A/B/L constant replaces it, or (b) it provably affects no retained observable and is deleted from the model, or (c) it is absorbed in a combination (then the combination still counts as 1). Insensitivity alone is NOT derivation: label it "fixed by convention" and keep it in the count as 0.5-evidence, reported separately.
- Known broken, do not use as evidence: vacuum_decay_cdl (floored action), tadpole_cancellation_certified flag, local proofs/BuscherRules.lean.
- Commit your files on loop/zero-param with a clear message ending with the line: Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
`

const GROUND_LEAN = {
  type: 'object',
  properties: {
    leanmaster_tag: { type: 'string' },
    constants: { type: 'array', items: { type: 'object', properties: {
      name: { type: 'string' }, file_line: { type: 'string' }, statement_verbatim: { type: 'string' },
      value: { type: 'string' }, candidate_param: { type: 'string', description: 'which simulator parameter it could fix, or none' },
      identification_tier: { type: 'string' }, caveat: { type: 'string' } },
      required: ['name', 'file_line', 'statement_verbatim', 'candidate_param', 'identification_tier'] } },
  },
  required: ['leanmaster_tag', 'constants'],
}
const GROUND_DATA = {
  type: 'object',
  properties: { datasets: { type: 'array', items: { type: 'object', properties: {
    id: { type: 'string' }, domain: { type: 'string' }, path: { type: 'string' }, url: { type: 'string' },
    sha256: { type: 'string' }, rows: { type: 'number' }, status: { type: 'string', description: 'FETCHED or ABSENT' }, note: { type: 'string' } },
    required: ['id', 'domain', 'status'] } } },
  required: ['datasets'],
}
const HARNESS = {
  type: 'object',
  properties: {
    script: { type: 'string' }, params: { type: 'array', items: { type: 'string' } },
    observables: { type: 'array', items: { type: 'string' } }, wired_pta: { type: 'boolean' },
    nominal_run_command: { type: 'string' }, nominal_output: { type: 'string' }, commit: { type: 'string' }, problems: { type: 'array', items: { type: 'string' } } },
  required: ['script', 'params', 'observables', 'nominal_run_command', 'problems'],
}
const CHECK = {
  type: 'object',
  properties: { usable: { type: 'boolean' }, reproduced_nominal: { type: 'boolean' },
    data_verified: { type: 'array', items: { type: 'string' } }, data_rejected: { type: 'array', items: { type: 'string' } },
    blocking_issues: { type: 'array', items: { type: 'string' } }, fixes_applied: { type: 'array', items: { type: 'string' } } },
  required: ['usable', 'blocking_issues'],
}
const FORWARD = {
  type: 'object',
  properties: {
    sweep_file: { type: 'string' }, n_points: { type: 'number' },
    chi2: { type: 'object', properties: { best_model: { type: 'number' }, lcdm: { type: 'number' }, dof: { type: 'number' }, datasets_used: { type: 'array', items: { type: 'string' } } } },
    effective_dimension: { type: 'number', description: 'rank of observable response (Fisher/PCA), i.e. how many parameter combinations the data can see' },
    tda: { type: 'object', properties: { tool: { type: 'string' }, betti_model: { type: 'string' }, betti_real: { type: 'string' },
      bottleneck_model_vs_real: { type: 'number' }, null_control: { type: 'string' }, known_answer_control: { type: 'string' }, files: { type: 'array', items: { type: 'string' } } } },
    proposed_reductions: { type: 'array', items: { type: 'object', properties: {
      param: { type: 'string' }, mechanism: { type: 'string', description: 'lean_constant | delete_unobservable | absorb_combination | convention' },
      replacement: { type: 'string' }, evidence: { type: 'string' }, tier: { type: 'string' } },
      required: ['param', 'mechanism', 'evidence', 'tier'] } },
    commit: { type: 'string' },
  },
  required: ['sweep_file', 'effective_dimension', 'tda', 'proposed_reductions'],
}
const REVERSE = {
  type: 'object',
  properties: {
    reduced_params: { type: 'array', items: { type: 'string' } },
    predicted_topology: { type: 'string' }, observed_topology: { type: 'string' }, topology_match: { type: 'boolean' },
    chi2_reduced: { type: 'number' }, chi2_full: { type: 'number' }, chi2_lcdm: { type: 'number' }, delta_aic_vs_full: { type: 'number' }, delta_aic_vs_lcdm: { type: 'number' },
    fit_degraded: { type: 'boolean' }, theory_feedback: { type: 'string' }, files: { type: 'array', items: { type: 'string' } }, commit: { type: 'string' } },
  required: ['reduced_params', 'topology_match', 'fit_degraded', 'theory_feedback'],
}
const VERDICT = {
  type: 'object',
  properties: { per_reduction: { type: 'array', items: { type: 'object', properties: {
      param: { type: 'string' }, accepted: { type: 'boolean' }, counts_as: { type: 'string', description: 'derived | deleted | absorbed | convention_only | rejected' }, reason: { type: 'string' }, command_rerun: { type: 'string' } },
      required: ['param', 'accepted', 'counts_as', 'reason'] } },
    tda_valid: { type: 'boolean' }, tda_issues: { type: 'array', items: { type: 'string' } }, fabricated_items: { type: 'array', items: { type: 'string' } } },
  required: ['per_reduction', 'tda_valid'],
}

// ---------------- Ground ----------------
phase('Ground')
const [lean, data, harness] = await parallel([
  () => agent(`${RULES}
TASK: inventory LeanMaster constants that could REPLACE a free parameter of the dual-scale simulator.
Simulator free parameters: a_pot, b_pot (quintessence potential, workshopcosmo.py ~L156), mu_sym, lambda_sym (symmetron), pta_suppression, c4_c0_ratio (NANOGrav hexadecapole), plus hard-coded h1=1/4, h2=5/4.
Read (grep theorem/def, then open) ${LM}/DualScaleCosmology/*.lean (DarkEnergyScale: rhoLambda, selfDual_lambda4_eq, omegaLambda; SelfDualCutoff; CKNBound; CosmicString: fStringGmu_selfDual; DualTower; ScaleFactorDuality), ${LM}/DualScaleValidation/Observables.lean, ${LM}/DualScaleStream2 (dual-scale bound, dualScale_eq_iff, r_bps_from_eot), ${LM}/DualScaleMoonshine (conformal weights if any: grep "1/4", "5/4", "weight"). Run: cd ${LM} && git describe --tags.
For each usable constant give the VERBATIM statement and file:line, the numeric value, which parameter it could fix, and the tier of the physical identification (usually C). Do not build Lean. Do not claim a theorem exists without quoting it.`,
    { label: 'ground:lean-constants', phase: 'Ground', model: 'haiku', schema: GROUND_LEAN }),
  () => agent(`${RULES}
TASK: fetch REAL public data into ${WT}/data/real/ (create it), one subfolder per domain, with a data/real/MANIFEST.json (id, url, sha256, rows, columns, licence/citation as printed by the source itself).
Domains wanted: (1) dark energy: DESI BAO distance measurements (try the CobayaSampler/bao_data GitHub repo, desi_* files) and Pantheon+ SN (github PantheonPlusSH0ES/DataRelease, Pantheon+SH0ES.dat); (2) pulsar timing: NANOGrav 15-yr Hellings-Downs binned correlation or free-spectrum products (Zenodo / github nanograv/15yr_stochastic_analysis) ; (3) fifth force / screening: only if a machine-readable public table exists, else ABSENT (do NOT type bounds from memory); (4) cosmic web point cloud for TDA on real data: a small public galaxy catalogue slice (e.g. SDSS/2MRS/6dF via VizieR TAP or astroquery-free plain HTTP CSV), up to ~5000 objects with RA, Dec, z.
Use curl -L with timeouts. After each download: sha256sum, wc -l, head -5, and sanity-check it is data not an HTML error page. Anything that fails after 3 URL attempts: status ABSENT with the URLs tried. Keep total under 200 MB. Commit the manifest and small files (<5 MB); add larger ones to .gitignore and keep them on disk.`,
    { label: 'ground:real-data', phase: 'Ground', model: 'haiku', schema: GROUND_DATA }),
  () => agent(`${RULES}
TASK: build a parametrized experiment harness ${WT}/scripts/param_loop_sim.py (do not break existing tests; you may refactor workshopcosmo.py minimally so a_pot and b_pot become function arguments with the same defaults).
CLI: --params '{"a_pot":1.0,...}' --out file.json . It must compute, from a parameter vector (a_pot, b_pot, mu_sym, lambda_sym, pta_suppression, c4_c0_ratio):
 - dark-energy observables: w(z) on a z-grid, w0, wa (CPL fit), H(z)/H0 and D_M/r_d-style distance ratios at z = 0.3,0.5,0.7,1.0,1.5,2.3 and distance moduli on a z-grid (so they can be compared with BAO/SN data later);
 - screening observable from run_symmetron_screening_simulation (screening_suppression_factor, phi_center_ratio);
 - PTA observable: WIRE pta_suppression and c4_c0_ratio into an angular-correlation prediction Gamma(theta) on 15 angular bins = Hellings-Downs + hexadecapole (l=4 Legendre) term, using what run_nanograv_hexadecapole_simulation actually does; if pta_suppression appears nowhere in any computation, report that fact verbatim (grep output) rather than inventing a role for it.
Exclude vacuum decay. Must run in < 20 s per point. Also add --selftest that checks: defaults reproduce the current workshopcosmo w0/wa, and that changing each parameter by x10 changes at least one observable (print which do NOT: those are structurally unobservable). Run it, paste the selftest output into your result, run pytest -q tests (must stay green), commit.`,
    { label: 'ground:harness', phase: 'Ground', model: 'sonnet', schema: HARNESS }),
])

phase('Check harness')
const check = await agent(`${RULES}
You are an INDEPENDENT checker; you did not write any of this. Inputs:
HARNESS REPORT: ${JSON.stringify(harness)}
DATA REPORT: ${JSON.stringify(data)}
LEAN REPORT: ${JSON.stringify(lean)}
1. Re-run the harness selftest and nominal command yourself; confirm the outputs. Read the harness code: is any observable hard-coded, clamped, or independent of its parameters by construction? 
2. For each dataset marked FETCHED: open it, confirm it is real tabular data matching its manifest (sha256sum), reject HTML/error pages or files whose content does not match the claimed source.
3. For 5 Lean constants (or all if fewer): grep the quoted statement in the cited file; list any that are not found verbatim.
4. Negative control: confirm that the selftest FAILS if you temporarily make one observable ignore its parameter (do it on a scratch copy in /tmp, not in the repo).
Fix small blocking problems yourself (and commit); otherwise list them. usable=false only if no dark-energy data at all or the harness cannot run.`,
  { label: 'check:harness+data', phase: 'Check harness', model: 'sonnet', schema: CHECK })

if (!check || !check.usable) {
  log('Harness or data unusable - stopping before any round')
  return { stopped: 'ground', lean, data, harness, check }
}

// ---------------- Rounds ----------------
let free = ['a_pot', 'b_pot', 'mu_sym', 'lambda_sym', 'pta_suppression', 'c4_c0_ratio']
const ledger = []
const rounds = []
for (let r = 1; r <= MAX_ROUNDS && free.length > 0; r++) {
  const state = `ROUND ${r}. Currently free parameters (${free.length}): ${free.join(', ')}.
Accepted so far: ${JSON.stringify(ledger)}
Verified datasets: ${JSON.stringify(check.data_verified || [])}; rejected: ${JSON.stringify(check.data_rejected || [])}.
Lean constants available: ${JSON.stringify(lean && lean.constants)}
Harness: ${harness.script}; observables: ${JSON.stringify(harness.observables)}; harness problems: ${JSON.stringify(harness.problems)}`

  const fwd = await agent(`${RULES}
${state}
FORWARD LOOP: hypothesis -> experiment -> TDA -> reduce. Write everything under ${WT}/audit/zero_param_loop/round${r}/ .
1. EXPERIMENT: sweep ONLY the currently free parameters over WIDE log ranges (at least 3 decades each where positive; Latin hypercube or Sobol, seed 42, 300-600 points, run in parallel with multiprocessing, 6 workers) with the harness; already-fixed parameters take their ledger values. Save sweep.csv (params + observables). Compute chi2 against the verified real data (BAO with its covariance if provided; SN with diagonal errors if covariance is too large; PTA bins if present). Also compute chi2 of flat LCDM (Om fitted) on the same data as the baseline. Report honestly if the model is simply LCDM-like everywhere in the sweep.
2. RESPONSE RANK: numerical Jacobian of the standardized observable vector w.r.t. log-params at the best-fit point and at 5 random points; singular values; effective_dimension = number of singular values > 1e-3 of the largest. List which parameter directions are null.
3. TDA with INRIA GUDHI (import gudhi; no home-made substitutes, not scripts/tda_mapper.py): 
   a. model side: persistent homology (Rips or alpha complex, dims 0-2) of the standardized sweep point cloud in OBSERVABLE space -> does the image look like a point, a curve, a surface? (this is the topological estimate of how many parameters the observables can see). 
   b. real-data side: persistent homology of a real point cloud from a verified dataset (galaxy catalogue slice in comoving coordinates if present; otherwise SN/BAO residual embedding), same pipeline. 
   c. controls: known-answer (noisy circle must give one long H1 bar), null (shuffled/Poisson cloud with same N and volume). Report bottleneck distances (gudhi.bottleneck_distance) model-vs-real and null-vs-real. Save persistence diagrams (json + png).
4. REDUCE: propose the NEXT step only (aim to remove up to 2 parameters this round, target ladder 6 -> 4 -> 2 -> 0). For each, give mechanism (lean_constant / delete_unobservable / absorb_combination / convention), the replacement value or formula, the evidence file, and the tier. Prefer lean_constant where a quoted LeanMaster statement supplies the number (e.g. the self-dual dark-energy scale rhoLambda(lP, L, Omega) for the potential height); say plainly when the identification is tier C.
Commit.`, { label: `r${r}:forward`, phase: 'Forward', model: 'sonnet', schema: FORWARD })
  if (!fwd) { log(`round ${r}: forward agent failed`); break }

  const rev = await agent(`${RULES}
${state}
FORWARD RESULT (from another agent): ${JSON.stringify(fwd)}
REVERSE LOOP: reduced parameters -> TDA -> experiment -> hypothesis. Work under ${WT}/audit/zero_param_loop/round${r}/reverse/ .
1. Build the REDUCED model: apply the proposed reductions (freeze / delete / substitute Lean constant) via harness params; list the parameters still free.
2. TDA first: from the reduced model, state BEFORE looking at data what the persistence diagram of its observable image should be (dimension = number of remaining free params; a zero-parameter model predicts a single point = only H0, one component, no H1). Then compute it with GUDHI on a fresh sweep of the reduced model (seed 43) and compare with the prediction and with the real-data diagram from the forward step (bottleneck distance).
3. EXPERIMENT: refit the reduced model to the same real data; compare chi2 and AIC/BIC with the full model and with LCDM (use the forward agent's sweep.csv for the full model; recompute LCDM yourself). fit_degraded = true if delta chi2 (reduced - full) > 2 per removed parameter.
4. HYPOTHESIS feedback: in plain words, what does this round say the theory must change? (e.g. "potential parameters are invisible to data because the field is frozen: either the quintessence sector is removed and dark energy IS the self-dual constant, or the potential must be steepened by a derived, not tuned, amount").
Commit.`, { label: `r${r}:reverse`, phase: 'Reverse', model: 'sonnet', schema: REVERSE })

  const verdict = await agent(`${RULES}
${state}
You are the SKEPTIC for round ${r}. You produced none of this. Default to rejecting what you cannot reproduce.
FORWARD: ${JSON.stringify(fwd)}
REVERSE: ${JSON.stringify(rev)}
For EACH proposed reduction: re-run one decisive command yourself (e.g. vary that parameter x10 / x0.1 at the best fit with the harness and show the observables and chi2; or grep the quoted Lean statement in ${LM}); decide accepted or not and what it counts as: derived (tier A/B/L number AND a stated identification), deleted (parameter provably reaches no retained observable - then it must be removed from the model, not frozen), absorbed (say which combination remains free), convention_only (insensitive but still in the equations: does NOT reduce the count), rejected.
TDA check: confirm gudhi was really imported and run (read the script, re-run the known-answer control), that the null control is distinguishable from the real data, and that no topological claim exceeds what the diagrams show. List anything fabricated (bounds, citations, numbers without a producing command).`,
    { label: `r${r}:skeptic`, phase: 'Skeptic', schema: VERDICT })
  if (!verdict) { log(`round ${r}: skeptic failed - nothing accepted`); rounds.push({ r, fwd, rev, verdict: null }); break }

  const removedNow = verdict.per_reduction.filter(p => p.accepted && ['derived', 'deleted'].includes(p.counts_as))
  const absorbed = verdict.per_reduction.filter(p => p.accepted && p.counts_as === 'absorbed')
  removedNow.forEach(p => ledger.push({ round: r, param: p.param, how: p.counts_as, reason: p.reason }))
  absorbed.forEach(p => ledger.push({ round: r, param: p.param, how: 'absorbed', reason: p.reason }))
  const gone = new Set(removedNow.map(p => p.param))
  // absorbed: a pair becoming one combination removes one name; keep conservative: remove only if skeptic said accepted
  absorbed.forEach(p => gone.add(p.param))
  const before = free.length
  free = free.filter(p => !gone.has(p))
  rounds.push({ r, fwd, rev, verdict, free_after: [...free] })
  log(`round ${r}: ${before} -> ${free.length} free (${free.join(', ') || 'none'}); conventions/rejections: ${verdict.per_reduction.filter(p => !gone.has(p.param)).map(p => p.param + ':' + p.counts_as).join(', ') || 'none'}`)
  if (free.length === before) { log(`round ${r}: no accepted reduction - loop is dry, stopping`); break }
}

// ---------------- Adapt theory ----------------
phase('Adapt theory')
const ctx = `Remaining free parameters: ${free.join(', ') || 'none'}.
Ledger: ${JSON.stringify(ledger)}
Round summaries: ${JSON.stringify(rounds.map(x => ({ r: x.r, eff_dim: x.fwd && x.fwd.effective_dimension, chi2: x.fwd && x.fwd.chi2, tda: x.fwd && x.fwd.tda, reverse: x.rev && { match: x.rev.topology_match, degraded: x.rev.fit_degraded, d_aic_lcdm: x.rev.delta_aic_vs_lcdm, feedback: x.rev.theory_feedback }, verdict: x.verdict && x.verdict.per_reduction })))}
Lean constants: ${JSON.stringify(lean && lean.constants)}`
const ANGLES = [
  ['minimal', 'MINIMAL: delete every sector the data cannot see. What is the smallest zero-parameter theory left (constants only from LeanMaster + measured lP, H0), what does it predict that LCDM does not, and how could it be falsified with the fetched data?'],
  ['derive', 'DERIVE: keep the sectors but replace each remaining tuned number by a derivation target. For each remaining parameter name the precise LeanMaster statement that would have to be proved (give a Lean-style statement), what already exists toward it (quote), and whether the derivation is realistic or wishful.'],
  ['data', 'DATA-FIRST: start from what the real data and the GUDHI diagrams actually show (including the DESI w0-wa preference if the fetched BAO data support it). Which adaptation of the hypothesis (e.g. dual-scale evolving vacuum tied to a(t) + 1/a(t), thawing field with derived slope) could move the model toward the data WITHOUT a new free parameter? Give one decisive experiment per idea.'],
]
const ideas = (await parallel(ANGLES.map(([k, p]) => () => agent(`${RULES}
${ctx}
BRAINSTORM, angle = ${p}
Be concrete and short (max 500 words). Mark every sentence that is conjecture as (C). Do not claim anything is proved that is not quoted above. Read-only: write no files.`,
  { label: `adapt:${k}`, phase: 'Adapt theory' })))).filter(Boolean)

phase('Synthesize')
const report = await agent(`${RULES}
${ctx}
THEORY ADAPTATION IDEAS: ${JSON.stringify(ideas)}
GROUND: data=${JSON.stringify(data)} check=${JSON.stringify(check)}
Write ${WT}/audit/zero_param_loop/REPORT.md (plain language, for a physicist in a hurry) and ${WT}/audit/zero_param_loop/ledger.json:
1. Parameter ladder table: round, free before -> after, what was removed, how (derived/deleted/absorbed), tier, skeptic verdict; separate line for "convention only" items that do NOT count.
2. Positive and negative results per round (chi2 vs LCDM on real data, response rank, GUDHI findings with controls).
3. What the reverse loop says about the hypothesis; the adapted hypothesis chosen (merge the three angles, say which ideas were dropped and why).
4. Next decisive experiment and next Lean statement to prove (one each).
5. Everything marked ABSENT or fabricated by a skeptic, listed verbatim.
Header line: "Generated by workflow zero-param-loop; exploratory (tier X) unless stated." Commit. Return the report markdown as your final text.`,
  { label: 'synthesize:report', phase: 'Synthesize' })

return { free_remaining: free, ledger, rounds: rounds.map(x => ({ r: x.r, free_after: x.free_after, eff_dim: x.fwd && x.fwd.effective_dimension, chi2: x.fwd && x.fwd.chi2, verdict: x.verdict })), check, report }
