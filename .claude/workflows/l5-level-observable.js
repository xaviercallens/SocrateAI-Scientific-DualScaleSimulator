export const meta = {
  name: 'l5-level-observable',
  description: 'Decision D3 (2026-09-28): finish spec item L5 -- does the K3 level reach a reported observable through the dynamics? Pre-registered addendum A10 committed before any data; log mode only; starts at y in {0.5, 0.8, 1.2} from the SAME point at N=7 and N=12; bounded Jacobian with positive and negative controls; two skeptics; verdict in {met, not met, inconclusive}',
  whenToUse: 'Milestone M2 in ROADMAP.md, after level-explicit-parameter is merged. Needs >= 12 GB free RAM (check free -g first) and 1-3 h. Adopts nothing: it produces the evidence for the M3 adoption gate.',
  phases: [
    { title: 'Register', detail: 'A10 addendum to audit/PRE_REGISTRATION.md, committed before data; refuses to continue otherwise' },
    { title: 'Integrator', detail: 'bounded Jacobian through the dV/dx kink; positive control (completes at both N) and negative control (byte-identical repeat)' },
    { title: 'Runs', detail: '2 levels x 3 starts x 4 seeds x 2 tolerances in log mode, 3 processes' },
    { title: 'Skeptic', detail: 'numerics/tolerance skeptic and registration/circularity skeptic, independent' },
    { title: 'Synthesize', detail: 'verdict + report text (orchestrator saves audit/l5_level_runs/REPORT.md)' },
  ],
}

const REPO = '/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator'
const WT = '/mnt/disks/disk-socrateai-local-1/dualscale-wt-l5'
const BRANCH = 'loop/l5-level'
const OUT = 'audit/l5_level_runs'
const PROCS = (args && args.procs) || 3

const RULES = `
GROUND RULES:
1. Work ONLY in git worktree ${WT} on branch ${BRANCH} (create with: cd ${REPO} && git worktree add ${WT} -b ${BRANCH} main; omit -b if the branch exists). Outputs under ${OUT}/. Never touch proofs/, main, other worktrees; never push.
2. Evidence-bound: every number you return comes from a script you ran; save script + JSON; every claim carries the command and a verbatim output line.
3. No literal answers in a computation path. The 2026-09-27 numbers (0.1116 vs 0.0747, etc.) may appear only in an "expected" field, labelled as such.
4. Stage only your own files by explicit path; never git add -A. Commit messages end with: Co-Authored-By: Claude Haiku 4.5 <noreply@anthropic.com>
5. Wording: never "proved"/"proof"/"theorem" for your own work. "Measured", "holds for every case run", "tier B".
6. Memory: before launching runs, check free -g; if available < 10 GB, STOP and return {aborted: "memory"} instead of running.
`

const DESIGN = `
DESIGN (fixed by decision D3, audit/DUALSCALE_TO_STREAMS_DECISIONS_2026_09_28.md; do not change without saying so):
- Potential: leanflow.core.modular_potential, mode="log" ONLY (ratio mode has measured zero force outside y in [0.36, 0.90]; audit/l5_level_runs/force_scan.json).
- Levels: N = 7 and N = 12. Same start point for both levels (the level must enter only through the dynamics).
- Starts: (x, y) in {(0.2, 0.5), (0.2, 0.8), (0.2, 1.2)}; 4 seeds = jitter of the start by uniform(-1e-3, 1e-3), seed 0 = no jitter.
- Tolerance per observable and per N = |obs(rtol=1e-8) - obs(rtol=1e-9)|; the run tolerance is the larger of the two levels'.
- Observables: w0, wa (CPL fit via workshopcosmo.run_observables_analysis), w_final, omega_phi_final, final tau. t_max = 70, 500 output points, Radau, atol = 1e-2 * rtol, max_step = 0.5 unless the integrator phase changes it (then say so).
- Criterion (A10): L5 is MET iff for at least one observable, at at least one start, |obs(12) - obs(7)| > tolerance for EVERY seed with the SAME sign. NOT MET iff no observable at any start does. INCONCLUSIVE iff fewer than 1 completed run per (N, start) pair, or the negative control fails.
- Base script: scripts/l5_level_experiment.py (already reduces x mod 1 and clamps y to [1e-8, 60] -- reuse _bounded/potential; its verdicts() implements the criterion per start; extend it to the three starts and the 1e-9 tolerance reference, do not rewrite it).
`

const REG_SCHEMA = { type: 'object', properties: {
  commit: { type: 'string' }, addendum_heading: { type: 'string' }, addendum_text: { type: 'string' },
  committed_before_any_run: { type: 'boolean' } }, required: ['commit', 'addendum_heading', 'addendum_text', 'committed_before_any_run'] }

const INTEG_SCHEMA = { type: 'object', properties: {
  commit: { type: 'string' }, method: { type: 'string' },
  positive_control: { type: 'object', properties: { command: { type: 'string' }, N7_completed: { type: 'boolean' }, N12_completed: { type: 'boolean' }, output_lines: { type: 'array', items: { type: 'string' } } }, required: ['command', 'N7_completed', 'N12_completed'] },
  negative_control: { type: 'object', properties: { command: { type: 'string' }, byte_identical: { type: 'boolean' }, sha256_a: { type: 'string' }, sha256_b: { type: 'string' } }, required: ['command', 'byte_identical'] },
  design_changes: { type: 'array', items: { type: 'string' } }, aborted: { type: 'string' } },
  required: ['commit', 'method', 'positive_control', 'negative_control', 'design_changes'] }

const RUNS_SCHEMA = { type: 'object', properties: {
  commit: { type: 'string' }, results_json: { type: 'string' }, command: { type: 'string' },
  runs_total: { type: 'integer' }, runs_completed: { type: 'integer' },
  per_pair_completed: { type: 'array', items: { type: 'object', properties: { N: { type: 'integer' }, start_y: { type: 'number' }, completed: { type: 'integer' } }, required: ['N', 'start_y', 'completed'] } },
  verdict_rows: { type: 'array', items: { type: 'object', properties: { start_y: { type: 'number' }, observable: { type: 'string' }, meets: { type: 'boolean' }, min_abs_diff: { type: 'number' }, max_tolerance: { type: 'number' }, sign_consistent: { type: 'boolean' } }, required: ['start_y', 'observable', 'meets'] } },
  wall_seconds: { type: 'number' }, aborted: { type: 'string' } },
  required: ['commit', 'results_json', 'command', 'runs_total', 'runs_completed', 'per_pair_completed', 'verdict_rows'] }

const SKEPTIC_SCHEMA = { type: 'object', properties: {
  lens: { type: 'string' }, verdict_stands: { type: 'boolean' },
  findings: { type: 'array', items: { type: 'object', properties: { severity: { type: 'string' }, file: { type: 'string' }, quote: { type: 'string' }, why: { type: 'string' }, fix: { type: 'string' } }, required: ['severity', 'quote', 'why'] } },
  recomputed: { type: 'array', items: { type: 'object', properties: { what: { type: 'string' }, command: { type: 'string' }, agrees: { type: 'boolean' } }, required: ['what', 'command', 'agrees'] } } },
  required: ['lens', 'verdict_stands', 'findings', 'recomputed'] }

const SYNTH_SCHEMA = { type: 'object', properties: {
  verdict: { type: 'string', enum: ['met', 'not_met', 'inconclusive'] },
  one_paragraph: { type: 'string' }, report_markdown: { type: 'string' },
  manuscript_sentence: { type: 'string' }, adoption_gate_input: { type: 'object', properties: { l5_met: { type: 'boolean' }, generic_start_reached_stable_state: { type: 'boolean' }, evidence: { type: 'string' } }, required: ['l5_met', 'generic_start_reached_stable_state', 'evidence'] } },
  required: ['verdict', 'one_paragraph', 'report_markdown', 'manuscript_sentence', 'adoption_gate_input'] }

phase('Register')
const reg = await agent(`${RULES}\n${DESIGN}
Write addendum "## Addendum 2026-09-28 (A10, appended; A1-A9 and the original rules above unchanged)" to audit/PRE_REGISTRATION.md in ${WT}. It states, in the style of A6-A9: the design above verbatim, the criterion, the reading of each outcome ("not met" => "on this potential the level does not reach observables through the dynamics; no level-selection claim may appear in any manuscript"; "met" => "a level-dependent observable exists on the experimental potential; this is NOT a K3 selection and NOT an adoption"), and that ratio mode is excluded with the measured reason. Commit it ALONE (git add audit/PRE_REGISTRATION.md). Do NOT run any simulation. Return the commit SHA and the text.`,
  { label: 'register-A10', phase: 'Register', schema: REG_SCHEMA, model: 'haiku' })
if (!reg || !reg.committed_before_any_run) return { aborted: 'registration not committed', reg }

phase('Integrator')
const integ = await agent(`${RULES}\n${DESIGN}
Registration is committed at ${reg.commit}; do not edit PRE_REGISTRATION.md.
The 2026-09-27 attempt stalled: Radau's numerical Jacobian enlarges its probe step without bound where dV/dx is identically 0.0 in float64 (audit/l5_level_runs/README.md, second and third passes). Fix the integrator in scripts/l5_level_experiment.py (extend, do not rewrite): pass an explicit jac= to solve_ivp built by central differences of rhs with the SAME _bounded clamp (step 1e-6 on x, y, u, v; analytic in a), OR reduce max_step through y in [0.25, 0.5]; state which in "method".
POSITIVE CONTROL: one run each at N=7 and N=12, mode log, start (0.2, 0.8), seed 0, rtol 1e-8, t_max 70 -- both must complete (sol.success True). NEGATIVE CONTROL: repeat the N=7 run; the two results JSON must be byte-identical (sha256). Commit the script change. If free -g shows < 10 GB available, return aborted="memory".`,
  { label: 'integrator', phase: 'Integrator', schema: INTEG_SCHEMA, model: 'haiku' })
if (!integ || integ.aborted) return { aborted: (integ && integ.aborted) || 'integrator died', reg, integ }
if (!integ.positive_control.N7_completed || !integ.positive_control.N12_completed || !integ.negative_control.byte_identical) {
  return { verdict: 'inconclusive', reason: 'integrator controls failed', reg, integ }
}

phase('Runs')
const runs = await agent(`${RULES}\n${DESIGN}
Registration ${reg.commit}; integrator ${integ.commit} (method: ${integ.method}). Run the full design: 2 levels x 3 starts x 4 seeds x 2 tolerances = 48 runs, mode log, --procs ${PROCS}. Write ${OUT}/results_l5.json by script (never by hand), including per-run success, nfev, seconds, clamp_events, and the verdict rows per start. Log progress lines. If free -g shows < 10 GB available before starting, return aborted="memory". Commit results_l5.json and any script edits. Report wall time.`,
  { label: 'runs-48', phase: 'Runs', schema: RUNS_SCHEMA, model: 'haiku' })
if (!runs || runs.aborted) return { aborted: (runs && runs.aborted) || 'runs died', reg, integ }

phase('Skeptic')
const skeptics = await parallel([
  () => agent(`${RULES}
You are the NUMERICS SKEPTIC. Registration ${reg.commit}, results ${runs.results_json} at commit ${runs.commit} in ${WT}. Try to REFUTE the verdict rows ${JSON.stringify(runs.verdict_rows)}:
 - recompute at least two verdict rows yourself from results_l5.json with your own script;
 - check the tolerance is |obs(1e-8) - obs(1e-9)| per N and the max was taken;
 - check that clamp_events == 0 on every completed run's trajectory (a clamped STATE, not an FD probe, invalidates that run);
 - check every completed run has n_fit_points >= 5 for the CPL fit;
 - check that the same start was used at both N (x0, y0 identical per seed).
verdict_stands=false if any finding is severity high.`, { label: 'skeptic-numerics', phase: 'Skeptic', schema: SKEPTIC_SCHEMA }),
  () => agent(`${RULES}
You are the REGISTRATION/CIRCULARITY SKEPTIC. Compare the A10 addendum text at ${reg.commit} (git -C ${WT} show ${reg.commit}:audit/PRE_REGISTRATION.md) with what was actually run (${runs.command}, script at ${runs.commit}). Try to REFUTE:
 - any design change between A10 and the run (starts, seeds, tolerances, t_max, mode, criterion) not declared in integ.design_changes=${JSON.stringify(integ.design_changes)};
 - any literal expected value inside a computation path (grep the scripts for 0.1116, 0.0747, 0.378, 0.289);
 - whether the registration commit precedes every results file commit (git log --format=%H --reverse);
 - whether the criterion in verdicts() matches A10 word for word in effect.
verdict_stands=false if any finding is severity high.`, { label: 'skeptic-registration', phase: 'Skeptic', schema: SKEPTIC_SCHEMA }),
])
const [numerics, registration] = skeptics

phase('Synthesize')
const synth = await agent(`${RULES}
Write the L5 report. Inputs: registration ${JSON.stringify(reg)}; integrator ${JSON.stringify(integ)}; runs ${JSON.stringify(runs)}; skeptics ${JSON.stringify(skeptics.filter(Boolean))}.
Verdict rule: "met" only if >= 1 verdict row meets AND both skeptics' verdict_stands are true; "inconclusive" if any (N, start) pair has 0 completed runs, or a skeptic does not stand; else "not_met". Negative results first. report_markdown is the full audit/l5_level_runs/REPORT.md: design (from A10), what completed, the table of verdict rows with min |diff| and max tolerance, skeptic findings verbatim, verdict, and the "reading" sentence from A10 that applies. manuscript_sentence is the one sentence a paper may carry. adoption_gate_input.generic_start_reached_stable_state: true only if, at both N, a run from a generic start has |d tau/dt| decreasing over the last 20% of t and omega_phi_final in (0, 1) -- compute it from results_l5.json, do not guess.`,
  { label: 'synthesize', phase: 'Synthesize', schema: SYNTH_SCHEMA })

return {
  verdict: synth ? synth.verdict : 'inconclusive',
  registration: reg, integrator: integ, runs, skeptics: skeptics.filter(Boolean), synthesis: synth,
  report: synth ? synth.report_markdown : null,
  next: `Orchestrator: write report to ${WT}/${OUT}/REPORT.md, commit on ${BRANCH}, read it, then take the M3 decision (ROADMAP.md).`,
}
