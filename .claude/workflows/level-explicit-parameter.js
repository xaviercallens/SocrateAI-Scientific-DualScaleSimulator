export const meta = {
  name: 'level-explicit-parameter',
  description: 'Decision D2 (2026-09-28): promote the K3 level N from a hidden module constant (FRICKE_Y = 1/sqrt(12)) to an explicit, reported argument of workshopcosmo, with a selectable potential (double_well | modular_log | modular_ratio); Haiku implements, an independent Haiku re-verifies, a Haiku anti-rigging review gates',
  whenToUse: 'Milestone M1 in ROADMAP.md. Run before l5-level-observable. Never touches proofs/ or main; one branch loop/level-explicit.',
  phases: [
    { title: 'Implement', detail: 'level= and potential= arguments, FRICKE_Y derived, tests, default output byte-identical' },
    { title: 'Verify', detail: 'independent agent re-runs every check from a clean state' },
    { title: 'Review', detail: 'anti-rigging: were the tests written to pass, or to detect?' },
  ],
}

const REPO = '/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator'
const WT = '/mnt/disks/disk-socrateai-local-1/dualscale-wt-level'
const BRANCH = 'loop/level-explicit'
const MAX_ATTEMPTS = (args && args.maxAttempts) || 2

const RULES = `
GROUND RULES:
1. Work ONLY in git worktree ${WT} on branch ${BRANCH}. If the worktree does not exist: cd ${REPO} && git worktree add ${WT} -b ${BRANCH} main (if the branch exists already, omit -b). Never touch proofs/, main, other worktrees; never push.
2. Evidence-bound: every claim you return carries the command you ran and a verbatim line of its output.
3. Stage only your own files by explicit path; check git diff --cached --name-only before committing. Never git add -A.
4. Commit message ends with: Co-Authored-By: Claude Haiku 4.5 <noreply@anthropic.com>
5. Do not weaken or delete an existing test to make it pass; if a test must change, say which and why in your return value.
`

const SPEC = `
WHAT TO BUILD (specs/LEANFLOW_ARCHITECTURE.md item L5 precondition; audit/K3_SELECTION.md section 6):
- workshopcosmo.py: add module-level DEFAULT_LEVEL = 12 and a function fricke_y(level) = 1/sqrt(level). FRICKE_Y stays as a name for backwards compatibility but MUST equal fricke_y(DEFAULT_LEVEL) and must no longer be a hand-typed float.
- compute_potential(x, y, a_pot, b_pot, level=DEFAULT_LEVEL, potential="double_well"): potential in {"double_well", "modular_log", "modular_ratio"}. "double_well" = today's formula with fricke_y(level) in place of FRICKE_Y. The modular ones call leanflow.core.modular_potential.modular_potential_xy(x mod 1, clamp(y, 1e-8, 60), N=level, mode=...) and return a central-difference gradient (h = 1e-7). Reduce x mod 1 (exact: T is in Gamma_0(N)) and clamp y (see scripts/l5_level_experiment.py, functions _bounded/potential -- reuse, do not re-derive).
- cosmology_rhs and run_quintessence_simulation take and forward level= and potential=; the returned dict gains "level" and "potential"; simulation_results.json and specs/SIMULATION_PROOF_REPORT.md print both.
- CLI: --level INT and --potential NAME on workshopcosmo.main(), defaults unchanged.
- Tests in tests/test_level_parameter.py:
  (a) DEFAULT BYTE-IDENTITY: with defaults, run_quintessence_simulation(t_max=20, num_points=100) equals the pre-change output to 1e-12 on every array (compute the reference ONCE from git stash / a checkout of main into a temp dir, save as tests/fixtures/quintessence_ref_t20.json, and compare).
  (b) LEVEL SENSITIVITY, modular_log: compute_potential(0.2, 1.2, level=7, potential="modular_log")[2] != compute_potential(0.2, 1.2, level=12, potential="modular_log")[2] -- the 2026-09-27 force scan gives 0.1116 vs 0.0747; assert the difference > 1e-3.
  (c) LEVEL INSENSITIVITY, double_well (documents the defect, does not hide it): for the double well, changing level changes ONLY fricke_y; assert that compute_potential at a point with cos(pi x)=0 (x=0.5) is identical for level 7 and 12, and that at x=0 it differs only through (y - fricke_y)^2. Name the test test_double_well_level_reaches_only_the_constant.
  (d) NEGATIVE CONTROL: potential="nope" raises ValueError.
- Do NOT change any default. Do NOT regenerate telemetry. Do NOT edit proofs/.
`

const IMPL_SCHEMA = {
  type: 'object',
  properties: {
    branch: { type: 'string' }, commit: { type: 'string' },
    files_changed: { type: 'array', items: { type: 'string' } },
    checks: { type: 'array', items: { type: 'object', properties: {
      name: { type: 'string' }, command: { type: 'string' }, output_line: { type: 'string' }, passed: { type: 'boolean' } },
      required: ['name', 'command', 'output_line', 'passed'] } },
    tests_changed_or_removed: { type: 'array', items: { type: 'string' } },
    notes: { type: 'string' },
  },
  required: ['branch', 'commit', 'files_changed', 'checks'],
}

const VERIFY_SCHEMA = {
  type: 'object',
  properties: {
    ok: { type: 'boolean' },
    checks: { type: 'array', items: { type: 'object', properties: {
      name: { type: 'string' }, command: { type: 'string' }, output_line: { type: 'string' }, passed: { type: 'boolean' } },
      required: ['name', 'command', 'output_line', 'passed'] } },
    refuted_claims: { type: 'array', items: { type: 'string' } },
    notes: { type: 'string' },
  },
  required: ['ok', 'checks', 'refuted_claims'],
}

const REVIEW_SCHEMA = {
  type: 'object',
  properties: {
    rigged: { type: 'boolean' },
    findings: { type: 'array', items: { type: 'object', properties: {
      file: { type: 'string' }, line: { type: 'integer' }, quote: { type: 'string' }, why: { type: 'string' } },
      required: ['file', 'quote', 'why'] } },
    must_fix: { type: 'array', items: { type: 'string' } },
  },
  required: ['rigged', 'findings', 'must_fix'],
}

let last = null
for (let attempt = 1; attempt <= MAX_ATTEMPTS; attempt++) {
  log(`attempt ${attempt}/${MAX_ATTEMPTS}`)
  phase('Implement')
  const feedback = last ? `\nPREVIOUS ATTEMPT WAS REJECTED. Verifier said: ${JSON.stringify(last.verify && last.verify.refuted_claims)}. Reviewer must_fix: ${JSON.stringify(last.review && last.review.must_fix)}. Fix these; do not argue with them.\n` : ''
  const impl = await agent(`${RULES}\n${SPEC}\n${feedback}
Implement it. Run: python3 -m pytest tests/test_level_parameter.py tests/test_workshopcosmo.py -q, and python3 workshopcosmo.py --quintessence --level 7 --potential modular_log (must not crash; report its printed final tau line). Commit on ${BRANCH}. Return the structured result; every check carries a verbatim output line.`,
    { label: `implement#${attempt}`, phase: 'Implement', schema: IMPL_SCHEMA, model: 'haiku' })
  if (!impl) { log('implementer died'); continue }

  phase('Verify')
  const verify = await agent(`${RULES}
You are the INDEPENDENT VERIFIER. You did not write this code. Another agent claims commit ${impl.commit} on ${BRANCH} in ${WT} satisfies this spec:
${SPEC}
Its claimed checks: ${JSON.stringify(impl.checks)}.
From a clean state (git -C ${WT} status must be clean; git -C ${WT} rev-parse HEAD must equal ${impl.commit}), re-run EVERY check yourself and add these of your own:
 - grep -n "FRICKE_Y *= *1.0 */ *math.sqrt(12" workshopcosmo.py must return nothing (the constant must be derived);
 - the fixture tests/fixtures/quintessence_ref_t20.json must have been generated from main, not from the new code: check its git log and regenerate it from a temp checkout of main (git -C ${REPO} worktree add /tmp/ref-main main; python3 -c ...), diff the two;
 - python3 workshopcosmo.py --quintessence (defaults) must print the same "Final tau" line as on main.
Try to REFUTE the claim. Default ok=false if any check cannot be reproduced. Clean up any temp worktree you add.`,
    { label: `verify#${attempt}`, phase: 'Verify', schema: VERIFY_SCHEMA, model: 'haiku' })

  phase('Review')
  const review = await agent(`${RULES}
You are the ANTI-RIGGING REVIEWER for commit ${impl.commit} on ${BRANCH} in ${WT}. Read the diff (git -C ${WT} diff main...${BRANCH}) and the new tests. Questions, each answered with a file:line quote:
 1. Does any test assert a value it computed from the code under test (a tautology)?
 2. Was any pre-existing test weakened, skipped or deleted? (git diff main...${BRANCH} -- tests/)
 3. Does test_double_well_level_reaches_only_the_constant actually document the defect (level reaches only fricke_y), or does it merely pass?
 4. Is any default changed, any telemetry regenerated, anything under proofs/ touched? (git diff --stat main...${BRANCH})
 5. Is x reduced mod 1 and y clamped BEFORE mpmath is called, for both modular modes?
Return rigged=true if any of 1, 2 or 4 is a yes. must_fix lists what the implementer must change.`,
    { label: `review#${attempt}`, phase: 'Review', schema: REVIEW_SCHEMA, model: 'haiku' })

  last = { impl, verify, review, attempt }
  const accepted = !!(verify && verify.ok && review && !review.rigged)
  log(`attempt ${attempt}: verify.ok=${verify && verify.ok} review.rigged=${review && review.rigged} -> ${accepted ? 'ACCEPTED' : 'rejected'}`)
  if (accepted) break
}

return {
  accepted: !!(last && last.verify && last.verify.ok && last.review && !last.review.rigged),
  branch: BRANCH, worktree: WT,
  commit: last && last.impl ? last.impl.commit : null,
  attempts: last ? last.attempt : 0,
  implementation: last && last.impl, verification: last && last.verify, review: last && last.review,
  next: 'If accepted: merge loop/level-explicit into main after reading the diff yourself; then run l5-level-observable. If not: read review.must_fix and verification.refuted_claims; do not merge.',
}
