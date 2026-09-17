export const meta = {
  name: 'dualscale-improve',
  description: 'Haiku implement → independent Haiku re-verify → Haiku anti-rigging review, one git branch per fix (numerics, tests, reproducibility, gates, manuscript honesty). Never touches proofs/ or main.',
  whenToUse: 'After physics-claims-audit / technical-audit. lean-foundation-bridge (new lean_foundation/ project on LeanMaster v4.33.1) runs only with allowFoundationBuild: true. Args: {only?: string[] item ids, allowFoundationBuild?: boolean, reviewModel?: "haiku"|"sonnet", maxAttempts?: number (default 2), batchSize?: number (default 3)}',
  phases: [
    { title: 'Implement', detail: 'Haiku in a git worktree, commits to branch haiku/improve-<id>', model: 'haiku' },
    { title: 'Verify', detail: 'different Haiku re-runs required checks on the branch (producer never verifies itself)', model: 'haiku' },
    { title: 'Review', detail: 'diff review: statement adequacy, no rigged PASS, no fabricated numbers, tier wording', model: 'haiku' },
  ],
}

const REPO = '/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator'
const LM = '/home/callensxavier_gmail_com/SocrateAI-Scientific-Agora-LeanMaster'
const LMPIN = '/mnt/disks/disk-socrateai-local-1/leanmaster-v2.2.0'
const A = args || {}
const REVIEW_MODEL = A.reviewModel || 'haiku'
const MAX_ATTEMPTS = A.maxAttempts || 2
const BATCH = A.batchSize || 3

const RULES = `
HARD RULES:
- Work only inside your git worktree (cwd). Never modify proofs/**, the ROOT lakefile.lean / lean-toolchain, or anything in ${LM} (read-only; another session is working there). Never push. Never touch branch main.
- Before any lake command: \`[ -e .lake ] || ln -s ${REPO}/.lake .lake\` (shares the built Mathlib; never run lake update or lake build of Mathlib). Rust: \`export CARGO_TARGET_DIR=${REPO}/rust_simulator/target\`.
- No fabricated numbers: every number you write into code comments, JSON, docs or LaTeX must come from a command you ran in this session (quote it in your report). Never hard-code a PASS/true status; statuses must be computed.
- Every new checker/test needs a positive control AND a negative control that fails on purpose.
- Do not weaken a test or assertion to make it pass. If blocked after 3 honest attempts, stop and report blocked=true with the first error.
- Tier wording (Mathesis): "proved in Lean" only for kernel-checked, standard-axioms-only, adequate statements; exact arithmetic = tier B; literature = L; physical interpretation = C ("we conjecture"/"if"); floats/simulations = X. Benchmarks without an execution log: "pending hardware verification".
- No git identity is configured and git config must not be changed: commit with \`git -c user.name="Xavier Callens" -c user.email=callensxavier@gmail.com commit ...\`.
- Commit on branch BRANCH with a conventional message ending with the line: Co-Authored-By: Claude Haiku 4.5 <noreply@anthropic.com>`

const ITEMS = [
  {
    id: 'num-cdl-shooting', kind: 'numerics',
    spec: `rust_simulator/src/vacuum_decay_cdl.rs — fix the bounce solver.
Known defects (verify each yourself first): (1) phi_true is hard-coded 2.05 but V'(phi)=m2*phi-3*kappa*phi^2+4*lambda*phi^3 with defaults m2=1,kappa=1.2,lambda=0.4 has non-zero roots at (3*kappa ± sqrt(9*kappa^2-16*lambda*m2))/(8*lambda) — compute them with python3 and confirm; (2) phi is clamped to >= phi_false, so the overshoot branch \`final_phi < phi_false\` can never fire and bisection only moves phi0_low; (3) euclidean_bounce_action uses .max(1.0), a floor that hides failure; (4) it is a flat-space Coleman bounce (no gravity) but is labelled Coleman-De Luccia.
Required: compute phi_true as the root with V''>0 (return an error/None if no second minimum); detect overshoot (phi crosses below phi_false) BEFORE any clamp and stop that trajectory; detect undershoot (dphi changes sign to positive while phi > phi_false); bisect properly; report converged: bool and the final |phi(rho_max) - phi_false|; remove the action floor; add doc comment "flat-space Coleman bounce; gravitational (CDL) corrections not included"; keep the central-charge fields but document c(N)=100(N+1) as a model assumption, not a derived c-theorem.
Tests (#[cfg(test)] in the same file): root has |V'|<1e-9 and V''>0; negative control: parameters with no second minimum yield error; bisection converges for defaults (converged==true); action is finite and > 0 without floor; a deliberately broken potential_deriv sign makes convergence fail (negative control).
Then regenerate: (cd rust_simulator && cargo run --release --bin vacuum_decay_cdl) and copy the regenerated vacuum_decay_cdl_summary.json to the repo root copy too. Do NOT edit manuscript numbers (another item does that).`,
    checks: ['cargo test --release', 'cargo run --release --bin vacuum_decay_cdl', 'python3 -m pytest -q'],
  },
  {
    id: 'num-swampland-integrator', kind: 'numerics',
    spec: `rust_simulator/src/swampland_geodesic.rs — make the geodesic flow honest.
Known defects (verify): the update labelled "Stiff symplectic / Velocity-Verlet" is an explicit Taylor step using start-of-step acceleration only; bound_satisfied uses a 1.5x slack factor; alpha_sdc = 1/sqrt(2) is hard-coded rather than measured; mode_mass_sq only includes the (n,m) lattice of one torus sector.
Required: replace the update with classical RK4 on the 8-dimensional first-order system; compute bound_satisfied strictly (slack configurable via a new config field bound_slack defaulting to 1.0, and reported in the summary); add alpha_effective = -ln(final_mass_gap/initial_mass_gap)/total_geodesic_distance to the summary (document as tier X numerical estimate); document the truncated mode set.
Tests: hyperbolic-metric speed sqrt((u1^2+u2^2)/y^2) is conserved along the flow to < 1e-6 relative drift at ds=0.01 over s_max=5 (RK4); negative control: the old explicit scheme (keep as a private test helper) drifts more than 10x the RK4 drift; alpha_effective is finite; strict bound check reports false for a synthetic record list violating the bound (negative control).
Regenerate swampland_geodesic_summary.json (cargo run --release --bin swampland_distance) in rust_simulator/ and the repo-root copy.`,
    checks: ['cargo test --release', 'cargo run --release --bin swampland_distance', 'python3 -m pytest -q'],
  },
  {
    id: 'num-stochastic-seeding', kind: 'numerics',
    spec: `rust_simulator/src/tachyon_condensation.rs and rust_simulator/src/kummer_langevin.rs — reproducibility and computed statuses.
Required: add a seed: u64 config field (default 42) and use rand::rngs::StdRng::seed_from_u64 everywhere randomness is drawn (replace thread_rng); record the seed in the summary JSON. Find how k_theory_charge_conserved, rr_charge_match, grothendieck_defect_rank and symmetry_broken are set: if any is a literal true/constant, replace it with a computed quantity — e.g. net kink charge Q = (sign(phi(+L)) - sign(phi(-L)))/2 evaluated at t=0 and t=final, conserved iff equal — and document in a doc comment what is computed (tier X numerical diagnostic, not K-theory).
Tests: same seed → bit-identical summaries; different seed → at least one field differs (negative control); the computed charge diagnostic returns false on a hand-constructed field whose boundary signs change (negative control).
Regenerate tachyon_condensation_summary.json and kummer_langevin_summary.json via the rust_simulator binaries (see src/main.rs / src/bin) and update root copies.`,
    checks: ['cargo test --release', 'python3 -m pytest -q'],
  },
  {
    id: 'repro-paths-sundials', kind: 'reproducibility',
    spec: `workshopcosmo.py (and scripts/*.py) — remove machine-specific paths and phantom integrations.
Required: remove every hard-coded /home/xavkal/... path and ../SocrateAI-Lean-Lib path (grep -rn first; list them in your report). rusty-SUNDIALS CSV loading must use an env var RUSTY_SUNDIALS_DIR; when unset or file missing, results must say "sundials_integrated": false and "sundials_status": "not run: RUSTY_SUNDIALS_DIR unset" and "independent_engines" must list only engines that actually ran in that invocation (computed, not a literal list). Same honest-status rule for any other external engine.
Tests: add tests in tests/test_workshopcosmo.py: with RUSTY_SUNDIALS_DIR unset, sundials_integrated is False and "rusty-SUNDIALS" not in independent_engines; with it pointing to a tmp dir containing a minimal fixture CSV, integrated is True (positive control). grep -rn "xavkal" returns nothing outside git history.
Also: running pytest currently OVERWRITES committed results (observed 2026-09-17: tachyon_condensation_summary.json float drift ~1e-15, tda_mapper_skeleton.json ~2.5k changed lines, tda_mapper_graph.png) — make tests write outputs to tmp_path (pytest fixture) or an OUTPUT_DIR env var, and seed the TDA Mapper (clustering/lens randomness) so skeleton output is deterministic; add a test that runs the Mapper twice with the same seed and asserts identical node/edge counts. Verify with: git status --porcelain after pytest shows no modified tracked files.`,
    checks: ['python3 -m pytest -q', 'grep -rn xavkal --include=*.py . ; test $? -eq 1', 'test -z "$(git status --porcelain --untracked-files=no)"'],
  },
  {
    id: 'gate-lean-axiom-audit', kind: 'gate',
    spec: `Add an enforced Lean axiom gate WITHOUT modifying proofs/ and WITHOUT copying LeanMaster code. Load the \`lean-proof-gate\` skill first (Skill tool) if available.
Create scripts/lean_gate.sh that: builds a temporary symlinked root (lakefile.lean, lean-toolchain, lake-manifest.json, .lake, every proofs/*.lean flattened to the root, proofs/LeanscratchDB/*.lean under LeanscratchDB/) — this is required because ${LM}/tools/axiom_audit.py derives module names from paths relative to LEAN_PROJECT_ROOT and this repo uses srcDir "proofs"; runs \`LEAN_PROJECT_ROOT=<tmp> python3 ${LM}/tools/axiom_audit.py <module files>\` over all library modules (exclude ExportAxioms.lean and non-target FTheoryCosmology.lean, but list exclusions in the output); FAILS (exit 2) if the tool prints AUDIT ERROR or if the number audited is 0 (the tool exits 0 in that case — a silent false pass); parses OK/FAIL lines into audit/lean_axiom_report.json {module, theorem, axioms, status}; exits 1 if any FAIL.
Wire workshopcosmo.py: zero_sorry / lean_certification_verified / lean_voa_verified must come from audit/lean_axiom_report.json for the relevant module (true only if every theorem in that module is OK); add "axiom_footprint"; add "statement_adequacy_reviewed": false (a passing axiom audit is not adequacy — never set it true in code).
Tests (tests/test_lean_gate.py, skip if lake missing): positive control — a temp root with one extra module \`theorem t : 1 + 1 = 2 := rfl\` → OK; negative control — a module declaring \`axiom bad : False\` and \`theorem u : False := bad\` → FAIL and exit 1; an empty module list → exit 2 (guards the silent-pass pitfall).
Do not hard-code today's counts anywhere.`,
    checks: ['bash scripts/lean_gate.sh; test $? -le 1', 'python3 -m pytest -q tests/test_lean_gate.py', 'python3 -m pytest -q'],
  },
  {
    id: 'lean-foundation-bridge', kind: 'lean', requiresArg: 'allowFoundationBuild',
    spec: `Create a NEW downstream Lean project lean_foundation/ that builds on LeanMaster's kernel-verified libraries and restates, non-vacuously, the content that proofs/*.lean only pretended to prove. Do not modify proofs/ or the root lakefile. Load skills \`leanmaster-onboard\`, \`string-theory-foundation\`, \`leanmaster-theorem-search\`, \`lean-tiered-proving\`, \`lean-proof-gate\` (Skill tool) first.
Setup: copy ${LM}/examples/consumer_demo/{lakefile.lean,lean-toolchain} into lean_foundation/ (toolchain leanprover/lean4:v4.33.1). Rename package to «dualscale-foundation», lib «DualScaleFoundation». Depend on the PINNED, already-built clone of LeanMaster at release tag v2.2.0 — NOT the live ${LM} checkout (another session edits it): \`require «SocrateAI-Scientific-Agora-LeanMaster» from "${LMPIN}"\` and keep \`packagesDir := "/mnt/disks/disk-socrateai-local-1/leanmaster/lake/packages"\`. Before building verify \`git -C ${LMPIN} describe --tags\` prints v2.2.0 and \`tail -1 ${LMPIN}/build_v2.2.0.log\` shows EXIT 0; if not, stop with blocked=true. Never run lake update. Import specific modules only (e.g. import DualScaleStream2.DFT.GeneralizedMetric), not whole roots. In docs cite release v2.2.0 (certificate docs/VERIFIED_FOUNDATION.md pinned to v2.1.0; tree differs only in comments per its author — say so).
Modules (one per replaced local file; header comment "-- Supersedes proofs/<file>.lean (see audit/lean_replacement_map.md)"):
  DualScaleFoundation/TDuality.lean — Buscher/T-duality content via DualScaleStream2.TDuality.ODD and DFT.GeneralizedMetric (e.g. genMetric, etaR_genMetric_sq, tduality_inverts_metric, isODD_mul, thetaShift_isODD). Add an \`example\` showing the old claim "T-duality applied twice is the identity" in the form LeanMaster supports, or document why it is not stated.
  DualScaleFoundation/K3Lattice.lean — Mukai/K3 via DualScaleStream2.Lattice (mukaiPair_*, sigK3_eq, sigMukai_eq, reflection_isometry).
  DualScaleFoundation/Tadpole.lean — via DualScaleStream2.Flux (k3k3_anomaly, tadpole_budget, flux_half_selfIntersection_integral). Do NOT restate the O7/T^4/Z2 charge bookkeeping — it is unresolved (T0 question).
  DualScaleFoundation/Moonshine.lean — via DualScaleStream2.Moonshine.EOT (eotA; massive multiplicities are 2*A_n per the skill's convention). State the 77/60 ratio as an exact arithmetic identity derived FROM eotA values (not re-typed constants), with a docstring: tier A arithmetic; its reading as a bispectrum ratio is tier C.
Every new theorem: statement must mention at least one LeanMaster declaration or be a genuine consequence of one — never a def written to make it true. Use the lean-tiered-proving order (exact?/simp/norm_num/omega/decide → Haiku reasoning; three strikes → leave \`sorry\` and list it as an open goal; never weaken a statement silently).
Gates (all must be run and reported): \`cd lean_foundation && lake build\`; \`grep -rn "sorry" DualScaleFoundation | grep -v "^.*--"\` count; \`LEAN_PROJECT_ROOT=$PWD python3 ${LM}/tools/axiom_audit.py DualScaleFoundation\` (library layout here is standard, so this works — but also assert the "N theorems audited" line has N>0). Do NOT run statement_lock --update (statement review G2 belongs to the project lead): instead write lean_foundation/STATEMENTS_FOR_REVIEW.md listing each theorem's verbatim #check output, the LeanMaster declarations it rests on, its tier, and which proofs/*.lean theorem it supersedes. Add lean_foundation/.lake to .gitignore coverage if not already matched.`,
    checks: ['lake build', 'axiom_audit.py DualScaleFoundation', 'STATEMENTS_FOR_REVIEW'],
  },
  {
    id: 'gate-cross-consistency', kind: 'gate',
    spec: `Add a Rule-8 cross-consistency checker.
Required: scripts/cross_consistency_check.py with a register audit/parameters.json: for each quantity (mapper nodes/edges/beta_1; bounce action; bubble radius; SDC final mass gap; soliton width; string count; symmetron screening factor; |gamma-1|; w0; wa; Delta chi2; ln Bayes factor; Lean theorem count), a canonical source (JSON file + key path, or the output of scripts/lean_axiom_audit.py if present) and a list of occurrence locators (file + regex with one capture group) found by grep in papers/T-dulaity alone/T_duality_Alone.tex, zenodo_bundle/T_duality_Alone.tex, submission cover letters, *_summary.json copies and simulation_results.json. Exact comparison for integers; relative tolerance 1e-9 for floats unless the manuscript rounds (then compare at the printed precision). Output a markdown table audit/cross_consistency.md and exit 1 on any mismatch. The register itself must be built from grep output, not memory — include in your report the grep commands used.
Tests (tests/test_cross_consistency.py) on tmp fixtures: consistent fixture → exit 0 (positive control); one mutated number → exit 1 naming the quantity (negative control). Running on the real repo is expected to FAIL today — record that output in your report; do not "fix" numbers here.`,
    checks: ['python3 -m pytest -q tests/test_cross_consistency.py', 'python3 scripts/cross_consistency_check.py; test $? -le 1'],
  },
  {
    id: 'ci-pipeline', kind: 'gate',
    spec: `Add .github/workflows/ci.yml (new file, branch only — not pushed): jobs (1) rust: cargo test --release in rust_simulator; (2) python: pip install -e .[dev] numpy scipy scikit-learn matplotlib, pytest -q; (3) lean: install elan, use lean-toolchain, cache .lake keyed on lake-manifest.json + lean-toolchain, lake build, then python3 scripts/lean_axiom_audit.py in report-only mode (continue-on-error: true, upload audit/lean_axiom_report.json as artifact) — the gate is made blocking later by a T0 decision; (4) consistency: python3 scripts/cross_consistency_check.py with continue-on-error: true, upload audit/cross_consistency.md. Reference scripts may not exist on main yet: guard with if [ -f ... ]. Validate YAML with python3 -c "import yaml,sys; yaml.safe_load(open('.github/workflows/ci.yml'))" (pip install pyyaml --user if needed). Also add pytest-cov to [project.optional-dependencies].dev in pyproject.toml.`,
    checks: ["python3 -c \"import yaml; yaml.safe_load(open('.github/workflows/ci.yml'))\"", 'python3 -m pytest -q'],
  },
  {
    id: 'manuscript-tier-honesty', kind: 'manuscript',
    spec: `papers/T-dulaity alone/T_duality_Alone.tex AND the public deposit metadata .zenodo.json (root; its "description" repeats the same claims: Lean proofs of Fourier–Mukai/Kummer, 187/557/376, 12–18 ns) — tier-honest WORDING pass only (sentence-level edits; no new results). Do NOT change any simulation number in this item: numbers are reconciled after the numerics branches merge, by scripts/cross_consistency_check.py; where a number is known to be inconsistent or unsourced, wrap it as \\textcolor{red}{[NUM-RECONCILE: <quantity>]} next to the current value (LaTeX) or append " [to be reconciled]" (JSON description) so the post-merge pass can find it with grep NUM-RECONCILE.
If ${REPO}/audit/physics_claims_ledger.md exists, use its "Confirmed defects" rows touching the manuscript as the edit list. Otherwise apply exactly these, after verifying each against the repo yourself:
(a) Abstract/intro/conclusion: replace "formally verify/certify/kernel-certified/down to the proof kernel/zero sorry axioms" for anything not tier A. Verified facts you must re-check with commands: proofs/BuscherRules.lean declares 5 global axioms, and \`inv_inv\` instantiated at Nat with x=2 derives False (reproduce: scratch file importing BuscherRules with \`theorem t : False := by have h := @SocrateAI.StringTheory.Buscher.inv_inv Nat _ _ 2; simp at h\`, run \`lake env lean\`); LeanscratchDB.HoloAlg gravitino_anomaly_cancellation is an axiom restated; SwamplandDistance/FluxVacuumDecay/TachyonCondensation/KummerTDA/DbraneBoundaryStates theorems hold by construction of their definitions. Disclose these in a new subsection "Scope, tiers and known limitations" (tier table A/B/L/C/X per headline claim).
(b) R_BPS = 77/60: "exact arithmetic identity on EOT multiplicities (tier B)"; any reading as a bispectrum / non-Gaussianity ratio → "we conjecture" (tier C).
(c) Benchmarks (1,520x, 7.1x, $ costs, 99.1%): no execution log exists in the repo (verify with grep -rn) → replace with "pending hardware verification" and remove the numbers from abstract/conclusion; keep a sentence saying benchmark scripts will be released.
(d) "Coleman-De Luccia" → "flat-space Coleman bounce (gravitational CDL corrections not included)" where it describes rust_simulator/src/vacuum_decay_cdl.rs.
(e) Orientifold charge bookkeeping: do NOT state which O-plane dimension or charges are correct (no physics from memory). Record only the internal inconsistency you can quote: proofs/TadpoleCancellation.lean and KummerTDAAnomalyCertification.lean cite Gimon-Polchinski (hep-th/9601038) for the 16 fixed points of T^4/Z2 while assigning O7- planes of charge -4 and 32 D7-branes of charge +2 — quote those lines, state in the limitations subsection that the construction cited and the O-plane type assigned have not been reconciled against the cited source (tier C, under revision), and add a t0_question asking the project lead to rule on the correct construction with a quoted source.
(e2) Where the manuscript needs a kernel-checked statement (T-duality/O(d,d), generalized metric, K3/Mukai signature, tadpole arithmetic, EOT multiplicities), cite LeanMaster instead of proofs/*.lean, using the exact form from the \`string-theory-foundation\` skill ("Kernel-checked in Lean 4 (v4.33.1, Mathlib v4.33.1) in DualScaleStream2.<Module>.<name>, repository xaviercallens/SocrateAI-Scientific-Agora-LeanMaster, release v2.1.0; depends only on Lean's three standard axioms"), and only for declarations you have seen verbatim in ${LM}/docs/VERIFIED_FOUNDATION.md §3. Carry its scope limits (one-direction results; dualScale_one is attainment not uniqueness; physical interpretation is tier L/C). Never write "zero axioms", "100% verified", "proves string theory".
(f) Mapper counts (187/557/376), bounce action, bubble radius and every other simulation number: tag with NUM-RECONCILE as above; do not replace.
Apply the same edits to papers/T-dulaity alone/submission/COVER_LETTER_*.tex where the same claims appear. Leave zenodo_bundle/ untouched (it is regenerated from the corrected sources at release). Also keep .zenodo.json valid JSON (python3 -m json.tool).
Checks: forbidden-phrase count before vs after (grep -ciE "formally verif|certif|kernel-certif|zero .sorry|1,520|7\\.1" on the .tex) must strictly decrease; LaTeX sanity: \\begin/\\end counts balanced (python3 script), and if pdflatex exists run it twice in a tmp copy of the directory and report errors.`,
    checks: ['python3 -m pytest -q', 'json.tool .zenodo.json'],
  },
]

const IMPL = {
  type: 'object',
  required: ['branch', 'commit', 'blocked', 'summary', 'checks', 'negative_controls', 'numbers_written'],
  properties: {
    branch: { type: 'string' },
    commit: { type: 'string', description: 'full SHA of the final commit on the branch' },
    blocked: { type: 'boolean' },
    blocked_reason: { type: 'string' },
    summary: { type: 'string' },
    files_changed: { type: 'array', items: { type: 'string' } },
    checks: { type: 'array', items: { type: 'object', required: ['command', 'exit_code', 'output_tail'], properties: { command: { type: 'string' }, exit_code: { type: 'integer' }, output_tail: { type: 'string' } } } },
    negative_controls: { type: 'array', items: { type: 'string' }, description: 'test names that fail on purpose and how' },
    numbers_written: { type: 'array', items: { type: 'string' }, description: 'each number written + the command output it came from' },
    t0_questions: { type: 'array', items: { type: 'string' }, description: 'decisions only the project lead should make' },
  },
}

const VERIFY = {
  type: 'object',
  required: ['pass', 'results'],
  properties: {
    pass: { type: 'boolean' },
    results: { type: 'array', items: { type: 'object', required: ['command', 'exit_code', 'output_tail'], properties: { command: { type: 'string' }, exit_code: { type: 'integer' }, output_tail: { type: 'string' } } } },
    negative_controls_really_fail: { type: 'boolean', description: 'you temporarily broke the code/fixture and saw the named negative-control test fail, then restored' },
    notes: { type: 'string' },
  },
}

const REVIEW = {
  type: 'object',
  required: ['approve', 'problems'],
  properties: {
    approve: { type: 'boolean' },
    problems: { type: 'array', items: { type: 'object', required: ['severity', 'where', 'issue'], properties: { severity: { type: 'string', enum: ['blocking', 'major', 'minor'] }, where: { type: 'string' }, issue: { type: 'string' } } } },
    t0_questions: { type: 'array', items: { type: 'string' } },
  },
}

const gateOk = (impl, item) => {
  if (!impl || impl.blocked || !impl.commit) return { ok: false, why: impl ? (impl.blocked_reason || 'no commit') : 'implementer died' }
  const failing = impl.checks.filter(c => c.exit_code !== 0)
  if (failing.length) return { ok: false, why: `failing checks: ${failing.map(c => c.command).join('; ')}` }
  const token = req => (req.match(/cargo (test|run)|pytest|lean_gate|cross_consistency_check|xavkal|git status|yaml|json.tool|lake build|axiom_audit|STATEMENTS_FOR_REVIEW/) || [req.split(' ')[0]])[0]
  const missing = item.checks.filter(req => !impl.checks.some(c => c.command.includes(token(req))))
  if (missing.length) return { ok: false, why: `required checks not run: ${missing.join('; ')}` }
  if (!impl.negative_controls.length && item.kind !== 'manuscript') return { ok: false, why: 'no negative control' }
  return { ok: true }
}

async function attempt(item, n, prev) {
  const branch = `haiku/improve-${item.id}${n > 1 ? `-r${n}` : ''}`
  const base = prev ? `Start from the previous attempt: \`git checkout -b ${branch} ${prev.commit}\`. Previous attempt problems to fix:\n${prev.feedback}` : `Create the branch: \`git checkout -b ${branch}\` (from main).`
  const impl = await agent(`Implement improvement "${item.id}" in the DualScaleSimulator repo (your cwd is an isolated git worktree of ${REPO}).
${base}
SPEC:
${item.spec}
Required checks to run at the end from the worktree root (report each with exit code and last ~30 lines): ${item.checks.map(c => `\`${c}\``).join(', ')}.
${RULES.replace('BRANCH', branch)}`, { label: `impl:${item.id}#${n}`, phase: 'Implement', schema: IMPL, model: 'haiku', isolation: 'worktree' })

  const g = gateOk(impl, item)
  if (!g.ok) return { item: item.id, attempt: n, branch, impl, stage: 'gate', ok: false, feedback: g.why }

  const [ver, rev] = await parallel([
    () => agent(`Independent verifier (you did not write this code). In your cwd (a fresh worktree of ${REPO}) run: \`git checkout --detach ${impl.commit}\`, then \`[ -e .lake ] || ln -s ${REPO}/.lake .lake\` and \`export CARGO_TARGET_DIR=${REPO}/rust_simulator/target\`.
Re-run exactly: ${item.checks.map(c => `\`${c}\``).join(', ')}. pass=true only if all exit 0 (or satisfy their own \`test\` clause).
Then pick ONE negative control claimed by the implementer (${JSON.stringify(impl.negative_controls)}), temporarily break the code/fixture it guards, confirm that test fails, restore with git checkout -- ., and set negative_controls_really_fail accordingly. Do not commit anything.`,
      { label: `verify:${item.id}#${n}`, phase: 'Verify', schema: VERIFY, model: 'haiku', isolation: 'worktree' }),
    () => agent(`Anti-rigging review of branch commit ${impl.commit} for "${item.id}". Run: \`git -C ${REPO} diff main...${impl.commit}\` and \`git -C ${REPO} show --stat ${impl.commit}\` (read-only; do not check out or modify anything).
Spec it should satisfy:
${item.spec}
Implementer's number provenance: ${JSON.stringify(impl.numbers_written)}
Blocking problems to look for: hard-coded PASS/true statuses or observables; numbers without provenance; assertions weakened, tests skipped or tolerance inflated to pass; negative controls that cannot fail; changes to proofs/**, lakefile.lean, lean-toolchain, main, or ${LM}; physics wording above its tier (e.g. "proves", "certifies", "CDL" for a flat-space bounce); scope creep beyond the spec. Major: missing doc of what a diagnostic actually computes; regenerated JSON not matching code output.
approve=false if any blocking problem. Put decisions that belong to the project lead in t0_questions.`,
      { label: `review:${item.id}#${n}`, phase: 'Review', schema: REVIEW, model: REVIEW_MODEL }),
  ])

  const verOk = ver && ver.pass && (item.kind === 'manuscript' || ver.negative_controls_really_fail !== false)
  const revOk = rev && rev.approve
  const feedback = [
    !verOk ? `Verifier: ${ver ? ver.notes || JSON.stringify(ver.results.filter(r => r.exit_code !== 0)) : 'verifier died'}${ver && ver.negative_controls_really_fail === false ? ' | negative control did not fail when code was broken' : ''}` : '',
    !revOk ? `Reviewer: ${rev ? rev.problems.map(p => `[${p.severity}] ${p.where}: ${p.issue}`).join(' | ') : 'reviewer died'}` : '',
  ].filter(Boolean).join('\n')
  return { item: item.id, attempt: n, branch, commit: impl.commit, impl, verify: ver, review: rev, ok: verOk && revOk, stage: 'final', feedback }
}

async function runItem(item) {
  let prev = null, last = null
  for (let n = 1; n <= MAX_ATTEMPTS; n++) {
    last = await attempt(item, n, prev)
    if (last.ok) break
    if (!last.impl || !last.impl.commit) { prev = null } else { prev = { commit: last.impl.commit, feedback: last.feedback } }
    log(`${item.id}: attempt ${n} not accepted (${last.stage}) — ${String(last.feedback).slice(0, 160)}`)
  }
  const t0 = [...((last.impl && last.impl.t0_questions) || []), ...((last.review && last.review.t0_questions) || [])]
  return {
    id: item.id, kind: item.kind, accepted: last.ok, branch: last.branch, commit: last.commit || null, attempts: last.attempt,
    status: last.ok ? 'ready_for_T0_review' : (last.stage === 'gate' ? 'blocked_at_gate' : 'rejected_escalate_to_T1'),
    summary: last.impl && last.impl.summary, feedback: last.ok ? '' : last.feedback, t0_questions: t0,
  }
}

const gated = ITEMS.filter(i => i.requiresArg && !A[i.requiresArg])
if (gated.length) log(`Skipped (needs args.${gated.map(i => i.requiresArg).join('/')}=true): ${gated.map(i => i.id).join(', ')}`)
const eligible = ITEMS.filter(i => !i.requiresArg || A[i.requiresArg])
const items = A.only ? eligible.filter(i => A.only.includes(i.id)) : eligible
if (A.only) log(`Restricted to: ${items.map(i => i.id).join(', ')}`)
log(`Running ${items.length} item(s) in batches of ${BATCH} (compile load cap), up to ${MAX_ATTEMPTS} attempts each; review model: ${REVIEW_MODEL}`)

const LEAN_TOUCHING = new Set(['gate-lean-axiom-audit', 'lean-foundation-bridge', 'manuscript-tier-honesty'])
const batches = []
const light = items.filter(i => !LEAN_TOUCHING.has(i.id)), heavy = items.filter(i => LEAN_TOUCHING.has(i.id))
for (let i = 0; i < light.length; i += BATCH) batches.push(light.slice(i, i + BATCH))
for (const h of heavy) batches.push([h]) // serialise lake env lean on the shared .lake
const results = []
for (let bi = 0; bi < batches.length; bi++) {
  const batch = batches[bi]
  const out = await pipeline(batch, item => runItem(item))
  results.push(...out.map((r, k) => r || { id: batch[k].id, accepted: false, status: 'crashed' }))
  log(`Batch ${bi + 1}/${batches.length}: ${out.filter(r => r && r.accepted).length}/${batch.length} accepted`)
}

return {
  accepted: results.filter(r => r.accepted).map(r => `${r.id} → ${r.branch} @ ${r.commit}`),
  not_accepted: results.filter(r => !r.accepted).map(r => `${r.id}: ${r.status}`),
  t0_questions: results.flatMap(r => r.t0_questions || []),
  results,
  note: 'Nothing merged. Review branches with: git log --oneline main..<branch>; git diff main...<branch>',
}
