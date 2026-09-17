export const meta = {
  name: 'technical-audit',
  description: 'Haiku fleet: evidence-bound engineering audit (build, tests, reproducibility, provenance, duplication, security), re-checks the earlier ungrounded Haiku report, writes audit/technical_audit.md',
  whenToUse: 'Periodic engineering health check, or before a release/Zenodo deposit. Caller writes returned `markdown` → audit/technical_audit.md. Args: {only?: string[] dimension ids}',
  phases: [
    { title: 'Probe', detail: 'one Haiku agent per dimension, must run commands', model: 'haiku' },
    { title: 'Verify', detail: '2 Haiku skeptics per finding: re-run evidence command, refute', model: 'haiku' },
  ],
}

const REPO = '/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator'
const LM = '/home/callensxavier_gmail_com/SocrateAI-Scientific-Agora-LeanMaster'
const A = args || {}

const GROUNDING = `
GROUNDING RULES (mandatory):
- You MUST run shell commands to gather evidence. A previous audit of this repo made zero tool calls and invented findings (e.g. "shell injection in publish_to_zenodo.py", "tokens in .zenodo.json", "zero type hints") — do not repeat that failure.
- Every finding must carry: evidence_command (exact command that demonstrates it, runnable from ${REPO}), evidence_output (short verbatim excerpt of its output), and file/line where applicable (line 0 if repo-wide).
- Report only what the command output shows. Also record notable things that are GOOD (kind "strength").
Repo root: ${REPO}. Toolchains available: lake/lean (elan), cargo, python3, pytest.`

const DIMENSIONS = [
  { id: 'build-lean', prompt: `Lean gate health using the reusable LeanMaster tools (skill \`lean-proof-gate\`; load it with the Skill tool if available). Do NOT run \`lake build\` and never write inside ${'${LM}'}.
Known tool pitfalls on THIS repo's layout (verified 2026-09-17): libraries use srcDir "proofs" with single-file modules, so \`LEAN_PROJECT_ROOT=${'${REPO}'} python3 ${'${LM}'}/tools/axiom_audit.py BuscherRules\` audits 0 theorems and exits 0 (a silent false pass), and passing proofs/X.lean fails with "unknown module prefix 'proofs'". Use the symlinked-root recipe instead:
  R=$(mktemp -d); for f in lakefile.lean lean-toolchain lake-manifest.json .lake; do ln -s ${'${REPO}'}/$f $R/$f; done; for f in ${'${REPO}'}/proofs/*.lean; do ln -s $f $R/; done; mkdir $R/LeanscratchDB; ln -s ${'${REPO}'}/proofs/LeanscratchDB/*.lean $R/LeanscratchDB/; cd $R && LEAN_PROJECT_ROOT=$R timeout 900 python3 ${'${LM}'}/tools/axiom_audit.py $(cd ${'${REPO}'}/proofs && ls *.lean LeanscratchDB/*.lean | grep -v -e ExportAxioms -e FTheoryCosmology)
Report: number audited and failing (a baseline run on 2026-09-17 gave 83 audited / 5 failing — re-measure, do not copy), each FAIL with its non-standard axioms, and a finding that "0 theorems audited" must be treated as failure. Then state explicitly that passing the axiom audit is NOT statement adequacy: sample 5 OK theorems and check whether they are true by definition. Also: lean-toolchain (${'${REPO}'} v4.34.0-rc2, Mathlib unpinned) vs LeanMaster's required v4.33.1; which proofs/*.lean are not lakefile targets; whether this repo has any enforced gate (script/test/CI); how workshopcosmo.py computes zero_sorry / lean_certification_verified.` },
  { id: 'build-rust', prompt: `Rust simulator health. cd rust_simulator && cargo build --release 2>&1 | tail; cargo clippy --release 2>&1 | tail -40 (skip if clippy missing); count #[test] / #[cfg(test)] (grep -rn). Check use of unwrap/expect/panics, hard-coded output paths, whether binaries write JSON into the CWD (duplicated *_summary.json at repo root, rust_simulator/, zenodo_bundle/ — diff them and report mismatches).` },
  { id: 'tests-python', prompt: `Python tests. Run: timeout 900 python3 -m pytest -q 2>&1 | tail -20. Then read tests/*.py and classify what the assertions check: physical invariants with tolerances, vs only "returns dict/status PASS", vs skipped when external tools missing. List tests that pass trivially or are skipped and why. Check pytest.ini. Check whether tests compile Lean or merely grep for "sorry".` },
  { id: 'reproducibility', prompt: `Reproducibility & provenance. grep -rn for absolute paths (/home/, xavkal) and sibling-repo paths (../rusty-SUNDIALS, ../SocrateAI-Lean-Lib) in *.py *.rs *.sh *.tex. For each JSON/CSV result file committed at repo root, find which script produces it (grep the filename) and whether its numbers match the manuscript (grep numbers in papers/T-dulaity alone/T_duality_Alone.tex). Is there a random seed set in the Langevin/MCMC code (grep seed)? Are hard-coded timestamps in results? Are there benchmark harnesses for the 1,520x / 7.1x claims (grep -rn speedup, timeit, criterion, benchmark)? Are large binaries committed (git ls-files | xargs du -ch | sort -h | tail; zenodo_deposit_bundle.zip, *.tar.gz, *.pt, *.pdf, 25k-line CSV)?` },
  { id: 'duplication-packaging', prompt: `Duplication & packaging. diff -rq leanflow deployment/huggingface/space/leanflow; zenodo_bundle/T_duality_Alone.tex vs papers/T-dulaity alone/T_duality_Alone.tex (diff | head); tarballs in zenodo_bundle vs sources (tar tzf | head). pyproject.toml references README.md — does it exist? Does pip install -e . --dry-run or python3 -m build work (use a temp venv only if fast; otherwise just inspect)? LaTeX build artefacts (*.aux, *.out, *.toc) tracked in git? Directory names with typos/spaces ("T-dulaity alone", "Simulatot Phase 2.md").` },
  { id: 'code-quality-security', prompt: `Code quality & security, grounded. grep -rn "shell=True\\|os.system\\|eval(\\|exec(\\|pickle.load\\|torch.load" --include=*.py; inspect scripts/publish_to_zenodo.py and scripts/deploy_leanflow_serverless.sh for how tokens are read (env var vs file vs hard-coded) — quote the lines; check .zenodo.json content for secrets; check leanflow/bridge/lean_ipc.py for subprocess usage and exception handling (grep -n "except"). Measure type-hint coverage roughly: python3 - <<'EOF' using ast to count functions with/without annotations in leanflow/ and workshopcosmo.py EOF. Report bare excepts with line numbers.` },
  { id: 'cross-consistency', prompt: `Rule-8 cross-consistency gate: every numeric value must be identical across code, JSON, LaTeX, Lean and README. Build a table of the key quantities — Mapper nodes/edges/beta_1, bounce action S_E, bubble radius, final mass gap, SDC alpha, soliton width, string counts, symmetron screening factor, |gamma-1|, w0/wa, Delta chi2, Bayes factor, speedups/costs, theorem counts ("14 Lean 4 formal proofs", theorems_verified_count) — and for each list every occurrence with file:line (grep -rn across *.json *.tex *.py *.rs *.lean *.md, excluding .lake/ and target/) and whether values agree. Also diff the duplicated summary JSONs (repo root vs rust_simulator/ vs zenodo_bundle/) and papers/…/T_duality_Alone.tex vs zenodo_bundle/T_duality_Alone.tex. One finding per mismatching quantity; one "strength" finding per quantity that is fully consistent.` },
  { id: 'prior-report-recheck', prompt: `An earlier audit report (produced without reading any files) made these claims. For EACH, run commands to CONFIRM or REFUTE and emit one finding with kind "prior_claim_confirmed" or "prior_claim_refuted":
1. No top-level README.md although pyproject.toml references it.
2. Zero type hints across leanflow/.
3. Bare "except:" in leanflow/bridge/lean_ipc.py.
4. Potential shell injection in scripts/publish_to_zenodo.py.
5. .zenodo.json likely contains API tokens.
6. deployment/huggingface/space/leanflow duplicates leanflow/.
7. No CI/CD (.github/workflows missing).
8. Rust simulator has no tests.
9. No pytest-cov / coverage config.
10. Tests status unknown (run them).
11. papers/"T-dulaity alone" directory name typo.
12. HuggingFace app lacks auth/CORS (check whether that is even applicable to a Gradio Space).
Emit exactly 12 findings, one per claim, even when the answer is "not applicable" (use kind prior_claim_refuted with the reason).` },
]

const FINDINGS = {
  type: 'object',
  required: ['findings', 'commands_run'],
  properties: {
    commands_run: { type: 'array', items: { type: 'string' } },
    findings: {
      type: 'array',
      items: {
        type: 'object',
        required: ['kind', 'severity', 'title', 'file', 'line', 'evidence_command', 'evidence_output', 'recommendation'],
        properties: {
          kind: { type: 'string', enum: ['defect', 'risk', 'strength', 'prior_claim_confirmed', 'prior_claim_refuted'] },
          severity: { type: 'string', enum: ['critical', 'high', 'medium', 'low', 'info'] },
          title: { type: 'string' },
          file: { type: 'string' },
          line: { type: 'integer' },
          evidence_command: { type: 'string' },
          evidence_output: { type: 'string' },
          recommendation: { type: 'string' },
        },
      },
    },
  },
}

const VERDICT = {
  type: 'object',
  required: ['refuted', 'evidence'],
  properties: {
    refuted: { type: 'boolean', description: 'the FACT in the finding is false (not: the recommendation is unwise)' },
    recommendation_problem: { type: 'string', description: 'if the fact is true but the recommendation would break something, explain here; else empty' },
    evidence: { type: 'string' },
  },
}

const dims = A.only ? DIMENSIONS.filter(d => A.only.includes(d.id)) : DIMENSIONS

const results = await pipeline(
  dims,
  d => agent(`${d.prompt}\n${GROUNDING}`, { label: `probe:${d.id}`, phase: 'Probe', schema: FINDINGS, model: 'haiku' }),
  (res, d) => {
    if (!res) return []
    const valid = res.findings.filter(f => f.evidence_command && f.evidence_output)
    if (valid.length < res.findings.length) log(`${d.id}: dropped ${res.findings.length - valid.length} finding(s) without evidence`)
    return parallel(valid.map((f, i) => async () => {
      if (f.kind === 'strength' || f.severity === 'info') return { dim: d.id, ...f, status: 'unverified_info' }
      const [rerun, refute] = await parallel([
        () => agent(`Re-run this evidence command from ${REPO} (it is read-only; if it would modify files, do not run it and set refuted=true): ${f.evidence_command}
Does the output still support the CONCLUSION of: "${f.title}"? (Claimed excerpt: ${JSON.stringify(f.evidence_output.slice(0, 600))}.) Judge the conclusion, not incidental counts — a different match count from a slightly different regex is not a refutation if the conclusion still holds. refuted=true only if the conclusion is not supported.`,
          { label: `rerun:${d.id}#${i}`, phase: 'Verify', schema: VERDICT, model: 'haiku', effort: 'low' }),
        () => agent(`Try to REFUTE this engineering finding in ${REPO}. Investigate independently with your own commands.
[${f.kind}/${f.severity}] ${f.title} (${f.file}:${f.line})
Recommendation given: ${f.recommendation}
Two separate questions: (1) refuted=true ONLY if the factual claim is false (or severity is off by 2+ levels) — default refuted=true if you cannot find supporting evidence yourself; (2) if the fact holds but the recommendation would break something (hard-coded paths, deployment layout, published artefacts), put that in recommendation_problem and keep refuted=false.`,
          { label: `refute:${d.id}#${i}`, phase: 'Verify', schema: VERDICT, model: 'haiku' }),
      ])
      const status = !rerun || rerun.refuted ? 'rejected' : (refute && refute.refuted ? 'disputed' : 'confirmed')
      const recommendation_problem = [rerun && rerun.recommendation_problem, refute && refute.recommendation_problem].filter(Boolean).join(' | ')
      return { dim: d.id, ...f, status, recommendation_problem, verify: { rerun, refute } }
    }))
  },
)

const all = results.filter(Boolean).flat().filter(Boolean)
const prior = all.filter(f => f.dim === 'prior-report-recheck').length
if (dims.some(d => d.id === 'prior-report-recheck') && prior < 12) log(`prior-report-recheck returned ${prior}/12 claim findings — some claims unaddressed`)
const sevOrder = { critical: 0, high: 1, medium: 2, low: 3, info: 4 }
const sorted = xs => xs.sort((a, b) => sevOrder[a.severity] - sevOrder[b.severity])
const esc = s => String(s || '').replace(/\|/g, '\\|').replace(/\n/g, ' ').slice(0, 400)
const row = f => `| ${f.severity} | ${f.kind} | ${esc(f.title)} | \`${esc(f.file)}:${f.line}\` | \`${esc(f.evidence_command)}\` | ${esc(f.recommendation)}${f.recommendation_problem ? ' ⚠ ' + esc(f.recommendation_problem) : ''} |`
const table = xs => xs.length ? ['| sev | kind | finding | location | evidence | recommendation |', '|---|---|---|---|---|---|', ...xs.map(row)].join('\n') : '_none_'

const confirmed = sorted(all.filter(f => f.status === 'confirmed'))
const disputed = sorted(all.filter(f => f.status === 'disputed'))
const rejected = all.filter(f => f.status === 'rejected')
const info = all.filter(f => f.status === 'unverified_info')

const footer = `\n---\nGenerated-by: technical-audit workflow (Haiku) | Verified-by: Haiku re-run + refute | Reviewed-by: T0 N\n`
const md = `# Technical Audit (evidence-bound)

Generated by the \`technical-audit\` workflow (Haiku probes + Haiku re-run/refute verification).
${all.length} findings: ${confirmed.length} confirmed, ${disputed.length} disputed, ${rejected.length} rejected, ${info.length} strengths/info.

## Confirmed
${table(confirmed)}

## Disputed
${table(disputed)}

## Strengths / info (not adversarially verified)
${table(info)}

## Rejected
${rejected.map(f => `- ${esc(f.title)} — ${esc(f.verify.rerun && f.verify.rerun.evidence)}`).join('\n') || '_none_'}
${footer}`

// Files are written by the caller from the returned markdown (avoids lossy copy through an agent).
return { total: all.length, confirmed: confirmed.length, disputed: disputed.length, rejected: rejected.length, report: 'audit/technical_audit.md', markdown: md, findings: all }
