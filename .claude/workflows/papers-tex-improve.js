export const meta = {
  name: 'papers-tex-improve',
  description: 'Tiered audit + rewrite of every LaTeX paper: Haiku audits and verifies each section, Sonnet rewrites it against audit/PAPER_FACTS.md, Haiku gates the rewrite (Sonnet repairs), Haiku splices and compiles, Haiku checks cross-paper consistency, Sonnet edits abstracts/correction notices',
  whenToUse: 'Run in the papers worktree after PAPER_FACTS.md is up to date. Args: {worktree?: string, only?: string[] unit ids, papers?: string[] paper ids}',
  phases: [
    { title: 'Audit', detail: 'Haiku: claims, tiers, forbidden phrases, numbers per section', model: 'haiku' },
    { title: 'Verify', detail: 'Haiku: quote fidelity of every finding in the section', model: 'haiku' },
    { title: 'Rewrite', detail: 'Sonnet (T1): tier-honest rewrite of the section into .paperwork/', model: 'sonnet' },
    { title: 'Gate', detail: 'Haiku mechanical gate; Sonnet repair on failure', model: 'haiku' },
    { title: 'Assemble', detail: 'Haiku: splice sections, pdflatex x2, fix-or-report', model: 'haiku' },
    { title: 'Consistency', detail: 'Haiku: same claim, same tier, same numbers across all papers + .zenodo.json', model: 'haiku' },
    { title: 'Editorial', detail: 'Sonnet: abstract + correction notice coherence per paper', model: 'sonnet' },
  ],
}

const A = args || {}
const WT = A.worktree || '/mnt/disks/disk-socrateai-local-1/dualscale-wt-papers'
const FACTS = `${WT}/audit/PAPER_FACTS.md`
const LEDGER = `${WT}/audit/physics_claims_ledger.md`

const PAPERS = [
  {
    id: 'tduality', file: 'papers/T-dulaity alone/T_duality_Alone.tex', pdfName: 'T_duality_Alone',
    units: [
      { id: 'td-preamble', range: [1, 45], role: 'preamble' },
      { id: 'td-abstract', range: [46, 57], role: 'abstract' },
      { id: 'td-intro', range: [58, 68], role: 'section' },
      { id: 'td-buscher', range: [69, 129], role: 'section' },
      { id: 'td-mukai', range: [130, 193], role: 'section' },
      { id: 'td-kummer', range: [194, 310], role: 'section' },
      { id: 'td-architecture', range: [311, 394], role: 'section' },
      { id: 'td-frontier', range: [395, 609], role: 'section' },
      { id: 'td-flops', range: [610, 637], role: 'section' },
      { id: 'td-gysin', range: [638, 708], role: 'section' },
      { id: 'td-cosmology', range: [709, 873], role: 'section' },
      { id: 'td-phase', range: [874, 1050], role: 'section' },
      { id: 'td-serverless', range: [1051, 1069], role: 'section' },
      { id: 'td-benchmark', range: [1070, 1122], role: 'section' },
      { id: 'td-conclusion', range: [1123, 1139], role: 'conclusion' },
      { id: 'td-biblio', range: [1140, 1279], role: 'bibliography' },
    ],
  },
  {
    id: 'leanflow', file: 'papers/leanflow-engine/leanflow_engine.tex', pdfName: 'leanflow_engine',
    units: [
      { id: 'lf-preamble', range: [1, 108], role: 'preamble' },
      { id: 'lf-abstract', range: [109, 133], role: 'abstract' },
      { id: 'lf-intro', range: [134, 171], role: 'section' },
      { id: 'lf-architecture', range: [172, 298], role: 'section' },
      { id: 'lf-dsl', range: [299, 334], role: 'section' },
      { id: 'lf-datasets', range: [335, 390], role: 'section' },
      { id: 'lf-toolkits', range: [391, 468], role: 'section' },
      { id: 'lf-benchmarks', range: [469, 505], role: 'section' },
      { id: 'lf-ecosystem', range: [506, 531], role: 'section' },
      { id: 'lf-conclusion', range: [532, 553], role: 'conclusion' },
      { id: 'lf-biblio', range: [554, 641], role: 'bibliography' },
    ],
  },
  { id: 'cover-cpc', file: 'papers/T-dulaity alone/submission/COVER_LETTER_CPC.tex', pdfName: 'COVER_LETTER_CPC', units: [{ id: 'cl-cpc', range: [1, 58], role: 'cover_letter' }] },
  { id: 'cover-jhep', file: 'papers/T-dulaity alone/submission/COVER_LETTER_JHEP.tex', pdfName: 'COVER_LETTER_JHEP', units: [{ id: 'cl-jhep', range: [1, 58], role: 'cover_letter' }] },
  { id: 'cover-scipost', file: 'papers/T-dulaity alone/submission/COVER_LETTER_SCIPOST.tex', pdfName: 'COVER_LETTER_SCIPOST', units: [{ id: 'cl-scipost', range: [1, 64], role: 'cover_letter' }] },
]

const COMMON = `
Worktree (your only writable area): ${WT}. Never edit main, proofs/, or any other repository.
Facts sheet — the ONLY admissible source for verification status, Lean results, simulation defects, benchmark status and citations: ${FACTS} (read it fully first).
Physics claims ledger (read if it exists; use its confirmed rows for your lines): ${LEDGER}.
Tiers: A kernel-checked+adequate · B exact arithmetic w/ negative control · L literature · C conjecture/interpretation · X numerics.
Forbidden below tier A: "formally verified", "certified", "certifies", "kernel-certified", "proves" (for physics), "zero axioms", "zero sorry" as a quality claim, "100%", "guarantee(s)", "first fully", "mechanized string theory". Tier C sentences need "we conjecture"/"if". Benchmarks per F4 → "pending hardware verification".
Never invent numbers, citations, DOIs or theorem names. Only theorem names that appear verbatim in the facts sheet may be cited.`

const pristine = (p, u) => `git -C "${WT}" show "HEAD:${p.file}" | sed -n '${u.range[0]},${u.range[1]}p' | cat -n   (displayed line + ${u.range[0] - 1} = true file line)`
const workPath = (p, u) => `${WT}/.paperwork/${p.id}/${u.id}.tex`

const AUDIT = {
  type: 'object', required: ['findings', 'summary'],
  properties: {
    summary: { type: 'string', description: 'one paragraph: what the section claims and at which tier it can honestly sit' },
    findings: { type: 'array', items: { type: 'object', required: ['line', 'quote', 'issue', 'category', 'justified_tier', 'fix'], properties: {
      line: { type: 'integer' }, quote: { type: 'string' }, issue: { type: 'string' },
      category: { type: 'string', enum: ['overclaimed_verification', 'withdrawn_lean_result', 'physics_error', 'unsupported_physical_claim', 'unreproducible_number', 'inconsistent_number', 'forbidden_phrase', 'mislabelled_method', 'citation_problem', 'latex_problem', 'sound'] },
      justified_tier: { type: 'string', enum: ['A', 'B', 'L', 'C', 'X', 'n/a'] },
      facts_ref: { type: 'string', description: 'F1..F6 item that supports the finding, or "source file <path>:<line>"' },
      fix: { type: 'string' } } } },
  },
}
const FIDELITY = { type: 'object', required: ['kept_indices', 'dropped'], properties: { kept_indices: { type: 'array', items: { type: 'integer' } }, dropped: { type: 'array', items: { type: 'string' } } } }
const REWRITE = { type: 'object', required: ['written_path', 'changes', 'numbers_kept', 't0_questions'], properties: {
  written_path: { type: 'string' }, changes: { type: 'array', items: { type: 'string' } },
  numbers_kept: { type: 'array', items: { type: 'string' }, description: 'each simulation/benchmark number still present and how it is wrapped' },
  t0_questions: { type: 'array', items: { type: 'string' } } } }
const GATE = { type: 'object', required: ['pass', 'problems'], properties: {
  pass: { type: 'boolean' },
  problems: { type: 'array', items: { type: 'object', required: ['severity', 'issue'], properties: { severity: { type: 'string', enum: ['blocking', 'minor'] }, issue: { type: 'string' } } } } } }

function roleRules(u) {
  switch (u.role) {
    case 'preamble': return 'Preamble: keep packages and macros; ensure \\newcommand{\\numreconcile}[2]{#1} is defined once (add after the last \\usepackage). Do NOT change the \\title or author block; if the title overclaims (e.g. "Lean 4 Certification"), put a proposed honest title in t0_questions.'
    case 'abstract': return 'Abstract: rewrite to tier-honest claims only (F1–F5). Immediately after \\end{abstract} add a short unnumbered paragraph "\\paragraph{Correction notice.}" stating, from F1/F3/F4/F5 only, what earlier versions claimed and what is withdrawn or downgraded, and that kernel-checked statements now come from LeanMaster v2.2.0 via lean_foundation/. Remove benchmark numbers.'
    case 'conclusion': return 'Conclusion: same tier discipline; no benchmark numbers; list open items (orientifold bookkeeping unresolved, numerics regenerated, statement review pending) as limitations.'
    case 'bibliography': return 'Bibliography: keep all existing \\bibitem entries unchanged unless a key is cited nowhere (check with grep on the pristine full file) — then leave it and note it in changes. Add \\bibitem{leanmaster2026} for "X. Callens et al., SocrateAI-Scientific-Agora-LeanMaster, release v2.2.0, https://github.com/xaviercallens/SocrateAI-Scientific-Agora-LeanMaster" if absent. Do not add any other reference.'
    case 'cover_letter': return 'Cover letter: rewrite so it describes the corrected manuscript honestly (tiers, withdrawn Buscher certification, benchmarks pending hardware verification, Zenodo new version per F5). Keep journal, editor salutation and signature. Add to t0_questions: "Should this submission letter still be sent, and to this journal?"'
    default: return 'Section: rewrite only what the confirmed findings require plus any forbidden phrase; preserve structure, labels (\\label), equations that are correct, figure/table environments and \\cite keys. Where the section describes a Lean result from proofs/*.lean, replace the claim with the LeanMaster/lean_foundation statement from F2 if one exists (cite with the F2 citation form and \\cite{leanmaster2026}), otherwise downgrade tier and say it is not kernel-checked. Where it describes the orientifold charges, mark as unresolved (F1). Wrap kept simulation numbers with \\numreconcile{value}{quantity}.'
  }
}

async function processUnit(p, u) {
  const audit = await agent(`Audit ONE section of a physics manuscript, line by line.
Paper: ${p.file} — unit ${u.id} (${u.role}), true lines ${u.range[0]}–${u.range[1]}.
Read it with: ${pristine(p, u)}
${COMMON}
For each problematic sentence give the TRUE line, a verbatim one-line quote, category, justified tier, the facts-sheet item (F1..F6) or source-file evidence, and a concrete fix. For the leanflow paper, claims about the software must be checked against the code in ${WT}/leanflow and ${WT}/scripts with grep (quote file:line in facts_ref). Record sound claims too (category "sound").`,
    { label: `audit:${u.id}`, phase: 'Audit', schema: AUDIT, model: 'haiku' })
  if (!audit) return { unit: u.id, paper: p.id, status: 'audit_failed' }

  const numbered = audit.findings.map((f, i) => `${i}. L${f.line}: ${JSON.stringify(f.quote)} — ${f.issue}`).join('\n')
  const fid = audit.findings.length ? await agent(`Quote-fidelity check. For each numbered finding, confirm the quote appears (whitespace-insensitive) within ±3 lines of the stated TRUE line of ${p.file} at git HEAD in ${WT}: use \`git -C "${WT}" show "HEAD:${p.file}" | sed -n 'A,Bp'\`. Also drop findings whose issue misreads the text. Return indices kept and a reason for each dropped one.\n${numbered}`,
    { label: `fidelity:${u.id}`, phase: 'Verify', schema: FIDELITY, model: 'haiku', effort: 'low' }) : { kept_indices: [], dropped: [] }
  const kept = (fid ? fid.kept_indices : []).map(i => audit.findings[i]).filter(Boolean)
  const actionable = kept.filter(f => f.category !== 'sound')

  const rewritePrompt = extra => `You are the T1 scientific editor. Rewrite unit ${u.id} of ${p.file} (true lines ${u.range[0]}–${u.range[1]}).
Source: ${pristine(p, u)}
${COMMON}
${roleRules(u)}
Verified findings for this unit (fix all of them; keep sound claims as they are):
${JSON.stringify(actionable, null, 1)}
Section summary from audit: ${audit.summary}
Output: write the COMPLETE replacement LaTeX for exactly these lines (no line numbers) to ${workPath(p, u)} (mkdir -p its directory). The file must replace the range 1:1 when spliced, so keep the unit's opening \\section/\\begin line and do not include text from neighbouring units. Minimal edits: do not restyle sound prose.
${extra || ''}`
  let rw = await agent(rewritePrompt(), { label: `rewrite:${u.id}`, phase: 'Rewrite', schema: REWRITE, model: 'sonnet' })
  if (!rw) return { unit: u.id, paper: p.id, status: 'rewrite_failed', findings: kept.length }

  const gatePrompt = () => `Mechanical gate for a rewritten manuscript unit. Compare original (${pristine(p, u)}) with the rewrite at ${workPath(p, u)}.
${COMMON}
Blocking if ANY: (1) a forbidden phrase remains for a claim below tier A (grep -niE "formally verif|certif|kernel-certif|zero axioms|zero .?sorry|100\\\\?%|guarantee|first fully" on the rewrite, then judge each hit — a sentence explicitly saying something is NOT certified is fine); (2) a number appears in the rewrite that is not in the original and not in ${FACTS} (extract numbers with grep -oE "[0-9]+([.,][0-9]+)?" from both and compare); (3) a benchmark number from F4 remains; (4) a theorem name is cited that is not in ${FACTS}; (5) \\begin/\\end environments or braces unbalanced (python3 count); (6) the unit's opening \\section/\\begin line or any \\label present in the original is missing; (7) text from outside the unit's range was added; (8) the rewrite claims something is kernel-checked that F2 does not list. Minor: style.
pass=true only if no blocking problem.`
  let gate = await agent(gatePrompt(), { label: `gate:${u.id}`, phase: 'Gate', schema: GATE, model: 'haiku' })
  let attempts = 1
  while (gate && !gate.pass && attempts < 2) {
    rw = await agent(rewritePrompt(`REPAIR: a gate rejected your previous file at ${workPath(p, u)}. Fix exactly these blocking problems, then overwrite the file:\n${gate.problems.filter(x => x.severity === 'blocking').map(x => '- ' + x.issue).join('\n')}`),
      { label: `repair:${u.id}`, phase: 'Gate', schema: REWRITE, model: 'sonnet' })
    gate = await agent(gatePrompt(), { label: `gate2:${u.id}`, phase: 'Gate', schema: GATE, model: 'haiku' })
    attempts++
  }
  return { unit: u.id, paper: p.id, range: u.range, status: gate && gate.pass ? 'ready' : 'gate_failed_escalate_T0', findings: kept.length, dropped: fid ? fid.dropped.length : 0,
    changes: rw ? rw.changes : [], numbers_kept: rw ? rw.numbers_kept : [], t0_questions: rw ? rw.t0_questions : [], gate_problems: gate ? gate.problems : [] }
}

const papers = (A.papers ? PAPERS.filter(p => A.papers.includes(p.id)) : PAPERS)
  .map(p => ({ ...p, units: A.only ? p.units.filter(u => A.only.includes(u.id)) : p.units }))
  .filter(p => p.units.length)
log(`${papers.length} paper(s), ${papers.reduce((n, p) => n + p.units.length, 0)} unit(s). zenodo_bundle/T_duality_Alone.tex is regenerated by copy at release, not edited.`)

const perPaper = await pipeline(
  papers,
  p => parallel(p.units.map(u => () => processUnit(p, u))),
  async (units, p) => {
    const ok = (units || []).filter(Boolean)
    const failed = ok.filter(r => r.status !== 'ready')
    const allUnits = PAPERS.find(x => x.id === p.id).units
    const partial = p.units.length !== allUnits.length
    const asm = await agent(`Assemble ${p.file} in ${WT}.
Ready units (their replacement files exist): ${JSON.stringify(ok.filter(r => r.status === 'ready').map(r => ({ id: r.unit, range: r.range, path: `${WT}/.paperwork/${p.id}/${r.unit}.tex` })))}
Units NOT ready (keep original lines): ${JSON.stringify(failed.map(r => r.unit))}${partial ? ' (partial run: all other units keep original lines)' : ''}
Steps: (1) write a python3 script that reads the file from \`git -C "${WT}" show "HEAD:${p.file}"\`, replaces each ready unit's TRUE line range with the replacement file content, processing ranges in descending order, and writes the result to ${WT}/${p.file}. (2) In a temp copy of the file's directory run \`pdflatex -interaction=nonstopmode\` twice (bibtex only if the file uses \\bibliography{}), and report the first 20 lines matching "^!" or "Undefined control sequence" or "Citation .* undefined" from the log. (3) If there are LaTeX errors caused by the splice (not pre-existing — check by compiling the HEAD version too), fix them minimally in ${WT}/${p.file} and recompile. (4) Copy the compiled PDF next to the .tex only if the original repo tracks a PDF for it (git -C "${WT}" ls-files). Return a short report: errors before/after, pages, and any unresolved error.`,
      { label: `assemble:${p.id}`, phase: 'Assemble', model: 'haiku' })
    return { paper: p.id, file: p.file, units: ok, assembleReport: asm }
  },
)

const results = perPaper.filter(Boolean)

phase('Consistency')
const consistency = await agent(`Cross-document consistency check in ${WT} after the rewrite of: ${results.map(r => r.file).join(', ')} and the deposit metadata ${WT}/.zenodo.json.
${COMMON}
For each headline claim (Buscher/T-duality verification status; Mukai/K3 lattice results; tadpole/orientifold status; 77/60 ratio tier; swampland/CDL/tachyon numerics status; TDA Mapper counts and meaning of beta_1; benchmark status; Zenodo correction notice) list how each document states it (file:line quote) and flag any document where tier, wording strength, or numbers differ. Also flag remaining forbidden phrases anywhere in these files (grep -niE). Do not edit files. Return markdown.`,
  { label: 'cross-paper-consistency', phase: 'Consistency', model: 'haiku' })

phase('Editorial')
const editorial = await parallel(results.filter(r => ['tduality', 'leanflow'].includes(r.paper)).map(r => () => agent(`T1 editorial pass on ${WT}/${r.file} (already rewritten and compiled).
${COMMON}
Read the abstract, the correction notice, the introduction and the conclusion together, plus this cross-document report:
${consistency}
Make them coherent with each other and with the body (same tiers, same scope, no contradictions, no leftover overclaim), editing ${WT}/${r.file} minimally. Do not touch simulation numbers except to keep \\numreconcile wrappers. Recompile with pdflatex twice in a temp copy and confirm no new errors. Return a bullet list of edits and any remaining T0 questions.`,
  { label: `editorial:${r.paper}`, phase: 'Editorial', model: 'sonnet' })))

return {
  worktree: WT,
  units: results.flatMap(r => r.units.map(u => ({ paper: r.paper, unit: u.unit, status: u.status, findings: u.findings, changes: u.changes, numbers_kept: u.numbers_kept, gate_problems: u.gate_problems }))),
  escalations: results.flatMap(r => r.units.filter(u => u.status !== 'ready').map(u => `${r.paper}/${u.unit}: ${u.status}`)),
  t0_questions: [...new Set(results.flatMap(r => r.units.flatMap(u => u.t0_questions || [])))],
  assemble: results.map(r => ({ paper: r.paper, report: r.assembleReport })),
  consistency, editorial,
  note: 'Changes are uncommitted in the worktree; orchestrator (T0) reviews git diff, then commits.',
}
