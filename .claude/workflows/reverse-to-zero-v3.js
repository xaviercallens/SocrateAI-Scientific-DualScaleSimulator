export const meta = {
  name: 'reverse-to-zero-v3',
  description: 'Reverse loop round 3 (round 2 = reverse-to-zero-v2, run 2026-09-19): uses the staged data in /mnt/disks/disk-socrateai-local-1/dualscale-data-r3 (Planck 2018 plik-lite likelihood, 1.34 M SDSS spectroscopic galaxies over 100 tiles, DESI/SDSS BAO repository) and the round-2 outcomes. commit-frozen registration before any data is loaded, then the four ranked follow-ups from round 1 (CMB TDA with a spectrum-matched null, cosmic-web TDA with redshift-space mocks and a powered gate, NANOGrav 15-yr c4 bound, DESI DR2 + CMB + SN); verbatim-only LeanMaster quotes at HEAD; two skeptics',
  whenToUse: 'Next round after reverse-to-zero-v2; read its REPORT and ERRATA first (audit/reverse_zero_r2/). Tests M0 (a hypothesis change, not a K3xT2 derivation) and tries to bound the untested c4_pta_product.',
  phases: [
    { title: 'Ground', detail: 'LeanMaster HEAD constraints, each quote grep-asserted' },
    { title: 'Register', detail: 'freeze statistic, null, seeds and decision rule per test in a commit, before any data are loaded' },
    { title: 'Experiment', detail: 'X1 CMB, X2 cosmic web, X3 PTA c4, X4 BAO+CMB+SN (memory-capped)' },
    { title: 'Skeptic', detail: 'statistics/reproducibility lens + claims/framing lens' },
    { title: 'Synthesize', detail: 'ledger, registered verdicts, errata, proposals; orchestrator saves REPORT.md' },
  ],
}

const WT = '/mnt/disks/disk-socrateai-local-1/dualscale-wt-reverse'
const OUT = `${WT}/audit/reverse_zero_r3`
const R1 = `${WT}/audit/reverse_zero`
const R2 = `${WT}/audit/reverse_zero_r2`
const PY = '/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python'
const LM = '/home/callensxavier_gmail_com/SocrateAI-Scientific-Agora-LeanMaster'
const PREREG = '/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/audit/PRE_REGISTRATION.md'
const DATA = `${WT}/data/real2`
const DR3 = '/mnt/disks/disk-socrateai-local-1/dualscale-data-r3'   // staged 2026-09-19, sha256 in MANIFEST.json

const RULES = `
GROUND RULES (strict; items 7-13 are the round-1 errata, see ${R1}/ERRATA.md):
1. Work ONLY in git worktree ${WT} (branch loop/reverse-zero). Outputs go under ${OUT}/<your-dir>/. Never touch proofs/, main or other worktrees, and never push. LeanMaster (${LM}) and ${PREREG} are READ-ONLY.
2. Python: ${PY} (gudhi 3.13, healpy, camb, astropy, numpy, scipy, pandas, sklearn). Seed every random call and state the seeds.
3. Evidence-bound: every number comes from a script you ran. Save the script and its JSON, and state the exact command with all arguments. Never invent a bound, DOI or dataset. If something cannot be fetched after 3 attempts, mark it ABSENT and list the URLs tried. Anything you did not try is NOT ATTEMPTED.
4. Tiers: A kernel-checked (a LeanMaster statement as written) / B exact arithmetic + negative control / L literature / C conjecture / X exploratory numerics. Data fits are X unless stated.
5. FRAMING RULE (non-negotiable): the zero-parameter model M0 (frozen flat LCDM with the imported Planck18 Omega_Lambda = 0.68885, GR tensor sector c4_pta_product = 0, no symmetron sector) is a HYPOTHESIS CHANGE, NOT a derivation from K3 x T2. It carries profiled nuisances (r_d*h, the SN offset); H0 is not frozen (PRE_REGISTRATION R1). Never write that K3 x T2 predicts or fixes mu_sym, c4, Omega_Lambda or any TDA outcome. LeanMaster Stream 8 records "Observables: none".
6. MEMORY SAFETY: run every python command as prlimit --as=10737418240 -- ${PY} ... in the FOREGROUND with timeout (no nohup or background jobs). Build one large complex or ensemble at a time, and free it before the next. The machine has 29 GB, no swap, and other agents running.
7. QUOTES ARE VERBATIM OR ABSENT. Every LeanMaster or PRE_REGISTRATION quote must be a fixed-string match. Assert it with git -C ${LM} show <commit>:<path> | grep -F -- '<quote>' (or grep -F on ${PREREG}) and record the command. Keep qualifiers such as "(the Tier C dual-scale identification of Stream 3)". Round 1 put invented text in LeanMaster's mouth; never paraphrase inside quotation marks.
8. TWO DIFFERENT P1s: LeanMaster Stream 6 P1 is the extra-dimension radius; PRE_REGISTRATION P1 is the cosmological constant. Name which one you mean.
9. "PRE-REGISTERED" means written in ${PREREG}, or in ${OUT}/registration/registration.json at a commit made BEFORE the data were loaded. The Register phase makes that commit. A threshold typed into a script is not pre-registered. The PRE_REGISTRATION threshold table (lines 34-36) is a table, not a rule for a comparison it does not name.
10. APPLY RULES TO THEIR DATASETS: P1's falsification rule names DESI DR3/final BAO + CMB + SN, Euclid or Rubin. DR2-based results are a pipeline action or an informational result, never a P1 verdict.
11. WORDING: "not rejected" / "consistent with" never becomes "confirms", "cleared decisively" or "the data reward". Convert Delta chi2 to sigma with scipy (chi2.sf, then norm.isf(p/2)), never by eye.
12. REPRODUCIBILITY: results JSON is written by scripts. Paths are relative to the repo root, found with pathlib from the script location. State commands as: cd <repo-relative dir> && ${PY} <script> <args>.
13. COMMITS: stage only your own files by explicit path, check git diff --cached --name-only, and never use git add -A. Message ends with: Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
`

const LEANC = { type: 'object', properties: { head: { type: 'string' }, describe: { type: 'string' }, quotes: { type: 'array', items: { type: 'object', properties: {
  topic: { type: 'string' }, quote: { type: 'string' }, file_line: { type: 'string' }, grep_command: { type: 'string' }, grep_matched: { type: 'boolean' }, consequence_for_simulator: { type: 'string' } },
  required: ['topic', 'quote', 'file_line', 'grep_command', 'grep_matched'] } },
  changed_since_eb791e7: { type: 'array', items: { type: 'string' } } }, required: ['head', 'quotes'] }
const REG = { type: 'object', properties: { commit: { type: 'string' }, file: { type: 'string' }, tests: { type: 'array', items: { type: 'object', properties: {
  id: { type: 'string' }, statistic: { type: 'string' }, null_model: { type: 'string' }, seeds: { type: 'string' }, decision_rule: { type: 'string' }, effective_tests: { type: 'string' }, data_files: { type: 'array', items: { type: 'string' } } },
  required: ['id', 'statistic', 'null_model', 'decision_rule'] } } }, required: ['commit', 'file', 'tests'] }
const EXP = { type: 'object', properties: {
  id: { type: 'string' }, question: { type: 'string' }, result: { type: 'string' }, numbers: { type: 'object', additionalProperties: true },
  registered_rule: { type: 'string', description: 'verbatim from registration.json or PRE_REGISTRATION.md, or "none"' },
  registration_commit_precedes_data_load: { type: 'boolean' },
  verdict: { type: 'string' }, parameter_effect: { type: 'string' },
  files: { type: 'array', items: { type: 'string' } }, caveats: { type: 'array', items: { type: 'string' } }, absent_or_not_attempted: { type: 'array', items: { type: 'string' } }, commit: { type: 'string' } },
  required: ['id', 'question', 'result', 'verdict', 'caveats'] }
const AUDIT = { type: 'object', properties: { lens: { type: 'string' },
  per_item: { type: 'array', items: { type: 'object', properties: { id: { type: 'string' }, reproduced: { type: 'boolean' }, stands: { type: 'boolean' }, reason: { type: 'string' } }, required: ['id', 'stands', 'reason'] } },
  framing_violations: { type: 'array', items: { type: 'string' } }, commit: { type: 'string' } }, required: ['lens', 'per_item', 'framing_violations'] }

phase('Ground')
const leanc = await agent(`${RULES}
LEANMASTER CONSTRAINTS AT HEAD (read-only). Run git -C ${LM} describe --tags and rev-parse HEAD. Round 1 read eb791e7 (v3.25.0). List the Stream 6/7/8 docs and Lean files changed since then (git log eb791e7..HEAD --stat).
Collect VERBATIM quotes (rule 7) for: the Stream 6 verdict on the extra-dimension radius and the kappa = 1 statement with its Tier C qualifier; the Stream 7 C-A/C-B verdicts, and the sentence on what a new hypothesis must drop; the Stream 8 "Observables" lines and any newer Stream 8 results; DarkEnergyScale.lean on where Omega_Lambda comes from; CosmicString.lean planckGmuBound; Tadpole.lean on tadpole_budget not being a construction. Every quote gets its grep command and a matched flag; drop any that do not match.
Also report whether ANY statement at HEAD mentions mu_sym, c4_pta_product, a unit-bearing observable, or a TDA prediction. Grep for them and quote what you find, or say none.`,
  { label: 'ground:leanmaster-head', phase: 'Ground', model: 'haiku', schema: LEANC })

const ctx = `ROUND 1: ${R1}/REPORT.md and ${R1}/ERRATA.md. ROUND 2: ${R2}/ (its REPORT.md if it was saved, else its registration, ledger and per-experiment JSON). Read them; reuse their code and fix the defects they name. Do NOT repeat a test whose round-2 outcome was conclusive; re-run only what round 2 left INCONCLUSIVE or NOT ATTEMPTED, on the better data below.
LEANMASTER AT HEAD: ${JSON.stringify(leanc)}
STAGED DATA (read ${DR3}/MANIFEST.json first; sha256 for every file; note its caveats and its not_found list): ${DR3}/sdss_spec_bulk/sdss_spec_ra100_250_dec-10_50.csv (1.34 M spectroscopic galaxies, z 0.05-0.6, heterogeneous selection, no random catalogue), ${DR3}/planck_2018 (plik-lite v22 TT and TTTEEE: binned data, covariance, bin weights, best-fit theory), ${DR3}/desi_sdss_bao (Cobaya bao_data: DESI DR1/DR2, eBOSS DR16, MGS, DR12), ${DR3}/{input,output}.json (NANOGrav per-frequency strain only; NOT angular HD).
DATA ON DISK (sha256-verified in round 1): ${DATA}/{dark_energy/desi_dr2, dark_energy/Pantheon+SH0ES_STAT+SYS.cov, cosmic_web/sdss_dr17_*.csv, cmb/wmap_ilc_9yr_v5.fits, cmb/wmap_temperature_kq85_analysis_mask_r9_9yr_v5.fits, pulsar_timing/NANOGrav15yr_PulsarTiming_v2.1.0}.`

phase('Register')
const reg = await agent(`${RULES}
${ctx}
REGISTRATION (you load NO data: no FITS, CSV, covariance or timing file may be opened, only file names and sizes). For each test X1-X4 below, write ${OUT}/registration/registration.json with: the statistic (exact definition), the null model and how it is calibrated, the seeds, the number of simulations or mocks, the effective number of tests and the multiplicity correction, the decision rule (threshold and what each outcome means for M0), and the data files. Commit it ALONE, before any other agent loads data, and return the commit hash.
X1 CMB TDA: the Betti/Euler curves (sublevel and superlevel) of the WMAP ILC map, and Planck SMICA if X1 can fetch it. The null is Gaussian sims whose MASKED PSEUDO-C_ell MATCHES THE DATA: lmax >= 3*nside, pixwin applied, a single mask pass, and C_ell from MASTER deconvolution or an iterated correction. Include a spectrum-match gate (max |z| per ell band) that must pass BEFORE the p-values are computed; if it fails, no p-values are reported. Fix the effective number of tests now, given that euler_chi = b0 - b1 and that sub/superlevel come from one map. The cosmic-string sensitivity is calibrated on the same statistic set as the test.
X2 cosmic-web TDA: SDSS DR17 against lognormal mocks in REDSHIFT SPACE (Kaiser plus a small-scale velocity dispersion, both stated), a footprint mask at nside >= 256, and the bias fitted on xi(r) over one r range. The GATE is a full-covariance chi2 on a DISJOINT r range, and its power must be shown first: it must reject xi = 0 and xi x 3. The topology statistic is restricted to the scale range the gate validates. Enlarge N, or use several non-overlapping subvolumes, so that range carries information (round 1 found all curves at [1,0,0] for r in [20,40] Mpc/h at N = 25000).
X3 PTA: bound c4_pta_product from binned Gamma(theta) with covariance, derived from NANOGrav 15-yr timing data (the enterprise/enterprise_extensions optimal statistic). The decision rule is stated as a bound, e.g. |c4| < x at 95%.
X4 BAO+CMB+SN: DESI DR2 BAO + a compressed CMB likelihood (distance priors: R, l_A, omega_b; cite the paper and the table) + Pantheon+ full covariance with ONE SN cut used everywhere (1580 cosmology-only primary, 1590 disclosed). Report CPL vs Lambda. This is informational for P1: PRE_REGISTRATION's P1 rule names DR3/final data.
Return the registration summary.`,
  { label: 'register:freeze', phase: 'Register', model: 'sonnet', schema: REG })
log(`registration commit: ${reg ? reg.commit : 'MISSING'}`)
if (!reg || !reg.commit) return { stopped: 'no registration commit; experiments not run', leanc, reg }

const regctx = `REGISTRATION (frozen at commit ${reg.commit}, file ${reg.file}): ${JSON.stringify(reg.tests)}
Apply the registered rule exactly. If you must deviate (a bug, an infeasible step), do not edit registration.json: record the deviation, its reason and its effect in your JSON, and report the registered-rule outcome as well, if it can be computed.`

phase('Experiment')
const exps = (await parallel([
  () => agent(`${RULES}
${ctx}
${regctx}
X1 - CMB TDA WITH A SPECTRUM-MATCHED NULL. Reuse ${R1}/E5-cmb-tda/cmb_tda.py. Fix the three round-1 causes of the pseudo-C_ell mismatch: lmax 256 at nside 128, pixwin off, and the mask applied twice. Try to fetch Planck PR3 SMICA and its common mask from the Planck Legacy Archive or IRSA; if that fails 3 times, mark it ABSENT. Run the registered gate first; only if it passes compute the p-values with N >= 100 sims (split into chunks if memory requires). Then the string-injection sensitivity on the registered statistic set, compared (not equated) with planckGmuBound. Save under ${OUT}/X1-cmb/ and commit.`,
    { label: 'exp:X1-cmb-matched-null', phase: 'Experiment', model: 'sonnet', schema: EXP }),
  () => agent(`${RULES}
${ctx}
${regctx}
X2 - COSMIC-WEB TDA IN REDSHIFT SPACE WITH A POWERED GATE, ON THE 1.34 M-GALAXY CATALOGUE (${DR3}/sdss_spec_bulk). Round 1 used 193 k galaxies (RA 140-220, Dec 0-50, z 0.02-0.12) and could not get information in the validated r range. The staged catalogue is ~7x larger and reaches z = 0.6. Use several NON-OVERLAPPING volume-limited subvolumes (state the cuts; the selection changes from MGS to LRG-like at z ~ 0.15, so define the volume-limited cuts from the n(z) you measure), each with N large enough that the gate-validated range carries topological information. Use the subvolumes as independent replicates for a variance estimate. Selection functions: build the angular mask from the data footprint at nside >= 256 and the radial selection from n(z), and state that no survey random catalogue exists. The registration governs the design; if it names the 193 k catalogue, run that as registered AND report the larger catalogue as a disclosed extension. Reuse ${R1}/E5-cosmic-web-tda-scaled/ (the e5b lognormal pipeline, including its round-1 bug fixes). Changes: redshift-space displacements in the mocks; a finer footprint mask; a Poisson null matched shell by shell (round 1 failed a KS test on r, p = 0.009); the full-covariance gate on a range disjoint from the bias fit, with its power demonstrated before use; and the statistic restricted to the validated range, with enough points there to carry information. Known-answer controls (planted voids, circle) at the same N. Save under ${OUT}/X2-cosmic-web/ and commit.`,
    { label: 'exp:X2-cosmic-web-redshift-space', phase: 'Experiment', model: 'sonnet', schema: EXP }),
  () => agent(`${RULES}
${ctx}
${regctx}
X3 - BOUND c4_pta_product WITH NANOGRAV 15-YR DATA. The timing data are at ${DATA}/pulsar_timing/NANOGrav15yr_PulsarTiming_v2.1.0; there is no angular Hellings-Downs table (round 1: ABSENT). Install enterprise, enterprise_extensions and their dependencies into a SEPARATE venv under /mnt/disks/disk-socrateai-local-1/venv-pta (do not modify ${PY}'s venv). If the install fails 3 times, record ABSENT with the errors and stop. Otherwise compute the optimal-statistic cross-correlations per pulsar pair (fixed CURN spectrum, as the registration says), bin them in angle with their covariance, and fit Gamma(theta) = A * [HD(theta) + c4 * the P4 term exactly as param_loop_sim.compute_pta_observable defines it]. Report c4 with its 68/95% interval, and a monopole/dipole negative control. Keep runtime under 60 min: use the pulsar subset and noise model the registration states. Save under ${OUT}/X3-pta/ and commit.`,
    { label: 'exp:X3-pta-c4-bound', phase: 'Experiment', model: 'sonnet', schema: EXP }),
  () => agent(`${RULES}
${ctx}
${regctx}
X4 - DESI DR2 BAO + PLANCK 2018 PLIK-LITE LIKELIHOOD + PANTHEON+. Reuse ${R1}/E1-desi-dr2/e1_desi_dr2.py and ${R1}/E3-M0/m0_model.py and unify them into one script with ONE SN cut (rule: 1580 primary, 1590 disclosed). Add the Planck 2018 plik-lite likelihood from ${DR3}/planck_2018: the binned C_ell (cl_cmb_plik_v22.dat), covariance (c_matrix_plik_v22.dat), bin edges and weights (blmin/blmax/bweight), and the Planck calibration nuisance A_planck. Compute theory C_ell with CAMB for each parameter point (w0, wa, Omega_m, H0, omega_b, omega_c, A_s, n_s, tau are all needed for a real CMB fit: state which are fitted, which are fixed and why, since a CMB likelihood cannot be used with Omega_m alone). CHECK the implementation first: reproduce the file's own best-fit theory (bf_lite_plikTTTEEE_v22b_lowl_simall.minimum.theory_cl) chi2 against the binned data, and check chi2/dof is about 1 for LCDM at Planck's best fit, before using it. plik-lite has no low-ell part (l < 30): say so and add a tau prior from a cited paper, or state that the result omits it. Report the 3.1 sigma comparison from the published DESI DR2 analysis only as a comparator. Report CPL vs Lambda with Omega_m free (2 dof); M0 vs fitted LCDM; and whether the published ~3.1 sigma is approached. For P1 this is informational (rule 10). Save under ${OUT}/X4-bao-cmb-sn/ and commit.`,
    { label: 'exp:X4-bao-cmb-sn', phase: 'Experiment', model: 'sonnet', schema: EXP }),
])).filter(Boolean)

phase('Skeptic')
const all = `REGISTRATION: ${JSON.stringify(reg)}
EXPERIMENTS: ${JSON.stringify(exps)}
LEANMASTER: ${JSON.stringify(leanc)}`
const audits = (await parallel([
  () => agent(`${RULES}
SKEPTIC - STATISTICS AND REPRODUCIBILITY lens. You produced none of this; default to distrust.
${all}
For each experiment:
1. Check with git log that the registration commit ${reg.commit} precedes the first commit that loaded data.
2. Rerun the main script in a clean clone (git clone ${WT} /tmp/rz2_clone, run there, delete it afterwards) and confirm the headline numbers.
3. Check that the registered rule was applied exactly, and that every deviation is disclosed.
4. X1: the spectrum-match gate, per ell band.
5. X2: the gate's power (it rejects xi = 0), that the gate range is disjoint from the bias-fit range, and the redshift-space treatment.
6. X3: the optimal-statistic normalisation, the pair covariance, and the P4-term definition against param_loop_sim.
7. X4: the CMB prior values against the cited table, a single SN cut, dof, and the sigma conversion.
Decide stands / does not stand per item. Write ${OUT}/skeptic_statistics/verdict.json and commit.`,
    { label: 'skeptic:statistics', phase: 'Skeptic', schema: AUDIT }),
  () => agent(`${RULES}
SKEPTIC - CLAIMS AND FRAMING lens. You produced none of this; default to distrust.
${all}
1. Grep-assert EVERY quotation mark in the committed files of this round against its source (rule 7). List each non-matching quote with its file:line.
2. Check rules 5 and 8-11 everywhere: hypothesis change, not derivation; the two P1s; "pre-registered" only where rule 9 allows; P1 verdicts only on P1's named datasets; no confirmation language.
3. Check that every ABSENT and NOT ATTEMPTED item is reported.
4. Check that nothing attributes to LeanMaster a statement about mu_sym, c4 or a TDA prediction that its text does not contain.
List violations with file:line. Write ${OUT}/skeptic_framing/verdict.json and commit.`,
    { label: 'skeptic:framing', phase: 'Skeptic', schema: AUDIT }),
])).filter(Boolean)

phase('Synthesize')
const report = await agent(`${RULES}
Write the REPORT as your final text (the harness blocks subagents from writing .md files; the orchestrator will save it). Also write ${OUT}/ledger_reverse_zero_r2.json with a deterministic builder script that reads only committed JSON, and commit both.
Inputs: ${all}
AUDITS: ${JSON.stringify(audits)}
Sections:
1. Bottom line: M0's parameter count (0 by hypothesis change, plus its profiled nuisances) and M2's (2). Say whether c4_pta_product moved from untested to bounded (X3).
2. Registered tests: for each of X1-X4, the registered rule verbatim, the registration commit, whether it preceded the data load, and the outcome under the rule. Deviations come first.
3. What changed relative to round 1 (REPORT 13756dc, ERRATA b77d11a).
4. LeanMaster at HEAD: verbatim quotes only, what changed since eb791e7, and what they close.
5. PROPOSED ADDENDUM text for PRE_REGISTRATION.md: dated, appended, never editing the earlier text. It records each X1-X4 registration and outcome, and any new bound on c4.
6. Ranked next steps (Tier C), each with the observable (with units), a freeze-before-compare plan and a realistic difficulty.
7. ERRATA table for this round: each skeptic finding, with file:line and the correct reading.
8. Everything ABSENT, NOT ATTEMPTED, disputed, or from memory.
Header: "Generated by workflow reverse-to-zero-v2; exploratory (tier X) unless stated; M0 is a hypothesis change, not a derivation."`,
  { label: 'synthesize:report', phase: 'Synthesize' })

return { leanc, reg, exps, audits, report }
