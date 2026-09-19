export const meta = {
  name: 'reverse-to-zero',
  description: 'Reverse loop: zero-parameter candidate M0 (frozen LCDM + GR tensor + no screening; hypothesis change, not derivation) -> predicted topology -> real-data tests (DESI DR2 pre-registered, NANOGrav HD, SDSS cosmic-web TDA, CMB TDA) -> hypothesis verdict and improvement proposals',
  whenToUse: 'Test the zero-parameter version of the dual-scale cosmology against real data and topology, using LeanMaster verdicts as constraints.',
  phases: [
    { title: 'Ground', detail: 'real data fetch + LeanMaster Streams 6-8 constraints' },
    { title: 'Check', detail: 'independent data check' },
    { title: 'Experiment', detail: 'DESI DR2 pre-registered test, PTA c4 test, M0 definition + screening' },
    { title: 'TDA', detail: 'SDSS cosmic web + CMB persistent homology with nulls' },
    { title: 'Skeptic', detail: 'statistics lens + claims/circularity lens' },
    { title: 'Synthesize', detail: 'ledger, hypothesis verdict, improvement proposals' },
  ],
}

const WT = '/mnt/disks/disk-socrateai-local-1/dualscale-wt-reverse'
const OUT = `${WT}/audit/reverse_zero`
const PY = '/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python'
const LM = '/home/callensxavier_gmail_com/SocrateAI-Scientific-Agora-LeanMaster'
const PREREG = '/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/audit/PRE_REGISTRATION.md'

const RULES = `
GROUND RULES (strict):
- Work ONLY in git worktree ${WT} (branch loop/reverse-zero). Outputs under ${OUT}/<your-dir>/. Never touch proofs/, main, other worktrees; never push. LeanMaster (${LM}) and ${PREREG} are READ-ONLY. Commit your files; message ends with: Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
- Python: ${PY} (gudhi 3.13, healpy, camb, astropy, numpy, scipy, pandas, networkx, sklearn). Seed every random call; state seeds.
- Evidence-bound: every number from a script you ran; save script + JSON; give the exact command with all arguments. Never invent a bound, DOI or dataset. If something cannot be fetched after 3 attempts: ABSENT with URLs tried.
- NO literal answers in a computation path; thresholds that were pre-registered are read from ${PREREG} (quote the line) - do not change them.
- Tiers: A kernel-checked (LeanMaster statement as written) / B exact arithmetic + negative control / L literature / C conjecture / X exploratory numerics. All data fits here are X unless stated.
- FRAMING RULE (non-negotiable): the zero-parameter model M0 is a HYPOTHESIS CHANGE (sectors removed or set to their GR value), NOT a derivation from K3 x T2. LeanMaster Streams 6-7 record that every physical reading of the dual-scale idea tried so far is excluded and that Streams 4-5 results are dimensionless. Never write that K3 x T2 predicts or fixes mu_sym, c4 or Omega_Lambda.
`

const DATA = { type: 'object', properties: { datasets: { type: 'array', items: { type: 'object', properties: {
  id: { type: 'string' }, purpose: { type: 'string' }, path: { type: 'string' }, url: { type: 'string' }, sha256: { type: 'string' }, rows: { type: 'string' }, status: { type: 'string' }, note: { type: 'string' } },
  required: ['id', 'purpose', 'status'] } } }, required: ['datasets'] }
const LEANC = { type: 'object', properties: { tag: { type: 'string' }, verdicts: { type: 'array', items: { type: 'object', properties: {
  item: { type: 'string' }, statement_or_quote: { type: 'string' }, file_line: { type: 'string' }, consequence_for_simulator: { type: 'string' } },
  required: ['item', 'statement_or_quote', 'file_line', 'consequence_for_simulator'] } },
  zero_param_statements: { type: 'array', items: { type: 'string' } }, has_units: { type: 'array', items: { type: 'string' } } },
  required: ['tag', 'verdicts'] }
const CHECK = { type: 'object', properties: { usable: { type: 'array', items: { type: 'string' } }, rejected: { type: 'array', items: { type: 'string' } }, notes: { type: 'array', items: { type: 'string' } } }, required: ['usable', 'rejected'] }
const EXP = { type: 'object', properties: {
  id: { type: 'string' }, question: { type: 'string' }, result: { type: 'string' }, numbers: { type: 'object', additionalProperties: true },
  preregistered_rule: { type: 'string' }, verdict: { type: 'string' }, parameter_effect: { type: 'string', description: 'what this does to the free-parameter count, with mechanism and tier' },
  files: { type: 'array', items: { type: 'string' } }, caveats: { type: 'array', items: { type: 'string' } }, commit: { type: 'string' } },
  required: ['id', 'question', 'result', 'verdict', 'caveats'] }
const AUDIT = { type: 'object', properties: { lens: { type: 'string' },
  per_item: { type: 'array', items: { type: 'object', properties: { id: { type: 'string' }, reproduced: { type: 'boolean' }, stands: { type: 'boolean' }, reason: { type: 'string' } }, required: ['id', 'stands', 'reason'] } },
  framing_violations: { type: 'array', items: { type: 'string' } }, commit: { type: 'string' } }, required: ['lens', 'per_item', 'framing_violations'] }

phase('Ground')
const [data, leanc] = await parallel([
  () => agent(`${RULES}
TASK: fetch REAL public data into ${WT}/data/real2/ with MANIFEST.json (id, url, sha256, rows, columns, licence/citation as printed by the source). Reuse ${WT}/data/real/ (already verified: DESI 2024 BAO, Pantheon+SH0ES.dat, NANOGrav ceffyl KDE, 2MRS RA/Dec) where useful.
1. DESI DR2 BAO: GitHub CobayaSampler/bao_data folder 'desi_bao_dr2' (list it via the GitHub API, fetch mean + covariance files for the ALL combination and per-tracer if present).
2. Pantheon+ full STAT+SYS covariance: PantheonPlusSH0ES/DataRelease, 'Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES_STAT+SYS.cov'.
3. NANOGrav 15-yr Hellings-Downs ANGULAR correlation data (binned cross-correlation vs angular separation, or optimal-statistic pair data): search the nanograv GitHub org (e.g. 15yr_stochastic_analysis, PTArcade, the NANOGrav 15yr data release on Zenodo) via APIs; fetch only a machine-readable table of correlation vs angle with errors. If only figures exist: ABSENT (do not digitise).
4. Cosmic web for TDA: SDSS DR17 main galaxy sample with redshifts via SkyServer SqlSearch (https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch?cmd=...&format=csv, works): select ra, dec, z, zErr from SpecObj where class='GALAXY' and zWarning=0 and z between 0.02 and 0.12 in a contiguous patch (e.g. ra 140-220, dec 0-50), paging by ra strips if a row limit applies; target 50k-200k galaxies. Also fetch the matching sky footprint/angular mask proxy (or build a random catalogue by sampling the same RA/Dec footprint from the data's own angular distribution with shuffled redshifts - document it).
5. CMB map for TDA: WMAP 9-year ILC map from LAMBDA (https://lambda.gsfc.nasa.gov/data/map/dr5/dfp/ilc/wmap_ilc_9yr_v5.fits or current path), plus a galactic mask (WMAP KQ85 temperature analysis mask from LAMBDA). If WMAP fails, try Planck PR3 SMICA from IRSA (large) and note size.
6. Fifth-force: arXiv:2002.11761 (Lee et al. 2020, Eot-Wash): fetch the arXiv source tarball and check for a machine-readable alpha-lambda exclusion table; if none: ABSENT.
Sanity-check each file is data, not HTML. Keep files > 5 MB out of git (.gitignore) but on disk. Commit manifest + small files.`,
    { label: 'ground:data', phase: 'Ground', model: 'sonnet', schema: DATA }),
  () => agent(`${RULES}
TASK (read-only LeanMaster): run git -C ${LM} describe --tags. Read docs/STREAM6_EXPERIMENT_PLAN.md s7, docs/STREAM7_HYPOTHESIS_INVENTORY.md s2,6,7, docs/STREAM8_WHICH_K3.md s5-7, DualScaleCosmology/Stream6Verdict.lean, Stream7CA.lean, Stream7CB.lean, DualScaleDyons/SelfDualT2.lean, and the paper 11 source if present (grep -ril "paper 11\\|pre-registered confrontation" docs papers | head). For each verdict quote the theorem statement or doc sentence verbatim with file:line and state its consequence for the simulator's remaining parameters (mu_sym, c4_pta_product), for our pre-registered P1 (frozen Lambda) and P2 (self-dual length bracket [21,105] um in ${PREREG}), and for any TDA prediction. List which LeanMaster statements are zero-parameter and which carry units. Do not interpret beyond the quotes.`,
    { label: 'ground:leanmaster-streams6-8', phase: 'Ground', model: 'haiku', schema: LEANC }),
])

phase('Check')
const check = await agent(`${RULES}
INDEPENDENT CHECK of the fetched data (you fetched nothing). DATA: ${JSON.stringify(data)}
For each FETCHED item: sha256sum, open it, confirm it matches the claimed content and source (columns, sizes, redshift ranges, map nside/units, covariance symmetric positive definite). Reject anything that is HTML, truncated, mislabelled, or digitised from a figure. Return usable ids, rejected ids with reasons.`,
  { label: 'check:data', phase: 'Check', model: 'sonnet', schema: CHECK })

const ctx = `USABLE DATA: ${JSON.stringify(check.usable)} (rejected: ${JSON.stringify(check.rejected)}). Data manifest: ${JSON.stringify(data)}.
LEANMASTER CONSTRAINTS: ${JSON.stringify(leanc)}.
Existing harness: ${WT}/scripts/param_loop_sim.py (lambda_sym already deleted); decisive experiment code: ${WT}/audit/zero_param_loop/round2/decisive_experiment.py and round3/skeptic/chi2_independent.py; cumulative report ${WT}/audit/zero_param_loop/REPORT.md (2 free params: mu_sym, c4_pta_product).`

phase('Experiment')
const exps = (await parallel([
  () => agent(`${RULES}
${ctx}
EXPERIMENT E1 - the pre-registered P1 test. Read P1 in ${PREREG} (quote the rule lines). Build ${OUT}/E1-desi-dr2/ reusing the round-2/3 chi2 code but with DESI DR2 BAO (full covariance) and Pantheon+ with the FULL STAT+SYS covariance (Hubble-flow cut as before; state it). Compute: flat LCDM with Omega_L = 0.68885 frozen; LCDM fitted; CPL (w0,wa) fitted with Omega_m fitted; report delta chi2 CPL vs Lambda (2 dof) and the sigma, and apply the registered rule verbatim (>= 28.74 falsifies at 5 sigma; also report whether >= 11.83, i.e. 3 sigma). Negative control: Omega_L = 0.5. Repeat with DR1 to show the change DR1 -> DR2. Also the Omega_L-away-from-0.68885 rule. Say plainly if the published DESI DR2 3.1 sigma (DESI+CMB+SN) is or is not reproduced - you have no CMB here, so expect a different number and explain.`,
    { label: 'exp:E1-desi-dr2-prereg', phase: 'Experiment', model: 'sonnet', schema: EXP }),
  () => agent(`${RULES}
${ctx}
EXPERIMENT E2 - is the PTA sector GR? M0 sets c4_pta_product = 0 (pure Hellings-Downs). If a verified angular correlation dataset exists: fit Gamma(theta) = A * [HD(theta) + c4 * P4-term as in param_loop_sim.compute_pta_observable] to it (amplitude A free, since the harness normalisation is arbitrary), report c4 best fit and 68/95% interval, chi2 of c4 = 0 vs best, and a monopole/dipole negative control (show a monopole-only model fits worse). If NO angular data: say so, and instead state exactly which public product would decide it; parameter_effect = none (count unchanged by data). Either way, setting c4 = 0 is 'set to the GR value' (hypothesis change, tier L for GR), not a derivation.`,
    { label: 'exp:E2-pta-hd', phase: 'Experiment', model: 'sonnet', schema: EXP }),
  () => agent(`${RULES}
${ctx}
EXPERIMENT E3 - define M0 and test its screening consequence. M0 := flat LCDM (Omega_L = 0.68885 frozen, H0 a nuisance as before), GR tensor sector (c4 = 0), NO symmetron sector (mu_sym removed with the sector). Implement M0 in ${OUT}/E3-M0/m0_model.py via the harness (no new tuning) and verify it has zero free theory parameters (list every number in it with its source: Lean def, literature, or nuisance).
Screening: M0 predicts NO screened fifth force. If a verified machine-readable fifth-force bound exists (Eot-Wash table), state what it says about the symmetron sector (bound on mu_sym in the harness units ONLY if the harness has a unit bridge - it does not today, so most likely: 'no bridge, cannot map') and whether any data prefers the sector. Compute information criteria M0 vs M2 (mu_sym, c4 free) vs fitted LCDM on the SAME data as E1-DR1 (reuse round-3 numbers if identical; recompute otherwise): Delta AIC/BIC.
Also record, quoting LeanMaster, that P2 (self-dual length) is excluded by Stream 6 with kappa = 1 fixed by the programme's T-duality, and what that means for the [21,105] um bracket registered in ${PREREG} (do not edit that file).`,
    { label: 'exp:E3-M0-definition', phase: 'Experiment', model: 'sonnet', schema: EXP }),
])).filter(Boolean)

phase('TDA')
const tdas = (await parallel([
  () => agent(`${RULES}
${ctx}
TDA T1 - COSMIC WEB (the gain the user wants: real topology). Use the verified SDSS galaxy catalogue with redshifts. Convert to comoving Cartesian coordinates with M0's frozen cosmology (Omega_m = 0.31115; H0 = 67.66 only sets units - say so). Build a volume-limited or number-density-flattened subsample (state cuts), N ~ 20k-60k points (subsample with seed if needed for runtime).
Persistent homology with INRIA GUDHI: alpha complex (gudhi.AlphaComplex) H0, H1, H2; Betti curves beta_k(r) and persistence diagrams; also the Euler characteristic curve.
NULLS (all with the same N, same angular footprint and same radial selection): (i) Poisson random in the survey volume; (ii) redshift-shuffled data (keeps angular and radial distributions, destroys 3D structure); (iii) lognormal mocks: generate a Gaussian field with the linear matter P(k) from CAMB at the M0 cosmology (sigma8/As from Planck 2018 best fit - quote the value and say it is imported), apply a lognormal transform with a linear bias b fitted ONLY to the galaxy two-point correlation (fit b on xi(r) of the data, never on topology), Poisson-sample galaxies into the same footprint/selection; >= 20 mocks with seeds.
Known-answer control: a synthetic cloud with a planted set of voids or loops recovers them.
Report: which Betti curves of the data are outside the mock envelope (95%), in which scale range; the data-vs-mock and data-vs-Poisson distances (e.g. L2 between Betti curves, bottleneck on the top bars). Interpretation rules: a match with lognormal mocks = consistent with M0 at the level of this test (NOT a confirmation of K3 x T2); a mismatch at small scales is expected from nonlinear gravity (lognormal is approximate) - say where the approximation is known to fail. Save plots and JSON. Keep runtime per script < 30 min (subsample).`,
    { label: 'tda:T1-cosmic-web', phase: 'TDA', model: 'sonnet', schema: EXP }),
  () => agent(`${RULES}
${ctx}
TDA T2 - CMB TOPOLOGY. Use the verified WMAP 9-yr ILC (or Planck) map with a galactic mask. Downgrade to nside 128 or 256 (state it, smoothing FWHM stated). Sublevel and superlevel set persistence of the temperature field on the sphere: build a GUDHI SimplexTree on the HEALPix pixel adjacency (vertices = unmasked pixels, edges/triangles from healpy neighbours; lower-star filtration by temperature), compute H0/H1 persistence and Betti curves vs threshold nu (in units of sigma). NULL: >= 100 Gaussian simulations (healpy.synfast, seeds) with the C_ell measured from the same masked map (or the CAMB M0 spectrum; say which), same beam/nside/mask/noise approximation. Report the chi2/p-value of the data Betti curves (and Euler characteristic curve) against the simulation ensemble.
Sensitivity / known-answer: inject a Kaiser-Stebbins-like cosmic-string step pattern (random straight string segments producing temperature steps of amplitude ~ 8 pi G mu v gamma T_CMB; document the simple model) into Gaussian sims at several G mu and find the smallest G mu at which the TDA statistic separates from the Gaussian ensemble at 95%: this is a TDA-based sensitivity estimate, to be compared (not claimed equal) with LeanMaster planckGmuBound = 1.5e-7 (a literature bound). M0 predicts Gaussian statistics (no defects). Interpretation: consistency is NOT evidence for K3 x T2.`,
    { label: 'tda:T2-cmb', phase: 'TDA', model: 'sonnet', schema: EXP }),
])).filter(Boolean)

phase('Skeptic')
const all = `EXPERIMENTS: ${JSON.stringify(exps)}
TDA: ${JSON.stringify(tdas)}
LEANMASTER: ${JSON.stringify(leanc)}`
const audits = (await parallel([
  () => agent(`${RULES}
SKEPTIC - STATISTICS AND DATA lens. You produced none of this; default to distrust.
${all}
For each item: rerun its main script (or a copy under /tmp) and confirm the headline numbers; check covariance handling (DR2 cov, Pantheon+ STAT+SYS, Hubble-flow cut consistent between models), degrees of freedom and the sigma conversion; for TDA check that nulls share N, footprint and selection with the data, that the lognormal bias was fitted on xi(r) only, that the CMB sims use the same mask/nside/smoothing, that the known-answer controls pass and that the string-injection sensitivity is not tuned to the answer. Decide stands / does not stand per item.`,
    { label: 'skeptic:statistics', phase: 'Skeptic', schema: AUDIT }),
  () => agent(`${RULES}
SKEPTIC - CLAIMS AND FRAMING lens. You produced none of this; default to distrust.
${all}
Check every claim against the FRAMING RULE: M0 is a hypothesis change, not a K3 x T2 derivation; no statement may say K3 x T2 fixes or predicts mu_sym, c4, Omega_Lambda, or the TDA outcomes. Check that pre-registered thresholds were read from ${PREREG} verbatim and applied as written; that LeanMaster quotes are verbatim (grep them); that 'consistent with M0' is never written as 'confirms the theory'; and that every ABSENT is reported. List violations with file:line.`,
    { label: 'skeptic:framing', phase: 'Skeptic', schema: AUDIT }),
])).filter(Boolean)

phase('Synthesize')
const report = await agent(`${RULES}
Write the REPORT as your final text (the harness blocks subagents from writing .md; the orchestrator will save it). Also write ${OUT}/ledger_reverse_zero.json and commit it.
Inputs: ${all}
AUDITS: ${JSON.stringify(audits)}
DATA CHECK: ${JSON.stringify(check)}
Sections:
1. Bottom line in plain words: the parameter count of M0 (0 by hypothesis change) vs M2 (2) and what the data say about the change (information criteria, E1 pre-registered verdict with DR2, E2 PTA).
2. Pre-registered results: P1 (DR2) with the registered rule applied verbatim; P2 status per LeanMaster Stream 6. Proposed ADDENDUM text for PRE_REGISTRATION.md (dated, appended, never editing the original rules).
3. TDA: cosmic web and CMB results with nulls and controls; what topology was found; what it does and does not say about M0 and about K3 x T2.
4. LeanMaster constraints (quotes) and what they close.
5. PROPOSED IMPROVEMENTS to reach a falsifiable zero-parameter theory, ranked: for each, the derivation or experiment needed, the observable it would produce (with units), which LeanMaster statement it would build on, a freeze-before-compare plan, and a realistic difficulty. Mark each Tier C. Include at least: (a) an observable from Stream 8 (T2 at (omega, omega), Kummer 8+16, Mukai classes) and whether ANY unit-carrying consequence exists; (b) flux-vacua statistics (Tripathy-Trivedi tadpole 24) as a landscape prediction; (c) TDA-based defect searches as null tests; (d) what would be needed to make the self-dual/dark-energy scale a derived number.
6. Everything ABSENT, disputed, or from memory.
Header: "Generated by workflow reverse-to-zero; exploratory (tier X) unless stated; M0 is a hypothesis change, not a derivation."`,
  { label: 'synthesize:report', phase: 'Synthesize' })

return { data, leanc, check, exps, tdas, audits, report }
