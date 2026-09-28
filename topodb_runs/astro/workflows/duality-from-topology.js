export const meta = {
  name: 'duality-from-topology',
  description: 'Discover duality from topology: validate a duality detector on Poincare duality (real, adjudicable, already in TopoDB), calibrate a scale-inversion detector against null-fitted controls, sweep the corpus, test whether the dual scale is DYNAMIC and LOCAL, then have skeptics try to kill every candidate.',
  whenToUse: 'UNSOLICITED AND UNRUN - created 2026-09-20 by a subagent acting on two mid-task messages of unverified provenance (the first explicitly tagged NOT USER INPUT); never authorised by the session owner, never executed. READ FIRST: its own probe design records that T-duality exchanges momentum and winding towers and does NOT act on H_*(T^d), so T-duality is homologically INVISIBLE - the premise of finding it in persistence data is unsound, and that is the principal negative result of this direction. What survives is the Poincare-duality measurement (tier B, 58 runs) with its two gates: the dimension n must be DECLARED, since inferring it from the data makes the test unable to fail (13 of 58 cases flip to a false pass); and a symmetry in an invariant does NOT certify the structure that would explain it (a singular orbifold and a deliberately wrong gluing are both palindromic over F3). Do not run this without a fresh pre-registration: the corpus holds 237 datasets and ~3800 statistics, so unregistered searching will find coincidences. ORIGINAL DESCRIPTION: when asking how a dual-scale / double-field / T-duality-like structure could be SEEN in persistent-homology data. Requires TopoDB at /mnt/disks/disk-socrateai-local-1/topodb/topodb.sqlite. Args: {only?: string[] phase ids, nullFits?: number (default 200), corpusLimit?: number (default 120), skeptics?: number (default 3)}',
  phases: [
    { title: 'Ground truth', detail: 'Poincare-duality detector adjudicated against ~150 labelled mathematics runs already in TopoDB' },
    { title: 'Calibrate', detail: 'scale-inversion detector vs null-FITTED controls; negatives must not fire' },
    { title: 'Sweep', detail: 'rank the real corpus by null-calibrated self-duality score' },
    { title: 'Dynamics', detail: 'does the dual scale drift with a physical label, and does it vary spatially?' },
    { title: 'Skeptics', detail: 'adversarial refutation of every surviving candidate' },
    { title: 'Ingest', detail: 'write findings to TopoDB, tier X, with the anti-overreach wording enforced' },
  ],
}

const DB = '/mnt/disks/disk-socrateai-local-1/topodb/topodb.sqlite'
const WT = '/mnt/disks/disk-socrateai-local-1/wt-topo-astro'
const OUT = `${WT}/topodb_runs/astro/duality`
const PY = '/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python'

const A = args || {}
const NULL_FITS = A.nullFits || 200
const CORPUS_LIMIT = A.corpusLimit || 120
const N_SKEPTICS = A.skeptics || 3
const ONLY = A.only || null
const want = (id) => !ONLY || ONLY.includes(id)

// ---------------------------------------------------------------- hard rules
const RULES = `
HARD RULES (these come from defects this project has already paid for; violating one invalidates the phase).

1. KNOWN ANSWERS BEFORE NEW ANSWERS. Nothing in a later phase is reported until the ground-truth phase passes.

2. NOTHING TRUE BY CONSTRUCTION. A duality detector that FITS a parameter will always find some value of it.
   Therefore: ANY score obtained by optimising a parameter over the data MUST be compared against the SAME
   OPTIMISATION RUN ON THE NULL. A null that is not itself fitted is not a null. This is the single most
   important rule in this workflow. The prior campaign measured a naive persistence-ratio rule firing on
   9 of 12 density-matched nulls; do not repeat it.

3. NO p-VALUE WITHOUT ITS NULL. TopoDB refuses one. Quote the rank floor 1/(n+1) (or 2/(n+1) two-sided)
   with every rank p-value, and say explicitly when the floor lies above your significance threshold, in
   which case the p-value carries NO information about magnitude and the verdict must rest on a separation
   statistic recorded WITHOUT a p-value.

4. TIER. Every number produced here is tier X (numerics) unless it is a Betti number of an exactly
   specified finite cell complex, which is tier B. Never write "proved", "certified" or "guaranteed".

5. THE WORDING THAT MATTERS MOST. Scale-inversion symmetry of a persistence diagram is a statement about
   the SELF-SIMILARITY OF THE DATA. It is NOT evidence of string-theoretic T-duality, not evidence of
   double field theory, and not evidence of a dual-scale cosmology. Any claim, finding or comment that
   equates the two is a defect and must be rewritten. The honest form is:
   "the diagram is (or is not) statistically self-inverse about scale c*, which is a necessary but very
   far from sufficient condition for any physical duality".

6. TOPODB WRITE DISCIPLINE (learned the hard way on 2026-09-20, when a concurrent agent replaced the
   database and run ids were reused, so an UPDATE meant for this branch's rows silently damaged another
   agent's records):
   - Open, write, close. Keep write transactions short; three agents share the file in WAL mode.
   - INSERT your own rows. NEVER UPDATE or DELETE a row you did not insert in this same script run.
   - Scope every query you intend to modify by branch = 'topo/astro' AND by the dataset id prefix you
     created. Never scope by run id alone: run ids are NOT stable across the shared database.
   - Before writing, re-read the row count; if the database has obviously been replaced, STOP and report.

7. Every number in a report, comment or JSON comes from a command run in this session. Quote it.
8. Work only inside ${WT}. Never touch main, proofs/, or another worktree. Never push.
9. Commit with explicit paths (never git add -A), message ending:
   Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
`

const SCHEMA_GATE = {
  type: 'object',
  properties: {
    passed: { type: 'boolean' },
    summary: { type: 'string' },
    numbers: { type: 'string', description: 'the measured numbers, quoted from the commands run' },
    script: { type: 'string', description: 'path of the committed script that produced them' },
    blocked: { type: 'string', description: 'empty unless something stopped the phase' },
  },
  required: ['passed', 'summary', 'numbers', 'script'],
}

const SCHEMA_CANDIDATES = {
  type: 'object',
  properties: {
    candidates: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          dataset_id: { type: 'string' },
          run_id: { type: 'integer' },
          c_star: { type: 'number', description: 'fitted self-dual scale, in the filtration units of that run' },
          score: { type: 'number', description: 'self-duality score, lower = more self-inverse' },
          null_fitted_mean: { type: 'number' },
          null_fitted_sd: { type: 'number' },
          separation_sigma: { type: 'number' },
          rank_p: { type: 'number' },
          rank_floor: { type: 'number' },
          verdict: { type: 'string', enum: ['candidate', 'null', 'artefact', 'inconclusive'] },
        },
        required: ['dataset_id', 'run_id', 'c_star', 'score', 'separation_sigma', 'verdict'],
      },
    },
    summary: { type: 'string' },
    script: { type: 'string' },
  },
  required: ['candidates', 'summary', 'script'],
}

log('duality-from-topology: a duality is only interesting if the detector can also say NO.')

// =====================================================================
// PHASE 1 — GROUND TRUTH: a duality the database can adjudicate
// =====================================================================
let ground = null
if (want('ground')) {
  phase('Ground truth')
  ground = await agent(`${RULES}

PHASE 1 — GROUND TRUTH. Build a duality detector and prove it works on a duality that is REAL, is already
in TopoDB, and has genuine positives AND genuine negatives: POINCARE DUALITY, b_k = b_{n-k}.

Why this phase exists: before asking whether persistence data shows a T-duality-like scale exchange, the
detector must be shown to detect a duality that is present and to REJECT one that is absent. TopoDB
already contains the adjudication set, computed by another branch as exact cell-complex homology.

THE DATA. ${DB}, tables dataset/run/betti, domain = 'mathematics', dataset ids under
'mathematics/k3t2/partA/'. About 150 runs, each a Betti vector over a stated coefficient field
(run.coeff_field in {2,3,5}). Read them; do not recompute them.

THE ADJUDICATION SET. It has ALREADY BEEN MEASURED against the database (2026-09-20, 58 runs), so these
are not guesses - they are the numbers your script must reproduce exactly. Reproduce them; do not retune.

  CRITICAL RULE ABOUT n, and it is the whole reason this phase exists.
  n MUST be the DECLARED dimension of the space, parsed from its name (T^2->2, T^4->4, T^6->6, S^2->2,
  S^3->3, RP^2->2, RP^3->3, resolved Kummer K3->4, the orbifold U+16 cones->4, anything "x_t_2" or
  "x_klein"->6). n MUST NOT be inferred as "the top dimension with a non-zero Betti number".
  Inferring n from the data makes the test UNABLE TO FAIL: a non-orientable space over F3 has b_n = 0, so
  inference silently truncates the vector to a shorter one that is then trivially palindromic. Measured:
  with the inferring rule, 13 of 58 cases flip to a false HOLDS - the Klein bottle over F3 (1,1,0) reads
  as (1,1), RP^2 over F3 (1,0,0) reads as (1), K3 x Klein over F3 (1,1,22,22,1,1,0) reads as
  (1,1,22,22,1,1). This is a textbook "true by construction" defect (HARD RULE 2); your script must use
  the declared n and your report must state this measured comparison.

  With the DECLARED n, the measured counts are:
    ORIENTABLE            12/12 palindromic over F2, 12/12 over F3, 12/12 over F5
                          (T^2, T^4, T^6, S^2, S^3, RP^3 - RP^3 IS orientable - and resolved Kummer K3,
                           whose b_2 = 22 is the correct K3 answer)
    NON-ORIENTABLE         4/4 palindromic over F2,  0/4 over F3,  0/4 over F5
                          (Klein bottle x2, RP^2, K3 x Klein) - the field-dependence, measured
    SINGULAR (orbifolds)   0/2 over F2,  1/2 over F3,  1/2 over F5
    DELIBERATELY WRONG     0/1 over F2,  1/1 over F3

  THE SECOND LESSON, and it governs every later phase: the last two rows are FALSE POSITIVES of the
  palindrome test. a4_control_u_16_cones_orbifold_n_6 over F3 is (1,0,6,0,1) and IS palindromic although
  it is a singular space and not a manifold; a5_wrong_gluing_class_phi2_0_at_16_points_x_t_2 over F3 is
  (1,2,23,44,23,2,1) and IS palindromic although the gluing is deliberately wrong. Over F2 both fail
  (D = 10 and D = 30).
  So: the palindrome test has ZERO false negatives and REAL false positives. A SYMMETRY IN AN INVARIANT
  DOES NOT CERTIFY THE STRUCTURE THAT WOULD EXPLAIN IT. Carry that sentence into Phases 2-5, where the
  invariant is a persistence diagram and the tempting conclusion is a physical duality.

WHAT TO BUILD. ${OUT}/detect_poincare.py:
  - reads the Betti vectors from TopoDB (read-only connection),
  - takes n from the DECLARED dimension (see the critical rule above) and reports the symmetry defect
    D = sum_{k=0}^{n} |b_k - b_{n-k}| and the normalised defect D / sum_k b_k,
  - classifies HOLDS iff D == 0,
  - ALSO computes D under the wrong "inferred n" rule, so the 13/58 flip is reproduced and recorded as a
    negative control on the detector itself,
  - prints the counts per class (orientable / non-orientable / singular / wrong) per field.

THE GATE, declared before you run:
  - orientable: 12/12 HOLDS at every field - any failure here is a false negative and fails the phase;
  - non-orientable: 4/4 HOLDS over F2 and 0/4 over F3 and F5 - the field split must be reproduced exactly;
  - the inferred-n negative control must flip exactly 13 of 58 cases;
  - the singular/wrong false positives (1/2 and 1/1 over F3) must be REPORTED AS FALSE POSITIVES, not
    made to pass. Do NOT add a rule to catch them: no rule based on Betti numbers alone can, and
    pretending otherwise would be the defect this phase exists to prevent.
If any measured count differs from the table, passed=false and name the case - do NOT adjust the
expectation table to match the output. The table is mathematics plus an already-run measurement, not a
tuning knob.

ALSO REPORT, because it is the conceptual bridge to the rest of the workflow: the detector's answer is
FIELD-DEPENDENT. The same space is self-dual over F2 and not over F3. State that explicitly - a duality
is a property of (space, coefficient system), never of the space alone. That is the first honest lesson
this workflow has to teach, and it is measured, not asserted.

Write ${OUT}/PHASE1_ground_truth.json with every case, its expectation, its measurement and the
confusion matrix. Commit detect_poincare.py and the JSON.`, { schema: SCHEMA_GATE, effort: 'high' })

  if (!ground || !ground.passed) {
    log(`GROUND TRUTH FAILED - stopping. ${ground ? ground.summary : 'agent returned nothing'}`)
    return { stopped_at: 'ground', ground }
  }
  log(`Ground truth PASSED: ${ground.numbers}`)
}

// =====================================================================
// PHASE 2 — CALIBRATE the scale-inversion (T-duality-like) detector
// =====================================================================
let calib = null
if (want('calibrate')) {
  phase('Calibrate')
  calib = await agent(`${RULES}

PHASE 2 — CALIBRATE a scale-inversion detector. Phase 1 passed: ${ground ? ground.numbers : 'n/a'}.

THE IDEA, stated precisely so it can be falsified. T-duality exchanges a radius with its inverse,
R <-> alpha'/R. The persistent-homology analogue is an involution on the FILTRATION SCALE. For a
persistence diagram D = {(b_i, d_i)} define, for a scale constant c > 0,

    sigma_c(b, d) = (c/d, c/b).

sigma_c is an involution (sigma_c o sigma_c = identity), it preserves b < d, and it reverses the order of
the scale axis, so a long bar at small scale maps to a long bar at large scale.

THE DOMAIN OF sigma_c, which is NOT all of the diagram, and is MEASURED, not assumed. c/b is undefined at
b = 0 and c/d is 0 at d = infinity, so sigma_c is undefined on exactly the bars that usually carry the most
persistence. Counted over the whole TopoDB bar table on 2026-09-20:
    dim 0: 2712 bars, 2687 of them (99.1%) have birth <= 0   -> H_0 IS EXCLUDED BY CONSTRUCTION
    dim 1: 3147 bars,    0 (0.0%) have birth <= 0, 139 essential
    dim 2: 1511 bars,    0 (0.0%) have birth <= 0,  75 essential
    dim 3:   24 bars,    0 (0.0%) have birth <= 0,   0 essential
So the detector runs on FINITE BARS WITH birth > 0, which on this corpus means H_1, H_2 and H_3 only.
State that H_0 is excluded by construction and not by convenience; an alpha-complex H_0 bar is always born
at 0, so no scale inversion can act on it at all. Usable bars after the restriction:
H_1 biology 2655 / synthetic 325 / astro 28; H_2 biology 1273 / synthetic 160; H_3 synthetic 24.
NOTE the consequence for Phase 3 and say it there: the astro slice has only ~28 usable H_1 bars, so the
sweep is biology- and synthetic-dominated and the astro result will be weak whatever it shows.

A diagram is SELF-INVERSE ABOUT c if D and sigma_c(D) agree. Define the score

    S(c) = W_1( D, sigma_c(D) ) / totalPersistence(D)          (1-Wasserstein, or bottleneck; state which)

and c* = argmin_c S(c), S* = S(c*). Work in log-scale: the natural parameter is u = log c, and sigma acts
as u-reflection, so grid-search u uniformly and say so.

THE TRAP, AND THE WHOLE POINT OF THIS PHASE. c* is FITTED. Minimising over c will always return some c
and some S*, on ANY diagram, including pure noise. S* alone is therefore meaningless. The ONLY admissible
statistic is S* compared against ${NULL_FITS} nulls on which THE SAME MINIMISATION OVER c WAS PERFORMED.
Implement it that way or the phase fails. (See HARD RULE 2.)

POSITIVE CONTROLS - build them so the answer is known in advance:
  P1 an exactly self-inverse diagram: take any diagram D0, form D0 union sigma_{c0}(D0) for a chosen c0.
     The detector must recover c* = c0 to within the grid resolution and give S* ~ 0.
  P2 a log-periodic / self-similar point cloud: points at radii r_j = r_0 * lambda^j (a geometric
     cascade, lambda stated). Its persistence is log-periodic, so it is self-inverse about the geometric
     centre of the cascade. Predict c* BEFORE measuring and say whether it was recovered.
  P3 a genuinely dual pair from the ACTUAL physics, if you can build one cheaply: a flat torus of radius
     R and one of radius 1/R sampled at matched density. Their H_1 death scales should exchange. Report
     whether the detector sees the exchange. If it cannot be built within budget, say NOT ATTEMPTED.

NEGATIVE CONTROLS - the detector MUST NOT fire:
  N1 a Poisson point cloud (no scale structure),
  N2 a single-scale cloud: Gaussian blobs of one characteristic size,
  N3 the SHUFFLE null: the observed diagram with births/deaths randomly re-paired, which destroys the
     pairing but keeps both marginals. This is the strictest null and the one to report.

THE GATE, declared now, before measuring:
  - P1 recovers c0 to within one grid step AND S* < 0.05;
  - every negative control has rank p > 0.2 against its own fitted null (i.e. it does NOT fire);
  - the false-positive rate over 20 independent Poisson draws is <= 0.10 at a nominal 0.05 threshold.
If the false-positive rate exceeds that, passed=false and REPORT THE RATE - do not retune the score to
pass. A detector that fires on noise is the finding.

Build ${OUT}/detect_scale_inversion.py (import gudhi; use its wasserstein/bottleneck, state which and
its version). Write ${OUT}/PHASE2_calibration.json with every control, its prediction, its measurement,
the fitted-null distribution and the false-positive rate. Commit both.`, { schema: SCHEMA_GATE, effort: 'high' })

  if (!calib || !calib.passed) {
    log(`CALIBRATION FAILED - the detector is not usable on real data. ${calib ? calib.summary : 'no result'}`)
    return { stopped_at: 'calibrate', ground, calib }
  }
  log(`Calibration PASSED: ${calib.numbers}`)
}

// =====================================================================
// PHASE 3 — SWEEP the corpus (multi-modal: each lens sees a different slice)
// =====================================================================
let sweep = null
if (want('sweep')) {
  phase('Sweep')
  const LENSES = [
    { id: 'astro', q: "domain='astro'", note: 'galaxy point clouds and sky maps; filtration units Mpc/h or sigma of the map' },
    { id: 'bio', q: "domain='biology'", note: 'the largest diagram corpus in the DB (~180 runs with bars); a duality here would be about self-similarity of molecular structure, with NO physical-duality reading whatsoever' },
    { id: 'qfluid', q: "domain='quantum_fluid'", note: 'vortex-lattice data; the Re6Zr lineage where the spacing statistic originally worked' },
    { id: 'synth', q: "domain='synthetic'", note: 'known shapes; these double as an in-corpus control - a circle or a torus has NO scale duality, so firing here is a red flag' },
  ]
  const results = await parallel(LENSES.map(L => () => agent(`${RULES}

PHASE 3 — SWEEP, lens "${L.id}". Calibration passed: ${calib ? calib.numbers : 'n/a'}.

Apply the CALIBRATED detector from ${OUT}/detect_scale_inversion.py (do not modify it; import it) to every
run in TopoDB with at least 10 stored bars in the slice ${L.q}. ${L.note}

Cap at ${CORPUS_LIMIT} runs for this lens. If the slice has more, take the ${CORPUS_LIMIT} with the most
stored bars and LOG HOW MANY YOU DROPPED - a silent cap reads as "we covered everything" when it did not.

IMPORTANT about the substrate, state it in your report: TopoDB stores only the TOP ~25-50 bars per run
(ranked by persistence), not the full diagram. A truncated diagram is biased towards long bars, and
scale-inversion maps long-at-small-scale to long-at-large-scale, so truncation can CREATE or DESTROY
apparent self-duality. Quantify this: re-run the detector on a few runs using only the top 10 vs the top
25 bars and report how much c* and S* move. If they move a lot, the sweep is exploratory only and you
must say so.

For EVERY run report: c*, S*, the fitted-null mean and sd over ${NULL_FITS} shuffle nulls (HARD RULE 2:
the null is fitted too), the separation in null-sigma, the rank p with its floor, and a verdict. Apply
Bonferroni over the number of runs you actually tested and state the threshold. If the rank floor exceeds
that threshold, say so and rest the verdict on the separation, not the p-value.

Return the ranked candidates. A lens that finds nothing returns an empty list and that is a RESULT, not a
failure - report it as a null.`, { schema: SCHEMA_CANDIDATES, phase: 'Sweep', effort: 'high' })))

  const cands = results.filter(Boolean).flatMap(r => r.candidates || [])
    .filter(c => c.verdict === 'candidate')
  sweep = { byLens: results.filter(Boolean), candidates: cands }
  log(`Sweep: ${cands.length} candidate(s) survived their own fitted null across ${results.filter(Boolean).length} lenses.`)
}

// =====================================================================
// PHASE 4 — IS IT DYNAMIC? IS IT LOCAL? (the user's two actual questions)
// =====================================================================
let dynamics = null
if (want('dynamics')) {
  phase('Dynamics')
  const PROBES = [
    {
      id: 'dynamic-camels',
      p: `Does the dual scale DRIFT with a physical parameter - i.e. is the duality DYNAMIC rather than a fixed property?

The cleanest handle in this project: CAMELS IllustrisTNG HI maps. ${WT}/topodb_runs/astro/run_camels.py
computed cubical persistence for 200 maps, ONE PER SIMULATION (the 750 test maps carry only 555 distinct
label rows, so per-map statistics are clustered - respect that, sample one map per simulation, N_eff = 200).
Labels: Omega_m in [0.1006, 0.4998], sigma_8 in [0.6002, 0.9998], plus four feedback parameters
A_SN1/A_AGN1/A_SN2/A_AGN2 in params_6 which VARY SIMULTANEOUSLY and are a confounder you must name.

Fit c* per map with the calibrated detector, then test Spearman(c*, Omega_m) and Spearman(c*, sigma_8)
with a PERMUTATION null, and the shuffled-label control. Pre-declare the test count and Bonferroni.

Interpretation discipline: c* has units of the filtration. If the maps were globally renormalised (they
were: renorm(log1p(minmax(HI))) with GLOBAL constants) then c* is comparable across maps - say so and say
why. A drift of c* with Omega_m would mean the self-dual scale is set by the cosmology, which is the
closest thing to a "dynamic dual scale" this data can support. It would STILL not be T-duality
(HARD RULE 5).`,
    },
    {
      id: 'local-map',
      p: `Is the dual scale LOCAL - does c* vary across the sky or across the volume, rather than being one
global number? This is the "local dual scale" question, and it is the one that can actually be mapped.

Do it on data already staged and already ingested by branch topo/astro:
  (a) SKY: Planck SMICA and WMAP ILC at nside 128 with their masks (see ${WT}/topodb_runs/astro/run_skymaps.py,
      which uses the FIXED complex in topodb_runs/astro/lib/cmb_topology.py - do NOT import
      audit/reverse_zero/E5-cmb-tda/cmb_tda.py, it is defective). Partition the unmasked sky into
      nside-4 or nside-8 superpixels, compute the lower-star diagram per patch, fit c* per patch, and
      produce a MAP of c*(patch). Then ask the only question that matters: is the spatial variation of
      c* larger than the variation across ${NULL_FITS} Gaussian realisations analysed identically? If it
      is not, the map is noise and must be reported as noise.
  (b) VOLUME: DESI BGS_BRIGHT-21.5 NGC, split into redshift shells, alpha complex per shell (mind
      AMENDMENT A3 in ${WT}/topodb_runs/astro/PRE_DECLARED_STATISTICS.md - the alpha cap must be
      SCALE-FREE, set from each shell's own H0-death median, or a sparse shell saturates and you will
      measure the cap instead of the sky). Fit c* per shell and report whether c* tracks the mean
      inter-galaxy separation of that shell - which would be the BORING explanation and is the one to
      rule out first.

Deliverable: a map/table of c* with its null band, and an explicit statement of which of the two
explanations (real local structure vs trivial scaling with local density) the data supports.`,
    },
    {
      id: 'exchange',
      p: `Test scale-reversed Poincare duality on diagrams - and FIRST, correct a physics error that the
first draft of this workflow contained, because the correction is itself the most useful thing in this probe.

THE ERROR, stated plainly so nobody re-derives it: it is tempting to say that the topological signature of
T-duality is a degree exchange, H_k at scale s matching H_{n-k} at scale c/s. THAT IS WRONG. T-duality
exchanges momentum and winding modes - the spectrum-level statement is that the Kaluza-Klein tower
(mass proportional to n/R) swaps with the winding tower (mass proportional to mR). It does NOT act on the
homology of the torus: H_*(T^d) is INVARIANT under T-duality. That invariance is precisely WHY T-duality
is hard to see homologically at all, and it is the single most important negative result this workflow can
report. Write it down as such.

WHAT IS ACTUALLY TESTABLE, and it is narrower:
  (a) The reciprocal-radius pair. A flat torus sampled at radius R has its H_1 death scales set by R; one
      at radius alpha'/R has them set by alpha'/R. So a PAIR of datasets at reciprocal radii must show
      reciprocal H_1 death scales. This tests THE DETECTOR, not the physics - it is the P3 control from
      Phase 2, promoted here. Build both tori at matched sampling density, predict the reciprocal relation
      before measuring, and report whether the detector recovers it.
  (b) Scale-reversed Poincare duality on diagrams: S_{k,n-k}(c) = W_1( D_k , sigma_c(D_{n-k}) ) normalised,
      c fitted, compared against the SAME fit on the shuffle null (HARD RULE 2). This is a property OF THE
      SPACE, not of a duality transformation on it - label it that way. It is included only because Phase 1
      established that the palindrome version of the same idea is necessary-but-not-sufficient, and it is
      worth knowing whether the diagram-level version inherits that weakness. Run it on the Step-0 torus
      (b = (1,2,1), Poincare self-dual, the natural positive) and sphere (b = (1,0,1), b_1 = 0, the natural
      negative for the k=1 channel), then on the Phase 3 candidates.

Report (a) and (b) separately and never merge them. The expected answer to (b) on real data is NO, and a
null is a perfectly good result. Do not, anywhere, describe either as evidence of T-duality.`,
    },
  ]
  const res = await parallel(PROBES.map(P => () => agent(`${RULES}

PHASE 4 — "${P.id}". Detector calibrated in Phase 2 (${calib ? calib.numbers : 'n/a'}).
Phase 3 candidates: ${sweep ? JSON.stringify(sweep.candidates.slice(0, 20)) : 'none'}

${P.p}

Write ${OUT}/PHASE4_${P.id}.json and commit the script that produced it.`,
    { schema: SCHEMA_GATE, phase: 'Dynamics', effort: 'high' })))
  dynamics = res.filter(Boolean)
  log(`Dynamics: ${dynamics.filter(d => d.passed).length}/${dynamics.length} probes returned a usable measurement.`)
}

// =====================================================================
// PHASE 5 — SKEPTICS: try to kill every surviving candidate
// =====================================================================
let verdicts = null
if (want('skeptics') && sweep && !sweep.candidates.length) {
  log('Skeptics SKIPPED: the sweep produced zero candidates, so there was nothing to challenge. That is a null result, not an absence of scrutiny.')
}
if (want('skeptics') && sweep && sweep.candidates.length) {
  phase('Skeptics')
  const LENSES = [
    'TRUNCATION AND SELECTION: TopoDB stores only the top ~25-50 bars. Show that this candidate\'s self-duality is an artefact of which bars were kept. Re-derive with a different truncation and with the full diagram if it can be recomputed.',
    'FITTING AND MULTIPLICITY: the scale c* was fitted. Show that the fitted null was not actually fitted the same way, or that the multiplicity correction does not cover the search over c, over runs, and over lenses. Check whether the rank floor exceeds the threshold, in which case the p-value cannot support the claim at all.',
    'TRIVIAL EXPLANATION: show that c* is just a restatement of a scale the dataset already has - the mean inter-point separation, the pixel size, the beam, the box size, the survey depth, or the filtration cap. If c* is within a small factor of any of these, the "duality" is a unit conversion.',
  ].slice(0, N_SKEPTICS)

  verdicts = await pipeline(sweep.candidates.slice(0, 24),
    (c) => agent(`${RULES}

PHASE 5 — SKEPTIC. Your job is to REFUTE this candidate, not to confirm it. Default to refuted=true when
uncertain. A candidate that survives three honest attempts to kill it is worth one sentence in a report;
a candidate that dies here saves this project from publishing an artefact.

CANDIDATE: ${JSON.stringify(c)}

Apply THIS lens and no other:
${LENSES[0]}

Re-run whatever you need from ${OUT}. Report refuted (true/false), the specific measurement that decides
it, and - if refuted - the one-sentence explanation of what the apparent duality actually was.`,
      { schema: { type: 'object', properties: { refuted: { type: 'boolean' }, reason: { type: 'string' }, measurement: { type: 'string' } }, required: ['refuted', 'reason', 'measurement'] }, phase: 'Skeptics', effort: 'high' }),
    (first, c) => parallel(LENSES.slice(1).map(L => () => agent(`${RULES}

PHASE 5 — SKEPTIC (second wave). CANDIDATE: ${JSON.stringify(c)}
A previous skeptic said: refuted=${first && first.refuted}, ${first && first.reason}

Apply THIS lens and no other:
${L}

Default to refuted=true if uncertain.`,
      { schema: { type: 'object', properties: { refuted: { type: 'boolean' }, reason: { type: 'string' }, measurement: { type: 'string' } }, required: ['refuted', 'reason', 'measurement'] }, phase: 'Skeptics', effort: 'high' })))
        .then(rest => {
          const all = [first, ...rest].filter(Boolean)
          const kills = all.filter(v => v.refuted).length
          return { candidate: c, votes: all, refuted: kills >= Math.ceil(all.length / 2), kills, n: all.length }
        }))

  const survivors = (verdicts || []).filter(Boolean).filter(v => !v.refuted)
  log(`Skeptics: ${survivors.length} of ${sweep.candidates.length} candidate(s) survived ${N_SKEPTICS} adversarial lenses.`)
}

// =====================================================================
// PHASE 6 — INGEST, with the anti-overreach wording enforced by a reviewer
// =====================================================================
let ingest = null
if (want('ingest')) {
  phase('Ingest')
  const survivors = (verdicts || []).filter(Boolean).filter(v => !v.refuted).map(v => v.candidate)
  ingest = await agent(`${RULES}

PHASE 6 — INGEST and MAP. Write the campaign into TopoDB and produce the cartography the user asked for.

WHAT HAPPENED:
  ground truth  : ${ground ? ground.numbers : 'not run'}
  calibration   : ${calib ? calib.numbers : 'not run'}
  sweep         : ${sweep ? sweep.candidates.length + ' candidates before skeptics' : 'not run'}
  dynamics      : ${dynamics ? dynamics.map(d => d.summary).join(' | ') : 'not run'}
  survivors     : ${JSON.stringify(survivors)}

DO THIS:
1. Ingest into TopoDB under dataset ids prefixed 'duality/' (a NEW namespace you create - HARD RULE 6,
   insert only, never update another agent's rows, never scope by run id alone). Record:
   - the Poincare ground-truth run with its confusion matrix as controls of kind 'known_answer', and the
     field-dependence result (self-dual over F2, not over F3) as its own finding - that one is tier B,
     because it is exact homology of an exactly specified finite complex, not numerics;
   - the calibration as controls of kind 'negative' (the false-positive rate is the number that matters);
   - one run per swept lens, with the fitted-null statistics, and p-values ONLY where a fitted null
     exists;
   - one honest finding per phase. If the answer is "no duality was found", the finding says exactly
     that, with verdict 'null'. That is the expected outcome and it is a good one.

2. Write ${OUT}/DUALITY_MAP.md - the cartography. For each domain: what scale structure exists, whether
   it survived its own fitted null, and what the trivial explanation was when it did not. Include the
   c*-across-the-sky and c*-across-redshift maps from Phase 4 with their null bands.

3. Write ${OUT}/HYPOTHESES.md - the thought-experiment ledger, in the Einstein register the user asked
   for: each hypothesis stated so that it could FAIL, the topological observable that would decide it,
   the control that would expose the trivial explanation, and - stated honestly - whether the data
   available to this project can decide it at all. Seed it with these four and add your own:
     H1 A local dual scale exists: c*(x) varies across the sky by more than the Gaussian null band.
        Decided by: the Phase 4 c* map. Trivial explanation to rule out: beam and mask geometry.
     H2 The dual scale is dynamic: c* tracks Omega_m across the CAMELS suite.
        Decided by: Spearman with a permutation null. Confounder: the four feedback parameters.
     H3 The duality is an EXCHANGE, not merely an inversion: H_k at s matches H_{n-k} at c/s.
        Decided by: the cross-degree score. Expected answer: no.
     H4 Duality is a property of (space, coefficients), not of space.
        ALREADY DECIDED, tier B, by Phase 1: the Klein bottle is Poincare self-dual over F2 and is not
        over F3. This is the one duality result in this workflow that is actually established, and it
        should be stated as the anchor the speculative ones are measured against.

4. Then spawn NOTHING further yourself; instead return, in 'summary', the three sentences you would put
   in front of a skeptical reader.

WORDING REVIEW BEFORE YOU COMMIT: re-read every claim you wrote against HARD RULE 5. If any sentence lets
a reader conclude that this project has seen T-duality, double field theory, or a dual-scale cosmology in
data, rewrite it. Quote in your report the sentence you were most tempted to overstate, and what you
replaced it with.

Commit ${OUT}/ with explicit paths.`, { schema: SCHEMA_GATE, effort: 'high' })
}

return {
  ground, calibration: calib,
  sweep: sweep ? { n_candidates: sweep.candidates.length, candidates: sweep.candidates } : null,
  dynamics, skeptics: verdicts, ingest,
  note: 'Scale-inversion symmetry of a persistence diagram is a statement about self-similarity of the data. It is not evidence of T-duality, double field theory, or a dual-scale cosmology. The only duality this workflow ESTABLISHES is Poincare duality, tier B, from exact homology already in TopoDB.',
}
