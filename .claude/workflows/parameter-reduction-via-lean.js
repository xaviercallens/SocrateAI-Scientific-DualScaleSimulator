
export const meta = {
  name: 'parameter-reduction-via-lean',
  description: 'Rerun simulations with falsification analysis; identify Lean-constrained parameters toward no-free-parameter theory',
  phases: [
    { title: 'Simulate', detail: 'Parameter sweeps with sensitivity analysis' },
    { title: 'Falsify', detail: 'Compute bounds vs experiment; identify excluded regions' },
    { title: 'Analyze Lean', detail: 'Extract Lean-derived constraints from LeanMaster' },
    { title: 'Synthesize', detail: 'Propose parameter reductions and derive new constraints' },
    { title: 'Report', detail: 'Generate improvement roadmap and verify claims' },
  ],
}

// Schema for simulation sweep results
const SIM_SWEEP_SCHEMA = {
  type: 'object',
  properties: {
    parameter_sweeps: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          param_name: { type: 'string' },
          nominal_value: { type: 'number' },
          tested_range: { type: 'string', description: '[min, max]' },
          sensitivity_metric: { type: 'number', description: 'output variance per 1% param change' },
          most_constraining_outcome: { type: 'string' },
        },
        required: ['param_name', 'nominal_value', 'sensitivity_metric'],
      },
    },
    key_findings: { type: 'array', items: { type: 'string' } },
  },
  required: ['parameter_sweeps', 'key_findings'],
}

// Schema for falsification analysis
const FALSIFY_SCHEMA = {
  type: 'object',
  properties: {
    observables_analyzed: { type: 'array', items: { type: 'string' } },
    falsification_results: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          observable: { type: 'string' },
          experimental_bound: { type: 'string', description: 'cited source + value ± error' },
          sim_prediction: { type: 'number' },
          sigma_from_bound: { type: 'number', description: '# standard deviations' },
          viable: { type: 'boolean', description: 'within 2-sigma?' },
          excluded_param_region: { type: 'string', description: 'if falsified, what param values ruled out' },
        },
        required: ['observable', 'experimental_bound', 'sim_prediction'],
      },
    },
    summary: {
      type: 'object',
      properties: {
        total_observables: { type: 'number' },
        viable_count: { type: 'number' },
        falsified_count: { type: 'number' },
        tightest_constraint: { type: 'string' },
      },
      required: ['total_observables', 'viable_count', 'falsified_count'],
    },
  },
  required: ['observables_analyzed', 'falsification_results', 'summary'],
}

// Schema for Lean constraint analysis
const LEAN_CONSTRAINTS_SCHEMA = {
  type: 'object',
  properties: {
    lean_derived_constants: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          constant_name: { type: 'string' },
          lean_module: { type: 'string', description: 'DualScaleStream2 or LeanMaster module path' },
          theorem_name: { type: 'string' },
          derived_value: { type: 'string', description: 'exact rational or computed value' },
          tier: { type: 'string', enum: ['A', 'B', 'L'] },
          fixes_current_free_param: { type: 'boolean' },
          current_code_value: { type: 'string' },
          mismatch: { type: 'string', description: 'if different from code' },
        },
        required: ['constant_name', 'lean_module', 'theorem_name', 'derived_value', 'tier'],
      },
    },
    newly_derivable: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          free_param_name: { type: 'string' },
          why_currently_free: { type: 'string' },
          lean_path_to_derive: { type: 'string', description: 'which theorems chain together' },
          expected_value: { type: 'string', description: 'prediction if derived' },
          confidence: { type: 'string', enum: ['high', 'medium', 'low'] },
        },
        required: ['free_param_name', 'why_currently_free', 'lean_path_to_derive'],
      },
    },
    external_conventions_reviewed: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          param_name: { type: 'string' },
          source: { type: 'string', description: 'Planck / DESI / PTA / other' },
          citation: { type: 'string', description: 'DOI or paper reference' },
          can_derive_from_first_principles: { type: 'boolean' },
        },
        required: ['param_name', 'source'],
      },
    },
    free_param_reduction_count: {
      type: 'object',
      properties: {
        before: { type: 'number' },
        after_lean_constraints: { type: 'number' },
        reduction_via_falsification: { type: 'number' },
        net_reduction: { type: 'number' },
      },
      required: ['before', 'after_lean_constraints', 'net_reduction'],
    },
  },
  required: ['lean_derived_constants', 'free_param_reduction_count'],
}

phase('Simulate')
log('Rerunning cosmological simulations with parameter sensitivity analysis...')

const simSweep = await agent(
  `Rerun workshopcosmo.py with parameter sweep. For each free parameter (a_pot, b_pot, mu_sym, lambda_sym, pta_suppression, c4_c0_ratio):
   1. Run simulation at nominal value
   2. Run at +10% and -10%
   3. Measure output sensitivity: how much do observables (w0, wa, screening_factor, boom_action) change per 1% parameter variation?
   4. Identify which parameters most constrain: (a) Cassini bound, (b) DESI consistency, (c) axiom audit passing rate
   
   Return parameter names, nominal values, tested ranges, sensitivity metrics, and the outcome each parameter constrains most.
   Be specific: "screening_factor sensitivity = 0.15 per 1% mu_sym change" not "screening_factor is sensitive".`,
  {phase: 'Simulate', schema: SIM_SWEEP_SCHEMA}
)

log(`Sensitivity analysis complete: ${simSweep.parameter_sweeps.length} parameters analyzed`)

phase('Falsify')
log('Computing falsification bounds against experimental constraints...')

const falsifyAnalysis = await agent(
  `Analyze falsification of current simulation against experiment:
   
   1. Extract observables from simulation_results.json: w0, wa, screening_factor, gamma-1 (Cassini), mapper_nodes, bounce_action
   2. For each observable, find experimental bound:
      - w0/wa: DESI 2024 DR1 (cite DOI 10.1038/s41586-024-07314-8 or similar)
      - screening: fifth-force bounds from Eöt-Wash / atomic equivalence principle tests
      - gamma-1: Cassini satellite bound (cite paper)
      - bounce_action: Coleman-De Luccia bounds (cite Kodama 1984 or similar)
      - TDA mapper: no direct experimental bound (theoretical, not falsifiable)
   
   3. For each: compute sigma distance = (sim - experiment) / experimental_error
   4. Mark "viable" if |sigma| < 2
   5. For falsified predictions, identify which parameter values are ruled out
   
   Return: observable name, experimental bound (with source), sim value, sigma distance, viable yes/no, excluded param region if falsified.`,
  {phase: 'Falsify', schema: FALSIFY_SCHEMA}
)

log(`Falsification analysis: ${falsifyAnalysis.summary.viable_count}/${falsifyAnalysis.summary.total_observables} viable`)
if (falsifyAnalysis.summary.falsified_count > 0) {
  log(`⚠ ${falsifyAnalysis.summary.falsified_count} observables falsified by experiment`)
}

phase('Analyze Lean')
log('Extracting Lean-derived constraints from LeanMaster v2.2.0...')

const leanConstraints = await agent(
  `Analyze which currently-free parameters can be derived from Lean 4 formalization:
   
   1. Query LeanMaster DualScaleStream2 (v2.2.0) and StringTheoryFormalization for theorems on:
      - Kummer K3 geometry (b2, chi, signature, Mukai pairing) → constrains lattice-derived ratios
      - Buscher T-duality (involutions on metric/B-field/dilaton) → if fixed, constrains field transformations
      - Tadpole cancellation (flux + brane charge) → constrains D7 D-brane count
      - Moonshine/M24 (R_BPS = 77/60, character decompositions) → constrains worldsheet VOA
   
   2. For each current free parameter (a_pot, b_pot, mu_sym, lambda_sym, etc.):
      - Is it already constrained by a Lean theorem? (Then move from tier C to tier A/B)
      - Could it be derived by extending current Lean modules? (Check LeanMaster downstream)
      - Does it conflict with a known Lean result?
   
   3. For external-convention parameters (Planck, DESI, PTA):
      - Can any be derived from first principles in Lean? (e.g., derive primordial power spectrum normalization from CFT)
   
   4. Count free parameters: before Lean constraints, after, after falsification, net reduction.
   
   Be specific: "mu_sym could be derived from cosmological perturbation theory × Lean screening theorem chain: [module A → B → C]" not "might be derived somehow".`,
  {phase: 'Analyze Lean', schema: LEAN_CONSTRAINTS_SCHEMA}
)

log(`Free parameter reduction: ${leanConstraints.free_param_reduction_count.before} → ${leanConstraints.free_param_reduction_count.after_lean_constraints}`)
log(`Lean-derived constants found: ${leanConstraints.lean_derived_constants.length}`)
if (leanConstraints.newly_derivable && leanConstraints.newly_derivable.length > 0) {
  log(`Candidates for derivation: ${leanConstraints.newly_derivable.map(d => d.free_param_name).join(', ')}`)
}

phase('Synthesize')
log('Synthesizing parameter reduction proposal and improvement roadmap...')

const synthesis = await agent(
  `Create a formal proposal to reduce free parameters toward "no free parameter" theory:
   
   **Input:** (1) sensitivity analysis showing which params most constrain observables, (2) falsification showing what experiment rules out, (3) Lean constraint analysis showing what's derivable
   
   **Output:** A structured proposal with:
   
   1. **Immediate reductions (no new theory needed):**
      - Parameters already Lean-constrained but coded as free (e.g., c112=77/60): fix in code to match Lean value
      - Parameters ruled out by falsification: propose new bounds and cite experiment
   
   2. **Short-term derivations (extend LeanMaster):**
      - Rank newly_derivable parameters by: (a) how much they constrain observables (high sensitivity), (b) confidence of derivation path
      - For top-3: write pseudo-code for the Lean theorem chain (don't implement, just describe the proof strategy)
   
   3. **Long-term reductions (new physics):**
      - Which parameters are pure model tuning (pta_suppression, c4_c0_ratio)? Propose experiments to measure them or reformulate the model to absorb them
   
   4. **No-free-parameter roadmap:**
      - Current count: N free parameters
      - After (1): N - M1 free
      - After (2): N - M1 - M2 free (achievable in 1-2 quarters if LeanMaster extensions merge)
      - After (3): 0 free (physics milestone; timeline TBD)
   
      Show the delta for each step and cite what changes in the paper/code.
   
   5. **Falsification-driven refinement:**
      - Where did experiment nearly falsify the theory (sigma closest to 2)? Propose tightening those parameters
      - Where is the theory most robust (sigma >> 2)? Those are the candidates for absorption into other parameters
   
   Format as a structured roadmap ready for presentation to the user.`,
  {phase: 'Synthesize'}
)

log('Proposal synthesized. Preparing verification report...')

phase('Report')
log('Generating final report and verification plan...')

const reportContent = await agent(
  `Write a comprehensive report: "Parameter Reduction Roadmap: Path to No-Free-Parameter String Cosmology"
   
   Structure:
   1. **Executive Summary:** Current free parameter count, proposed reductions per Lean constraints + falsification, net reduction, timeline
   2. **Sensitivity Analysis:** Which parameters most constrain which observables? (table: param → observable → sensitivity_metric)
   3. **Falsification Results:** Which predictions pass/fail experiment? (table: observable → bound → sim → sigma → viable?)
   4. **Lean Constraint Inventory:** (table: param → module → theorem → derived_value → code_value → mismatch?)
   5. **Derivation Candidates:** (table: param → reason_free → Lean_path → confidence)
   6. **Roadmap:** Three-phase timeline with deliverables per phase
   7. **Next Steps for Each Phase:**
      - Phase 1 (code fixes): PRs needed to align hard-coded values with Lean theorems
      - Phase 2 (LeanMaster extensions): Description of 3 new theorems to implement, priority ranking by constraint gain
      - Phase 3 (experimental): What experiments would most constrain remaining free parameters?
   
   Include caveats: any parameters where Lean derivation is blocked by missing formalization in LeanMaster, any experimental bounds with low confidence.
   
   Cite all experiments and Lean theorems by module/file/theorem name.`,
  {phase: 'Report'}
)

return {
  simulation_sweep: simSweep,
  falsification_analysis: falsifyAnalysis,
  lean_constraints: leanConstraints,
  proposal_synthesis: synthesis,
  final_report: reportContent,
  summary: {
    initial_free_params: leanConstraints.free_param_reduction_count.before,
    after_lean: leanConstraints.free_param_reduction_count.after_lean_constraints,
    after_falsification: leanConstraints.free_param_reduction_count.after_lean_constraints - falsifyAnalysis.summary.falsified_count,
    total_reduction: leanConstraints.free_param_reduction_count.net_reduction,
    tier_a_theorems_used: leanConstraints.lean_derived_constants.filter(c => c.tier === 'A').length,
    next_action: 'Review proposal and select phase for implementation (code fix / Lean extension / experimental)',
  }
}
