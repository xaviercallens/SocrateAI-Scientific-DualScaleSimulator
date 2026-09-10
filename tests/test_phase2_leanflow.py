"""
=============================================================================
LeanFlow Phase 2 Test Suite
=============================================================================
Validates all Phase 2 components:
- String theory dataset ingestion (CICY, Kreuzer-Skarke)
- Declarative Invariant Specification DSL (@guardrail)
- Stiff solver with modular domain folding & PyTorch interoperability
- Working Group 1: cymetric Kahler cone positive-definiteness & Monge-Ampère
- Working Group 2: Swampland MCMC filter & Lean 4 verification gate
- Working Group 3: Topological defect extraction & Ramond-Ramond neutrality
=============================================================================
"""

import pytest
import numpy as np

try:
    import torch
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False

import leanflow as lf
from leanflow.datasets.cicy import parse_cicy_matrix, CICYConfiguration
from leanflow.datasets.kreuzer_skarke import parse_ks_polytope, KreuzerSkarkePolytope
from leanflow.core.projections import (
    metric_positivity_projection,
    modular_domain_fold,
    weak_energy_condition_lock,
    kahler_cone_positivity_projection,
    tadpole_budget_check,
    FRICKE_Y
)
from leanflow.bridge.lean_ipc import LeanVerificationClient


class TestStringDatasets:
    def test_load_all_canonical_cicy(self):
        cicy_names = ["quintic", "cicy_7887", "tian_yau", "bicubic", "k3_x_t2"]
        for name in cicy_names:
            cicy = lf.load_cicy(name)
            assert isinstance(cicy, CICYConfiguration)
            assert cicy.is_calabi_yau is True
            assert cicy.euler_characteristic == 2 * (cicy.h11 - cicy.h21)
            assert cicy.tadpole_bound == abs(cicy.euler_characteristic) / 24.0

    def test_parse_custom_cicy(self):
        # Tian-Yau 3-fold: P3 x P3 with degrees [[3, 0, 1], [0, 3, 1]]
        matrix = np.array([
            [3, 0, 1],
            [0, 3, 1]
        ])
        ambient = [3, 3]
        cicy = parse_cicy_matrix(matrix, ambient, name="TianYau_Test", h11=14, h21=23)
        assert cicy.is_calabi_yau is True
        assert cicy.ambient_dimension == 6
        assert cicy.manifold_dimension == 3
        assert cicy.euler_characteristic == -18

    def test_load_all_canonical_ks(self):
        ks_names = ["ks_quintic", "ks_bicubic", "ks_k3_surface"]
        for name in ks_names:
            poly = lf.load_ks(name)
            assert isinstance(poly, KreuzerSkarkePolytope)
            assert poly.is_reflexive is True
            assert poly.euler_characteristic == 2 * (poly.h11 - poly.h21)
            assert poly.vertices.shape[1] == 4

    def test_parse_custom_ks_polytope(self):
        verts = np.array([
            [-1, -1, -1, -1],
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ])
        poly = parse_ks_polytope(verts, h11=1, h21=101, ks_id="Test_Reflexive")
        assert poly.is_reflexive is True
        assert poly.euler_characteristic == -200


class TestCoreProjectionsAndDSL:
    def test_modular_domain_fold(self):
        # Test point outside fundamental domain: tau = 0.2 + 0.4i (|tau|^2 = 0.2 < 1)
        x_in, y_in = 0.2, 0.4
        x_out, y_out, folds = modular_domain_fold(x_in, y_in)
        assert folds > 0
        # In fundamental domain: |x| <= 0.5, |tau| >= 1.0 - 1e-6, y >= Fricke
        assert abs(x_out) <= 0.5001
        assert (x_out**2 + y_out**2) >= 0.999
        assert y_out >= FRICKE_Y

    def test_weak_energy_condition_lock(self):
        # Violating WEC: rho = 1.0, p = -2.0 (rho + p = -1.0 < 0, w = -2.0)
        rho, p_proj, w_proj = weak_energy_condition_lock(1.0, -2.0, min_sum=0.0)
        assert rho + p_proj >= 0.0
        assert w_proj >= -1.0

    def test_guardrail_decorator_telemetry(self):
        @lf.guardrail(
            theorem="SocrateAI.Cosmology.wec_kinetic_identity",
            invariants=["metric_positivity", "wec_bound"]
        )
        def dummy_eom(t, state):
            return np.array([state[2], state[3], -0.1, -0.2])

        # Test state with dangerously small y (dilaton / scale)
        state = np.array([0.0, 1e-7, 1.0, 0.5])
        out = dummy_eom(0.0, state)
        assert dummy_eom.telemetry.total_calls >= 1
        assert dummy_eom.telemetry.projections_applied >= 1
        assert dummy_eom.telemetry.invariants_checked >= 2
        assert len(out) == 4

    def test_invariant_registry(self):
        thm = lf.InvariantRegistry.get_theorem("SocrateAI.Cosmology.wec_kinetic_identity")
        assert thm is not None
        assert "wec_bound" in thm["invariants"]


class TestStiffSolver:
    def test_solver_bdf_moduli(self):
        def simple_linear_moduli(t, y):
            # Damped oscillator in moduli space
            return np.array([y[2], y[3], -0.5 * y[0] - 0.1 * y[2], -0.5 * (y[1] - 1.0) - 0.1 * y[3]])

        solver = lf.Solver(method="BDF")
        y0 = [0.1, 0.5, 0.0, 0.0]
        sol = solver.solve(simple_linear_moduli, y0, (0.0, 5.0), t_eval=np.linspace(0.0, 5.0, 50))
        assert sol.success is True
        assert sol.y.shape == (4, 50)
        assert sol.telemetry.total_steps == 50
        assert sol.telemetry.step_latency_ms >= 0.0
        assert np.all(sol.y[1, :] >= FRICKE_Y)

    @pytest.mark.skipif(not TORCH_AVAILABLE, reason="PyTorch is not installed")
    def test_solver_torch_interoperability(self):
        def torch_rhs(t, y):
            return -0.5 * y

        solver = lf.Solver(method="BDF")
        y0 = torch.tensor([1.0, 2.0], dtype=torch.float64)
        sol = solver.solve(torch_rhs, y0, (0.0, 2.0), t_eval=np.linspace(0.0, 2.0, 20))
        assert sol.success is True
        assert isinstance(sol.tensor_y, torch.Tensor)
        assert sol.torch().shape == (2, 20)


class TestWorkingGroups:
    def test_cymetric_kahler_projector_numpy(self):
        # Indefinite Hermitian matrix (has negative eigenvalue)
        bad_metric = np.array([
            [1.0, 0.5j],
            [-0.5j, -0.2]
        ])
        min_eval_before = np.min(np.linalg.eigvalsh(bad_metric))
        assert min_eval_before < 0

        good_metric = lf.project_kahler_metric(bad_metric, min_eigenval=1e-3)
        evals_after = np.linalg.eigvalsh(good_metric)
        assert np.all(evals_after >= 1e-3 - 1e-7)

        # Monge-Ampere loss
        loss = lf.compute_monge_ampere_loss(good_metric, np.array([1.0]))
        assert loss >= 0.0

    @pytest.mark.skipif(not TORCH_AVAILABLE, reason="PyTorch is not installed")
    def test_cymetric_kahler_projector_torch(self):
        bad_metric = torch.tensor([
            [2.0, 0.0],
            [0.0, -1.0]
        ], dtype=torch.float64)
        projected = lf.project_kahler_metric(bad_metric, min_eigenval=0.01)
        evals = torch.linalg.eigvalsh(projected)
        assert torch.all(evals >= 0.01 - 1e-7)

    def test_swampland_mcmc_filter(self):
        mcmc_filter = lf.SwamplandMCMCFilter(euler_characteristic=24, delta_d_cutoff=6.0)
        batch = [
            {"flux_charge": 0.5, "delta_d": 2.0, "tau_im": 1.0},   # Acceptable
            {"flux_charge": 5.0, "delta_d": 2.0, "tau_im": 1.0},   # Tadpole overflow (> 24/24 = 1.0)
            {"flux_charge": 0.5, "delta_d": 10.0, "tau_im": 1.0},  # Distance breakdown (> 6.0)
            {"flux_charge": 0.5, "delta_d": 2.0, "tau_im": -0.5},  # Metric negative
        ]
        results = mcmc_filter.filter_batch(batch, async_certify=True)
        assert results["total_evaluated"] == 4
        assert results["accepted_count"] == 1
        assert results["rejection_breakdown"]["tadpole"] == 1
        assert results["rejection_breakdown"]["swampland_distance"] == 1
        assert results["rejection_breakdown"]["metric_positivity"] == 1
        assert results["lean_certification"] is not None
        assert results["lean_certification"]["is_satisfied"] is True

    def test_topological_defect_extractor(self):
        extractor = lf.TopologicalDefectExtractor(energy_threshold=0.5)
        # Synthetic 16x16 grid
        grid = np.zeros((16, 16, 2))
        grid[8:, :, 0] = 1.0  # Domain wall
        summary = lf.extract_topological_defects(grid)
        assert summary.num_nodes > 0
        assert summary.num_edges > 0
        assert summary.is_tadpole_neutral is True
        assert summary.lean_certified is True

    def test_lean_verification_client(self):
        client = LeanVerificationClient()
        # Neutral charges
        res_neutral = client.verify_tadpole_cancellation([1, -1, 2, -2])
        assert res_neutral["is_neutral"] is True
        assert res_neutral["lean_verified"] is True
        assert res_neutral["certificate"] is not None

        # Non-neutral charges
        res_charged = client.verify_tadpole_cancellation([1, 2, 3])
        assert res_charged["is_neutral"] is False
        assert res_charged["certificate"] is None

        # Tadpole budget
        res_budget = client.verify_tadpole_budget(total_flux=5.0, euler_char=-200)
        assert res_budget["is_satisfied"] is True
        assert res_budget["max_allowed"] == 200.0 / 24.0

        # K-theory Grothendieck class conservation
        res_k = client.verify_k_theory_conservation(
            brane_e=(2, 1, 0),
            antibrane_f=(1, 1, 0)
        )
        assert res_k["is_conserved"] is True
        assert res_k["soliton_defect_class"] == (1, 0, 0)
