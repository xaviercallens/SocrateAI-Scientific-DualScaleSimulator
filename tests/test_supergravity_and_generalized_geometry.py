"""
=============================================================================
Unit Test Suite: Supergravity Moduli, Generalized Geometry, and D-Brane Inflow
=============================================================================
Validates:
1. 4D N=1 Supergravity F-term scalar potential, GVW flux superpotentials,
   Kähler metrics, and Weil-Petersson geodesic flows.
2. O(D,D; Z) Generalized Metric, O(D,D) metric condition, Courant brackets,
   and non-geometric flux chains (H -> f -> Q -> R).
3. Microscopic D-brane boundary states, Neumann/Dirichlet reflection tensors,
   and Callan-Harvey anomaly inflow cancellation.
=============================================================================
"""

import math
import pytest
import numpy as np

import leanflow as lf
from leanflow.core.supergravity import (
    SupergravityParameters,
    SupergravityModuli,
    evaluate_sugra_landscape_grid,
)
from leanflow.core.generalized_geometry import (
    GeneralizedMetric,
    FluxChain,
    courant_bracket,
    dorfman_bracket,
    compute_t_fold_monodromy,
)
from leanflow.core.dbrane_inflow import DpBraneBoundaryState, CallanHarveyInflow


class TestSupergravityModuli:
    def setup_method(self):
        self.params = SupergravityParameters(
            a0=2.0,
            b0=1.0,
            a1=-1.0,
            b1=2.0,
            A_np=0.8,
            a_np=0.2 * math.pi,
            c_uplift=0.05,
            m_pl=1.0
        )
        self.sugra = SupergravityModuli(self.params)

    def test_kahler_potential_values(self):
        tau_im = 1.5
        u_im = 1.0
        sigma = 2.0
        k = self.sugra.kahler_potential(tau_im, u_im, sigma)
        assert np.isfinite(k)
        # K = -ln(2*1.5) - 3*ln(2*1.0) - 3*ln(2*2.0)
        expected = -math.log(3.0) - 3.0 * math.log(2.0) - 3.0 * math.log(4.0)
        np.testing.assert_allclose(k, expected, rtol=1e-10)

    def test_kahler_potential_unphysical_boundary(self):
        # Non-positive values return infinity (repulsive barrier)
        assert math.isinf(self.sugra.kahler_potential(-0.5, 1.0, 1.0))
        assert math.isinf(self.sugra.kahler_potential(1.0, -0.2, 1.0))
        assert math.isinf(self.sugra.kahler_potential(1.0, 1.0, 0.0))

    def test_superpotential_values(self):
        tau = complex(0.0, 1.0)
        u = complex(0.0, 1.0)
        t = complex(3.0, 0.0)
        w = self.sugra.superpotential(tau, u, t)
        # W = (a0 - tau*b0) + (a1 - tau*b1)*u + A*exp(-a*t)
        # = (2.0 - 1.0j) + (-1.0 - 2.0j)*1.0j + 0.8*exp(-0.2*pi*3)
        # (-1 - 2j)*j = -j - 2*(-1) = 2 - j
        # (2 - 1j) + (2 - 1j) = 4 - 2j
        expected = complex(4.0, -2.0) + 0.8 * math.exp(-0.2 * math.pi * 3.0)
        np.testing.assert_allclose(w.real, expected.real, rtol=1e-10)
        np.testing.assert_allclose(w.imag, expected.imag, rtol=1e-10)

    def test_covariant_derivatives_structure(self):
        tau = complex(0.5, 1.2)
        u = complex(0.0, 1.0)
        t = complex(4.0, 0.0)
        w, d_tau, d_u, d_t = self.sugra.covariant_derivatives(tau, u, t)
        assert np.isfinite(abs(w))
        assert np.isfinite(abs(d_tau))
        assert np.isfinite(abs(d_u))
        assert np.isfinite(abs(d_t))

    def test_f_term_potential_positivity_and_uplift(self):
        tau = complex(0.1, 1.2)
        u = complex(0.0, 1.0)
        t = complex(4.0, 0.0)
        v_with_up = self.sugra.f_term_potential(tau, u, t, include_uplift=True)
        v_no_up = self.sugra.f_term_potential(tau, u, t, include_uplift=False)

        assert np.isfinite(v_with_up)
        assert np.isfinite(v_no_up)
        assert v_with_up > v_no_up
        # Uplift difference should be c_uplift / (2*sigma)^2
        expected_diff = self.params.c_uplift / ((2.0 * t.real) ** 2)
        np.testing.assert_allclose(v_with_up - v_no_up, expected_diff, rtol=1e-10)

    def test_potential_derivatives_tau(self):
        tau = complex(0.2, 1.5)
        u = complex(0.0, 1.0)
        t = complex(3.5, 0.0)
        v0, dv1, dv2 = self.sugra.potential_derivatives_tau(tau, u, t, eps=1e-5)

        assert np.isfinite(v0)
        assert np.isfinite(dv1)
        assert np.isfinite(dv2)

    def test_weil_petersson_christoffel_tau(self):
        tau_im = 2.0
        gammas = self.sugra.weil_petersson_christoffel_tau(tau_im)
        assert gammas[(0, 0, 1)] == -0.5
        assert gammas[(0, 1, 0)] == -0.5
        assert gammas[(1, 0, 0)] == 0.5
        assert gammas[(1, 1, 1)] == -0.5

    def test_geodesic_step_tau_trajectory(self):
        tau0 = complex(0.0, 2.0)
        v_tau0 = complex(0.05, -0.02)
        u = complex(0.0, 1.0)
        t = complex(3.0, 0.0)

        tau_next, v_next = self.sugra.geodesic_step_tau(tau0, v_tau0, u, t, dt=0.01)
        assert tau_next.imag > 0.0
        assert np.isfinite(tau_next.real)
        assert np.isfinite(tau_next.imag)

    def test_landscape_grid_evaluation(self):
        grid_data = evaluate_sugra_landscape_grid(n_points=10)
        assert grid_data["potential_grid"].shape == (10, 10)
        assert np.all(np.isfinite(grid_data["potential_grid"]))
        assert grid_data["min_potential"] <= grid_data["max_potential"]


class TestGeneralizedGeometry:
    def test_generalized_metric_block_structure_and_odd_condition(self):
        d = 3
        rng = np.random.default_rng(42)
        a = rng.standard_normal((d, d))
        g = a.T @ a + np.eye(d)
        b = rng.standard_normal((d, d))
        b = 0.5 * (b - b.T)

        gen_metric = GeneralizedMetric(g, b)
        assert gen_metric.d == d
        assert gen_metric.h.shape == (2 * d, 2 * d)

        # Check O(D,D) metric condition: H^T eta H = eta
        assert gen_metric.verify_odd_invariance(tol=1e-8) is True

        # Check inverse property: H^{-1} = eta H eta
        assert gen_metric.verify_inverse_relation(tol=1e-8) is True

    def test_generalized_metric_t_duality_transformation(self):
        d = 2
        g = np.diag([2.0, 1.5])
        b = np.array([[0.0, 0.4], [-0.4, 0.0]])
        gen_metric = GeneralizedMetric(g, b)

        gen_dual = gen_metric.t_duality_transform(direction=0)
        assert gen_dual.verify_odd_invariance(tol=1e-8) is True
        assert gen_dual.verify_inverse_relation(tol=1e-8) is True

        # Test T-duality radius inversion when B=0: R -> 1/R
        gen_metric_nob = GeneralizedMetric(np.diag([2.0, 1.0]), np.zeros((2, 2)))
        gen_dual_nob = gen_metric_nob.t_duality_transform(direction=0)
        np.testing.assert_allclose(gen_dual_nob.g[0, 0], 0.5, rtol=1e-10)

    def test_courant_and_dorfman_brackets(self):
        d = 3
        x = np.array([1.0, 0.0, 0.0])
        xi = np.array([0.0, 1.0, 0.0])
        u = np.concatenate([x, xi])

        y = np.array([0.0, 1.0, 0.0])
        eta = np.array([0.0, 0.0, 1.0])
        v = np.concatenate([y, eta])

        # Antisymmetry of Courant bracket: [u, u]_C = 0
        c_uu = courant_bracket(u, u)
        np.testing.assert_allclose(c_uu, np.zeros(2 * d), atol=1e-10)

        # Skew-symmetry: [u, v]_C = -[v, u]_C
        c_uv = courant_bracket(u, v)
        c_vu = courant_bracket(v, u)
        np.testing.assert_allclose(c_uv, -c_vu, atol=1e-10)

    def test_flux_chain_progression(self):
        chain = FluxChain.from_standard_three_torus_h_flux(h_val=1.0)
        assert chain.active_flux_type() == "Standard Geometric (H-flux)"

        c1 = chain.apply_t_duality(cycle=0)
        assert c1.active_flux_type() == "Geometric Twisted Torus (f-flux)"
        assert abs(c1.f[0, 1, 2]) > 0.5

        c2 = c1.apply_t_duality(cycle=1)
        assert c2.active_flux_type() == "Non-Geometric Locally-Geometric (Q-flux T-Fold)"
        assert abs(c2.q[0, 1, 2]) > 0.5

        c3 = c2.apply_t_duality(cycle=2)
        assert c3.active_flux_type() == "Non-Geometric Non-Local (R-flux)"
        assert abs(c3.r[0, 1, 2]) > 0.5

    def test_t_fold_monodromy_odd_invariance(self):
        q_val = 1.0
        m = compute_t_fold_monodromy(q_val)
        assert m.shape == (4, 4)
        
        # O(2, 2) metric eta = [[0, I], [I, 0]]
        i2 = np.eye(2)
        z2 = np.zeros((2, 2))
        eta = np.block([[z2, i2], [i2, z2]])

        # Verify M^T eta M = eta
        m_eta_m = m.T @ eta @ m
        np.testing.assert_allclose(m_eta_m, eta, atol=1e-10)


class TestDbraneBoundaryStatesAndInflow:
    def test_dp_brane_boundary_state_properties(self):
        for p in [1, 2, 3, 5, 7, 9]:
            state = DpBraneBoundaryState(p=p)
            assert state.neumann_count == p + 1
            assert state.dirichlet_count == 9 - p

            s_mat = state.s_tensor
            # S^2 = I
            np.testing.assert_allclose(s_mat @ s_mat, np.eye(10), atol=1e-10)

            # Trace should be (p+1) - (9-p) = 2p - 8
            expected_trace = 2 * p - 8
            np.testing.assert_allclose(np.trace(s_mat), expected_trace, atol=1e-10)

            # Test oscillator annihilation
            alpha_tilde = np.ones(10)
            # alpha = -S * alpha_tilde
            alpha = -s_mat @ alpha_tilde
            assert state.verify_oscillator_annihilation(alpha, alpha_tilde) is True

    def test_callan_harvey_inflow_cancellation(self):
        inflow = CallanHarveyInflow(p1_tangent=-48.0, f_flux=1.0)
        res = inflow.compute_anomaly_descent(gauge_variation=1.0)

        assert res["is_anomaly_cancelled"] is True
        assert abs(res["delta_s_total"]) < 1e-12
        assert res["delta_s_defect"] == -res["delta_s_bulk"]
        assert "callan_harvey_exact_anomaly_cancellation" in res["lean4_certificate"]

    def test_leanflow_module_exports(self):
        assert hasattr(lf, "SupergravityModuli")
        assert hasattr(lf, "SupergravityParameters")
        assert hasattr(lf, "GeneralizedMetric")
        assert hasattr(lf, "FluxChain")
        assert hasattr(lf, "DpBraneBoundaryState")
        assert hasattr(lf, "CallanHarveyInflow")
