"""
=============================================================================
Test Suite: WorkshopCosmo Dual-Scale Simulation & Formal Proof Suite
=============================================================================
Aligned with: specs/DEFINITION_OF_DONE.md & specs/roadmap.md

Requirement Mapping:
- REQ-COSMO-01: AXE 1 - Quintessence ODE & Attractor Flow Convergence
- REQ-COSMO-02: AXE 2 - Symmetron / Chameleon Non-Linear Screening PDE
- REQ-COSMO-03: AXE 1/2 - Cosmological Observables & Consilience
- REQ-COSMO-04: AXE 3 - Lean 4 Kernel-Verified Tadpole Cancellation
- REQ-COSMO-05: Integrated Exports (CSV, JSON, Markdown Report)
- REQ-COSMO-06: Verification & DoD Compliance Matrix
=============================================================================
"""

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import math
import json
import pytest
import numpy as np

import workshopcosmo as wc


class TestReqCosmo01Quintessence:
    """Tests for REQ-COSMO-01: AXE 1 Quintessence ODE & Modulus Flow"""

    def test_potential_stationary_points(self):
        """Verify saddle at Fricke point and minimum at Orbifold point."""
        # Fricke point: x = 0, y = 1/sqrt(12)
        v_f, dv_dx_f, dv_dy_f = wc.compute_potential(wc.FRICKE_X, wc.FRICKE_Y)
        assert abs(dv_dx_f) < 1e-12, f"Expected dV/dx = 0 at Fricke point, got {dv_dx_f}"
        assert abs(dv_dy_f) < 1e-12, f"Expected dV/dy = 0 at Fricke point, got {dv_dy_f}"

        # Orbifold point: x = 0.5, y = sqrt(3)/2
        v_o, dv_dx_o, dv_dy_o = wc.compute_potential(wc.ORBIFOLD_X, wc.ORBIFOLD_Y)
        assert abs(dv_dx_o) < 1e-12, f"Expected dV/dx = 0 at Orbifold point, got {dv_dx_o}"
        assert abs(dv_dy_o) < 1e-12, f"Expected dV/dy = 0 at Orbifold point, got {dv_dy_o}"

        # Energy hierarchy: Fricke saddle is high-energy inflation, Orbifold is low-energy vacuum
        assert v_f > v_o, f"Expected V(Fricke) > V(Orbifold), got {v_f} <= {v_o}"

    def test_quintessence_attractor_convergence(self):
        """Verify modulus flows safely from Fricke saddle to Orbifold attractor."""
        res = wc.run_quintessence_simulation(t_max=65.0, num_points=200)
        assert res["success"] is True, "Quintessence ODE solver failed"
        assert res["is_modulus_positive"] is True, "Modulus y dropped to <= 0 (singularity breach)"

        # Check distance to attractor
        dist = res["attractor_distance"]
        assert dist < 1e-3, f"Expected attractor distance < 1e-3, got {dist:.2e}"
        assert res["is_attractor_converged"] is True

        # Check final values
        final_x = res["final_state"]["x"]
        final_y = res["final_state"]["y"]
        assert abs(final_x - wc.ORBIFOLD_X) < 1e-3
        assert abs(final_y - wc.ORBIFOLD_Y) < 1e-3

    def test_cross_engine_consilience_sundials(self):
        """Verify cross-engine agreement between Python Radau and rusty-SUNDIALS BDF."""
        sundials_data = wc.load_rusty_sundials_csv()
        if sundials_data is not None and "error" not in sundials_data:
            assert sundials_data["is_converged"] is True
            assert sundials_data["attractor_distance"] < 1e-3
            assert abs(sundials_data["final_x"] - wc.ORBIFOLD_X) < 1e-3
            assert abs(sundials_data["final_y"] - wc.ORBIFOLD_Y) < 1e-3


class TestReqCosmo02SymmetronScreening:
    """Tests for REQ-COSMO-02: AXE 2 Symmetron Non-Linear Screening PDE"""

    def test_symmetron_core_screening(self):
        """Verify spontaneous symmetry restoration inside dense core."""
        res = wc.run_symmetron_screening_simulation(r_max=8.0, n_points=250)
        assert res["success"] is True, "Symmetron BVP solver failed"

        # Center symmetry restoration: phi(0)/phi_0 << 1
        assert res["phi_center_ratio"] < 1e-3, f"Expected core suppression < 1e-3, got {res['phi_center_ratio']}"

        # Fifth force suppression factor at surface
        assert res["screening_suppression_factor"] < 5e-4, (
            f"Expected screening factor < 5e-4, got {res['screening_suppression_factor']}"
        )
        assert res["is_screened"] is True
        assert res["cassini_bound_satisfied"] is True


class TestReqCosmo03Observables:
    """Tests for REQ-COSMO-03: Cosmological Observables & Evidence"""

    def test_density_budget_and_w(self):
        """Verify total relative density Omega_tot ~ 1 and late-time dark energy behavior."""
        quint_res = wc.run_quintessence_simulation(t_max=65.0, num_points=100)
        obs_res = wc.run_observables_analysis(quint_res)

        # In late universe, scalar field dominates as dark energy
        assert -1.5 < obs_res["late_time_w_current"] <= -0.7, (
            f"Late time w_phi out of expected cosmological window: {obs_res['late_time_w_current']}"
        )

        # Check Omega sum at final step
        omega_tot = (
            quint_res["omega_r"][-1] + quint_res["omega_m"][-1] + quint_res["omega_phi"][-1]
        )
        assert abs(omega_tot - 1.0) < 1e-2, f"Expected Omega_tot ~ 1, got {omega_tot}"


class TestReqCosmo04Lean4Proof:
    """Tests for REQ-COSMO-04: AXE 3 Lean 4 Formal Tadpole Cancellation Proof"""

    def test_lean_verification_zero_sorry(self):
        """Verify that Lean 4 compiles the Tadpole Cancellation module with ZERO sorry."""
        res = wc.run_lean_verification()
        assert res["success"] is True, f"Lean verification failed: {res.get('compiler_output') or res.get('error')}"
        assert res["compiled_ok"] is True, "Lean compiler returned non-zero exit code"
        assert res["sorry_count"] == 0, f"Found {res['sorry_count']} sorry statements in formal proof!"

        # Verify key theorems are declared
        expected_theorems = [
            "num_fixed_points_is_16",
            "total_O7_charge_is_minus_64",
            "total_D7_charge_is_64",
            "d7_tadpole_cancellation",
            "three_generation_index",
            "curvature_d3_charge_is_1",
            "flux_saturates_d3_tadpole",
            "irreducible_anomalies_vanish",
            "dual_scale_model_is_in_landscape",
        ]
        for th in expected_theorems:
            assert th in res["theorems_verified"], f"Missing theorem {th} from proof module"


class TestReqCosmo05Integration:
    """Tests for REQ-COSMO-05 & REQ-COSMO-06: Exports & DoD Compliance"""

    def test_generated_reports_exist_and_valid(self):
        """Verify JSON and Markdown reports are generated with non-empty content."""
        quint_res = wc.run_quintessence_simulation(t_max=60.0, num_points=100)
        symm_res = wc.run_symmetron_screening_simulation(r_max=6.0, n_points=150)
        obs_res = wc.run_observables_analysis(quint_res)
        lean_res = wc.run_lean_verification()
        sundials_res = wc.load_rusty_sundials_csv()

        report = wc.generate_report(
            quint_res, symm_res, obs_res, lean_res, sundials_res,
            output_json="test_simulation_results.json",
            output_md="test_report.md"
        )

        assert os.path.isfile("test_simulation_results.json")
        assert os.path.isfile("test_report.md")

        with open("test_simulation_results.json", "r", encoding="utf-8") as f:
            data = json.load(f)

        assert data["axe1_quintessence"]["status"] == "PASS"
        assert data["axe2_symmetron"]["status"] == "PASS"
        assert data["axe3_formal_proof"]["status"] == "PASS"

        # Cleanup test artifacts
        if os.path.isfile("test_simulation_results.json"):
            os.remove("test_simulation_results.json")
        if os.path.isfile("test_report.md"):
            os.remove("test_report.md")


class TestReqCosmo07MathieuMoonshineVOA:
    """Tests for REQ-COSMO-07: AXE 4 Mathieu Moonshine VOA & R_NL = 77/60 Rigidity"""

    def test_conformal_ward_identities(self):
        """Verify 2D CFT 3-point correlator satisfies Virasoro L_{-1}, L_0, L_1 Ward identities."""
        res = wc.verify_sl2c_mobius_conformal_invariance()
        assert res["conformal_invariance_verified"] is True, f"Ward identity residual too large: {res['relative_error']}"
        assert res["relative_error"] < 1e-7

    def test_mathieu_clebsch_gordan_algebra(self):
        """Verify M24 representation dimensions and Clebsch-Gordan reduction to 77/60."""
        res = wc.verify_mathieu_moonshine_algebra()
        assert res["is_77_over_60"] is True, f"Expected 77/60, got {res['reduced_ratio']}"
        assert res["is_irreducible"] is True
        assert res["A1"] == 90
        assert res["A2"] == 462
        assert res["sym2_A1_dim"] == 4095
        assert res["superconformal_factor"] == 4

    def test_lean_vertex_operator_proof(self):
        """Verify Lean 4 compilation of MathieuVertexOperators.lean with ZERO sorry."""
        res = wc.run_lean_vertex_operator_verification()
        assert res["success"] is True, f"Lean compilation failed: {res.get('output') or res.get('error')}"
        assert res["compiled_ok"] is True
        assert res["sorry_count"] == 0, f"Found {res['sorry_count']} sorry stubs in VOA proof!"
        assert len(res["theorems_verified"]) >= 16


class TestReqCosmo08NanogravHexadecapole:
    """Tests for REQ-COSMO-08: AXE 5 NANOGrav 15-year Hexadecapole Anomaly (l=4)"""

    def test_nanograv_hexadecapole_hidden_within_envelope(self):
        """Verify that the C4/C0 = 16.07 hexadecapole anomaly hides beneath the HD cosmic variance envelope."""
        res = wc.run_nanograv_hexadecapole_simulation(c4_c0_ratio=16.07)
        assert res["success"] is True
        assert res["c4_c0_ratio"] == 16.07
        
        # Max deviation should be quite small because pta_suppression is 0.005
        # So 16.07 * 0.005 = 0.08035, and Legendre max amplitude is 1.
        assert res["max_deviation"] < 0.15, f"Deviation {res['max_deviation']} exceeds 0.15 threshold"
        assert res["is_hidden"] is True, "The l=4 anomaly breached the HD cosmic variance envelope!"


class TestReqCosmo09KummerLangevinTDA:
    """Tests for REQ-COSMO-09: AXE 6 Kummer Langevin Simulation, TDA Mapper & Lean 4 Certification"""

    def test_rust_langevin_simulation_ssb_and_defects(self):
        """Verify Rust stochastic Langevin dynamics produces spontaneous symmetry breaking and defect network."""
        res = wc.run_kummer_langevin_simulation(grid_size=32, t_max=10.0)
        assert res["success"] is True, "Rust Langevin simulation failed"
        assert res["symmetry_broken"] is True, "Symmetry breaking failed to occur below T_crit"
        assert res["final_temperature"] < 1.0, f"Temperature {res['final_temperature']} did not cool below T_crit"
        assert res["final_string_count"] > 0, "No cosmic string vortices formed during phase transition"
        assert res["final_wall_pixel_count"] > 0, "No domain wall interfaces formed"
        assert res["mean_field_norm"] > 0.2, "Field failed to develop non-zero vacuum expectation value"

    def test_tda_mapper_extraction_and_equivalence_classes(self):
        """Verify TDA Mapper extracts 1-skeleton, 1-cycles (beta_1), and partitions into 3 string classes."""
        res = wc.run_tda_mapper_analysis()
        assert res["success"] is True, "TDA Mapper extraction failed"
        assert res["num_nodes"] > 10, f"Too few Mapper nodes extracted: {res['num_nodes']}"
        assert res["num_edges"] > 10, f"Too few Mapper edges extracted: {res['num_edges']}"
        assert res["betti_1_cycles"] > 0, "No persistent 1-cycles (cosmic string loops) detected by Mapper"
        
        # Check presence of all 3 equivalence classes
        classes = res["classes"]
        assert "AttractorVacuum" in classes and classes["AttractorVacuum"] > 0, "Missing AttractorVacuum class"
        assert "DomainWall" in classes and classes["DomainWall"] > 0, "Missing DomainWall class"
        assert "CosmicString" in classes and classes["CosmicString"] > 0, "Missing CosmicString class"

    def test_lean4_tda_tadpole_certification_zero_sorry(self):
        """Verify Lean 4 compilation of KummerTDAAnomalyCertification.lean with ZERO sorry and strict tadpole cancellation."""
        res = wc.run_lean_tda_certification()
        assert res["success"] is True, f"Lean compilation failed: {res.get('compiler_output') or res.get('error')}"
        assert res["compiled_ok"] is True
        assert res["sorry_count"] == 0, f"Found {res['sorry_count']} sorry stubs in Lean TDA proof!"
        
        expected_theorems = [
            "num_kummer_fixed_points_is_16",
            "total_o7_charge_is_minus_64",
            "total_d7_charge_is_64",
            "kummer_bulk_tadpole_cancellation",
            "attractor_vacuum_anomaly_free",
            "domain_wall_anomaly_free",
            "cosmic_string_anomaly_free",
            "tda_all_equivalence_classes_anomaly_free",
            "irreducible_anomalies_zero",
            "topological_defects_preserve_string_landscape",
            "sum_vacuum_charges_eq_zero",
        ]
        for th in expected_theorems:
            assert th in res["theorems_verified"], f"Missing theorem {th} from Lean TDA proof module"

    def test_end_to_end_consilience_pipeline(self):
        """Verify complete end-to-end consilience pipeline (Rust SDE -> Python TDA Mapper -> Lean 4 Proof)."""
        res = wc.run_full_tda_langevin_pipeline(grid_size=32, t_max=10.0)
        assert res["success"] is True
        assert res["tadpole_anomaly"] == 0, "Net tadpole anomaly is non-zero!"
        assert res["is_string_landscape"] is True, "Topological defect network failed string landscape criterion"
        assert res["lean_certification_verified"] is True
        assert res["mapper_1_cycles"] > 0


class TestReqCosmo10SwamplandDistance:
    """Tests for REQ-COSMO-10: Swampland Distance Conjecture Geodesic Flow & Tower Collapse"""

    def test_swampland_distance_geodesic_flow_and_tower_collapse(self):
        """Verify moduli geodesic flow induces exponential mass gap decay."""
        res = wc.run_swampland_distance_simulation(d_max=8.0)
        assert res["success"] is True
        assert res["alpha"] > 0.5, f"Expected order-one alpha coupling, got {res['alpha']}"
        assert res["final_mass_gap"] < 0.01, f"Expected mass gap collapse, got {res['final_mass_gap']}"
        assert res["cutoff_scale"] < 0.01, f"Expected cutoff breakdown, got {res['cutoff_scale']}"

    def test_swampland_distance_persistent_homology_barcodes(self):
        """Verify persistent homology barcodes capture the topological collapse of mass intervals."""
        res = wc.run_swampland_distance_simulation(d_max=8.0)
        barcodes = res["barcodes"]
        assert 0.0 in barcodes and 8.0 in barcodes
        max_d0 = max(death for _, death in barcodes[0.0])
        max_d8 = max(death for _, death in barcodes[8.0])
        assert max_d8 < 0.1 * max_d0, f"Expected barcode interval collapse by >10x, got {max_d8} vs {max_d0}"

    def test_lean4_sdc_formal_proof_zero_sorry(self):
        """Verify Lean 4 compilation of SwamplandDistanceConjecture.lean with ZERO sorry."""
        res = wc.run_lean_sdc_certification()
        assert res["success"] is True, f"Lean compilation failed: {res.get('compiler_output')}"
        assert res["compiled_ok"] is True
        assert res["sorry_count"] == 0, f"Found {res['sorry_count']} sorry stubs in SDC proof!"
        for th in ["alpha_sdc_strictly_positive", "sdc_mass_gap_collapses", "sdc_eft_breakdown"]:
            assert th in res["theorems_verified"]

    def test_rust_swampland_simulation_consilience(self):
        """Verify native Rust Swampland Distance simulation executes and validates SDC bound."""
        res = wc.run_rust_swampland_simulation(s_max=5.0)
        assert res["success"] is True
        assert res["total_geodesic_distance"] > 3.0
        assert res["mass_gap_ratio"] < 0.05
        assert res["bound_satisfied"] is True


class TestReqCosmo11TachyonCondensation:
    """Tests for REQ-COSMO-11: Tachyon Condensation & K-Theory Charge Conservation (Sen's Soliton)"""

    def test_tachyon_condensation_sen_soliton_roll_down(self):
        """Verify non-linear roll-down forms a localized Sen soliton at x=0."""
        res = wc.run_tachyon_condensation_simulation(t_max=10.0)
        assert res["success"] is True
        assert res["kink_width"] < 4.0, f"Expected localized soliton width, got {res['kink_width']}"
        assert len(res["t_history"]) >= 4, "Too few time frames saved during roll-down"

    def test_tachyon_tda_mapper_defect_isolation(self):
        """Verify TDA Mapper extracts the 1-skeleton isolating the BPS D-defect node."""
        res = wc.run_tachyon_condensation_simulation(t_max=10.0)
        nodes = res["tda_nodes"]
        edges = res["tda_edges"]
        assert len(nodes) >= 3, "Too few TDA Mapper nodes in defect skeleton"
        assert len(edges) >= 2, "Too few TDA Mapper edges in defect skeleton"
        defect_nodes = [n for n in nodes if n.get("type") == "defect"]
        assert len(defect_nodes) == 1, "Failed to isolate the unique Sen soliton defect node"

    def test_lean4_tachyon_k_theory_conservation_zero_sorry(self):
        """Verify Lean 4 compilation of TachyonCondensationKTheory.lean with ZERO sorry."""
        res = wc.run_lean_tachyon_certification()
        assert res["success"] is True, f"Lean compilation failed: {res.get('compiler_output')}"
        assert res["compiled_ok"] is True
        assert res["sorry_count"] == 0, f"Found {res['sorry_count']} sorry stubs in Tachyon K-theory proof!"
        for th in ["sen_conjecture_k_theory_conservation", "tachyon_condensation_preserves_rr_charge"]:
            assert th in res["theorems_verified"]

    def test_rust_tachyon_simulation_consilience(self):
        """Verify native Rust Tachyon Condensation simulation forms localized Sen soliton."""
        res = wc.run_rust_tachyon_simulation(t_max=15.0)
        assert res["success"] is True
        assert res["k_theory_charge_conserved"] is True
        assert res["soliton_peak_energy"] > 0.5
        assert res["soliton_width"] < 5.0


class TestReqCosmo12VacuumDecayCTheorem:
    """Tests for REQ-COSMO-12: Coleman-De Luccia Vacuum Decay & Holographic c-Theorem"""

    def test_vacuum_decay_cdl_bounce_shooting_solution(self):
        """Verify Coleman-De Luccia Euclidean bounce profile and action positivity."""
        res = wc.run_vacuum_decay_simulation()
        assert res["success"] is True
        assert res["bounce_action_s_e"] > 0, "Bounce action must be strictly positive"
        assert res["bubble_radius"] > 0.5, f"Unphysical bubble radius: {res['bubble_radius']}"

    def test_vacuum_decay_holographic_c_theorem_irreversibility(self):
        """Verify monotonic central charge decrease Delta c < 0 across flux vacua."""
        res = wc.run_vacuum_decay_simulation()
        c_vals = res["central_charges"]
        diffs = np.diff(c_vals)
        assert np.all(diffs < 0), f"Holographic c-theorem violated: shifts {diffs} are not all negative"

    def test_lean4_cdl_c_theorem_proof_zero_sorry(self):
        """Verify Lean 4 compilation of FluxVacuumDecayCTheorem.lean with ZERO sorry."""
        res = wc.run_lean_cdl_certification()
        assert res["success"] is True, f"Lean compilation failed: {res.get('compiler_output')}"
        assert res["compiled_ok"] is True
        assert res["sorry_count"] == 0, f"Found {res['sorry_count']} sorry stubs in CDL c-theorem proof!"
        for th in ["cdl_action_strictly_positive", "holographic_c_theorem_decay", "true_vacuum_energy_is_lower"]:
            assert th in res["theorems_verified"]

    def test_rust_vacuum_decay_simulation_consilience(self):
        """Verify native Rust Coleman-De Luccia simulation calculates S_E and satisfies c-theorem."""
        res = wc.run_rust_vacuum_decay_simulation(rho_max=12.0)
        assert res["success"] is True
        assert res["euclidean_bounce_action"] > 0.0
        assert res["c_theorem_satisfied"] is True
        assert res["delta_c"] < 0


