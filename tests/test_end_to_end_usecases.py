"""
=============================================================================
End-to-End Use Cases & Hugging Face Deployment Test Suite
=============================================================================
Tests:
1. Pretrained PyTorch QuinticMetricPINN model loading, forward pass & autograd.
2. Serialized datasets integrity (CICY, Kreuzer-Skarke, 1,000 MCMC vacua).
3. Interactive Gradio app construction, callbacks, and plotting pipelines.
4. Publication articles (PDFs) and Google Colab notebook schema validity.
=============================================================================
"""

import os
import json
import pytest
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import torch
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False

try:
    import gradio as gr
    GRADIO_AVAILABLE = True
except ImportError:
    GRADIO_AVAILABLE = False

import leanflow as lf


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
HF_ROOT = os.path.join(REPO_ROOT, "deployment", "huggingface")


@pytest.mark.skipif(not TORCH_AVAILABLE, reason="PyTorch required for PINN tests")
class TestPretrainedPINNModel:
    def test_checkpoint_exists_and_loads(self):
        ckpt_path = os.path.join(HF_ROOT, "models", "quintic_metric_pinn.pt")
        assert os.path.exists(ckpt_path), f"Checkpoint missing at {ckpt_path}"
        assert os.path.getsize(ckpt_path) > 10000, "Checkpoint file abnormally small"
        
        # Import model class
        import sys
        sys.path.insert(0, os.path.join(HF_ROOT, "models"))
        from quintic_metric_pinn import QuinticMetricPINN
        
        model = QuinticMetricPINN(input_dim=6, hidden_dim=64)
        ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=True)
        state_dict = ckpt["model_state_dict"] if "model_state_dict" in ckpt else ckpt
        model.load_state_dict(state_dict)
        model.eval()
        
        # Test forward pass with real patch coordinates (dim=6 for 3 complex coords)
        x = torch.randn(4, 6, dtype=torch.float64)
        raw_metric = model(x)
        assert raw_metric.shape == (4, 3, 3)
        assert raw_metric.dtype == torch.complex128
        
        # Verify Hermiticity
        diff = raw_metric - raw_metric.mH
        assert torch.max(torch.abs(diff)).item() < 1e-6
        
    def test_forward_guarded_and_autograd(self):
        import sys
        sys.path.insert(0, os.path.join(HF_ROOT, "models"))
        from quintic_metric_pinn import QuinticMetricPINN
        
        model = QuinticMetricPINN(input_dim=6, hidden_dim=64)
        x = torch.randn(2, 6, dtype=torch.float64)
        
        # Guarded pass
        guarded_metric, min_evals = model.forward_guarded(x, min_eigenval=0.01)
        assert guarded_metric.shape == (2, 3, 3)
        assert torch.min(min_evals).item() >= 0.0099, f"Kähler cone violated: min_eval = {torch.min(min_evals).item()}"
        
        # Test loss and backward pass
        omega_sq = torch.ones(2, dtype=torch.float64)
        loss, _ = model.compute_monge_ampere_loss(x, omega_sq)
        assert loss.item() >= 0.0
        loss.backward()
        
        # Verify gradients exist in model parameters
        has_grads = any(p.grad is not None and torch.norm(p.grad) > 0 for p in model.parameters())
        assert has_grads, "Autograd backward pass failed to compute gradients"


class TestDeploymentDatasets:
    def test_cicy_canonical_dataset(self):
        path = os.path.join(HF_ROOT, "datasets", "cicy_canonical.json")
        assert os.path.exists(path), f"CICY dataset missing at {path}"
        
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
            
        assert "quintic" in data
        assert "k3_x_t2" in data
        quintic = data["quintic"]
        assert quintic["chi"] == -200
        assert quintic["tadpole_bound"] == pytest.approx(200.0 / 24.0, abs=1e-3)
        assert quintic["h11"] == 1
        assert quintic["h21"] == 101

    def test_kreuzer_skarke_dataset(self):
        path = os.path.join(HF_ROOT, "datasets", "kreuzer_skarke_canonical.json")
        assert os.path.exists(path), f"KS dataset missing at {path}"
        
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
            
        assert "ks_quintic" in data
        assert "ks_k3_surface" in data
        poly_q = data["ks_quintic"]
        assert poly_q["num_vertices"] == 5
        assert poly_q["chi"] == -200

    def test_mcmc_vacuum_samples(self):
        path = os.path.join(HF_ROOT, "datasets", "mcmc_vacuum_samples.json")
        assert os.path.exists(path), f"MCMC samples missing at {path}"
        
        with open(path, "r", encoding="utf-8") as f:
            payload = json.load(f)
            
        assert payload["metadata"]["total_samples"] == 1000
        assert len(payload["samples"]) > 0
        assert payload["metadata"]["acceptance_rate"] > 0.0
        assert payload["metadata"]["lean4_certified"] is True
        
        # Verify schema of individual samples
        sample0 = payload["samples"][0]
        assert "flux_charge" in sample0
        assert "delta_d" in sample0
        assert "tau_im" in sample0


@pytest.mark.skipif(not GRADIO_AVAILABLE, reason="Gradio required for Web App tests")
class TestHuggingFaceGradioApp:
    def test_callbacks_execution(self):
        import sys
        sys.path.insert(0, os.path.join(HF_ROOT, "space"))
        from app import evaluate_cymetric_pinn, run_swampland_explorer, run_defect_extractor
        
        # 1. Calabi-Yau metric callback
        status, guaranteed, fig = evaluate_cymetric_pinn(-0.8, 0.001)
        assert "UNPHYSICAL DRIFT" in status
        assert "LeanFlow Certified" in guaranteed
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
        
        # 2. Swampland MCMC explorer callback
        metrics, cert, fig_mcmc = run_swampland_explorer("quintic", 150, 7.0)
        assert "Evaluated Candidates" in metrics
        assert "150" in metrics
        assert "SocrateAI.StringTheory.AtiyahSingerK3" in cert
        assert isinstance(fig_mcmc, plt.Figure)
        plt.close(fig_mcmc)
        
        # 3. Cosmic defect extractor callback
        summary, fig_defect = run_defect_extractor(16, 0.35)
        assert "Simplicial Nodes" in summary
        assert "Tadpole Neutral" in summary
        assert isinstance(fig_defect, plt.Figure)
        plt.close(fig_defect)

    def test_build_app_structure(self):
        import sys
        sys.path.insert(0, os.path.join(HF_ROOT, "space"))
        from app import build_app
        
        demo = build_app()
        assert isinstance(demo, gr.Blocks)
        assert demo.title == "LeanFlow: Neuro-Symbolic String Theory & SciML Engine"


class TestPublicationArticlesAndColabNotebook:
    def test_preprints_exist(self):
        t_duality_pdf = os.path.join(HF_ROOT, "articles", "T_duality_Alone.pdf")
        engine_pdf = os.path.join(HF_ROOT, "articles", "leanflow_engine.pdf")
        
        assert os.path.exists(t_duality_pdf), f"Missing {t_duality_pdf}"
        assert os.path.getsize(t_duality_pdf) > 1000000, "T_duality_Alone.pdf too small (< 1MB)"
        
        assert os.path.exists(engine_pdf), f"Missing {engine_pdf}"
        assert os.path.getsize(engine_pdf) > 500000, "leanflow_engine.pdf too small (< 500KB)"

    def test_colab_notebook_schema(self):
        nb_path = os.path.join(REPO_ROOT, "notebooks", "leanflow_community_showcase.ipynb")
        assert os.path.exists(nb_path), f"Missing Colab notebook at {nb_path}"
        
        with open(nb_path, "r", encoding="utf-8") as f:
            nb = json.load(f)
            
        assert nb["nbformat"] == 4
        assert len(nb["cells"]) >= 10
        
        # Check that key topics are present in markdown / code
        all_text = " ".join("".join(cell.get("source", [])) for cell in nb["cells"])
        assert "Kähler Cone" in all_text or "Kahler Cone" in all_text
        assert "Swampland" in all_text
        assert "Atiyah-Singer" in all_text or "AtiyahSinger" in all_text
        assert "TopologicalDefectExtractor" in all_text
        assert "Gradio" in all_text
