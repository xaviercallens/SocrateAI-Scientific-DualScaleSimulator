# Hugging Face Deployment Bundle for LeanFlow 2.0

This directory contains the release artifacts and deployment packages for the **LeanFlow Neuro-Symbolic String Theory & SciML Suite** ready for publication on Hugging Face:

```
deployment/huggingface/
├── space/               # Interactive Gradio Web App for Hugging Face Spaces
│   ├── app.py           # Multi-tab interactive application
│   ├── requirements.txt # Python dependencies
│   ├── README.md        # HF Space metadata & documentation
│   └── leanflow/        # Bundled core neuro-symbolic engine
├── models/              # Pretrained Neural Network Surrogates
│   ├── quintic_metric_pinn.py # PyTorch complex Hermitian metric PINN
│   └── quintic_metric_pinn.pt # Pretrained weights for Calabi-Yau metric learning
├── datasets/            # Topological & Geometric String Theory Datasets
│   ├── cicy_canonical.json           # Canonical CICY Calabi-Yau threefold database
│   ├── kreuzer_skarke_canonical.json # Kreuzer-Skarke 4D reflexive polyhedra benchmarks
│   └── mcmc_vacuum_samples.json      # 1,000 screened MCMC string vacua with certificates
└── articles/            # Scientific Preprints (PDFs)
    ├── T_duality_Alone.pdf  # 31-page flagship physics & cosmology paper
    └── leanflow_engine.pdf  # 12-page computational architecture paper
```

---

## 🚀 Deployment Instructions

### 1. Deploy Space to Hugging Face Spaces
```bash
# Create Space repo on Hugging Face (SDK: Gradio)
git clone https://huggingface.co/spaces/<YOUR_USERNAME>/leanflow-engine
cp -r deployment/huggingface/space/* <YOUR_USERNAME>/leanflow-engine/
cd <YOUR_USERNAME>/leanflow-engine
git add .
git commit -m "feat: deploy LeanFlow 2.0 interactive space"
git push
```

### 2. Upload Model to Hugging Face Hub
```python
from huggingface_hub import HfApi
api = HfApi()
api.upload_file(
    path_or_fileobj="deployment/huggingface/models/quintic_metric_pinn.pt",
    path_in_repo="quintic_metric_pinn.pt",
    repo_id="<YOUR_USERNAME>/leanflow-quintic-pinn",
    repo_type="model",
)
```

### 3. Upload Datasets to Hugging Face Hub
```python
api.upload_folder(
    folder_path="deployment/huggingface/datasets",
    repo_id="<YOUR_USERNAME>/leanflow-calabi-yau-vacua",
    repo_type="dataset",
)
```
