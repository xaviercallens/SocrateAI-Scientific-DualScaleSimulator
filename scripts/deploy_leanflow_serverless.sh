#!/bin/bash
set -e

# Prompt for GCP Project ID if not provided
if [ -z "$1" ]; then
  echo "Usage: $0 <GCP_PROJECT_ID>"
  exit 1
fi

PROJECT_ID=$1
IMAGE_NAME="gcr.io/${PROJECT_ID}/leanflow-sim:latest"
JOB_NAME="leanflow-weatherbench-sim"
REGION="us-central1"

echo "============================================================"
echo "🚀 Deploying LeanFlow Serverless Simulation Pipeline"
echo "============================================================"

# Move to the parent directory to include both repositories in the context
cd /home/xavkal/xdev

echo "[1/3] Building the unified Docker container..."
docker build -t ${IMAGE_NAME} -f Dockerfile.leanflow .

echo "[2/3] Pushing to Google Container Registry..."
# Requires `gcloud auth configure-docker`
docker push ${IMAGE_NAME}

echo "[3/3] Executing Cloud Run Job on Preemptible Hardware (Min=0)..."
gcloud run jobs create ${JOB_NAME} \
  --image ${IMAGE_NAME} \
  --tasks 1 \
  --max-retries 3 \
  --cpu 8 \
  --memory 32Gi \
  --set-env-vars="LEAN_PATH=/opt/lean/lib,SUNDIALS_TOLERANCE=1e-6" \
  --region ${REGION} \
  --execute-now

echo "============================================================"
echo "✅ Deployment successful. Monitor logs in the GCP Console."
echo "============================================================"
