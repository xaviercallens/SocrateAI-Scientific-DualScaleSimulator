#!/bin/bash
DEADLINE=$(( $(date +%s) + 1200 ))
D_TDA=/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/D-tda/exports.json
A_GENUS=/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/A-genus/exports.json
while [ "$(date +%s)" -lt "$DEADLINE" ]; do
  if [ -f "$D_TDA" ]; then echo "FOUND_D_TDA"; exit 0; fi
  if [ -f "$A_GENUS" ]; then echo "FOUND_A_GENUS"; exit 0; fi
  sleep 15
done
echo "TIMEOUT_NEITHER_FOUND"
exit 1
