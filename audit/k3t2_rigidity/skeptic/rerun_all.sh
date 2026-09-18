#!/bin/bash
# Skeptic rerun of every forward-track and reverse script, in place.
# Byte-comparison is done afterwards against the committed JSON with git diff.
PY=/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python
R=/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity
L=$R/skeptic/logs
run() { d=$1; s=$2; ( cd $R/$d && /usr/bin/time -f "%e s" $PY $s > $L/${d}__${s}.log 2>&1; echo "exit=$? $d/$s" >> $L/${d}__${s}.log ); }
( run genus-moonshine run_track_a.py ) &
( run dyons theta_forms.py; run dyons hurwitz_class_numbers.py; run dyons part1_euler_numbers.py; run dyons part2_psi_m.py; run dyons part3_polar.py; for s in p1_p2_dmz516_m4.py p3_p4_polar_m4_m5.py p5_trace_bound_d7_10.py; do run reverse $s; done ) &
( for s in 01_e8.py 02_k3_mukai_gamma.py 03_oddz_checks.py 04_dual_scale_bound.py 05_tadpole_arithmetic.py 06_rigidity_a_lattice_enum.py 07_rigidity_b_tduality_radius.py; do run lattices-duality $s; done ) &
( for s in 01_controls.py 02_invariance_and_quotient.py 03_resolution_hybrid.py 04_signature.py 05_kunneth_k3xt2.py 07_rigidity_scan.py 06_pointcloud_clifford_torus.py; do run tda-gudhi $s; done ) &
wait
echo ALLDONE > $L/ALLDONE
