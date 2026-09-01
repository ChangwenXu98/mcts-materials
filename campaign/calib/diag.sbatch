#!/bin/bash
#SBATCH --job-name=sh-diag
#SBATCH --partition=venkvis-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=campaign/calib/diag-%j.out
#SBATCH --error=campaign/calib/diag-%j.err

set -euo pipefail
cd /nfs/turbo/coe-venkvis/changwex/projects/mcts-materials
export OMP_NUM_THREADS=16
export PYTHONPATH=/nfs/turbo/coe-venkvis/changwex/projects/mcts-materials/src
PY=/home/changwex/turbo/anaconda3/envs/mace-agent/bin/python
$PY campaign/calib/diag_prerelax.py campaign/calib/diag_result.txt
echo "Done: $(date)"
cat campaign/calib/diag_result.txt
