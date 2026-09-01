#!/bin/bash
#SBATCH --job-name=sh-mcts
#SBATCH --partition=venkvis-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=campaign/logs/slurm-%j.out
#SBATCH --error=campaign/logs/slurm-%j.err

set -euo pipefail
REPO=/nfs/turbo/coe-venkvis/changwex/projects/mcts-materials
cd "$REPO"

module load gcc/13.2.0 openmpi/4.1.6

# QE's own threading fights the MPI decomposition; QERunner pins
# OMP_NUM_THREADS=1 for every pw.x/pp.x/projwfc.x subprocess. The parent
# process keeps its threads for the MACE pre-relaxation.
export OMP_NUM_THREADS=4
export PYTHONPATH="$REPO/src"

# Cluster paths stay in the environment rather than in the shared config.
export QE_BIN_DIR=/nfs/turbo/coe-venkvis/changwex/qe-build/qe-7.3.1/bin
export ESPRESSO_PSEUDO=/nfs/turbo/coe-venkvis/changwex/projects/SuperScientist/opt/pseudo

PY=/home/changwex/turbo/anaconda3/envs/mace-agent/bin/python
OUT="$REPO/campaign/results"

# Wavefunctions and cubes go to scratch, never to turbo.
mkdir -p /scratch/venkvis_root/venkvis/changwex/superhydride-runs "$OUT" campaign/logs

echo "host $(hostname)  start $(date)"
$PY campaign/run_campaign.py campaign/config_campaign.yaml "$OUT" 5
echo "Done: $(date)"
