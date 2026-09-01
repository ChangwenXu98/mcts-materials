#!/bin/bash
#SBATCH --job-name=sh-valid
#SBATCH --partition=venkvis-debug
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=00:30:00
#SBATCH --output=campaign/calib/slurmval-%j.out
#SBATCH --error=campaign/calib/slurmval-%j.err

set -euo pipefail
cd /nfs/turbo/coe-venkvis/changwex/projects/mcts-materials

module load gcc/13.2.0 openmpi/4.1.6
export OMP_NUM_THREADS=1
export PYTHONPATH=/nfs/turbo/coe-venkvis/changwex/projects/mcts-materials/src

PY=/home/changwex/turbo/anaconda3/envs/mace-agent/bin/python
PSEUDO=/nfs/turbo/coe-venkvis/changwex/projects/SuperScientist/opt/pseudo
BIN=/nfs/turbo/coe-venkvis/changwex/qe-build/qe-7.3.1/bin
WORK=/scratch/venkvis_root/venkvis/changwex/superhydride-val

mkdir -p "$WORK"
cp campaign/LaBeH8_root.cif "$WORK/"

$PY campaign/calib/validate.py "$WORK" "$PSEUDO" "$BIN" 24
cp -r "$WORK/calib/results.json" campaign/calib/ 2>/dev/null || true
echo "Done: $(date)"
