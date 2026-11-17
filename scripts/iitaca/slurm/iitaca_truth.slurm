#!/bin/bash
# Monte Carlo truth of the scored events of iItaca runs (prepare_iitaca_truth.jl). Hyperion cluster.
#   sbatch --output=LOGDIR/%x.%j.out --error=LOGDIR/%x.%j.err iitaca_truth.slurm OUTDIR
# Reads the NEXUS files under NAUSICAA_DATA (read-only) and writes
# OUTDIR/truth_bi214.h5 and OUTDIR/truth_single_electrons.h5.
#SBATCH --job-name=iitaca_truth
#SBATCH --partition=general
#SBATCH --qos=serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --time=04:00:00
set +e
type module >/dev/null 2>&1 || source /etc/profile
set -eo pipefail
module load Julia/1.11.6-linux-x86_64
REPO=${NAUSICAA_REPO:-/scratch/jjgomez/NAUSICAA}
export JULIA_DEPOT_PATH=${NAUSICAA_JULIA_DEPOT:-/scratch/jjgomez/julia-depot}
export NAUSICAA_ROOT=$REPO
: "${NAUSICAA_DATA:?Set NAUSICAA_DATA}"
OUT=$1; SCRIPT=${NAUSICAA_TRUTH_SCRIPT:-$REPO/scripts/iitaca/prepare_iitaca_truth.jl}
NET=$REPO/results/network/II_L1000_D500/trainings/itrack_sparse_e100
echo "node $(hostname) commit $(cd "$REPO" && git rev-parse HEAD) script $SCRIPT"
julia --startup-file=no --project="$REPO" "$SCRIPT" "$REPO/config/iitaca/runs/II_L1000_D500_single_electrons.toml" "$NET/training.h5" "$OUT/truth_single_electrons.h5"
julia --startup-file=no --project="$REPO" "$SCRIPT" "$REPO/config/iitaca/runs/II_L1000_D500_bi214.toml" "$NET/scores_bi214.h5" "$OUT/truth_bi214.h5"
