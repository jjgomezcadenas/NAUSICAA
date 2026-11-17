#!/bin/bash
# One step of the eITACA analysis: a Julia script of the repository with its arguments. Hyperion cluster.
#   sbatch --output=LOGDIR/%x.%j.out --error=LOGDIR/%x.%j.err eitaca_analysis.slurm SCRIPT [ARGS...]
# SCRIPT is relative to the repository NAUSICAA_REPO, for instance
#   scripts/eitaca/prepare_eitaca_features.jl results/eitaca/IE_L1500_D400_bi214      (features of every batch, statistics)
#   scripts/eitaca/build_eitaca_network_events.jl results/eitaca/IE_L1500_D400_bb0nu   (training events)
#   scripts/eitaca/prepare_eitaca_truth.jl CARD SCORES OUTPUT.h5                       (needs NAUSICAA_DATA, read-only)
#SBATCH --job-name=eitaca_analysis
#SBATCH --partition=general
#SBATCH --exclude=hyperion-219   # this node fails to resolve the user: tasks there are requeued and held ("user env retrieval failed", 2026-11-17)
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
cd "${NAUSICAA_REPO:-/scratch/jjgomez/NAUSICAA}"
export JULIA_DEPOT_PATH=${NAUSICAA_JULIA_DEPOT:-/scratch/jjgomez/julia-depot}
echo "node $(hostname) commit $(git rev-parse HEAD) script $1"
julia --startup-file=no --project=. "$@"
