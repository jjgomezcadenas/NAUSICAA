#!/bin/bash
# Chain the hit files of batches FIRST..LAST of a run card. Hyperion cluster.
#   sbatch --dependency=afterok:<array job> --output=... itaca_chain.slurm RUN_CARD FIRST LAST
#SBATCH --job-name=itaca_chain
#SBATCH --partition=general
#SBATCH --qos=serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --time=02:00:00
set +e
type module >/dev/null 2>&1 || source /etc/profile
set -eo pipefail
module load Julia/1.11.6-linux-x86_64
REPO=${NAUSICAA_REPO:-/scratch/jjgomez/NAUSICAA}
export JULIA_DEPOT_PATH=${NAUSICAA_JULIA_DEPOT:-/scratch/jjgomez/julia-depot}
: "${NAUSICAA_DATA:?Set NAUSICAA_DATA}"
cd "$REPO"
echo "node $(hostname) commit $(git rev-parse HEAD) card $1 batches $2..$3"
julia --startup-file=no --project=. scripts/production/run_batches.jl "$1" "$2" "$3"
