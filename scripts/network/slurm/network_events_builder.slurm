#!/bin/bash
# Training file of one run for the classification network. Hyperion cluster.
#   sbatch --output=LOGDIR/%x.%j.out --error=LOGDIR/%x.%j.err \
#          network_events_builder.slurm RUN_DIR OUTPUT.h5 [--roi=COMPARISON.h5]
# Paths are relative to the repository. RUN_DIR is results/readout/<run_label>.
#SBATCH --job-name=network_events
#SBATCH --partition=general
#SBATCH --qos=serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --time=01:00:00
set +e
type module >/dev/null 2>&1 || source /etc/profile
set -eo pipefail
module load Julia/1.11.6-linux-x86_64
export JULIA_DEPOT_PATH=${NAUSICAA_JULIA_DEPOT:-/scratch/jjgomez/julia-depot}
cd "${NAUSICAA_REPO:-/scratch/jjgomez/NAUSICAA}"
echo "node $(hostname) commit $(git rev-parse HEAD)"
julia --startup-file=no --project=. scripts/network/network_events_builder.jl "$@"
