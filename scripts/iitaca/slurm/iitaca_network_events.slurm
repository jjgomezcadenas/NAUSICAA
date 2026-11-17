#!/bin/bash
# Statistics (HDF5 only) and training events of the network of one iItaca run. Hyperion cluster.
#   sbatch --output=LOGDIR/%x.%j.out --error=LOGDIR/%x.%j.err iitaca_network_events.slurm RUN_DIR
# RUN_DIR (results/iitaca/<run_label>, relative to the repository) holds the
# feature tables and hit files of the batches. Writes RUN_DIR/statistics/events.h5
# (prepare_iitaca_statistics.jl, read by prepare_iitaca_truth.jl) if absent and
# results/network/<variant>/events/<sample>/events.h5 (build_iitaca_network_events.jl).
#SBATCH --job-name=iitaca_network_events
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
cd "${NAUSICAA_REPO:-/scratch/jjgomez/NAUSICAA}"
export JULIA_DEPOT_PATH=${NAUSICAA_JULIA_DEPOT:-/scratch/jjgomez/julia-depot}
RUN=$1
echo "node $(hostname) commit $(git rev-parse HEAD) run $RUN"
[ -e "$RUN/statistics/events.h5" ] || julia --startup-file=no --project=. scripts/iitaca/prepare_iitaca_statistics.jl "$RUN"
julia --startup-file=no --project=. scripts/iitaca/build_iitaca_network_events.jl "$RUN"
