#!/bin/bash
# One input batch of an eITACA run card per array task: reduced hit table and
# electron file. Hyperion cluster.
#   sbatch --array=0-59 --output=LOGDIR/%x.%A_%a.out --error=LOGDIR/%x.%A_%a.err \
#          eitaca_batches.slurm RUN_CARD NEXUS_FILE_PATTERN
# RUN_CARD is relative to the repository. NEXUS_FILE_PATTERN is the NEXUS file
# of a batch with {batch} for the batch tag; its reduced hit table is written
# beside it if absent. NAUSICAA_DATA must be exported to the job (--export).
# Output: results/eitaca/<run_label>/batches/<tag>/electrons.h5 in the working
# copy NAUSICAA_REPO; a batch whose output exists is skipped.
#SBATCH --job-name=eitaca_batches
#SBATCH --partition=general
#SBATCH --exclude=hyperion-219   # this node fails to resolve the user: tasks there are requeued and held ("user env retrieval failed", 2026-11-17)
#SBATCH --qos=serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --time=04:00:00
set +e
type module >/dev/null 2>&1 || source /etc/profile
set -eo pipefail
module load Julia/1.11.6-linux-x86_64
REPO=${NAUSICAA_REPO:-/scratch/jjgomez/NAUSICAA}
export JULIA_DEPOT_PATH=${NAUSICAA_JULIA_DEPOT:-/scratch/jjgomez/julia-depot}
: "${NAUSICAA_DATA:?Set NAUSICAA_DATA}"
CARD=$1; PATTERN=$2; B=$SLURM_ARRAY_TASK_ID; TAG=$(printf "b%05d" "$B")
NEXUS_FILE=${PATTERN//\{batch\}/$TAG}
cd "$REPO"
echo "node $(hostname) batch $TAG commit $(git rev-parse HEAD) card $CARD"
[ -e "${NEXUS_FILE%.next.h5}_hitdf.h5" ] || julia --startup-file=no --project=. scripts/events/extract_hits.jl "$NEXUS_FILE"
BATCH_DIR=results/eitaca/$(sed -n 's/^run_label *= *"\(.*\)".*/\1/p' "$CARD")/batches/$TAG
[ -e "$BATCH_DIR/electrons.h5" ] || julia --startup-file=no --project=. scripts/eitaca/run_eitaca_batch.jl "$CARD" --batch="$B"
