#!/bin/bash
# Run one NEXUS batch of the ITACA gas-cylinder simulation.
#   run_nexus_batch.sh SAMPLE BATCH EVENTS OUTDIR [SEED_OFFSET]
# SAMPLE: bb0nu | electron.  BATCH: integer >= 0.
# Event numbers start at BATCH*EVENTS; the seed is SEED_BASE + SEED_OFFSET +
# BATCH, with SEED_BASE 100000 (bb0nu) or 200000 (electron). SEED_OFFSET
# (default 0) makes a production independent of another: the 10 x 100 pilot
# used 0. Output:
#   OUTDIR/SAMPLE_bBBBBB.next.h5, its two macros and a log.
set -eo pipefail
SAMPLE=$1; BATCH=$2; EVENTS=$3; OUTDIR=$4; SEED_OFFSET=${5:-0}
NEXUS=${NEXUS_DIR:-/scratch/jjgomez/nexus}
BUILD=${NEXUS_BUILD:-$NEXUS/build-icelake}
case $SAMPLE in bb0nu) SEED_BASE=100000;; electron) SEED_BASE=200000;; *) echo "Unknown sample $SAMPLE"; exit 2;; esac
SEED=$((SEED_BASE + SEED_OFFSET + BATCH)); START=$((BATCH * EVENTS))
TAG=$(printf "%s_b%05d" "$SAMPLE" "$BATCH")
mkdir -p "$OUTDIR"
cp -n "$NEXUS/macros/geometries/ItacaChamber.mac" "$OUTDIR/" 2>/dev/null || true
OUT="$OUTDIR/$TAG.next"
[ -e "$OUT.h5" ] && { echo "Output exists: $OUT.h5"; exit 3; }

# Batch shells do not define the module command. The profile may print an
# error on nodes that cannot resolve the user name; that must not stop the job.
set +e
type module >/dev/null 2>&1 || source /etc/profile
set -e
type module >/dev/null 2>&1 || { echo "module command unavailable on $(hostname)"; exit 4; }
module purge
module load Geant4/11.1.2-GCC-11.3.0 HDF5/1.12.2-gompi-2022a GSL/2.7-GCC-11.3.0
cd "$(dirname "$(which geant4-config)")" && source geant4.sh
export LD_LIBRARY_PATH="$BUILD:${LD_LIBRARY_PATH:-}"
set -u

sed -e "s|@SEED@|$SEED|" -e "s|@START_ID@|$START|" -e "s|@OUTPUT@|$OUT|" \
    "$NEXUS/macros/ItacaChamber.$SAMPLE.config.mac" > "$OUTDIR/$TAG.config.mac"
sed -e "s|@CONFIG@|$OUTDIR/$TAG.config.mac|" \
    "$NEXUS/macros/ItacaChamber.$SAMPLE.init.mac" > "$OUTDIR/$TAG.init.mac"

cd "$NEXUS"   # macros use paths relative to the NEXUS directory
echo "nexus commit $(git -C "$NEXUS" rev-parse HEAD) node $(hostname) seed $SEED start_id $START events $EVENTS"
echo "ItacaChamber.cc sha256 $(sha256sum source/geometries/ItacaChamber.cc | cut -c1-64)"
SECONDS=0
"$BUILD/nexus" -b -n "$EVENTS" "$OUTDIR/$TAG.init.mac" > "$OUTDIR/$TAG.log" 2>&1
echo "TIMING sample=$SAMPLE batch=$BATCH events=$EVENTS wall_s=$SECONDS bytes=$(stat -c %s "$OUT.h5")"
