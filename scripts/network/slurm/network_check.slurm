#!/bin/bash
# Check of the network software on one GPU: the tests, a short training of each
# first-stage option and the time per training step. Hyperion cluster.
#   sbatch --qos=test --time=00:10:00 --output=LOGDIR/%x.%j.out network_check.slurm SIGNAL.h5 BACKGROUND.h5 SCRATCH_DIR
#SBATCH --job-name=network_check
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --gres=gpu:1
set +e
type module >/dev/null 2>&1 || source /etc/profile
set -eo pipefail
module purge
module load PyTorch/2.1.2-foss-2021a-CUDA-12.1.1
# The environment adds h5py and spconv on top of the module's Python (requirements-network-hyperion.txt)
PY=${NAUSICAA_NETWORK_PYTHON:-/scratch/jjgomez/network-venv/bin/python}
cd "${NAUSICAA_REPO:-/scratch/jjgomez/NAUSICAA}"
echo "node $(hostname) commit $(git rev-parse HEAD) gpu $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | head -1)"
$PY tests/network/test_network.py
$PY - <<'P'
import sys,time,torch
sys.path.insert(0,'scripts/network')
from network_classifier import Classifier
for arms,shape in ((('itrack',),{'itrack':(56,56,56)}),(('etrack',),{'etrack':(48,48,64)})):
    for stem in ('stride2','full'):
        m=Classifier(arms,stem=stem).cuda();opt=torch.optim.AdamW(m.parameters());g={a:torch.rand(64,1,*s,device='cuda') for a,s in shape.items()};y=torch.ones(64,device='cuda');ts=[]
        for n in range(12):
            torch.cuda.synchronize();t=time.time();opt.zero_grad();l=torch.nn.functional.binary_cross_entropy_with_logits(m(g),y);l.backward();opt.step();torch.cuda.synchronize();ts.append(time.time()-t)
        print(f'TIMING {arms[0]} {stem}: {1000*sorted(ts)[len(ts)//2]:.0f} ms per step of 64 events, peak memory {torch.cuda.max_memory_allocated()/2**30:.1f} GiB',flush=True)
        del m,opt,g;torch.cuda.empty_cache();torch.cuda.reset_peak_memory_stats()
P
$PY scripts/network/network_trainer.py "$1" "$2" "$3" --input=both --epochs=1 --fraction=0.01 --validation-batches=5 --test-batches=0
