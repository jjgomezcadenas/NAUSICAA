#!/usr/bin/env python3
"""Apply a trained classification network to other samples.

  python scripts/network/network_scorer.py TRAINING_DIR SAMPLE.h5 [SAMPLE.h5 ...]
  python scripts/network/network_scorer.py TRAINING_DIR SAMPLE.h5 --other-selection=TAG

TRAINING_DIR holds training.h5 and model.pt of network_trainer.py; each SAMPLE.h5
is a file of network_events_builder.jl or build_iitaca_network_events.jl,
selected as the training files were (the same selection text). For each sample:

 1. Rebuild the network of the training (input, model, channels, first stage) and load the selected weights.
 2. Score every event of the sample, without symmetry operations, as the test events of the training.
 3. Take the thresholds that keep 90, 80, 70, 60 and 50 % of the signal test events of the training.
 4. Give the fraction of the sample above each threshold, with its binomial error.

Writes TRAINING_DIR/scores_<species>.h5 (schema itaca.network_scores.v1): the
score of every event, the thresholds, the accepted fractions, and the identity
of the training and of the sample. An existing file is preserved.

--other-selection=TAG scores samples selected differently from the training
files (for instance on an energy window instead of the ROI): the thresholds
still come from the training's test events. The output is
TRAINING_DIR/scores_<species>_<TAG>.h5 and records both selection texts.
"""
import argparse,hashlib,sys
from pathlib import Path
import h5py,numpy as np,torch
sys.path.insert(0,str(Path(__file__).parent))
from network_events import Sample,Transform,merge,text
from network_classifier import Classifier,separation,sparse_available

def load_training(directory):
    with h5py.File(directory/'training.h5') as f:
        s=f['settings']
        return dict(input=text(f['input'][()]),model=text(s['model'][()]) if 'model' in s else 'dense',stem=text(s['stem'][()]),
                    channels=tuple(int(c) for c in s['channels'][:]) if 'channels' in s else (16,32,64,128),
                    selection=text(f['signal_selection'][()]),signal=text(f['signal_species'][()]),background=text(f['background_species'][()]),
                    commit=text(f['commit'][()]),score=f['test/score'][:].astype(float),label=f['test/label'][:]>0)

def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('training');p.add_argument('samples',nargs='+');p.add_argument('--batch-size',type=int,default=256)
    p.add_argument('--other-selection',metavar='TAG',help='accept samples selected differently from the training files; write scores_<species>_<TAG>.h5')
    p.add_argument('--device',default='cuda' if torch.cuda.is_available() else 'mps' if torch.backends.mps.is_available() else 'cpu')
    a=p.parse_args();directory=Path(a.training);t=load_training(directory);sparse=t['model']=='sparse'
    if sparse and not sparse_available(): raise SystemExit('The sparse model needs spconv and a CUDA GPU; it runs on Hyperion (requirements-network-hyperion.txt).')
    if sparse: from network_sparse import sparse_batch
    else: from network_dense import dense_batch
    arms=('itrack','etrack') if t['input']=='both' else (t['input'],);device=torch.device(a.device)
    model=Classifier(arms,t['channels'],'full' if sparse else t['stem'],sparse).to(device)
    model.load_state_dict(torch.load(directory/'model.pt',map_location=device));model.eval()
    weights=hashlib.sha256((directory/'model.pt').read_bytes()).hexdigest()
    reference=separation(t['score'],t['label']);targets=reference['target_signal_efficiency'];thresholds=reference['threshold']
    print(f'training {directory} ({t["signal"]} against {t["background"]}, input {t["input"]}, model {t["model"]}); device {device}',flush=True)
    for path in a.samples:
        sample=Sample(path)
        if sample.selection!=t['selection'] and not a.other_selection: raise SystemExit(f'{path} was selected differently from the training files (use --other-selection):\n  {sample.selection}\n  {t["selection"]}')
        target=directory/(f'scores_{sample.species}_{a.other_selection}.h5' if a.other_selection else f'scores_{sample.species}.h5')
        if target.exists(): raise SystemExit(f'Preserve existing scores: {target}')
        data,events,_=merge([sample]);missing=[arm for arm in arms if arm not in data]
        if missing: raise SystemExit(f'{path} holds no {", ".join(missing)}')
        ids=np.arange(sample.n);scores=np.empty(sample.n,np.float32)
        with torch.no_grad():
            for i in range(0,sample.n,a.batch_size):
                part=ids[i:i+a.batch_size];tr=Transform(len(part))
                z=model({arm:(sparse_batch if sparse else dense_batch)(data[arm],part,tr,device) for arm in arms});scores[i:i+len(part)]=z.float().cpu().numpy()
        accepted=np.array([(scores>h).sum() for h in thresholds]);fraction=accepted/sample.n;error=np.sqrt(fraction*(1-fraction)/sample.n)
        tmp=target.with_suffix('.h5.tmp')
        with h5py.File(tmp,'w') as f:
            f['schema']='itaca.network_scores.v1';f['species']=sample.species;f['run_label']=sample.run_label;f['selection']=sample.selection;f['training_selection']=t['selection']
            f['sample_sha256']=hashlib.sha256(Path(path).read_bytes()).hexdigest();f['model_sha256']=weights;f['training_commit']=t['commit']
            f['events/event_number']=events['event_number'];f['events/score']=scores
            for k in ('ion_drift_cm','itrack_E2_keV','energy_keV'):
                if k in events: f[f'events/{k}']=events[k]
            f['target_signal_efficiency']=targets;f['threshold']=thresholds;f['signal_efficiency']=reference['signal_efficiency']
            f['accepted_events']=accepted;f['acceptance']=fraction;f['acceptance_error']=error
        tmp.replace(target)
        print(f'{sample.species}: {sample.n} events; '+'; '.join(f'{100*e:.0f} % signal: {100*x:.2f} ± {100*s:.2f} % ({k})' for e,x,s,k in zip(targets,fraction,error,accepted)),flush=True)
        print(target,flush=True)

if __name__=='__main__': main()
