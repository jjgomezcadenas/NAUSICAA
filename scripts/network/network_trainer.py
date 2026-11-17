#!/usr/bin/env python3
"""Train and test the classification network.

  python scripts/network/network_trainer.py SIGNAL.h5 BACKGROUND.h5 OUTPUT_DIR --input=itrack|etrack|both
         [--epochs=20] [--batch-size=64] [--learning-rate=1e-3] [--weight-decay=1e-4] [--fraction=1.0]
         [--validation-batches=8] [--test-batches=9] [--seed=1] [--stem=stride2|full] [--device=mps|cuda|cpu]
         [--model=dense|sparse] [--channels=16,32,64,128]

SIGNAL and BACKGROUND are the files of network_events_builder.jl. Events are
split by production batch: the validation and test batches are held out and
the others train. --fraction keeps that part of the training events (learning
curve). --stem selects the resolution of the first stage (network_dense.Trunk).
--model=sparse computes only on the cells with charge, at full resolution
(network_sparse.SparseTrunk; needs spconv and a CUDA GPU, so it runs on
Hyperion). --channels gives the number of feature maps of each of the four
stages.
Each epoch trains with random symmetry operations and shifts, then
scores the validation events without them. The weights with the lowest
validation loss score the test events; the files are rewritten after every
epoch, so an interrupted run keeps its selected weights and their test scores.

Writes in OUTPUT_DIR: training.h5 (settings, the history per epoch, and for
each test event its number, label, score, drift, reconstructed energies and
the classical E2 of each track) and model.pt (the selected weights).
"""
import argparse,hashlib,subprocess,sys,time
from pathlib import Path
import h5py,numpy as np,torch
sys.path.insert(0,str(Path(__file__).parent))
from network_events import ARMS,ENERGY_SCALE_KEV,Sample,Transform,merge
from network_dense import dense_batch
from network_classifier import Classifier,separation,sparse_available

def batches_of(text): return [int(v) for v in str(text).split(',') if v!='']

def commit():
    """Commit of the repository that holds this script, with '+changes' when tracked files differ from it."""
    def git(*a): return subprocess.run(('git','-C',str(Path(__file__).parent),*a),capture_output=True,text=True,check=True).stdout.strip()
    try: return git('rev-parse','HEAD')+('+changes' if git('status','--porcelain','--untracked-files=no') else '')
    except (OSError,subprocess.CalledProcessError): return 'unknown'

def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('signal');p.add_argument('background');p.add_argument('output');p.add_argument('--input',required=True,choices=('itrack','etrack','both'))
    p.add_argument('--epochs',type=int,default=20);p.add_argument('--batch-size',type=int,default=64);p.add_argument('--learning-rate',type=float,default=1e-3)
    p.add_argument('--weight-decay',type=float,default=1e-4);p.add_argument('--fraction',type=float,default=1.0);p.add_argument('--validation-batches',default='8')
    p.add_argument('--test-batches',default='9');p.add_argument('--seed',type=int,default=1);p.add_argument('--stem',default='stride2',choices=('stride2','full'))
    p.add_argument('--model',default='dense',choices=('dense','sparse'));p.add_argument('--channels',default='16,32,64,128')
    p.add_argument('--device',default='cuda' if torch.cuda.is_available() else 'mps' if torch.backends.mps.is_available() else 'cpu')
    a=p.parse_args();out=Path(a.output);target=out/'training.h5'
    if target.exists(): raise SystemExit(f'Preserve existing training: {target}')
    channels=tuple(batches_of(a.channels));sparse=a.model=='sparse'
    if len(channels)!=4 or min(channels)<1: raise SystemExit('--channels takes four positive numbers')
    if sparse and not sparse_available(): raise SystemExit('The sparse model needs spconv and a CUDA GPU; it runs on Hyperion (requirements-network-hyperion.txt).')
    if sparse: a.stem='full';from network_sparse import sparse_batch
    out.mkdir(parents=True,exist_ok=True);arms=ARMS if a.input=='both' else (a.input,);device=torch.device(a.device)
    torch.manual_seed(a.seed);rng=np.random.default_rng(a.seed)
    samples=[Sample(a.signal),Sample(a.background)]
    if samples[0].species==samples[1].species: raise SystemExit('Signal and background are the same species')
    if samples[0].selection!=samples[1].selection: raise SystemExit(f'The two files were selected differently:\n  {samples[0].selection}\n  {samples[1].selection}')
    data,events,which=merge(samples);batch=events['batch'];label=(which==0).astype(np.float32)   # label 1 for the signal sample
    missing=[arm for arm in arms if arm not in data]
    if missing: raise SystemExit(f'The files hold no {", ".join(missing)}: choose --input among {", ".join(data)}')
    val_b,test_b=batches_of(a.validation_batches),batches_of(a.test_batches)
    if set(val_b)&set(test_b): raise SystemExit('Validation and test batches overlap')
    is_val=np.isin(batch,val_b);is_test=np.isin(batch,test_b);train=np.flatnonzero(~is_val&~is_test);val=np.flatnonzero(is_val);test=np.flatnonzero(is_test)
    if not (len(train) and len(val) and len(test)): raise SystemExit('An empty set: check the batch numbers')
    if a.fraction<1: train=np.sort(rng.choice(train,int(round(a.fraction*len(train))),replace=False))
    print(f'device {device}; model {a.model}; channels {channels}; input {a.input}; train {len(train)} validation {len(val)} test {len(test)} events',flush=True)

    def grids(ids,augment):
        """Grids of the events `ids` (dense, or sparse tensors); the same symmetry operation is applied to both tracks of an event."""
        t=Transform(len(ids),rng if augment else None);return {arm:(sparse_batch if sparse else dense_batch)(data[arm],ids,t,device) for arm in arms}

    model=Classifier(arms,channels,a.stem,sparse).to(device);parameters=sum(v.numel() for v in model.parameters())
    steps=a.epochs*int(np.ceil(len(train)/a.batch_size));optimizer=torch.optim.AdamW(model.parameters(),lr=a.learning_rate,weight_decay=a.weight_decay)
    schedule=torch.optim.lr_scheduler.OneCycleLR(optimizer,max_lr=a.learning_rate,total_steps=steps);loss_fn=torch.nn.BCEWithLogitsLoss()

    def score(ids):
        model.eval();scores=np.empty(len(ids),np.float32);total=0.0
        with torch.no_grad():
            for i in range(0,len(ids),4*a.batch_size):
                part=ids[i:i+4*a.batch_size];z=model(grids(part,False));total+=float(torch.nn.functional.binary_cross_entropy_with_logits(z,torch.as_tensor(label[part],device=device),reduction='sum'))
                scores[i:i+len(part)]=z.float().cpu().numpy()
        return scores,total/len(ids)

    def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    hashes=(sha(a.signal),sha(a.background));code=commit()
    def save(history,best_epoch,ts,tl,m):
        """Write training.h5 with the history so far and the test scores of the selected weights."""
        tmp=out/'training.h5.tmp'
        with h5py.File(tmp,'w') as f:
            f['schema']='itaca.network_training.v1';f['input']=a.input;f['signal_species']=samples[0].species;f['background_species']=samples[1].species
            f['signal_run_label']=samples[0].run_label;f['background_run_label']=samples[1].run_label;f['signal_selection']=samples[0].selection;f['background_selection']=samples[1].selection;f['commit']=code;f['signal_sha256']=hashes[0];f['background_sha256']=hashes[1]
            for key,v in dict(epochs=a.epochs,completed_epochs=len(history['train_loss']),batch_size=a.batch_size,learning_rate=a.learning_rate,weight_decay=a.weight_decay,fraction=a.fraction,seed=a.seed,
                              parameters=parameters,best_epoch=best_epoch+1,train_events=len(train),validation_events=len(val),test_events=len(test),time_s=time.time()-started).items(): f[f'settings/{key}']=v
            f['settings/validation_batches']=np.array(val_b);f['settings/test_batches']=np.array(test_b);f['settings/stem']=a.stem;f['settings/model']=a.model;f['settings/channels']=np.array(channels);f['settings/device']=str(device);f['settings/torch']=str(torch.__version__)
            for arm in arms: f[f'settings/energy_scale_keV/{arm}']=ENERGY_SCALE_KEV[arm]
            for key,v in history.items(): f[f'history/{key}']=np.array(v)
            f['test/event_number']=events['event_number'][test];f['test/batch']=batch[test];f['test/label']=label[test].astype(np.int8);f['test/score']=ts;f['test/loss']=tl
            for name in ('ion_drift_cm','electron_drift_cm','ion_energy_keV','electron_energy_keV','energy_keV','itrack_E2_keV','etrack_E2_keV','itrack_E1_keV','etrack_E1_keV','itrack_length_mm','etrack_length_mm'):
                if name in events: f[f'test/{name}']=events[name][test]   # the columns of the files; iItaca files have no etrack
            for key,v in m.items(): f[f'test/separation/{key}']=v
        tmp.replace(target)

    history={k:[] for k in ('train_loss','validation_loss','validation_area','validation_background_at_80','epoch_s')};best=np.inf;best_epoch=-1;started=time.time();result=None
    for epoch in range(a.epochs):
        t0=time.time();model.train();order=rng.permutation(train);running=0.0
        for i in range(0,len(order),a.batch_size):
            part=order[i:i+a.batch_size]
            if len(part)<2: continue
            optimizer.zero_grad();loss=loss_fn(model(grids(part,True)),torch.as_tensor(label[part],device=device));loss.backward();optimizer.step();schedule.step();running+=loss.item()*len(part)
        vs,vl=score(val);m=separation(vs,label[val]>0.5);k=int(np.argmin(np.abs(m['target_signal_efficiency']-0.8)))
        for key,v in (('train_loss',running/len(order)),('validation_loss',vl),('validation_area',m['area']),('validation_background_at_80',m['background_acceptance'][k]),('epoch_s',time.time()-t0)): history[key].append(v)
        print(f"epoch {epoch+1}/{a.epochs}: train loss {history['train_loss'][-1]:.4f}, validation loss {vl:.4f}, area {m['area']:.4f}, background at 80% {100*m['background_acceptance'][k]:.2f}%, {time.time()-t0:.0f} s",flush=True)
        if vl<best:      # these weights are the selection so far: keep them and their test scores
            best=vl;best_epoch=epoch;torch.save(model.state_dict(),out/'model.pt.tmp');(out/'model.pt.tmp').replace(out/'model.pt');ts,tl=score(test);result=(ts,tl,separation(ts,label[test]>0.5))
        save(history,best_epoch,*result)
    ts,tl,m=result;k=int(np.argmin(np.abs(m['target_signal_efficiency']-0.8)))
    print(f"test: area {m['area']:.4f}, background at 80% signal {100*m['background_acceptance'][k]:.2f} ± {100*m['background_acceptance_error'][k]:.2f}% ({m['signal_events']} signal, {m['background_events']} background events; weights of epoch {best_epoch+1})",flush=True)
    print(target)

if __name__=='__main__': main()
