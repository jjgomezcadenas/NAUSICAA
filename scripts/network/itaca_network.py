"""Data, network and measures for the bb0nu / single-electron classification.

Events come from the files of build_network_events.jl: for each track the cell
indices and energies of the saved hits inside the fixed frame. A batch is
filled into a dense grid on the training device. The network is a 3D residual
network; with two tracks, one trunk per track joined before the output.
"""
import h5py,numpy as np,torch
from torch import nn

ARMS=('itrack','etrack')
ENERGY_SCALE_KEV={'itrack':20.0,'etrack':10.0}   # typical hit energy; the grids hold energy / scale
def text(v): return v.decode() if isinstance(v,bytes) else v

class Sample:
    """Events of one species: the sparse hits of both tracks and the event table."""
    def __init__(self,path):
        with h5py.File(path) as f:
            if text(f['schema'][()])!='itaca.network_events.v1': raise ValueError(f'Unexpected schema in {path}')
            self.species=text(f['species'][()]);self.run_label=text(f['run_label'][()]);self.frame_side_mm=float(f['frame_side_mm'][()])
            self.events={k:f['events'][k][:] for k in f['events']}
            self.arms={a:{k:f[a][k][()] for k in f[a]} for a in ARMS}
        self.n=len(self.events['event_number'])

def merge(samples):
    """Events of several samples as one set: arms, event table and the index of the sample of each event."""
    arms={}
    for a in ARMS:
        parts=[s.arms[a] for s in samples]
        for q in parts[1:]:
            if not (np.array_equal(q['cells'],parts[0]['cells']) and np.allclose(q['pitch_mm'],parts[0]['pitch_mm'])): raise ValueError(f'The {a} lattices of the samples differ')
        shift=np.cumsum([0]+[len(q['energy_keV']) for q in parts[:-1]])
        arms[a]={k:np.concatenate([q[k] for q in parts]) for k in ('ix','iy','iz','energy_keV')}
        arms[a]['offsets']=np.concatenate([parts[0]['offsets']]+[q['offsets'][1:]+d for q,d in zip(parts[1:],shift[1:])])
        arms[a]['cells']=parts[0]['cells'];arms[a]['pitch_mm']=parts[0]['pitch_mm'];arms[a]['scale']=np.float32(ENERGY_SCALE_KEV[a])
    events={k:np.concatenate([s.events[k] for s in samples]) for k in samples[0].events}
    return arms,events,np.concatenate([np.full(s.n,i) for i,s in enumerate(samples)])

def padded(cells,multiple=8):
    """Grid size: the lattice cells rounded up to a multiple of 8."""
    return tuple(int(-(-int(c)//multiple)*multiple) for c in cells)

class Transform:
    """One operation per event on the cell indices: transverse rotation or reflection, Z reflection and a shift.

    swap exchanges X and Y; flip reverses an axis. Together they give the eight
    symmetries of the square in the transverse plane and the Z reflection. The
    shift places the lattice inside the padded grid; the default centres it.
    """
    def __init__(self,n,rng=None):
        self.swap=np.zeros(n,bool);self.flip=np.zeros((n,3),bool);self.position=None
        if rng is not None:
            self.swap=rng.random(n)<0.5;self.flip=rng.random((n,3))<0.5;self.position=rng.random((n,3))

def dense_batch(arm,rows,transform,device):
    """Dense grid (events, 1, nx, ny, nz) of the hit energies of the given events, divided by the energy scale."""
    off=arm['offsets'];cells=arm['cells'].astype(int);grid=padded(cells);start=off[rows];length=off[rows+1]-start;total=int(length.sum())
    first=np.cumsum(length)-length;b=np.repeat(np.arange(len(rows)),length);pos=np.arange(total)-np.repeat(first,length)+np.repeat(start,length)
    idx=[arm[k][pos].astype(np.int64) for k in ('ix','iy','iz')]
    if cells[0]==cells[1]:
        s=transform.swap[b];idx[0],idx[1]=np.where(s,idx[1],idx[0]),np.where(s,idx[0],idx[1])
    elif transform.swap.any(): raise ValueError('X and Y can be exchanged only on a square lattice')
    for i in range(3):
        idx[i]=np.where(transform.flip[b,i],cells[i]-1-idx[i],idx[i]);slack=grid[i]-cells[i]
        shift=np.full(len(rows),slack//2) if transform.position is None else np.minimum((transform.position[:,i]*(slack+1)).astype(np.int64),slack)
        idx[i]=idx[i]+shift[b]
    x=torch.zeros((len(rows),1,*grid),device=device)
    x[torch.as_tensor(b,device=device),0,torch.as_tensor(idx[0],device=device),torch.as_tensor(idx[1],device=device),torch.as_tensor(idx[2],device=device)]=torch.as_tensor(arm['energy_keV'][pos]/arm['scale'],device=device)
    return x

class Block(nn.Module):
    """Two 3x3x3 convolutions with a skip connection."""
    def __init__(self,c):
        super().__init__();self.a=nn.Conv3d(c,c,3,padding=1,bias=False);self.na=nn.BatchNorm3d(c);self.b=nn.Conv3d(c,c,3,padding=1,bias=False);self.nb=nn.BatchNorm3d(c)
    def forward(self,x): return torch.relu(x+self.nb(self.b(torch.relu(self.na(self.a(x))))))

class Trunk(nn.Module):
    """Four stages; each after the first halves the grid and doubles the channels. Output: one vector per event.

    stem 'stride2': the first convolution (5x5x5) already halves the grid, so the
    first stage works at half resolution. stem 'full': the first stage works on
    the input lattice (3x3x3 convolution), at about six times the computing time.
    """
    def __init__(self,channels=(16,32,64,128),stem='stride2'):
        if stem not in ('stride2','full'): raise ValueError(f'Unknown stem {stem}')
        first=nn.Conv3d(1,channels[0],5,stride=2,padding=2,bias=False) if stem=='stride2' else nn.Conv3d(1,channels[0],3,padding=1,bias=False)
        super().__init__();layers=[first,nn.BatchNorm3d(channels[0]),nn.ReLU(),Block(channels[0])]
        for a,b in zip(channels[:-1],channels[1:]): layers+=[nn.Conv3d(a,b,3,stride=2,padding=1,bias=False),nn.BatchNorm3d(b),nn.ReLU(),Block(b)]
        self.layers=nn.Sequential(*layers);self.features=channels[-1]
    def forward(self,x): return self.layers(x).mean(dim=(2,3,4))

class Classifier(nn.Module):
    """One trunk per input track; the vectors are concatenated and reduced to one output (logit of bb0nu)."""
    def __init__(self,arms,channels=(16,32,64,128),stem='stride2'):
        super().__init__();self.arms=tuple(arms);self.trunks=nn.ModuleDict({a:Trunk(channels,stem) for a in self.arms});n=channels[-1]*len(self.arms)
        self.head=nn.Linear(n,1) if len(self.arms)==1 else nn.Sequential(nn.Linear(n,64),nn.ReLU(),nn.Linear(64,1))
    def forward(self,grids): return self.head(torch.cat([self.trunks[a](grids[a]) for a in self.arms],dim=1)).squeeze(1)

def separation(score,label,targets=(0.9,0.8,0.7,0.6,0.5)):
    """Area under the curve and, for each signal efficiency, the threshold, the background acceptance and its binomial error.

    An event is accepted when score > threshold. The threshold of a target is
    the score quantile that keeps that fraction of the signal.
    """
    score=np.asarray(score,float);label=np.asarray(label,bool);s=score[label];b=score[~label];order=np.argsort(np.r_[s,b],kind='stable')
    rank=np.empty(len(order));rank[order]=np.arange(1,len(order)+1)
    # average ranks of tied scores
    allv=np.r_[s,b];u,inv,cnt=np.unique(allv,return_inverse=True,return_counts=True);mean_rank=np.bincount(inv,rank)/cnt;rank=mean_rank[inv]
    area=(rank[:len(s)].sum()-len(s)*(len(s)+1)/2)/(len(s)*len(b))
    thr=np.array([np.quantile(s,1-t) for t in targets]);eff=np.array([(s>t).mean() for t in thr]);acc=np.array([(b>t).mean() for t in thr])
    return dict(area=float(area),target_signal_efficiency=np.array(targets),threshold=thr,signal_efficiency=eff,background_acceptance=acc,
                background_acceptance_error=np.sqrt(acc*(1-acc)/len(b)),signal_events=len(s),background_events=len(b))
