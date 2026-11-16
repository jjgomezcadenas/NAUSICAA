"""The classifier and its measures of separation."""
import numpy as np,torch
from torch import nn
from network_dense import Trunk

def sparse_available():
    """Whether the sparse network can run here: spconv installed and a CUDA GPU present."""
    try: import spconv.pytorch
    except ImportError: return False
    return torch.cuda.is_available()

class Classifier(nn.Module):
    """One trunk per input track; the vectors are concatenated and reduced to one output (logit of bb0nu).

    sparse=True uses network_sparse.SparseTrunk (full resolution) and takes the tensors of sparse_batch.
    """
    def __init__(self,arms,channels=(16,32,64,128),stem='stride2',sparse=False):
        super().__init__();self.arms=tuple(arms)
        if sparse:
            from network_sparse import SparseTrunk
            self.trunks=nn.ModuleDict({a:SparseTrunk(channels) for a in self.arms})
        else: self.trunks=nn.ModuleDict({a:Trunk(channels,stem) for a in self.arms})
        n=channels[-1]*len(self.arms)
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
