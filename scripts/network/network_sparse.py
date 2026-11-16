"""Sparse 3D residual network: only the cells with charge are computed.

The stages are those of network_dense.Trunk with stem 'full'. Needs the library
spconv and a CUDA GPU (Hyperion, requirements-network-hyperion.txt); this is
the only file that imports spconv.
"""
import numpy as np,torch
from torch import nn
import spconv.pytorch as spconv
from network_events import batch_cells,padded

def sparse_batch(arm,rows,transform,device):
    """The same events as a sparse tensor: one row per hit, with its event index, cell indices and energy / scale."""
    b,idx,energy=batch_cells(arm,rows,transform);indices=torch.as_tensor(np.stack([b,*idx],axis=1).astype(np.int32),device=device)
    return spconv.SparseConvTensor(torch.as_tensor(energy,dtype=torch.float32,device=device).unsqueeze(1),indices,list(padded(arm['cells'])),len(rows))

class SparseBlock(nn.Module):
    """Block on the cells with charge: two 3x3x3 convolutions that keep the set of cells, with a skip connection."""
    def __init__(self,c,key):
            super().__init__();self.a=spconv.SubMConv3d(c,c,3,padding=1,bias=False,indice_key=key);self.na=nn.BatchNorm1d(c);self.b=spconv.SubMConv3d(c,c,3,padding=1,bias=False,indice_key=key);self.nb=nn.BatchNorm1d(c)
    def forward(self,x):
        y=self.a(x);y=self.b(y.replace_feature(torch.relu(self.na(y.features))));return y.replace_feature(torch.relu(x.features+self.nb(y.features)))

class Normalised(nn.Module):
    """Batch normalisation and ReLU of the features of a sparse tensor."""
    def __init__(self,c): super().__init__();self.n=nn.BatchNorm1d(c)
    def forward(self,x): return x.replace_feature(torch.relu(self.n(x.features)))

class SparseTrunk(nn.Module):
    """The stages of Trunk with stem 'full', computed only on cells with charge. Output: one vector per event.

    Inside a stage the convolutions keep the set of occupied cells (submanifold
    convolutions). The stride-2 convolution between stages occupies every
    coarse cell that has an occupied cell in its 3x3x3 input. The output is the
    mean over the occupied cells of the last stage.
    """
    def __init__(self,channels=(16,32,64,128)):
        super().__init__();layers=[spconv.SubMConv3d(1,channels[0],3,padding=1,bias=False,indice_key='stage0'),Normalised(channels[0]),SparseBlock(channels[0],'stage0')]
        for i,(a,b) in enumerate(zip(channels[:-1],channels[1:])): layers+=[spconv.SparseConv3d(a,b,3,stride=2,padding=1,bias=False),Normalised(b),SparseBlock(b,f'stage{i+1}')]
        self.layers=nn.Sequential(*layers);self.features=channels[-1]
    def forward(self,x):
        y=self.layers(x);event=y.indices[:,0].long();total=torch.zeros((y.batch_size,self.features),device=y.features.device,dtype=y.features.dtype).index_add_(0,event,y.features)
        return total/torch.bincount(event,minlength=y.batch_size).clamp(min=1).unsqueeze(1)
