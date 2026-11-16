"""Dense 3D residual network: every cell of the grid is computed. Runs on any device."""
import torch
from torch import nn
from network_events import batch_cells,padded

def dense_batch(arm,rows,transform,device):
    """Dense grid (events, 1, nx, ny, nz) of the hit energies of the given events, divided by the energy scale."""
    b,idx,energy=batch_cells(arm,rows,transform);x=torch.zeros((len(rows),1,*padded(arm['cells'])),device=device)
    x[torch.as_tensor(b,device=device),0,torch.as_tensor(idx[0],device=device),torch.as_tensor(idx[1],device=device),torch.as_tensor(idx[2],device=device)]=torch.as_tensor(energy,device=device)
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
