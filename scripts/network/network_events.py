"""Events for the classification network: reading, merging and symmetry operations.

Events come from the files of network_events_builder.jl: for each track the
cell indices and energies of the saved hits inside the fixed frame. A batch is
filled into a dense grid (network_dense.py) or a sparse tensor
(network_sparse.py) on the training device.
"""
import h5py,numpy as np

ARMS=('itrack','etrack')
ENERGY_SCALE_KEV={'itrack':20.0,'etrack':10.0}   # typical hit energy; the grids hold energy / scale
def text(v): return v.decode() if isinstance(v,bytes) else v

class Sample:
    """Events of one species: the sparse hits of the tracks in the file (itrack, etrack or both) and the event table."""
    def __init__(self,path):
        with h5py.File(path) as f:
            if text(f['schema'][()])!='itaca.network_events.v1': raise ValueError(f'Unexpected schema in {path}')
            self.species=text(f['species'][()]);self.run_label=text(f['run_label'][()]);self.frame_side_mm=float(f['frame_side_mm'][()]);self.selection=text(f['selection'][()])
            self.events={k:f['events'][k][:] for k in f['events']}
            self.arms={a:{k:f[a][k][()] for k in f[a]} for a in ARMS if a in f}
        self.n=len(self.events['event_number'])

def merge(samples):
    """Events of several samples as one set: arms, event table and the index of the sample of each event.

    The arms are those present in every sample; the event table keeps the columns common to all samples.
    """
    arms={}
    for a in [a for a in ARMS if all(a in s.arms for s in samples)]:
        parts=[s.arms[a] for s in samples]
        for q in parts[1:]:
            if not (np.array_equal(q['cells'],parts[0]['cells']) and np.allclose(q['pitch_mm'],parts[0]['pitch_mm'])): raise ValueError(f'The {a} lattices of the samples differ')
        shift=np.cumsum([0]+[len(q['energy_keV']) for q in parts[:-1]])
        arms[a]={k:np.concatenate([q[k] for q in parts]) for k in ('ix','iy','iz','energy_keV')}
        arms[a]['offsets']=np.concatenate([parts[0]['offsets']]+[q['offsets'][1:]+d for q,d in zip(parts[1:],shift[1:])])
        arms[a]['cells']=parts[0]['cells'];arms[a]['pitch_mm']=parts[0]['pitch_mm'];arms[a]['scale']=np.float32(ENERGY_SCALE_KEV[a])
    events={k:np.concatenate([s.events[k] for s in samples]) for k in samples[0].events if all(k in s.events for s in samples)}
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

def batch_cells(arm,rows,transform):
    """Hits of the given events after the transform: event index in the batch, the three cell indices in the padded grid, energy / scale."""
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
    return b,idx,arm['energy_keV'][pos]/arm['scale']
