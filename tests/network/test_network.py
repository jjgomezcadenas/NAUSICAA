"""Checks of the network data and measures: python tests/network/test_network.py"""
import sys
from pathlib import Path
import numpy as np,torch
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/network'))
from network_events import Transform,padded
from network_dense import dense_batch
from network_classifier import Classifier,separation,sparse_available

def arm(cells):
    # Three events with 2, 0 and 3 hits
    return dict(offsets=np.array([0,2,2,5]),ix=np.array([0,4,1,2,3],np.uint8),iy=np.array([0,1,0,5,2],np.uint8),iz=np.array([0,2,7,1,9],np.uint8),
                energy_keV=np.array([10,20,30,40,50],np.float32),cells=np.array(cells),scale=np.float32(10))

def test_dense():
    a=arm((6,6,10));rows=np.array([2,0,1]);x=dense_batch(a,rows,Transform(3),'cpu').numpy()
    assert x.shape==(3,1,8,8,16) and padded((6,6,10))==(8,8,16)
    s=(1,1,3)                                   # centred in the padding
    assert x[0,0,1+s[0],0+s[1],7+s[2]]==3 and x[0,0,2+s[0],5+s[1],1+s[2]]==4 and x[0,0,3+s[0],2+s[1],9+s[2]]==5
    assert x[1,0,s[0],s[1],s[2]]==1 and x[1,0,4+s[0],1+s[1],2+s[2]]==2 and x[2].sum()==0
    assert np.isclose(x.sum(),15) and (x>0).sum()==5

def test_symmetries():
    a=arm((6,6,10));rows=np.array([0]);t=Transform(1);t.swap[:]=True;t.flip[0]=(True,False,True);x=dense_batch(a,rows,t,'cpu').numpy()[0,0]
    # hit (0,0,0): swap -> (0,0,0); flip x,z -> (5,0,9); hit (4,1,2): swap -> (1,4,2); flip -> (4,4,7)
    assert x[5+1,0+1,9+3]==1 and x[4+1,4+1,7+3]==2 and (x>0).sum()==2
    # every operation keeps the hits and their energies
    rng=np.random.default_rng(0)
    for _ in range(50):
        y=dense_batch(a,np.array([0,1,2]),Transform(3,rng),'cpu').numpy();assert np.allclose(np.sort(y[y>0]),[1,2,3,4,5])
    rect=arm((6,5,10));t=Transform(1);t.swap[:]=True
    try: dense_batch(rect,rows,t,'cpu');raise AssertionError('swap on a rectangular lattice')
    except ValueError: pass

def test_network():
    torch.manual_seed(0);g={'itrack':torch.rand(2,1,16,16,16),'etrack':torch.rand(2,1,16,16,24)}
    for arms in (('itrack',),('etrack',),('itrack','etrack')):
        m=Classifier(arms,channels=(4,8,16,32));z=m(g);assert z.shape==(2,) and torch.isfinite(z).all()
        z.sum().backward();assert all(p.grad is not None for p in m.parameters())

def test_sparse():
    """The sparse tensor holds the hits of the dense grid; the sparse network runs. Needs spconv (cluster)."""
    if not sparse_available(): return 'skipped: needs spconv and a CUDA GPU (Hyperion)'
    from network_sparse import sparse_batch
    device='cuda';a=arm((6,6,10));rows=np.array([2,0,1]);rng=np.random.default_rng(3);state=rng.bit_generator.state
    x=sparse_batch(a,rows,Transform(3,rng),device);rng.bit_generator.state=state;d=dense_batch(a,rows,Transform(3,rng),device)
    assert x.batch_size==3 and tuple(x.spatial_shape)==(8,8,16) and torch.equal(x.dense(),d)
    torch.manual_seed(0);m=Classifier(('itrack',),channels=(4,8,16,32),sparse=True).to(device);z=m({'itrack':x})
    assert z.shape==(3,) and torch.isfinite(z).all()          # the event without hits gives a finite output
    z.sum().backward();assert all(p.grad is not None for p in m.parameters())
    assert sum(p.numel() for p in m.parameters())==sum(p.numel() for p in Classifier(('itrack',),channels=(4,8,16,32),stem='full').parameters())

def test_separation():
    rng=np.random.default_rng(1);s=rng.normal(1,1,20000);b=rng.normal(-1,1,20000);m=separation(np.r_[s,b],np.r_[np.ones(20000,bool),np.zeros(20000,bool)])
    assert abs(m['area']-0.9214)<0.005            # Phi(2/sqrt(2))
    k=list(m['target_signal_efficiency']).index(0.8);assert abs(m['signal_efficiency'][k]-0.8)<1e-3 and abs(m['background_acceptance'][k]-0.1232)<0.01   # 1-Phi(1.158)
    same=separation(np.r_[np.ones(10),np.ones(10)],np.r_[np.ones(10,bool),np.zeros(10,bool)]);assert same['area']==0.5

if __name__=='__main__':
    for name,f in list(globals().items()):
        if name.startswith('test_'): print(f() or 'passed',name)
