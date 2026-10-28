"""Verify full-cloud ideal pad collection and the interactive pad/time display."""
from pathlib import Path
import sys,tempfile
import h5py
import numpy as np
import plotly.graph_objects as go
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from event_viewer.server import Session

def check(session,p):
    assert p['active_tab']=='pads';info=p['pads'];t=p['transport']
    assert info['total_ions']==t['retained_ions']
    assert all(0<=v<info['pitch_mm'] for v in info['grid_offset_mm'])
    with h5py.File(session.work/'current.h5') as f:
        b=f['transport/bins'];g=f['pads'];h=g['history']
        assert h['ions'][:].sum()==info['total_ions']
        assert np.allclose(np.diff(g['time_edges_ms'][:]),20,atol=.01)
        assert np.allclose(np.diff(g['display/axes/x_edges_mm'][:]),2,atol=.001)
        assert np.allclose(np.diff(g['display/axes/y_edges_mm'][:]),2,atol=.001)
        ix,iy,it=(h[k][:]-1 for k in ('pad_ix','pad_iy','slice_index'))
        # Independent charge assignment from saved diffusion centres and grid axes.
        edges=[g['display/axes/'+a+'_edges_mm'][:] for a in 'xy']+[g['z_edges_mm'][:]]
        idx=[np.minimum(np.searchsorted(e,b[a][:],side='right')-1,len(e)-2) for e,a in zip(edges,'xyz')]
        if info['drift_direction']==-1:idx[2]=len(edges[2])-2-idx[2]
        expected={}
        for key,n in zip(zip(*idx),b['ions'][:]):expected[key]=expected.get(key,0)+int(n)
        assert expected==dict(zip(zip(ix,iy,it),map(int,h['ions'][:])))
        for plane in ('xy','xz','yz'):
            gdisplay=g['display/projections/'+plane+'_ions'][:]
            assert gdisplay.sum()==info['total_ions']
            plotted=np.nan_to_num(np.array(p['figures'][plane]['data'][0]['z'],dtype=float))
            np.testing.assert_allclose(plotted,gdisplay*info['charge_per_ion_C']*1e15,rtol=1e-14,atol=1e-15)
        signature=h['ions'][:].tobytes()+ix.tobytes()+iy.tobytes()+it.tobytes()
    for fig in p['figures'].values():go.Figure(fig)
    assert p['figures']['xz']['layout']['yaxis']['title']['text']=='t [ms]'
    return signature

with tempfile.TemporaryDirectory(prefix='pad-check-',dir=ROOT/'work') as work:
    s=Session(ROOT/'config/event_viewer.toml',work)
    try:
        p=s.navigate(dict(action='sample',sample='bb0nu_100um',tab='pads'));sig=check(s,p)
        seed=p['transport']['seed']
        p=s.navigate(dict(action='tab',tab='transport'));assert p['transport']['seed']==seed
        assert check(s,s.navigate(dict(action='tab',tab='pads')))==sig
        p=s.navigate(dict(action='id',event_id='2'));check(s,p)
        p=s.navigate(dict(action='sample',sample='electrons'));sig=check(s,p)
        p=s.navigate(dict(action='realisation'));assert check(s,p)!=sig
        p=s.navigate(dict(action='transport',drift_mode='fixed',fixed_cm=150));sig=check(s,p)
        assert p['transport']['drift_cm']==150
    finally:s.close()
    s=Session(ROOT/'config/event_viewer.toml',work)
    try:
        assert check(s,s.navigate({}))==sig
        print('Both species: full ion collection, independent pad assignment, exact charge projections, 2 mm pitch/20 ms slices, shared transport, realisations and restart passed.')
    finally:s.close()
