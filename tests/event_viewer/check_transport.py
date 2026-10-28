"""Practical transported-hit checks against exact Julia HDF5 counts and state."""
from pathlib import Path
import sys,tempfile
import h5py
import numpy as np
import plotly.graph_objects as go
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from event_viewer.server import Session

def check(s,p):
    t=p['transport'];assert t['status']=='complete'
    assert 10<=t['drift_cm']<=150
    assert t['sampled_ions']==t['retained_ions']+t['outside_ions']
    with h5py.File(s.work/'current.h5') as f:
        counts=f['transport/bins/ions'][:]
        assert counts.sum()==t['retained_ions']
        for plane in ('xy','xz','yz'):
            projection=f['transport/display/projections/'+plane+'_ions'][:]
            assert projection.sum()==t['retained_ions']
            plotted=np.array(p['figures'][plane]['data'][0]['z'],dtype=float)
            np.testing.assert_array_equal(np.nan_to_num(plotted),projection)
        signature=tuple(f['transport/bins/'+a][:].tobytes() for a in ('x','y','z','ions'))
    for fig in p['figures'].values():go.Figure(fig)
    return signature

with tempfile.TemporaryDirectory(prefix='transport-check-',dir=ROOT/'work') as work:
    s=Session(ROOT/'config/event_viewer.toml',work)
    try:
        p=s.navigate(dict(action='tab',tab='transport'));signature=check(s,p)
        t=p['transport'];print('First transport:',t['drift_cm'],'cm;',t['retained_ions'],'retained ions;',p['occupied_cells'],'cells',flush=True)
        s.navigate(dict(action='tab',tab='energy'))
        assert check(s,s.navigate(dict(action='tab',tab='transport')))==signature
        s.navigate(dict(action='next'))
        assert check(s,s.navigate(dict(action='previous')))==signature
        p=s.navigate(dict(action='realisation'));assert p['transport']['realisation']==1
        assert check(s,p)!=signature
        p=s.navigate(dict(action='transport',drift_mode='fixed',fixed_cm=10));check(s,p)
        sigma=p['transport']['sigma_t_mm']
        p=s.navigate(dict(action='transport',drift_mode='fixed',fixed_cm=150));signature=check(s,p)
        np.testing.assert_allclose(p['transport']['sigma_t_mm']/sigma,np.sqrt(15),rtol=2e-6)
        assert p['transport']['drift_cm']==150
        try:s.navigate(dict(action='transport',drift_mode='fixed',fixed_cm=-10))
        except ValueError:pass
        else:raise AssertionError('Invalid distance accepted')
        assert s.navigate({})['transport']['drift_cm']==150
        p=s.navigate(dict(action='sample',sample='electrons'));check(s,p)
        signature=check(s,p)
    finally:s.close()
    s=Session(ROOT/'config/event_viewer.toml',work)
    try:
        p=s.navigate({});assert p['active_tab']=='transport' and p['sample']=='electrons'
        assert p['transport']['drift_cm']==150 and check(s,p)==signature
        print('Both species, ion projections/ledger, reproducibility, new realisation, fixed widths, invalid-distance recovery and restart passed.')
    finally:s.close()
