"""Focused handoff and navigation checks with real events; run from the repo root."""
from pathlib import Path
import sys
import tempfile
import json
import h5py
import numpy as np
import plotly.graph_objects as go

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from event_viewer.server import Session


def verify_snapshot(session,payload):
    """Compare Julia output with original hits and independent NumPy histograms."""
    with h5py.File(session.work/'current.h5') as f:
        assert str(f['event/event_id'][()])==payload['event_id']
        source=f['context/source_path'].asstr()[()]
        with h5py.File(source) as src:
            row=int(f['event/position'][()])-1
            start=int(src['events/first_row'][row])-1
            stop=start+int(src['events/n_hits'][row])
            for col in ('x_mm','y_mm','z_mm','energy_MeV','time_ns','particle_id','hit_id'):
                assert np.array_equal(f['event/hits/'+col][:],src['hits/'+col][start:stop]),col
        e=f['event/hits/energy_MeV'][:]
        for plane in ('xy','xz','yz'):
            a,b=plane
            xe=f[f'display/axes/{a}_edges_mm'][:];ye=f[f'display/axes/{b}_edges_mm'][:]
            expected,_x,_y=np.histogram2d(f[f'event/hits/{a}_mm'][:],f[f'event/hits/{b}_mm'][:],bins=(xe,ye),weights=e.astype(float))
            actual=f[f'display/projections/{plane}_energy_MeV'][:]
            np.testing.assert_allclose(actual,expected.T,rtol=2e-6,atol=1e-9)
            np.testing.assert_allclose(actual.sum(dtype=float),e.sum(dtype=float),rtol=2e-6)
            assert int(f[f'display/projections/{plane}_hit_count'][:].sum())==len(e)
        for figure in payload['figures'].values():go.Figure(figure) # Validate Plotly schema.
        assert len(payload['figures']['xyz']['data'][0]['x'])==8*payload['occupied_cells']


with tempfile.TemporaryDirectory(prefix='viewer-check-',dir=ROOT/'work') as work:
    s=Session(ROOT/'config/event_viewer.toml',work)
    try:
        a=s.navigate({});assert a['sample']=='bb0nu' and a['position']==1
        verify_snapshot(s,a)
        b=s.navigate(dict(action='next'));assert b['position']==2
        # Event 2 has distant clusters: ~38 × 88 × 960 mm. Its full scene must
        # use the panel area with independent axis scales, retaining mm bounds.
        far=s.navigate(dict(action='id',event_id='2'));verify_snapshot(s,far)
        scene=far['figures']['xyz']['layout']['scene']
        spans=np.array([np.diff(scene[a+'axis']['range'])[0] for a in 'xyz'])
        assert spans.max()/spans.min()>20
        assert scene['aspectmode']=='cube'
        for plane in ('xy','xz','yz'):
            layout=far['figures'][plane]['layout']
            for axis,coordinate in zip(('xaxis','yaxis'),plane):
                assert 'scaleanchor' not in layout[axis]
                assert layout[axis]['range']==scene[coordinate+'axis']['range']
        s.navigate(dict(action='id',event_id=b['event_id']))
        c=s.navigate(dict(action='previous'));assert c['event_id']==a['event_id']
        before=(s.work/'current.h5').read_bytes()
        try:s.navigate(dict(action='id',event_id='999999999999999'))
        except ValueError:pass
        else:raise AssertionError('Missing event should fail')
        assert (s.work/'current.h5').read_bytes()==before
        d=s.navigate(dict(action='sample',sample='electrons'));verify_snapshot(s,d)
        e=s.navigate(dict(action='bins',nbins=200));verify_snapshot(s,e)
        assert e['event_id']==d['event_id'] and e['energy_MeV']==d['energy_MeV']
        view={'camera':{'eye':{'x':2,'y':1,'z':1}},
              'ranges':{'xy':{'xaxis':[0,10],'yaxis':[0,20]},
                        'xz':{'xaxis':[-40,40],'yaxis':[100,500]},
                        'yz':{'xaxis':[-20,60],'yaxis':[0,1000]}}}
        s.save_view(dict(sample=e['sample'],event_id=e['event_id'],view=view,opacity=.5))
        expected=(e['sample'],e['event_id'],e['nbins'])
        print('Both samples, navigation, failed-ID preservation, energy/count maps and Plotly schemas passed.')
    finally:s.close()
    restored=Session(ROOT/'config/event_viewer.toml',work)
    try:
        p=restored.navigate({})
        assert (p['sample'],p['event_id'],p['nbins'])==expected
        assert p['view']==view and p['opacity']==.5
        print('Session restart restores sample/event, binning, camera and opacity.')
    finally:restored.close()
