"""Focused handoff and navigation checks with real events; run from the repo root."""
from pathlib import Path
import sys
import tempfile
import json
import os
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


def verify_graph(session,payload):
    """Check saved graph numbering, energy accounting and plotted path geometry."""
    graph=payload['graph'];assert graph['status']=='complete'
    with h5py.File(session.work/'current.h5') as f:
        assert graph['component_count']==len(f['reconstruction/components'])
        total=float(f['event/summary/total_energy_MeV'][()])
        energies=[]
        for c in graph['components']:
            g=f['reconstruction/components/'+str(c['component_id'])]
            n=c['node_count'];e=g['nodes/energy_MeV'][:]
            np.testing.assert_allclose(c['energy_MeV'],e.sum(dtype=float),rtol=2e-6)
            np.testing.assert_allclose(c['energy_fraction'],c['energy_MeV']/total,rtol=2e-6)
            energies.append(c['energy_MeV'])
            for name in ('edges/src','edges/dst','backbone/path_indices','backbone/endpoint_indices'):
                values=g[name][:]
                assert g[name].attrs['index_base']==1
                assert np.all((values>=1)&(values<=n))
            assert len(g['edges/src'])==len(g['edges/dst'])==c['edge_count']
            path=g['backbone/path_indices'][:]-1
            for plane in ('xy','xz','yz','xyz'):
                overlays=[t for t in graph['overlays'][plane] if t['meta']['component_id']==c['component_id']]
                backbone=next(t for t in overlays if t['meta']['layer']=='backbone')
                for target,coordinate in zip('xyz',plane):
                    np.testing.assert_array_equal(backbone[target],g['nodes/'+coordinate+'_mm'][:][path])
        assert energies==sorted(energies,reverse=True)
        np.testing.assert_allclose(sum(energies),graph['retained_energy_MeV'],rtol=2e-6)
        np.testing.assert_allclose(graph['retained_energy_MeV']+graph['below_threshold_energy_MeV'],total,rtol=2e-6)
    for plane,traces in graph['overlays'].items():
        go.Figure(data=payload['figures'][plane]['data']+traces,layout=payload['figures'][plane]['layout'])


with tempfile.TemporaryDirectory(prefix='viewer-check-',dir=ROOT/'work') as work:
    s=Session(ROOT/'config/event_viewer.toml',work)
    try:
        a=s.navigate({});assert a['sample']=='bb0nu' and a['position']==1
        verify_snapshot(s,a)
        b=s.navigate(dict(action='next'));assert b['position']==2
        # Event 2 has distant clusters: ~38 × 88 × 960 mm. Its full scene must
        # use the panel area with independent axis scales, retaining mm bounds.
        far=s.navigate(dict(action='id',event_id='2'));verify_snapshot(s,far)
        assert far['graph']['status']=='pending'
        far=s.navigate(dict(action='tab',tab='tracks'));verify_graph(s,far)
        assert far['active_tab']=='tracks' and far['graph']['component_count']>1
        build=far['graph']['build_request_id']
        rebin=s.navigate(dict(action='bins',nbins=200));verify_graph(s,rebin)
        assert rebin['graph']['build_request_id']==build
        far=s.navigate(dict(action='bins',nbins=100))
        assert s.navigate(dict(action='tab',tab='energy'))['graph']['build_request_id']==build
        assert s.navigate(dict(action='tab',tab='tracks'))['graph']['build_request_id']==build
        other=s.navigate(dict(action='next'));verify_graph(s,other)
        assert other['event_id']!=far['event_id'] and other['graph']['build_request_id']!=build
        far=s.navigate(dict(action='id',event_id='2'))
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
        d=s.navigate(dict(action='sample',sample='electrons'));verify_snapshot(s,d);verify_graph(s,d)
        e=s.navigate(dict(action='bins',nbins=200));verify_snapshot(s,e)
        assert e['event_id']==d['event_id'] and e['energy_MeV']==d['energy_MeV']
        view={'camera':{'eye':{'x':2,'y':1,'z':1}},'selected_component':1,'layers':{'edges':False},
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
        assert p['active_tab']=='tracks';verify_graph(restored,p)
        following=restored.navigate(dict(action='next'))
        assert following['view']=={'layers':view['layers']}
        assert following['active_tab']=='tracks'
        print('Session restart restores sample/event, binning, camera and opacity.')
    finally:restored.close()
    card=Path(work)/'empty_graph.toml'
    card.write_text((ROOT/'config/event_viewer.toml').read_text().replace('energy_threshold_keV = 1.0','energy_threshold_keV = 100000.0'))
    empty=Session(card,Path(work)/'empty')
    try:
        p=empty.navigate(dict(action='tab',tab='tracks'));verify_graph(empty,p)
        assert p['graph']['component_count']==0 and p['graph']['components']==[]
        assert all(not traces for traces in p['graph']['overlays'].values())
        assert p['hit_count']>0 and p['energy_MeV']>0
        print('MC graph HDF5/plot parity, cache reuse, event changes and empty reconstruction passed.')
    finally:empty.close()

# File opening: conversion, reuse, filenames with spaces, failure recovery and restore.
with tempfile.TemporaryDirectory(prefix='viewer-files-',dir=ROOT/'work') as work:
    directory=Path(work)/'event files';directory.mkdir()
    source=directory/'test event μ.next.h5'
    with h5py.File(Path(os.environ.get('NAUSICAA_RAW_DIRECTORY',str(Path.home()/'NAUSICAA')))/'0nubb.next.h5') as src,h5py.File(source,'w') as dst:
        rows=src['MC/hits'][:4096];event_id=rows['event_id'][0]
        assert rows['event_id'][-1]!=event_id
        dst['MC/hits']=rows[rows['event_id']==event_id]
        particles=src['MC/particles'][:]
        dst['MC/particles']=particles[particles['event_id']==event_id]
    card=Path(work)/'files.toml'
    card.write_text((ROOT/'config/event_viewer.toml').read_text().replace(
        'file_directory = "events"','file_directory = '+json.dumps(str(directory))))
    session_work=Path(work)/'session'
    s=Session(card,session_work)
    try:
        listing=s.list_files()
        assert listing['directory']==str(directory) and listing['entries'][0]['name']==source.name
        loaded=s.navigate(dict(action='file',path=source.name))
        assert loaded['event_count']==1 and loaded['sample'].startswith('file:')
        verify_snapshot(s,loaded)
        reduced=directory/'test event μ_hitdf.h5';assert reduced.is_file()
        assert loaded['file_path']==str(reduced)
        mtime=reduced.stat().st_mtime_ns
        assert s.navigate(dict(action='file',path=str(source)))['sample']==loaded['sample']
        assert reduced.stat().st_mtime_ns==mtime
        assert s.navigate(dict(action='file',path=str(reduced)))['event_id']==loaded['event_id']
        snapshot=(s.work/'current.h5').read_bytes()
        bad=directory/'unrelated.h5'
        with h5py.File(bad,'w') as f:f['x']=[1,2]
        for filename in ('absent.h5',str(bad)):
            try:s.navigate(dict(action='file',path=filename))
            except ValueError:pass
            else:raise AssertionError('Invalid event file should fail')
            assert (s.work/'current.h5').read_bytes()==snapshot
            assert s.navigate({})['sample']==loaded['sample']
        with h5py.File(source,'r+') as f:f.attrs['changed_for_test']=True
        try:s.navigate(dict(action='file',path=str(source)))
        except ValueError as exc:assert 'different source contents' in str(exc)
        else:raise AssertionError('Stale reduced table should fail')
        assert reduced.stat().st_mtime_ns==mtime
        assert (s.work/'current.h5').read_bytes()==snapshot
    finally:s.close()
    restored=Session(card,session_work)
    try:
        p=restored.navigate({});verify_snapshot(restored,p)
        assert p['sample']==loaded['sample'] and p['file_path']==str(reduced)
        original=restored.navigate(dict(action='sample',sample='bb0nu'))
        assert original['sample']=='bb0nu'
        assert restored.navigate(dict(action='sample',sample=loaded['sample']))['event_count']==1
        print('File browser, NEXUS conversion/reuse, arbitrary paths, failure preservation and file-session restore passed.')
    finally:restored.close()
