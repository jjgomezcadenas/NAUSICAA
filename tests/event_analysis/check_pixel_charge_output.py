"""Independently validate saved batch pixel/frame reductions with NumPy."""
from pathlib import Path
import h5py,numpy as np,tomllib,argparse
ROOT=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser()
parser.add_argument('--results-dir',type=Path,default=ROOT/'results/pixel_charge')
args=parser.parse_args()
for sample in ('bb0nu','electrons'):
    with h5py.File(args.results_dir/(sample+'.h5')) as f:
        e=f['events'];c=f['cells'];nmin=int(f['settings/n_min'][()])
        ids=c['event_id'][:];ions=c['n_ions'][:];keep=ions>0
        config=tomllib.loads(f['provenance/card_toml'].asstr()[()])
        with h5py.File(f['provenance/source_path'].asstr()[()]) as source:
            assert len(e['event_id'])==min(config['run']['n_events'],len(source['events/event_id']))
        np.testing.assert_array_equal(e['produced_ions'][:],e['retained_ions'][:]+e['outside_ions'][:])
        np.testing.assert_array_equal(e['collected_ions'][:],e['retained_ions'][:])
        starts=np.r_[0,np.flatnonzero(ids[1:]!=ids[:-1])+1];stops=np.r_[starts[1:],len(ids)]
        for i,(start,stop) in enumerate(zip(starts,stops)):
            assert e['event_id'][i]==ids[start]
            v=ions[start:stop];above=v[v>0]
            assert v.sum()==e['collected_ions'][i]
            assert len(above)==e['pulse_count'][i]
            assert bool(len(above))==bool(e['has_pulse'][i])
            if len(above):
                assert above.max()==e['maximum_ions'][i] and above.min()==e['minimum_ions'][i]
                assert above.sum()==e['pulse_ion_sum'][i] and len(above)==e['pulse_count'][i]
                np.testing.assert_allclose(above.mean(),e['mean_ions'][i],rtol=2e-7)
            assert np.allclose(np.diff(f['event_axes/'+str(ids[start])+'/time_edges_ms'][:]),f["settings/frame_ms"][()],atol=.03)
        has=e['has_pulse'][:].astype(bool)
        minimum_mask=e['minimum_histogram_accepted'][:].astype(bool)
        np.testing.assert_array_equal(minimum_mask,has & (e['minimum_ions'][:]>nmin))
        for name,values in (('inclusive',ions[keep]),('mean',e['mean_ions'][:][has]),('maximum',e['maximum_ions'][:][has]),('minimum',e['minimum_ions'][:][minimum_mask])):
            h=f['histograms/'+name];expected,_=np.histogram(values,h['edges_ions'][:])
            np.testing.assert_array_equal(expected,h['counts'][:])
            assert expected.sum()==len(values)
        for name in ('maximum','minimum','mean'):
            values=e[name+'_ions'][:][minimum_mask if name=='minimum' else has].astype(float)
            stats=f['distribution_statistics/'+name]
            assert stats['entries'][()]==len(values)
            if len(values):
                np.testing.assert_allclose(stats['mean_ions'][()],values.mean(),rtol=3e-7)
                np.testing.assert_allclose(stats['median_ions'][()],np.median(values),rtol=1e-7)
            else:assert np.isnan(stats['mean_ions'][()]) and np.isnan(stats['median_ions'][()])
        produced=e['produced_ions'][:]
        for name,values in (('all_pixels',e['collected_ions'][:]),):
            retention=f['ion_retention/'+name];ratio=values/produced
            np.testing.assert_allclose(retention['aggregate_retention'][()],values.sum()/produced.sum(),rtol=2e-7)
            np.testing.assert_allclose(retention['retention_per_event'][:],ratio,rtol=2e-7)
            np.testing.assert_allclose(retention['mean_event_retention'][()],ratio.mean(),rtol=2e-7)
            np.testing.assert_allclose(retention['median_event_retention'][()],np.median(ratio),rtol=2e-7)
        assert 'above_threshold' not in f['ion_retention']
        print(sample+': '+str(len(e['event_id']))+' events; uncut inclusive/mean/max, true event minimum selection, retention and all four histograms passed.')
