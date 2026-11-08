"""Projection conservation and signed response checks for the production viewer."""
import sys
from pathlib import Path
import unittest
import tempfile
import json
import h5py
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from cloud_viewer.production import ProductionSource
class Scalar:
    def __init__(self,x):self.x=x
    def __getitem__(self,_):return self.x
class ProjectionTest(unittest.TestCase):
    def test_signed_coordinates_and_charge(self):
        g={k:Scalar(v) for k,v in dict(pitch_mm=3,bucket_width_mm=3,drift_direction=1,reference_mm=np.array([10,20,30]),phase_mm=np.array([1,2]),event_id=7,ion_drift_cm=50,electron_drift_cm=100,pair_checksum=b'abc').items()}
        panel=ProductionSource.panel(g,[-1,0],[0,1],[10,11],np.array([-5.,12.]),30,'ion','ions',1,'Test')
        self.assertEqual(panel['edges']['x'],[8.,11.,14.])
        self.assertEqual(panel['edges']['z'],[60.,63.,66.])
        self.assertTrue(panel['signed'])
        for image in panel['projections'].values():self.assertEqual(np.sum(image),7)
        self.assertEqual(panel['cloud3d']['ions'],[-5.,12.])
        # The same physical cells arrive in reverse order at the other electrode.
        g['drift_direction']=Scalar(-1)
        other=ProductionSource.panel(g,[-1,0],[0,1],[-11,-12],np.array([-5.,12.]),3,'xe_he','PE',1,'Test')
        self.assertEqual(other['edges'],panel['edges'])
        self.assertEqual(other['projections'],panel['projections'])
        self.assertEqual(other['cloud3d']['z'],panel['cloud3d']['z'])
        self.assertEqual(other['cloud3d']['arrival_time'],[4.5,1.5])

class MatchedProductionTest(unittest.TestCase):
    """Read sparse, nonconsecutive MC IDs through a single-gas production."""
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.paths={k:self.root/f'{k}.h5' for k in ('ions','photons','sensor','asic')}
        for kind,path in self.paths.items():
            with h5py.File(path,'w') as f:
                f['metadata/complete']=True
                if kind=='ions':f['metadata/frame_s']=.04
                elif kind=='photons':f['metadata/photon_band']=np.bytes_('blue_tpb')
                else:
                    f['metadata/source_path']=np.bytes_(str(self.paths['photons']))
                    f['metadata/bucket_width_s']=.25e-6
                for eid in (7,101):
                    g=f.create_group(f'events/{eid}')
                    for k,v in dict(event_id=eid,ion_drift_cm=60.,electron_drift_cm=90.,
                                    pair_checksum=np.bytes_(f'pairs{eid}'),pitch_mm=4.,
                                    reference_mm=[0.,0.,0.],phase_mm=[.5,1.],
                                    bucket_width_mm=4. if kind=='ions' else .25,
                                    drift_direction=1 if kind=='ions' else -1).items():g[k]=v
                    if kind=='ions':
                        g['pixel_ix']=[0,1];g['pixel_iy']=[1,0];g['first_bucket']=3
                        g['counts']=[[12,0],[0,20]];g['recovered_ions']=[[-3.,0],[0,18.]]
                    else:
                        g['pixel_ix']=[0,1];g['pixel_iy']=[1,0];g['bucket']=[3,5]
                        if kind=='photons':g['count']=[100,200]
                        elif kind=='sensor':
                            g['detected_pe']=[32,64];g['sensor_pe_equiv']=[31.,62.]
                        else:
                            g['readout_pe_equiv']=[-2.,60.]
                            g['saturated']=[False,True];g['clip_loss_pe_equiv']=1.25
        self.manifest=self.root/'manifest.toml'
        self.manifest.write_text('complete = true\n[[shards]]\nsample = "bb0nu"\n'
            'gas = "xe_he15"\nevent_ids = [7, 101]\n'
            f'sensor = {json.dumps(str(self.paths["sensor"]))}\n'
            f'asic = {json.dumps(str(self.paths["asic"]))}\n')
        self.card=self.root/'viewer.toml'
        self.card.write_text('sample = "bb0nu"\nlabel = "Matched ELM"\n'
            'default_gas = "xe_he"\ndefault_pitch_mm = 4\nchamber_length_cm = 150\n'
            f'sipm_manifest = {json.dumps(str(self.manifest))}\n[[ions]]\npitch_mm = 4\n'
            f'path = {json.dumps(str(self.paths["ions"]))}\n')

    def test_retired_inputs_have_a_current_viewer_pointer(self):
        self.card.write_text('historical_sources_retired = true\n'+self.card.read_text())
        self.assertEqual(ProductionSource(self.card).ids,[7,101])
        self.paths['ions'].unlink()
        with self.assertRaisesRegex(ValueError,'elm_matched_response_viewer.toml'):
            ProductionSource(self.card)

    def test_single_gas_defaults_and_file_position(self):
        source=ProductionSource(self.card)
        self.assertEqual(source.info()['gases'],['xe_he'])
        self.assertEqual(source.info()['pitches'],['4'])
        payload=source.event_payload(2)
        self.assertEqual(payload['panels'][0]['source_event_id'],101)
        self.assertEqual(payload['panels'][1]['source_event_id'],101)
        self.assertEqual(payload['event_count'],2)
        ions,electrons=payload['panels']
        self.assertEqual(ions['view_ranges'],electrons['view_ranges'])
        for axis in 'xyz':
            lo=min(p['edges'][axis][0] for p in (ions,electrons))
            hi=max(p['edges'][axis][-1] for p in (ions,electrons))
            view=ions['view_ranges'][axis]
            self.assertAlmostEqual(view[1]-view[0],1.25*(hi-lo))
        with self.assertRaisesRegex(ValueError,'unavailable'):source.event_payload(1,gas='pure_xe')

    def test_response_stages_clipping_and_conservation(self):
        source=ProductionSource(self.card)
        for stage,total in [('photons',300),('detected',96),('sensor',93),('readout',58)]:
            payload=source.event_payload(1,electron_stage=stage)
            ion,electron=payload['panels']
            self.assertEqual(ion['signal_count'],15)
            self.assertIn('40 ms',ion['title'])
            self.assertIn('4 mm / 250 ns',electron['title'])
            self.assertEqual(electron['clipped_buckets'],1)
            self.assertEqual(electron['clip_loss_pe'],1.25)
            self.assertEqual(electron['signal_count'],total)
            for projection in electron['projections'].values():self.assertEqual(np.sum(projection),total)
        self.assertEqual(source.event_payload(1,stage='input')['panels'][0]['signal_count'],32)

    def test_different_ion_placement_is_rejected(self):
        with h5py.File(self.paths['ions'],'r+') as f:f['events/7/ion_drift_cm'][...]=61.
        with self.assertRaisesRegex(ValueError,'Paired event mismatch: ion_drift_cm'):
            ProductionSource(self.card).event_payload(1)

    def test_electronics_placement_is_checked(self):
        with h5py.File(self.paths['asic'],'r+') as f:f['events/7/electron_drift_cm'][...]=89.
        with self.assertRaisesRegex(ValueError,'Electron history mismatch: electron_drift_cm'):
            ProductionSource(self.card).event_payload(1)

if __name__=='__main__':unittest.main()
