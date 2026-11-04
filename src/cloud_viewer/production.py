"""Read completed ion CSA and SiPM production on paired positive-input cells."""
from pathlib import Path
from functools import lru_cache
import tomllib
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
def resolve(path):
    """Resolve a repository-relative source path."""
    p = Path(path)
    return p if p.is_absolute() else ROOT / p

class ProductionSource:
    """Index one sample's ion pitches and electron response shards by event ID."""
    def __init__(self, card, operating_point="legacy"):
        self.card = tomllib.loads(resolve(card).read_text())
        if operating_point != "legacy":
            entry = next(x for x in self.card['operating_points'] if x['id']==operating_point)
            index = tomllib.loads(resolve(entry['index']).read_text())
            if not index['complete']: raise ValueError('Incomplete operating point')
            self.card['sipm_manifest'] = index['sipm_manifest']
        self.ions = {str(x['pitch_mm']): resolve(x['path']) for x in self.card['ions']}
        self.shards = {}
        manifest = tomllib.loads(resolve(self.card['sipm_manifest']).read_text())
        if not manifest['complete']:
            raise ValueError('Incomplete SiPM production')
        id_sets = []
        for path in self.ions.values():
            with h5py.File(path) as f:
                if not f['metadata/complete'][()]:
                    raise ValueError('Incomplete ion response')
                id_sets.append(set(map(int, f['events'])))
        for shard in manifest['shards']:
            if shard['sample'] != self.card['sample']:
                continue
            gas = 'xe_he' if shard['gas'] == 'xe_he15' else shard['gas']
            for event in shard['event_ids']:
                self.shards[gas, event] = shard
        for gas in ('pure_xe', 'xe_he'):
            id_sets.append({event for g,event in self.shards if g == gas})
        if any(ids != id_sets[0] for ids in id_sets):
            raise ValueError('Production branches have different event IDs')
        self.ids = sorted(id_sets[0])
        self.event_counts = {self.card['sample']: len(self.ids)}

    @staticmethod
    def panel(g, ix, iy, bucket, values, time_width, carrier, unit, ordinal, title):
        """Project signed cell amplitudes, preserving charge in each projection.

        Lattice indices are zero-based floors relative to reference + phase.
        Time is relative to the first occupied bucket, in ms for ions and us
        for electrons. Each view has independent axes.
        """
        pitch = float(g['pitch_mm'][()])
        base = np.asarray(g['reference_mm'][()], float)[:2] + g['phase_mm'][()]
        indices = [np.asarray(v, dtype=np.int64) for v in (ix,iy,bucket)]
        low = [int(v.min()) for v in indices]; high = [int(v.max()) for v in indices]
        edges = {'x': base[0]+np.arange(low[0], high[0]+2)*pitch,
                 'y': base[1]+np.arange(low[1], high[1]+2)*pitch,
                 'z': np.arange(high[2]-low[2]+2)*time_width}
        shifted = [v-lo for v,lo in zip(indices,low)]
        projections = {}
        for plane,a,b in [('xy',0,1),('xz',0,2),('yz',1,2)]:
            image = np.zeros((high[b]-low[b]+1, high[a]-low[a]+1))
            np.add.at(image,(shifted[b],shifted[a]),values)
            projections[plane] = image.tolist()
        points = {axis: ((edges[axis][v]+edges[axis][v+1])/2).tolist()
                  for axis,v in zip('xyz',shifted)}
        points.update(ions=np.asarray(values,float).tolist(),factor=1)
        return dict(carrier=carrier,title=title,event_number=ordinal,
                    source_event_id=int(g['event_id'][()]),
                    ion_drift_cm=float(g['ion_drift_cm'][()]),
                    electron_drift_cm=float(g['electron_drift_cm'][()]),
                    pair_checksum=g['pair_checksum'][()].decode(),
                    signal_count=float(np.sum(values)),signal_unit=unit,
                    signed=bool(np.any(values<0)),occupied_cells=len(values),
                    axis_units=dict(x='mm',y='mm',z='ms' if carrier=='ion' else 'us'),
                    axis_names=dict(x='x',y='y',z='relative arrival time'),
                    edges={k:v.tolist() for k,v in edges.items()},
                    projections=projections,cloud3d=points)

    @lru_cache(maxsize=6)
    def event_payload(self, event_number, sample='bb0nu', gas='pure_xe', pitch='3', stage='response', electron_stage='readout'):
        """Return paired ion/electron panels by one-based position in the file."""
        if sample not in self.event_counts or not 1 <= event_number <= len(self.ids):
            raise ValueError('Sample or event position is unavailable')
        if stage not in ('input','response') or electron_stage not in ('photons','detected','sensor','readout'):
            raise ValueError('Unknown response stage')
        event = self.ids[event_number-1];key=f'events/{event}'
        with h5py.File(self.ions[pitch]) as f:
            g=f[key];counts=g['counts'][()];pad,frame=np.nonzero(counts>0)
            values=counts[pad,frame] if stage=='input' else g['recovered_ions'][()][pad,frame]
            ions=self.panel(g,g['pixel_ix'][()][pad],g['pixel_iy'][()][pad],
                            frame+int(g['first_bucket'][()]),values,float(f['metadata/frame_s'][()])*1000,
                            'ion','ions',event_number,f'Ions · {pitch} mm / {int(float(pitch)*10)} ms · {stage}')
        shard=self.shards[gas,event]
        with h5py.File(resolve(shard['sensor'])) as sensor, h5py.File(resolve(shard['asic'])) as asic:
            with h5py.File(sensor['metadata/source_path'][()].decode()) as photons:
                g=photons[key];s=sensor[key];a=asic[key]
                for paired in (s,a):
                    if paired['pair_checksum'][()] != g['pair_checksum'][()]:
                        raise ValueError('Electron pair checksum mismatch')
                px,py,pb=(g[k][()] for k in ('pixel_ix','pixel_iy','bucket'))
                if electron_stage=='photons': values=g['count'][()]
                else:
                    group=a if electron_stage=='readout' else s
                    keys=zip(*(group[k][()] for k in ('pixel_ix','pixel_iy','bucket')))
                    lookup={k:i for i,k in enumerate(keys)}
                    positions=[lookup[k] for k in zip(px,py,pb)]
                    field={'detected':'detected_pe','sensor':'sensor_pe_equiv','readout':'readout_pe_equiv'}[electron_stage]
                    values=group[field][()][positions]
                electrons=self.panel(g,px,py,pb,values,
                    float(sensor['metadata/bucket_width_s'][()])*1e6,gas,
                    'photons' if electron_stage=='photons' else 'PE',event_number,
                    f'Electrons · {"pure Xe" if gas=="pure_xe" else "Xe + 15% He"} · 250 ns · {electron_stage}')
        for k in ('pair_checksum','ion_drift_cm','electron_drift_cm'):
            if ions[k]!=electrons[k]: raise ValueError(f'Paired event mismatch: {k}')
        return dict(event=event_number,event_count=len(self.ids),sample=sample,gas=gas,
                    panels=[ions,electrons],view='production')
