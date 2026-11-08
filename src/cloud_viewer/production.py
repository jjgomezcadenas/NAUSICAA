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
        if self.card.get('historical_sources_retired') and any(not p.is_file() for p in self.ions.values()):
            raise ValueError('Historical event responses were retired. Use config/specs/runs/elm_matched_response_viewer.toml for the current matched production; saved historical plots remain under results/response_overview.')
        self.label = self.card.get('label', 'Detector response')
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
            if gas not in ('pure_xe', 'xe_he'):
                raise ValueError(f'Unknown gas in production: {gas}')
            for event in shard['event_ids']:
                self.shards[gas, event] = shard
        self.gases = sorted({gas for gas, _ in self.shards})
        if not self.ions or not self.gases:
            raise ValueError('Production needs completed ion and electron branches')
        for gas in self.gases:
            id_sets.append({event for g,event in self.shards if g == gas})
        if any(ids != id_sets[0] for ids in id_sets):
            raise ValueError('Production branches have different event IDs')
        self.ids = sorted(id_sets[0])
        self.event_counts = {self.card['sample']: len(self.ids)}
        self.default_gas = self.card.get('default_gas', self.gases[0])
        self.default_pitch = str(self.card.get('default_pitch_mm',
                                              '3' if '3' in self.ions else next(iter(self.ions))))
        if self.default_gas not in self.gases or self.default_pitch not in self.ions:
            raise ValueError('Viewer defaults must select available production branches')

    def info(self):
        """Describe available branches and defaults for this completed production."""
        return dict(label=self.label, event_counts=self.event_counts,
                    samples=list(self.event_counts), gases=self.gases,
                    pitches=list(self.ions), default_gas=self.default_gas,
                    default_pitch=self.default_pitch)

    @staticmethod
    def panel(g, ix, iy, bucket, values, time_width, carrier, unit, ordinal, title):
        """Project signed cell amplitudes, preserving charge in each projection.

        Lattice indices are zero-based floors relative to reference + phase.
        Longitudinal cell edges are mapped back to source z in mm using the
        saved reference, bucket width and drift direction. Electron rows are
        reversed into increasing z. Relative arrival times remain attached to
        the 3D points for hover information.
        """
        pitch = float(g['pitch_mm'][()])
        base = np.asarray(g['reference_mm'][()], float)[:2] + g['phase_mm'][()]
        indices = [np.asarray(v, dtype=np.int64) for v in (ix,iy,bucket)]
        low = [int(v.min()) for v in indices]; high = [int(v.max()) for v in indices]
        edges = {'x': base[0]+np.arange(low[0], high[0]+2)*pitch,
                 'y': base[1]+np.arange(low[1], high[1]+2)*pitch,
                 'z': np.arange(high[2]-low[2]+2)*time_width}
        shifted = [v-lo for v,lo in zip(indices,low)]
        arrival_time = (shifted[2]+0.5)*time_width
        direction = int(g['drift_direction'][()])
        if direction not in (-1, 1):
            raise ValueError('Drift direction must be +1 or -1')
        width = float(g['bucket_width_mm'][()])
        edges['z'] = (float(g['reference_mm'][()][2]) + direction *
                      np.arange(low[2], high[2]+2)*width)
        if direction == -1:
            edges['z'] = edges['z'][::-1]
            shifted[2] = high[2]-indices[2]
        projections = {}
        for plane,a,b in [('xy',0,1),('xz',0,2),('yz',1,2)]:
            image = np.zeros((high[b]-low[b]+1, high[a]-low[a]+1))
            np.add.at(image,(shifted[b],shifted[a]),values)
            projections[plane] = image.tolist()
        points = {axis: ((edges[axis][v]+edges[axis][v+1])/2).tolist()
                  for axis,v in zip('xyz',shifted)}
        points.update(ions=np.asarray(values,float).tolist(),factor=1,
                      arrival_time=arrival_time.tolist())
        return dict(carrier=carrier,title=title,event_number=ordinal,
                    source_event_id=int(g['event_id'][()]),
                    ion_drift_cm=float(g['ion_drift_cm'][()]),
                    electron_drift_cm=float(g['electron_drift_cm'][()]),
                    pair_checksum=g['pair_checksum'][()].decode(),
                    signal_count=float(np.sum(values)),signal_unit=unit,
                    signed=bool(np.any(values<0)),occupied_cells=len(values),
                    axis_units=dict(x='mm',y='mm',z='mm'),
                    axis_names=dict(x='x',y='y',z='z'),
                    arrival_time_unit='ms' if carrier=='ion' else 'us',
                    edges={k:v.tolist() for k,v in edges.items()},
                    projections=projections,cloud3d=points)

    @lru_cache(maxsize=6)
    def event_payload(self, event_number, sample='bb0nu', gas=None, pitch=None, stage='response', electron_stage='readout'):
        """Return paired ion/electron panels by one-based position in the file."""
        if sample not in self.event_counts or not 1 <= event_number <= len(self.ids):
            raise ValueError('Sample or event position is unavailable')
        gas = self.default_gas if gas is None else gas
        pitch = self.default_pitch if pitch is None else str(pitch)
        if gas not in self.gases or pitch not in self.ions:
            raise ValueError('Gas or ion pitch is unavailable in this production')
        if stage not in ('input','response') or electron_stage not in ('photons','detected','sensor','readout'):
            raise ValueError('Unknown response stage')
        event = self.ids[event_number-1];key=f'events/{event}'
        with h5py.File(self.ions[pitch]) as f:
            g=f[key];counts=g['counts'][()];pad,frame=np.nonzero(counts>0)
            values=counts[pad,frame] if stage=='input' else g['recovered_ions'][()][pad,frame]
            frame_ms=float(f['metadata/frame_s'][()])*1000
            ions=self.panel(g,g['pixel_ix'][()][pad],g['pixel_iy'][()][pad],
                            frame+int(g['first_bucket'][()]),values,frame_ms,
                            'ion','ions',event_number,f'Ions · {pitch} mm / {frame_ms:g} ms · {stage}')
        shard=self.shards[gas,event]
        with h5py.File(resolve(shard['sensor'])) as sensor, h5py.File(resolve(shard['asic'])) as asic:
            with h5py.File(sensor['metadata/source_path'][()].decode()) as photons:
                g=photons[key];s=sensor[key];a=asic[key]
                for paired in (s,a):
                    for field in ('pair_checksum','ion_drift_cm','electron_drift_cm'):
                        if paired[field][()] != g[field][()]:
                            raise ValueError(f'Electron history mismatch: {field}')
                px,py,pb=(g[k][()] for k in ('pixel_ix','pixel_iy','bucket'))
                if electron_stage=='photons': values=g['count'][()]
                else:
                    group=a if electron_stage=='readout' else s
                    keys=zip(*(group[k][()] for k in ('pixel_ix','pixel_iy','bucket')))
                    lookup={k:i for i,k in enumerate(keys)}
                    positions=[lookup[k] for k in zip(px,py,pb)]
                    field={'detected':'detected_pe','sensor':'sensor_pe_equiv','readout':'readout_pe_equiv'}[electron_stage]
                    values=group[field][()][positions]
                bucket_us=float(sensor['metadata/bucket_width_s'][()])*1e6
                electrons=self.panel(g,px,py,pb,values,
                    bucket_us,gas,
                    'photons' if electron_stage=='photons' else 'PE',event_number,
                    f'Electrons · {"pure Xe" if gas=="pure_xe" else "Xe + 15% He"} · '
                    f'{float(g["pitch_mm"][()]):g} mm / {bucket_us*1000:g} ns · {electron_stage}')
                electrons['clipped_buckets'] = int(np.count_nonzero(a['saturated'][()]))
                electrons['clip_loss_pe'] = float(a['clip_loss_pe_equiv'][()])
                electrons['photon_band'] = (photons['metadata/photon_band'][()].decode()
                                             if 'metadata/photon_band' in photons else 'vuv')
        for k in ('pair_checksum','ion_drift_cm','electron_drift_cm'):
            if ions[k]!=electrons[k]: raise ValueError(f'Paired event mismatch: {k}')
        if 'chamber_length_cm' in self.card and not np.isclose(
                ions['ion_drift_cm']+ions['electron_drift_cm'],
                self.card['chamber_length_cm'], rtol=0, atol=2e-5):
            raise ValueError('Paired drift distances do not sum to the chamber length')
        # Shared initial framing: 25% larger than the combined occupied cells.
        ranges = {}
        for axis in 'xyz':
            lo = min(p['edges'][axis][0] for p in (ions,electrons))
            hi = max(p['edges'][axis][-1] for p in (ions,electrons))
            margin = (hi-lo)*0.125
            ranges[axis] = [lo-margin,hi+margin]
        for panel in (ions,electrons):
            panel['view_ranges'] = ranges
        return dict(event=event_number,event_count=len(self.ids),sample=sample,gas=gas,
                    label=self.label,panels=[ions,electrons],view='production')
