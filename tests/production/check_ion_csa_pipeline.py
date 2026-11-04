"""Independently verify published ion matrices, provenance and source-cell replay."""
from pathlib import Path
import argparse
import hashlib
import json
import tomllib
import h5py
import numpy as np


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def string(dataset):
    value = dataset[()]
    return value.decode() if isinstance(value, bytes) else str(value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    args = parser.parse_args()
    manifest = tomllib.loads(args.manifest.read_text())
    assert manifest['complete']
    assert manifest['schema'] == 'nausicaa.ion_response_run.v1'
    checks = 2
    events = 0
    paired = {}
    for node in manifest['responses']:
        path, parent = Path(node['path']), Path(node['parent_corrected_response'])
        assert sha(path) == node['sha256']
        checks += 1
        with h5py.File(path) as out, h5py.File(parent) as previous:
            meta = out['metadata']
            assert out.attrs['schema'] in (b'nausicaa.ion_csa_response.v1', 'nausicaa.ion_csa_response.v1')
            assert bool(meta['complete'][()])
            assert string(meta['channel_selection']) == 'pixels_with_positive_input_charge'
            assert string(meta['absent_channel_policy']) == 'no_input_charge_noise_not_simulated'
            assert not bool(meta['measured_zero_suppression'][()])
            assert not bool(meta['negative_estimate_clipping'][()])
            assert sha(parent) == string(meta['parent_response_sha256'])
            assert string(meta['numerical_payload']) == 'unchanged_from_corrected_parent'
            checks += 8
            # Independent exact comparison of every numerical payload dataset.
            for group_name in ('events', 'calibration', 'covariance_by_frames'):
                datasets = []
                out[group_name].visititems(lambda key, value: datasets.append(key) if isinstance(value, h5py.Dataset) else None)
                old = []
                previous[group_name].visititems(lambda key, value: old.append(key) if isinstance(value, h5py.Dataset) else None)
                assert datasets == old
                checks += 1
                for key in datasets:
                    assert np.array_equal(out[group_name][key][()], previous[group_name][key][()]), key
                    checks += 1
            source_path = string(meta['source_path'])
            assert sha(source_path) == string(meta['source_sha256'])
            ids = meta['selected_event_ids'][:].tolist()
            assert ids == node['selected_event_ids']
            assert sorted(map(int, out['events'])) == ids
            checks += 3
            with h5py.File(source_path) as source:
                for event_id in ids:
                    events += 1
                    event = out[f'events/{event_id}']
                    upstream = source[f'events/{event_id}']
                    count = event['counts'][:]  # HDF5/Python: pad, frame
                    reconstructed = event['recovered_ions'][:]
                    x, y = event['pixel_ix'][:], event['pixel_iy'][:]
                    assert count.shape == reconstructed.shape
                    assert np.all(count >= 0) and np.all(count.sum(axis=1) > 0)
                    assert count.sum() == upstream['initial_pairs'][()]
                    assert string(event['pair_checksum']) == string(upstream['pair_checksum'])
                    positions = dict(zip(zip(x.tolist(), y.tolist()), range(len(x))))
                    first = int(event['first_bucket'][()])
                    cells = list(zip(upstream['pixel_ix'][:], upstream['pixel_iy'][:], upstream['bucket'][:], upstream['count'][:]))
                    assert set(positions) == {(int(ix), int(iy)) for ix, iy, _, _ in cells}
                    assert np.count_nonzero(count) == len(cells)
                    for ix, iy, bucket, ions in cells:
                        assert count[positions[(int(ix), int(iy))], int(bucket)-first] == ions
                    sigma = event['noise_sigma_ions'][:]
                    cov = out[f'covariance_by_frames/{count.shape[1]}/covariance_e2'][:]
                    assert np.allclose(np.diag(cov), sigma**2, rtol=1e-10)
                    assert np.max(np.abs(reconstructed.astype(float)-event['noiseless_estimated_ions'][:]-event['noise_ions'][:])/sigma) < .001
                    checks += 8 + len(cells)
                    identity = (int(upstream['initial_pairs'][()]), string(event['pair_checksum']))
                    key = (node['sample'], event_id)
                    assert paired.setdefault(key, identity) == identity
                    checks += 1
    result = {'passed': True, 'responses': len(manifest['responses']), 'events': events, 'checks': checks,
              'payload': 'all event, calibration and covariance datasets exactly equal to corrected parent'}
    target = args.manifest.parent / 'independent_validation.json'
    target.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
