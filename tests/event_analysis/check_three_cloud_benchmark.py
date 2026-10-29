"""Check paired drift, count conservation and concurrency invariance."""
from pathlib import Path
from contextlib import ExitStack
import argparse
import hashlib
import h5py
import numpy as np


def digest_file(path, expected_events):
    digest = hashlib.sha256()
    with h5py.File(path) as f:
        assert f["metadata/complete"][()]
        assert f["metadata/event_count"][()] == expected_events
        carrier = path.stem.split("_", 1)[1]
        expected_helium = 0.0 if carrier == "pure_xe" else 0.15
        assert abs(float(f["metadata/helium_fraction"][()]) - expected_helium) < 1e-6
        events = f["events"]
        assert len(events) == expected_events
        for event_id in sorted(events, key=int):
            g = events[event_id]
            counts = g["ions"][()]
            assert counts.dtype == np.int64
            assert sum(map(int, counts)) == int(g["retained_ions"][()])
            assert int(g["sampled_ions"][()]) == (
                int(g["retained_ions"][()]) + int(g["outside_ions"][()]))
            for field in ("x", "y", "z", "x_edges_mm", "y_edges_mm", "z_edges_mm"):
                assert g[field].dtype == np.float32
            digest.update(event_id.encode())
            for field in ("x", "y", "z", "ions", "sampled_ions", "retained_ions",
                          "outside_ions", "ion_drift_cm", "electron_drift_cm",
                          "sigma_t_mm", "sigma_l_mm"):
                digest.update(np.asarray(g[field][()]).tobytes())
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--events", type=int, default=100)
    args = parser.parse_args()
    roots = ([args.directory] if (args.directory / "bb0nu_ion.h5").exists()
             else [args.directory / f"workers_{w}" for w in (1, 2, 4, 6)])
    for sample in ("bb0nu", "electrons"):
        baseline = None
        for root in roots:
            paths = [root / f"{sample}_{carrier}.h5"
                     for carrier in ("ion", "pure_xe", "xe_he")]
            current = tuple(digest_file(p, args.events) for p in paths)
            if baseline is None:
                baseline = current
            else:
                assert current == baseline, (sample, root, "clouds changed")
            with ExitStack() as stack:
                files = [stack.enter_context(h5py.File(p)) for p in paths]
                for event_id in files[0]["events"]:
                    groups = [f["events"][event_id] for f in files]
                    positions = [g["ion_drift_cm"][()] for g in groups]
                    electron_distances = [g["electron_drift_cm"][()] for g in groups]
                    assert positions[0] == positions[1] == positions[2]
                    assert electron_distances[0] == electron_distances[1] == electron_distances[2]
                    assert abs(float(positions[0]) + float(electron_distances[0]) - 150) < 1e-5
                    assert len({int(g["sampled_ions"][()]) for g in groups}) == 1
        print(f"{sample}: all {args.events} paired events conserve counts"
              + (" and match across worker counts" if len(roots) > 1 else ""))


if __name__ == "__main__":
    main()
