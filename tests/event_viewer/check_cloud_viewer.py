"""Check six-panel event payloads against the saved diffusion clouds."""
from pathlib import Path
import sys
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from cloud_viewer.server import CloudSource


def main():
    source = CloudSource(ROOT / "results/three_clouds_8000")
    assert source.event_count == 8000
    for event in (0, 2, 7999):
        payload = source.event_payload(event)
        assert len(payload["panels"]) == 6
        for panel in payload["panels"]:
            retained = panel["retained_ions"]
            assert sum(panel["cloud3d"]["ions"]) == retained
            assert all(sum(map(sum, panel["projections"][plane])) == retained
                       for plane in ("xy", "xz", "yz"))
            assert panel["sampled_ions"] == retained + panel["outside_ions"]
            assert abs(panel["ion_drift_cm"] + panel["electron_drift_cm"] - 150) < 1e-5
            with h5py.File(source.paths[panel["sample"], panel["carrier"]]) as file:
                g = file[f"events/{event}"]
                x, y, z = (g[axis][()] for axis in "xyz")
                counts = g["ions"][()]
                xy, _, _ = np.histogram2d(y, x, bins=(g["y_edges_mm"][()],
                                                         g["x_edges_mm"][()]), weights=counts)
                assert np.array_equal(xy, panel["projections"]["xy"])
        print(f"event {event}: six panels, complementary drift and all display sums verified")
    try:
        source.event_payload(8000)
    except ValueError:
        pass
    else:
        raise AssertionError("Out-of-range event accepted")


if __name__ == "__main__":
    main()
