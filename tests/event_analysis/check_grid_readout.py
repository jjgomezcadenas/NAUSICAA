"""Validate paired ion-pad and electron-SiPM time histories."""
from pathlib import Path
import argparse
import h5py
import numpy as np

STREAMS = ((sample, carrier) for sample in ("bb0nu", "electrons")
           for carrier in ("ion", "pure_xe", "xe_he"))


def text(value):
    return value.decode() if isinstance(value, bytes) else value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--events", type=int, required=True)
    args = parser.parse_args()
    for sample in ("bb0nu", "electrons"):
        paths = [args.directory / f"{sample}_{carrier}.h5"
                 for carrier in ("ion", "pure_xe", "xe_he")]
        files = [h5py.File(path) for path in paths]
        try:
            for file, carrier in zip(files, ("ion", "pure_xe", "xe_he")):
                assert text(file.attrs["schema"]) == "nausicaa.direct_grid_readout.v2"
                assert file["metadata/complete"][()]
                assert file["metadata/event_count"][()] == args.events
                assert len(file["events"]) == args.events
                assert text(file["metadata/signal_unit"][()]) == ("ions" if carrier == "ion" else "pe")
            total_electrons = [0, 0]
            total_pe = [0, 0]
            for event in range(args.events):
                groups = [file[f"events/{event}"] for file in files]
                assert len({float(g["ion_drift_cm"][()]) for g in groups}) == 1
                assert len({float(g["electron_drift_cm"][()]) for g in groups}) == 1
                for axis in "xyz":
                    assert len({float(g[f"reference_{axis}_mm"][()]) for g in groups}) == 1
                for axis in "xy":
                    assert groups[1][f"grid_offset_{axis}_mm"][()] == groups[2][f"grid_offset_{axis}_mm"][()]
                for index, group in enumerate(groups):
                    values = group["signal_count"][()]
                    assert int(group["sampled_carriers"][()]) == (
                        int(group["input_carriers"][()]) +
                        int(group["outside_carriers"][()]))
                    assert values.dtype == np.int64 and np.all(values > 0)
                    assert int(values.sum()) == int(group["total_signal"][()])
                    assert all(len(group[field]) == len(values)
                               for field in ("pad_ix", "pad_iy", "slice_index"))
                    assert "photons" not in group and "electron_count" not in group
                    edges = group["time_edges_ms" if index == 0 else "time_edges_us"][()]
                    assert len(edges) == len(group["z_edges_mm"])
                    assert np.all(np.diff(edges) > 0)
                    width = 20 if index == 0 else 2
                    assert abs(float(np.median(np.diff(edges))) - width) < 0.003
                    if index == 0:
                        assert group["total_signal"][()] == group["input_carriers"][()]
                    else:
                        total_electrons[index - 1] += int(group["input_carriers"][()])
                        total_pe[index - 1] += int(group["total_signal"][()])
            for index, file in enumerate(files[1:]):
                expectation = (total_electrons[index] *
                               float(file["metadata/el_photons_per_electron"][()]) *
                               float(file["metadata/pe_probability_per_photon"][()]))
                assert abs(total_pe[index] / expectation - 1) < 0.003
            print(f"{sample}: {args.events} events; ion conservation, paired grid phase, "
                  "time buckets and PE yield verified")
        finally:
            for file in files:
                file.close()


if __name__ == "__main__":
    main()
