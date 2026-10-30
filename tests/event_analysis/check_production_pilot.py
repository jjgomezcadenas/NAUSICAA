#!/usr/bin/env python3
"""Validate the four saved pilot files and shared pair identity per sample."""
from __future__ import annotations

import argparse
from pathlib import Path

import h5py


def value(group: h5py.Group, name: str):
    item = group[name][()]
    item = item.item() if hasattr(item, "item") else item
    return item.decode() if isinstance(item, bytes) else item


def check(summary: Path) -> None:
    with h5py.File(summary, "r") as results:
        assert results.attrs["schema"] == b"nausicaa.detector_pilot_summary.v1"
        assert value(results["metadata"], "complete")
        for sample in ("bb0nu", "single_electrons"):
            if sample not in results:
                continue
            group = results[sample]
            names = ("realisation", "ions", "pure_xe", "xe_he15")
            handles = {name: h5py.File(value(group, name + "_path"), "r") for name in names}
            try:
                assert all(value(handle["metadata"], "complete") for handle in handles.values())
                event_sets = [set(handle["events"].keys()) for handle in handles.values()]
                assert all(ids == event_sets[0] for ids in event_sets)
                assert len(event_sets[0]) == value(group, "events")
                for event_id in event_sets[0]:
                    rows = {name: handle["events"][event_id] for name, handle in handles.items()}
                    checksums = [value(row, "pair_checksum") for row in rows.values()]
                    assert len(set(checksums)) == 1
                    initial = value(rows["realisation"], "initial_pairs")
                    assert initial == value(rows["ions"], "initial_pairs")
                    assert initial == sum(rows["ions"]["count"][()])
                    for gas in ("pure_xe", "xe_he15"):
                        row = rows[gas]
                        assert value(row, "initial_pairs") == initial
                        assert value(row, "surviving_electrons") <= initial
                        assert sum(row["count"][()]) <= value(row, "emitted_photons")
                        keys = list(zip(row["pixel_ix"][()], row["pixel_iy"][()], row["bucket"][()]))
                        assert keys == sorted(set(keys))
            finally:
                for handle in handles.values():
                    handle.close()


def compare_common_event(first: Path, second: Path, sample: str, event_id: str) -> None:
    """Check that changing only the event limit leaves one event's streams intact."""
    with h5py.File(first, "r") as a, h5py.File(second, "r") as b:
        assert value(a[sample], "stream_identity") == value(b[sample], "stream_identity")
        for branch in ("realisation", "ions", "pure_xe", "xe_he15"):
            with h5py.File(value(a[sample], branch + "_path"), "r") as left, \
                 h5py.File(value(b[sample], branch + "_path"), "r") as right:
                x, y = left["events"][event_id], right["events"][event_id]
                assert value(x, "pair_checksum") == value(y, "pair_checksum")
                if branch != "realisation":
                    for name in ("pixel_ix", "pixel_iy", "bucket", "count"):
                        assert x[name][()].tolist() == y[name][()].tolist()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("--compare-summary", type=Path)
    parser.add_argument("--sample", default="single_electrons")
    parser.add_argument("--event-id", default="0")
    args = parser.parse_args()
    check(args.summary)
    if args.compare_summary:
        compare_common_event(args.summary, args.compare_summary, args.sample, args.event_id)
    print("Production pilot HDF5 identities and ledgers: passed")
