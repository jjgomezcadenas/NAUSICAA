"""Check selected ion/PE panels and one-based file positions against HDF5."""
from pathlib import Path
import sys
import h5py

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from cloud_viewer.server import GridSource


def main():
    source = GridSource(ROOT / "results/grid_readout_direct_8000")
    assert source.event_counts == {"bb0nu": 8000, "electrons": 8000}
    for sample, gas, event in (
        ("bb0nu", "pure_xe", 1),
        ("bb0nu", "xe_he", 100),
        ("electrons", "pure_xe", 3791),
        ("electrons", "xe_he", 8000),
    ):
        payload = source.event_payload(event, sample, gas)
        assert payload["view"] == "readout"
        assert payload["sample"] == sample and payload["gas"] == gas
        assert payload["event"] == event and payload["event_count"] == 8000
        assert [panel["carrier"] for panel in payload["panels"]] == ["ion", gas]
        for panel in payload["panels"]:
            assert panel["sample"] == sample and panel["event_number"] == event
            total = panel["signal_count"]
            assert sum(panel["cloud3d"]["ions"]) == total
            assert all(sum(map(sum, panel["projections"][plane])) == total
                       for plane in ("xy", "xz", "yz"))
            assert panel["axis_names"]["z"] == "t"
            assert panel["signal_unit"] == ("ions" if panel["carrier"] == "ion" else "PE")
            assert panel["axis_units"]["z"] == ("ms" if panel["carrier"] == "ion" else "us")
            if panel["carrier"] == "ion":
                assert total == panel["input_carriers"]
            with h5py.File(source.paths[panel["sample"], panel["carrier"]]) as f:
                g = f[f"events/{event - 1}"]
                assert total == int(g["signal_count"][()].sum())
        print(f"{sample}/{gas} event {event}: two panels and all readout projections verified")
    assert source.event_payload(100, "bb0nu", "pure_xe")["panels"][0] == source.read_panel("bb0nu", "ion", 99)
    for args in ((0, "bb0nu", "pure_xe"), (8001, "electrons", "xe_he"),
                 (1, "unknown", "pure_xe"), (1, "bb0nu", "unknown")):
        try:
            source.event_payload(*args)
        except ValueError:
            pass
        else:
            raise AssertionError(f"Invalid selection accepted: {args}")


if __name__ == "__main__":
    main()
