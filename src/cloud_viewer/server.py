"""Serve one sample's ion and selected-gas PE histories from direct HDF5 readout."""
from functools import lru_cache
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlsplit
import gzip
import json
import math
import webbrowser
import h5py
import numpy as np
from plotly.offline import get_plotlyjs

ROOT = Path(__file__).resolve().parents[2]
WEB = Path(__file__).parent / "web"
SAMPLES = ("bb0nu", "electrons")
CARRIERS = ("ion", "pure_xe", "xe_he")


class CloudSource:
    """Read one event from each of six completed HDF5 streams."""

    def __init__(self, directory):
        self.directory = Path(directory).resolve()
        self.paths = {(sample, carrier): self.directory / f"{sample}_{carrier}.h5"
                      for sample in SAMPLES for carrier in CARRIERS}
        counts = []
        for path in self.paths.values():
            with h5py.File(path) as file:
                if not bool(file["metadata/complete"][()]):
                    raise ValueError(f"Incomplete cloud file: {path}")
                counts.append(int(file["metadata/event_count"][()]))
        if len(set(counts)) != 1:
            raise ValueError("The six cloud files have different event counts")
        self.event_count = counts[0]

    @staticmethod
    def project(ix, iy, iz, counts, shape):
        """Sum occupied cell counts along the hidden axis for display."""
        nx, ny, nz = shape
        projections = {}
        for plane, row, col, dims in (
            ("xy", iy, ix, (ny, nx)),
            ("xz", iz, ix, (nz, nx)),
            ("yz", iz, iy, (nz, ny)),
        ):
            image = np.zeros(dims, dtype=np.int64)
            np.add.at(image, (row, col), counts)
            projections[plane] = image.tolist()
        return projections

    @staticmethod
    def coarse_cloud(ix, iy, iz, counts, edges, factor=4):
        """Aggregate 3D cells for a responsive plot, preserving all counts."""
        bins = [len(edges[axis]) - 1 for axis in "xyz"]
        shape = tuple(math.ceil(n / factor) for n in bins)
        coords = (ix // factor, iy // factor, iz // factor)
        flat = np.ravel_multi_index(coords, shape)
        aggregate = np.bincount(flat, weights=counts, minlength=np.prod(shape)).reshape(shape)
        occupied = np.nonzero(aggregate)
        values = np.rint(aggregate[occupied]).astype(np.int64)
        centres = {}
        for axis, index, n in zip("xyz", occupied, bins):
            lo = edges[axis][index * factor]
            hi = edges[axis][np.minimum((index + 1) * factor, n)]
            centres[axis] = ((lo.astype(np.float64) + hi) / 2).tolist()
        assert int(values.sum()) == int(counts.sum())
        return {**centres, "ions": values.tolist(), "factor": factor}

    def read_panel(self, sample, carrier, event):
        with h5py.File(self.paths[sample, carrier]) as file:
            key = f"events/{event}"
            if key not in file:
                raise KeyError(f"Event {event} is absent from {sample}/{carrier}")
            g = file[key]
            edges = {axis: g[f"{axis}_edges_mm"][()] for axis in "xyz"}
            coords = {axis: g[axis][()] for axis in "xyz"}
            counts = g["ions"][()]
            indices = [np.searchsorted(edges[axis], coords[axis], side="right") - 1
                       for axis in "xyz"]
            shape = tuple(len(edges[axis]) - 1 for axis in "xyz")
            if any(np.any((i < 0) | (i >= n)) for i, n in zip(indices, shape)):
                raise ValueError(f"Cloud cell outside its saved edges: {sample}/{carrier}/{event}")
            projections = self.project(*indices, counts, shape)
            cloud3d = self.coarse_cloud(*indices, counts, edges)
            retained = int(g["retained_ions"][()])
            assert retained == int(counts.sum())
            for plane in ("xy", "xz", "yz"):
                assert sum(map(sum, projections[plane])) == retained
            return dict(sample=sample, carrier=carrier, event_id=event,
                        ion_drift_cm=float(g["ion_drift_cm"][()]),
                        electron_drift_cm=float(g["electron_drift_cm"][()]),
                        sampled_ions=int(g["sampled_ions"][()]),
                        retained_ions=retained, outside_ions=int(g["outside_ions"][()]),
                        occupied_cells=len(counts),
                        signal_count=retained,
                        signal_unit="ions" if carrier == "ion" else "electrons",
                        axis_units={"x": "mm", "y": "mm", "z": "mm"},
                        axis_names={"x": "x", "y": "y", "z": "z"},
                        edges={a: value.tolist() for a, value in edges.items()},
                        projections=projections, cloud3d=cloud3d)

    @lru_cache(maxsize=3)
    def event_payload(self, event):
        if not 0 <= event < self.event_count:
            raise ValueError(f"Event number must be between 0 and {self.event_count - 1}")
        panels = [self.read_panel(sample, carrier, event)
                  for sample in SAMPLES for carrier in CARRIERS]
        for start in (0, 3):
            positions = {(p["ion_drift_cm"], p["electron_drift_cm"])
                         for p in panels[start:start + 3]}
            if len(positions) != 1:
                raise ValueError("Paired cloud distances disagree")
        return dict(event=event, event_count=self.event_count, view="cloud", panels=panels)


class GridSource:
    """Read paired ion pad and selected-gas SiPM histories by file position."""

    def __init__(self, directory):
        self.directory = Path(directory).resolve()
        self.paths = {(sample, carrier): self.directory / f"{sample}_{carrier}.h5"
                      for sample in SAMPLES for carrier in CARRIERS}
        self.event_counts = {}
        for (sample, carrier), path in self.paths.items():
            with h5py.File(path) as file:
                schema = file.attrs.get("schema")
                if isinstance(schema, bytes):
                    schema = schema.decode()
                if schema != "nausicaa.direct_grid_readout.v2":
                    raise ValueError(f"Expected direct readout schema v2: {path}")
                if not bool(file["metadata/complete"][()]):
                    raise ValueError(f"Incomplete grid readout: {path}")
                count = int(file["metadata/event_count"][()])
                if sample in self.event_counts and self.event_counts[sample] != count:
                    raise ValueError(f"Grid readout event counts disagree for {sample}")
                self.event_counts[sample] = count

    def read_panel(self, sample, carrier, event_index):
        with h5py.File(self.paths[sample, carrier]) as file:
            key = f"events/{event_index}"
            if key not in file:
                raise KeyError(f"Grid readout event {event_index} is absent from {sample}/{carrier}")
            g = file[key]
            ion = carrier == "ion"
            time_unit = "ms" if ion else "us"
            edges = {"x": g["x_edges_mm"][()], "y": g["y_edges_mm"][()],
                     "z": g[f"time_edges_{time_unit}"][()]}
            ix, iy, iz = (g[field][()].astype(np.int64) - 1
                          for field in ("pad_ix", "pad_iy", "slice_index"))
            counts = g["signal_count"][()]
            shape = tuple(len(edges[axis]) - 1 for axis in "xyz")
            if any(np.any((index < 0) | (index >= size))
                   for index, size in zip((ix, iy, iz), shape)):
                raise ValueError(f"Grid readout index outside saved edges: {key}")
            projections = CloudSource.project(ix, iy, iz, counts, shape)
            display3d = {axis: ((edges[axis][index].astype(np.float64)
                                  + edges[axis][index + 1]) / 2).tolist()
                         for axis, index in zip("xyz", (ix, iy, iz))}
            display3d.update(ions=counts.tolist(), factor=1)
            total = int(g["total_signal"][()])
            if total != int(counts.sum()):
                raise ValueError(f"Signal ledger mismatch: {key}")
            return dict(sample=sample, carrier=carrier, event_number=event_index + 1,
                        ion_drift_cm=float(g["ion_drift_cm"][()]),
                        electron_drift_cm=float(g["electron_drift_cm"][()]),
                        signal_count=total, signal_unit="ions" if ion else "PE",
                        input_carriers=int(g["input_carriers"][()]),
                        occupied_cells=len(counts),
                        axis_units={"x": "mm", "y": "mm", "z": time_unit},
                        axis_names={"x": "x", "y": "y", "z": "t"},
                        edges={axis: value.tolist() for axis, value in edges.items()},
                        projections=projections, cloud3d=display3d)

    @lru_cache(maxsize=12)
    def event_payload(self, event_number, sample="bb0nu", gas="pure_xe"):
        if sample not in SAMPLES:
            raise ValueError(f"Unknown sample: {sample}")
        if gas not in ("pure_xe", "xe_he"):
            raise ValueError(f"Unknown electron gas: {gas}")
        event_count = self.event_counts[sample]
        if not 1 <= event_number <= event_count:
            raise ValueError(f"Event position must be between 1 and {event_count}")
        event_index = event_number - 1
        panels = [self.read_panel(sample, carrier, event_index)
                  for carrier in ("ion", gas)]
        positions = {(p["ion_drift_cm"], p["electron_drift_cm"]) for p in panels}
        if len(positions) != 1:
            raise ValueError("Paired grid readout distances disagree")
        return dict(event=event_number, event_count=event_count, sample=sample,
                    gas=gas, view="readout", panels=panels)


def serve(readout_directory, host="127.0.0.1", port=8767, open_browser=False):
    readout = GridSource(readout_directory)
    bundle = get_plotlyjs().encode()

    class Handler(BaseHTTPRequestHandler):
        def send_data(self, data, content_type, code=200, compress=False):
            payload = gzip.compress(data, compresslevel=4) if compress else data
            self.send_response(code)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(payload)))
            self.send_header("Cache-Control", "no-store")
            if compress:
                self.send_header("Content-Encoding", "gzip")
            self.end_headers()
            self.wfile.write(payload)

        def do_GET(self):
            route = urlsplit(self.path)
            try:
                if route.path == "/api/event":
                    query = parse_qs(route.query)
                    event = int(query.get("event", ["1"])[0])
                    sample = query.get("sample", ["bb0nu"])[0]
                    gas = query.get("gas", ["pure_xe"])[0]
                    payload = readout.event_payload(event, sample, gas)
                    payload = json.dumps(payload, separators=(",", ":")).encode()
                    self.send_data(payload, "application/json", compress=True)
                elif route.path == "/api/info":
                    self.send_data(json.dumps(dict(event_counts=readout.event_counts,
                                                   samples=SAMPLES, gases=("pure_xe", "xe_he"))).encode(),
                                   "application/json")
                elif route.path == "/plotly.js":
                    self.send_data(bundle, "application/javascript", compress=True)
                elif route.path in ("/", "/index.html", "/app.js", "/style.css"):
                    name = "index.html" if route.path == "/" else route.path[1:]
                    content_type = {"index.html": "text/html", "app.js": "application/javascript",
                                    "style.css": "text/css"}[name]
                    self.send_data((WEB / name).read_bytes(), content_type)
                else:
                    self.send_data(b"Not found", "text/plain", 404)
            except (ValueError, KeyError, OSError) as error:
                self.send_data(json.dumps({"error": str(error)}).encode(), "application/json", 400)

    httpd = ThreadingHTTPServer((host, port), Handler)
    url = f"http://{host}:{httpd.server_port}/"
    counts = ", ".join(f"{sample}: {count}" for sample, count in readout.event_counts.items())
    print(f"Ion/PE viewer: {url} ({counts})", flush=True)
    if open_browser:
        webbrowser.open(url)
    try:
        httpd.serve_forever()
    finally:
        httpd.server_close()
