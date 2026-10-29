"""Serve six saved diffusion clouds for a paired MC event index."""
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
        return dict(event=event, event_count=self.event_count, panels=panels)


def serve(directory, host="127.0.0.1", port=8767, open_browser=False):
    source = CloudSource(directory)
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
                    event = int(parse_qs(route.query).get("event", ["0"])[0])
                    payload = json.dumps(source.event_payload(event), separators=(",", ":")).encode()
                    self.send_data(payload, "application/json", compress=True)
                elif route.path == "/api/info":
                    self.send_data(json.dumps(dict(event_count=source.event_count)).encode(), "application/json")
                elif route.path == "/plotly.js":
                    self.send_data(bundle, "application/javascript", compress=True)
                elif route.path in ("/", "/index.html", "/app.js", "/style.css"):
                    name = "index.html" if route.path == "/" else route.path[1:]
                    content_type = {"index.html": "text/html", "app.js": "application/javascript",
                                    "style.css": "text/css"}[name]
                    self.send_data((WEB / name).read_bytes(), content_type)
                else:
                    self.send_data(b"Not found", "text/plain", 404)
            except (ValueError, KeyError) as error:
                self.send_data(json.dumps({"error": str(error)}).encode(), "application/json", 400)

    httpd = ThreadingHTTPServer((host, port), Handler)
    url = f"http://{host}:{httpd.server_port}/"
    print(f"Cloud viewer: {url} ({source.event_count} paired events)", flush=True)
    if open_browser:
        webbrowser.open(url)
    try:
        httpd.serve_forever()
    finally:
        httpd.server_close()
