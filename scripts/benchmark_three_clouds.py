#!/usr/bin/env python3
"""Benchmark six Julia cloud streams at bounded process concurrency."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import argparse
import json
import os
import subprocess
import time
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
STREAMS = [(sample, carrier) for sample in ("bb0nu", "electrons")
           for carrier in ("ion", "pure_xe", "xe_he")]


def run_one(julia, sample, carrier, n_events, directory):
    output = directory / f"{sample}_{carrier}.h5"
    log = directory / f"{sample}_{carrier}.log"
    env = dict(os.environ, JULIA_NUM_THREADS="1")
    start = time.monotonic()
    with log.open("w") as handle:
        result = subprocess.run([julia, "--startup-file=no", f"--project={ROOT}",
                                 str(ROOT / "scripts/propagate_cloud_stream.jl"),
                                 sample, carrier, str(n_events), str(output)],
                                cwd=ROOT, env=env, stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"{sample}/{carrier} failed; see {log}")
    peak_rss_bytes = next(int(line.split("=", 1)[1]) for line in log.read_text().splitlines()
                          if line.startswith("peak_rss_bytes="))
    return {"sample": sample, "carrier": carrier,
            "seconds": time.monotonic() - start, "bytes": output.stat().st_size,
            "peak_rss_bytes": peak_rss_bytes}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--events", type=int, default=100)
    parser.add_argument("--workers", type=int, nargs="+", default=[1, 2, 4, 6])
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--output", default="results/three_cloud_benchmark")
    args = parser.parse_args()
    root = ROOT / args.output
    root.mkdir(parents=True, exist_ok=True)
    summary_path = root / "benchmark.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else []
    for workers in args.workers:
        directory = root / f"workers_{workers}"
        directory.mkdir(exist_ok=True)
        start = time.monotonic()
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(run_one, args.julia, s, c, args.events, directory)
                       for s, c in STREAMS]
            rows = [future.result() for future in as_completed(futures)]
        row = {"workers": workers, "events_per_stream": args.events,
               "wall_seconds": time.monotonic() - start,
               "total_bytes": sum(item["bytes"] for item in rows),
               "sum_worker_peak_rss_bytes": sum(item["peak_rss_bytes"] for item in rows),
               "streams": sorted(rows, key=lambda item: (item["sample"], item["carrier"]))}
        summary = sorted([item for item in summary if item["workers"] != workers] + [row],
                         key=lambda item: item["workers"])
        summary_path.write_text(json.dumps(summary, indent=2) + "\n")
        with h5py.File(root / "benchmark.h5", "w") as h5:
            h5.attrs["schema"] = "nausicaa.three_cloud_benchmark.v1"
            h5["workers"] = np.array([item["workers"] for item in summary], dtype="i8")
            h5["events_per_stream"] = np.array([item["events_per_stream"] for item in summary], dtype="i8")
            h5["wall_seconds"] = np.array([item["wall_seconds"] for item in summary], dtype="f8")
            h5["output_bytes"] = np.array([item["total_bytes"] for item in summary], dtype="i8")
            h5["sum_worker_peak_rss_bytes"] = np.array(
                [item.get("sum_worker_peak_rss_bytes", -1) for item in summary], dtype="i8")
        print(f"workers={workers}: {row['wall_seconds']:.2f} s, "
              f"{row['total_bytes']/1e6:.1f} MB", flush=True)


if __name__ == "__main__":
    main()
