# Historical-run retirement — 2026-11-08

The approved second cleanup removed **80.40 GiB**
in **244,812 generated files**. The repository now occupies
**111.0 GiB**, compared with approximately 233 GiB before both cleanups.

| Retired numerical outputs | Allocated GiB removed | Files |
|---|---:|---:|
| `results/ion_csa_0p5_enc_production_4mm` | 22.28 | 163,024 |
| `results/ion_csa_bb0nu_3mm_production` | 3.07 | 16,760 |
| `results/ion_csa_v3_pilot` | 4.67 | 34,020 |
| `results/ion_csa_v3_pilot_4mm` | 3.28 | 24,098 |
| `results/ion_csa_v3_pilot_2mm` | 1.17 | 6,148 |
| `c25084b9…/response/` | 45.93 | 762 |

The earlier bb0nu realization is
`results/simulations/production_pilot/bb0nu/c25084b9e721e407c0688bfdba837f8f7fe5e453d8dad73034681d2d994ad95c/`.
Its detector inputs and primary ledger remain; the deletion scope was its
`response/` subtree. Current matched productions, background inputs and earlier
pure-Xe detector inputs remain available.

## Preserved scientific record

Cards, source snapshots, decks/models, summaries, calibration files and generated
reports remain in their original locations. Most retired SPICE acquisitions,
per-case numerical responses, PWL inputs, raw simulator files, logs and obsolete
plans were removed. Their SHA-256 identities and allocated sizes were recorded
before each deletion. Response metadata and preparation-card datasets were saved
in bounded HDF5 metadata archives. Tracked scientific source artifacts were
excluded from the deletion list.

**48 complete reference case folders** remain, comprising **328 checksum-verified
files**. For pilot plans still available, the retained cases span low, median and
high predicted minimum output voltage within each sample. For historical studies
whose preparation plans had already been retired, cases were selected by sorted
case-key quantiles. These historical reference cases complement the complete,
unchanged current nonlinear calibration evidence.

The retained guard-boundary case is event **6625**, pixel **(-42, 221)**, including
its raw waveform, current input, circuit, model and numerical response. Its saved
predicted minimum is **0.61713055 V**. This is the case that
justified the **0.63 V** guard. It remains under
`results/ion_csa_v3_pilot_4mm/289ba0ab6997a94f/spice_cases/bb0nu_p4_6625_-42_221`;
the corresponding other pilot/production cases were also retained.

The full list and hashes are in
[reference_histories.json](../../results/storage_cleanup/2026-11-08_historical/reference_histories.json).
Each retired root has a `RETIRED.md` marker. The broad candidate pool was about
82 GiB; retained evidence and source records account for the smaller deletion total.

## Current production checks

- **37,840 protected files** retain their original sizes and modification times.
  These include current calibration inputs, reduced event tables and the retained
  detector inputs from the earlier realization.
- All **9,391 validation histories and 240 fitting waveforms** remain present.
  The current nonlinear calibration HDF5 has its original SHA-256 hash.
- All **328 files** in the 48 retained historical case folders match their saved
  SHA-256 hashes. Their raw time/current/voltage arrays passed dimension, ordering
  and finite-value checks.
- The matched viewer produced paired ion/electron event payloads at **2, 3 and
  4 mm**, for positions **1, 3791 and 8000**: nine checks passed.
- All three live services responded to `/api/info` and event requests: ports
  **8768**, **8769** and **8774**. The full services loaded event 8000 and the
  smoke service loaded event 1. No restart was needed during this retirement.
- **75 Python tests**, `make physics`, and **7 Julia overview checks** passed.
  `git diff --check` passed.

## Historical overview and handoffs

`config/specs/runs/response_overview.toml` is now explicitly historical. Its
reducer checks source availability before opening the saved HDF5 reduction for
replacement. Missing retired inputs produce instructions for the saved
presentation and current viewer. A regression test verifies that a saved
reduction survives this error unchanged. The viewer likewise gives a current
configuration pointer when the historical card references retired inputs.

The old overview's numerical HDF5 summary, figures, HTML and PDF remain.
**`make response-overview-present` regenerated the HTML and compiled PDF** from
the preserved summary. They are labelled as the historical 10–150 cm production.
Use **`make response-viewer`** for the current matched 5–145 cm configuration.
The older CSA/SiPM handoffs and 4 mm production note now link to this retirement
record and the current matched production.

## Exception and accounting

One old unfinished **283 MiB** `ion_response.h5.tmp` had an unreadable HDF5 object
header. The first pass stopped during its metadata inspection. The file was
confirmed to be an approved historical temporary; its original bytes were hashed,
the error recorded, and retirement resumed using the existing deletion ledger.
All other successfully readable response metadata were preserved. There are no
remaining temporary files from this cleanup.

Filesystem free space was **477.6 GiB** at the final check.
The observed free-space increase also includes concurrent activity and filesystem
accounting; the **80.40 GiB** removal figure is
based on the enumerated deleted files, not that global change.

The local full deletion inventory is
`results/storage_cleanup/2026-11-08_historical/deleted_sha256.jsonl.gz`.
Metadata archives and protected-file inventories remain local beside it.
[Completion totals](../../results/storage_cleanup/2026-11-08_historical/completed.json),
[verification results](../../results/storage_cleanup/2026-11-08_historical/verification.json)
and [disk measurements](../../results/storage_cleanup/2026-11-08_historical/disk_measurements.json)
are retained as compact records. Inspection, execution and verification scripts
and logs are under `work/historical_retirement/`.
