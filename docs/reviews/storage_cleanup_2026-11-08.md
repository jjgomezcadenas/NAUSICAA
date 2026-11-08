# Completed storage cleanup — 2026-11-08

The approved cleanup reduced the repository from **233.0 GiB to 191.3 GiB**.
Observed filesystem free space increased by **42.0 GiB**
to **385.0 GiB**. Filesystem free-space changes include
concurrent activity and APFS accounting; the itemized figures below are allocated-file measurements.

## Removed and preserved

| Action | Allocated space removed | Preservation |
|---|---:|---|
| Old CSA plans and PWL current inputs | 22.06 GiB | 177,666 file identities and SHA-256 hashes saved; 117,697 acquired HDF5 histories unchanged |
| Legacy binned cloud/readout files | 5.74 GiB | 18 large files retired; summaries, source metadata and hashes retained; 126 KiB event fixture extracted |
| Duplicate corrected ion responses | 12.77 GiB | Current 2/3/4 mm published responses retained; original metadata and hashes archived; verified canonical receipts installed |
| Obsolete review worktree | about 1.00 GiB | Both modified PDFs preserved and checksum-verified; review notes and logs retained |
| Old `work/julia_depot` cache | 0.22 GiB | Active `work/julia-depot` and viewer environment retained |

The old prepared-input selection was restricted to `plan.h5` and `.pwl` files
under `results/ion_csa_bb0nu_3mm_production` and
`results/ion_csa_0p5_enc_production_4mm`. Historical reruns must prepare those inputs
again. Raw acquisitions, model/deck records, jobs, cards, summaries and reports
remain available. The report presenter now reads its saved response card from
the response file when a prepared plan has been retired.

The retired binned caches were in `results/three_clouds_8000`,
`results/grid_readout_8000` and `results/three_cloud_sample_500`. Their summary
HDF5 files remain. The regression fixture is
[results/validation/direct_replay_fixture](../../results/validation/direct_replay_fixture/README.md).
It contains the exact original event arrays for single-electron events 0 and
3790. The replay test now reads that fixture.

## Response deduplication

All **601,989 non-metadata datasets** in the three
corrected parents were compared with their published copies, including shapes,
types, values and attributes. Every comparison passed. Original metadata were
saved separately. Published files and their hashes are unchanged.

Each removed `ion_response.h5` has an adjacent `ion_response.canonical.toml`
receipt. The receipt records the original hash, retained-file relative path and
retained-file hash. Julia readers and saved-study registration verify the receipt,
the published file checksum, its declared parent and its numerical-payload
contract. The Python presenter uses the same checks. Original summary files and
their recorded parent hashes remain unchanged. Publication reuses an already
validated production file, avoiding another numerical copy.

The original files were first renamed within their directories. Python and Julia
then read through the receipts, and Julia verified saved-study references and
production-node reuse. Only after those checks passed were the renamed duplicates
deleted. No `.retiring` files remain.

Keep the receipts alongside the retained output. Copies of all three receipts,
full paths and verification results are saved in
[the cleanup records](../../results/storage_cleanup/2026-11-08/verified_duplicates.json).

## Validation

- **74 Python tests passed**, including canonical-response hash mismatch and
  incomplete-output rejection; `make physics` passed.
- **44 Julia pipeline checks passed**, including original-hash preservation,
  missing-parent lookup, repeat registration and rejection of altered receipts.
- Legacy replay passed for events **0 and 3790**, including the extended event's
  direct assignment comparison.
- Matched viewer data loaded at **2, 3 and 4 mm**, for event positions
  **1, 3791 and 8000**, including the electron branch: nine checks passed.
- Current 2/3/4 mm corrected response presentations regenerated. The 2 and 3 mm
  matched production PDF reports compiled. The older completed 4 mm presentation
  also regenerated with its prepared plan absent.
- All **9,391 calibration validation waveform paths**
  remain present; all 117,697 acquisition files checked in the affected
  old-study directories retain their original sizes and modification times.
- All four original G4 files, reduced hit tables, current detector inputs,
  production outputs, primary ledgers, circuit calibration evidence, manuscripts,
  papers and reports remain.

The existing viewers on ports **8768, 8769 and 8774** initially closed HTTP
connections even for `/api/info`. Their output descriptors pointed to old session
pipes. They were restarted with the same inputs and ports, with stdout/stderr
redirected to `work/storage_cleanup/viewer_<port>.log`. All three then passed
`/api/info` and event-1 requests. The matched production remains at
[localhost:8769](http://127.0.0.1:8769/).

`make note` remains blocked by the previously reported build failure in the
older `notes/physics/nausica.tex` draft. Its log is in
`work/storage_cleanup/make_note.log`; author text was preserved. The generated
matched production PDFs build successfully.

## Worktree preservation and follow-up

The two modified review-worktree PDFs were saved under
[worktree_preserved](../../results/storage_cleanup/2026-11-08/worktree_preserved/manifest.json)
before removing that worktree with Git. Both differ from the main working-tree
PDFs, which were preserved. The obsolete cache had no open files at the cleanup
check; documented launcher examples now use the active hyphenated depot path.

The broader retirement of older electronics configurations, the additional
filtered-stage responses and `results/diffused_tracks` was **not performed**.
Those require separate scientific retention decisions or dependency checks.
Repository history was preserved; three previously tracked large pilot HDF5 files
now appear as working-tree deletions. A future commit will record their removal,
while their older Git versions remain in history. No commit or push was made.

Full local hash inventories, metadata snapshots, receipts, preservation checks
and disk measurements are under `results/storage_cleanup/2026-11-08/`.
Inspection and cleanup scripts and logs are in `work/storage_cleanup/`.
The [original audit](storage_audit_2026-11-08.html) records the approved scope.
