# Storage archive of 2026-11-17

Status at 13:50 CET: every folder below is archived and verified, and its
original removed. Mac: 566 files (production_pilot) and 209,587 files
(elm_matched_production) identical. Hyperion: `/scratch/jjgomez/NAUSICAA-runs/nexus`
is empty; for `prod_bi214_20M` the only difference was a `SHA256SUMS` list
present only in the `/data` copy.

Data moved off the working disks to `/data/jjgomez` on the DIPC Hyperion cluster
(`hyperion.sw.ehu.es`), on the user's instruction of 2026-11-17. Every file is
compared by SHA-256 on both sides before the original is removed; the checksum
lists are kept beside the archived folders (`<folder>.sha256`).

`/data` is read-only on the compute nodes: read archived files from the login
node, or copy them back to `/scratch` before a job uses them.

## From the Mac (`results/`, untracked in Git)

| Folder on the Mac | Archive on Hyperion | Size | Content |
|---|---|---|---|
| `results/simulations/production_pilot` | `/data/jjgomez/NAUSICAA-mac-archive/results/simulations/production_pilot` | 37 GB | cached detector-response realisations of the matched ELM productions (bb0nu 35 GB, single electrons 2.4 GB) and their pilot and scan summaries |
| `results/elm_matched_production` | `/data/jjgomez/NAUSICAA-mac-archive/results/elm_matched_production` | 16 GB | matched ELM production of 2026-11-08: analysis files, pages and `ion_csa/359613d487dfeee8`, the training study of the nonlinear CSA calibration v1 |

The nonlinear calibration itself, `results/ion_csa_nonlinear/v1/calibration.h5`,
stays on the Mac. Its card
`config/specs/response/ion_csa/nonlinear_calibration_v1.toml` names the training
study by its old path; a revalidation of the calibration needs the study copied
back to that path first.

To restore a folder on the Mac, from the repository root:

    rsync -a jjgomez@hyperion.sw.ehu.es:/data/jjgomez/NAUSICAA-mac-archive/results/<folder> results/<parent>/

Files that refer to these folders by path (cards, scripts and documents of the
historical productions; none is used by the ITACA or iItaca chains of
2026-11-17):

`results/simulations/production_pilot`:

- `config/specs/generated/19c062288425/run.toml`
- `config/specs/generated/4642ca9fdf8c/run.toml`
- `config/specs/generated/4725adea70b9/run.toml`
- `config/specs/generated/602b106e8f72/run.toml`
- `config/specs/generated/91f0e47f63fc/run.toml`
- `config/specs/generated/9d5f403be80c/run.toml`
- `config/specs/generated/ab59d0545f4b/run.toml`
- `config/specs/generated/df1414afcdc2/run.toml`
- `config/specs/runs/electron_250ns.toml`
- `config/specs/runs/electron_250ns_ep2.toml`
- `config/specs/runs/electron_250ns_ep2p5.toml`
- `config/specs/runs/electron_250ns_ep3.toml`
- `config/specs/runs/elm_matched_detector_full.toml`
- `config/specs/runs/elm_matched_detector_smoke.toml`
- `config/specs/runs/elm_matched_ions_pitch2.toml`
- `config/specs/runs/elm_matched_ions_pitch3.toml`
- `config/specs/runs/elm_matched_response_viewer.toml`
- `config/specs/runs/elm_matched_smoke_viewer.toml`
- `config/specs/runs/elm_photons_ep2.toml`
- `config/specs/runs/elm_photons_ep2_drift10cm.toml`
- `config/specs/runs/elm_photons_pitch4mm_ep2.toml`
- `config/specs/runs/elm_photons_pitch4mm_ep2_full.toml`
- `config/specs/runs/ion_csa_0p5_enc_benchmark.toml`
- `config/specs/runs/ion_csa_0p5_enc_production.toml`
- `config/specs/runs/ion_csa_0p5_enc_production_2mm.toml`
- `config/specs/runs/ion_csa_0p5_enc_production_4mm.toml`
- `config/specs/runs/ion_csa_0p5_enc_production_4mm_guard063.toml`
- `config/specs/runs/ion_csa_bb0nu_3mm_production.toml`
- `config/specs/runs/ion_csa_v3_integration_smoke.toml`
- `config/specs/runs/ion_csa_v3_pilot.toml`
- `config/specs/runs/ion_csa_v3_pilot_2mm.toml`
- `config/specs/runs/ion_csa_v3_pilot_4mm.toml`
- `config/specs/runs/ion_csa_v3_pilot_all.toml`
- `config/specs/runs/ion_csa_v3_production.toml`
- `config/specs/runs/production_pilot.toml`
- `config/specs/runs/response_overview.toml`
- `config/specs/runs/sipm_bb0nu_ep2p5.toml`
- `config/specs/runs/sipm_bb0nu_ep3.toml`
- `config/specs/runs/sipm_bb0nu_production.toml`
- `config/specs/runs/sipm_bb0nu_statistical_ep2_pilot100.toml`
- `config/specs/runs/sipm_bb0nu_statistical_ep2p5_pilot100.toml`
- `config/specs/runs/sipm_bb0nu_statistical_ep3_pilot100.toml`
- `config/specs/runs/sipm_bb0nu_statistical_pilot.toml`
- `config/specs/runs/sipm_response_ep2p5.toml`
- `config/specs/runs/sipm_response_ep3.toml`
- `config/specs/runs/sipm_response_pilot.toml`
- `config/specs/runs/sipm_response_statistical_ep2_pilot.toml`
- `config/specs/runs/sipm_response_statistical_ep2p5_pilot.toml`
- `config/specs/runs/sipm_response_statistical_ep3_pilot.toml`
- `config/specs/runs/sipm_response_statistical_pilot.toml`
- `docs/ion_csa_handoff.md`
- `docs/ion_csa_response.md`
- `docs/matched_ion_csa_pitch_scan.md`
- `docs/response_production_spec.html`
- `docs/response_production_spec.md`
- `docs/reviews/historical_retirement_2026-11-08.html`
- `docs/reviews/historical_retirement_2026-11-08.md`
- `docs/reviews/storage_audit_2026-11-08.html`
- `docs/reviews/storage_audit_2026-11-08.json`
- `docs/reviews/storage_audit_2026-11-08.md`
- `docs/sipm_electronics_handoff.md`
- `docs/sipm_production_handoff.md`
- `scripts/production/README.md`
- `scripts/production/run_electron_ep2_pilot.jl`
- `scripts/production/run_readout_pulses.py`
- `scripts/production/run_sipm_tpb_pilot.jl`

`results/elm_matched_production`:

- `config/specs/response/ion_csa/nonlinear_calibration_v1.toml`
- `config/specs/runs/elm_matched_csa_full.toml`
- `config/specs/runs/elm_matched_response_viewer.toml`
- `config/specs/runs/elm_matched_sipm_full.toml`
- `config/specs/runs/elm_matched_sipm_response_full.toml`
- `config/specs/runs/ion_csa_nonlinear_pilot.toml`
- `config/specs/runs/ion_csa_nonlinear_production_100.toml`
- `docs/elm_matched_production.md`
- `docs/generated_reports.json`
- `docs/reviews/independent_repository_review_handoff.md`
- `docs/reviews/repository_review_2026-11-05.html`
- `docs/reviews/repository_review_2026-11-05.md`
- `docs/reviews/storage_audit_2026-11-08.html`
- `docs/reviews/storage_audit_2026-11-08.md`
- `scripts/production/present_ion_csa_scaling.py`
- `scripts/production/present_matched_analysis.py`
- `scripts/production/reduce_ion_csa_nonlinear_report.jl`
- `scripts/production/reduce_matched_analysis.jl`
- `scripts/production/run_matched_ion_pitch.jl`

## On Hyperion, from `/scratch` to `/data`

Geant4 (NEXUS) files of the ITACA hybrid productions, Xe-15 % He, 16 bar, R = 1.8 m,
L = 1.3 m, from `/scratch/jjgomez/NAUSICAA-runs/nexus/itaca_XeHe15_16bar_290K_R180_L130/`
to `/data/jjgomez/NAUSICAA-runs/nexus/itaca_XeHe15_16bar_290K_R180_L130/`:

| Folder | Size | Note |
|---|---|---|
| `prod_bi214_20M` | 59 GB | already copied to `/data` earlier; the scratch copy is removed once identical |
| `prod_100k` | 30 GB | bb0nu and single electrons, 20 x 5,000 each |
| `prod_10x5000` | 6.4 GB | |
| `pilot_10x100` | 0.9 GB | |
| `pilot_bi214_10x10000` | 0.6 GB | |

A rerun of a hybrid readout batch needs `NAUSICAA_DATA` pointing to a copy of
these files on `/scratch`.

Geant4 files of the iItaca variant II_L1000_D500 (pilot and production, 212 GB),
`/scratch/jjgomez/NAUSICAA-runs/nexus/iitaca/II_L1000_D500/`, go to
`/data/jjgomez/NAUSICAA-runs/nexus/iitaca/II_L1000_D500/` once their readout is
complete. The hit files and feature tables stay in the working copy of the
`iItaca` branch, `/scratch/jjgomez/NAUSICAA-iitaca/results/iitaca/`.

Script: `/scratch/jjgomez/NAUSICAA-runs/bookkeeping/archive.sh SRC DEST` (rsync
without overwriting, SHA-256 of every file on both sides, removal of SRC only
when the lists are identical). Logs in the same folder.
