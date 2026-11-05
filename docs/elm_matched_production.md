# Matched ELM production

The matched production uses one common realization for the ion and electron arms of 8,000 bb0nu events. The ion reference drift is uniform between 5 and 145 cm; the electron reference drift is its complement in the 150 cm chamber. Both arms retain the same source event IDs and primary-pair checksums.

The primary-reference card anchors the primary sampling stream independently of the new placement. `realisation_identities` in `scripts/production/specs.jl` verifies the reference seed, source, ionization work and Fano factor. The new realization identity includes the new placement specification. `run_common_realisation.jl` verifies every primary count and checksum against the saved reference.

## Physical stages

1. `run_matched_detector.jl` creates the common ledger, collects ions into 4 mm pixels and 40 ms frames, and collects Xe–15% He photons into 4 mm EL channels and 250 ns buckets.
2. The ion CSA pipeline and statistical SiPM/ASIC pipeline run concurrently. Each uses four workers for its independent histories or event shards. Their input files and response manifests remain separate.
3. `assemble_matched_response.jl` checks event IDs, primary checksums and both drift distances across the common ledger, ion collection, photon collection, corrected CSA, SiPM charge and ASIC response. It publishes the combined manifest and viewer card after these checks pass.
4. Julia reduces the saved histories. Python produces the plots, HTML and self-contained LaTeX/PDF from those reductions.

The electron optical geometry is a 2.5 mm hole/guide diameter, 5 mm EL depth and a 0.5 mm radial xenon gap at E/P = 2 kV/cm/bar and 15 bar. The optical model includes analytical hole acceptance, TPB conversion 0.65, forward emission 0.5 and the agreed light-guide capture. Saved photons arrive at the SiPM; PDE is applied in the response stage. The S14160-3015PS model uses 21,808 illuminated cells within 39,984 physical cells and a provisional 12.5 ns recovery constant. The full dark rate remains 700 kcps. See the generated report and resolved cards for the complete assumptions.

The CSA uses the existing circuit calibration with the 0.5 ENC fast-model tolerance and a 0.63 V initial voltage guard. Selected nonlinear histories and validation controls receive transistor simulation. Charge quadrature and independent time-step refinement are checked before publication.

## Commands and files

```sh
make elm-matched-production PYTHON=work/event-viewer-venv/bin/python
make elm-matched-analysis PYTHON=work/event-viewer-venv/bin/python
make elm-matched-report PYTHON=work/event-viewer-venv/bin/python
```

The first command runs the complete workflow. The second regenerates numerical reductions and presentation from the saved response. The third regenerates presentation alone.

- Combined graph: `results/elm_matched_production/manifest.toml`.
- Pairing validation: `results/elm_matched_production/matching_validation.h5`.
- Numerical reductions: `matched_analysis.h5` and `electron_analysis.h5` in that directory.
- Viewer card: `config/specs/runs/elm_matched_response_viewer.toml`.
- HTML: `results/elm_matched_production/elm_matched_production.html`.
- TeX/PDF: `latex_generated/electron_readout/elm_matched_production/`.
- Full-run cards: `config/specs/runs/elm_matched_*full.toml`.
- Shared realization: `config/specs/realisations/elm_matched_full.toml`.

The `*_smoke.toml` cards run ten matched events. Their files are saved under `results/elm_matched_smoke` with a separate viewer card.

Pulse mean/max comparisons use identical positive-input channel–bucket cells at every response stage. Signed reconstructed charge is retained. Ion total charge includes all stored quiet/recovery frames. Electron energy widths correct known attachment survival and deposited energy, and report both the full ensemble and the subset with no clipped ASIC bucket.
