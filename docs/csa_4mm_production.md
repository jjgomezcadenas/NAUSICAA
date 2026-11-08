# Ion CSA production at 4 mm / 40 ms

Historical production. Its large event-response files and most acquired histories
were retired on 2026-11-08. Reports, summaries, cards and selected reference cases,
including the event-6625 guard-boundary case, remain. See
[the retirement record](reviews/historical_retirement_2026-11-08.md).
The current matched production is documented in [elm_matched_production.md](elm_matched_production.md).

The 8,000-event bb0nu study uses `config/specs/runs/ion_csa_0p5_enc_production_4mm_guard063.toml`, `config/specs/response/ion_csa/csa_0p5_enc_4mm.toml` and the established three-frame reconstruction card. The input is the saved 4 mm ion pad history, one positive Xe₂⁺ charge per counted ion. The circuit, weighting-field approximation, stationary-noise model and charge-integral checks are specified in the generated reports and the response card.

## Guard selection

The first 4 mm run expanded from 0.35 to 0.61 V and passed its sparse controls. The earlier 4 mm bb0nu pilot contains an independent direct-SPICE case at event 6625, pixel (−42, 221), with predicted minimum 0.61713055 V and filtered maximum error 0.50529 ENC. The same pad exists in the full plan, above the 0.61 V guard, and had no direct waveform in that first run. Its response writer was stopped before a completed response file was accepted. Its raw waveforms were retained for reuse.

The production card starts at 0.63 V and uses 30 controls per voltage stratum. This selects 29,339 direct-SPICE pixel histories among 715,442 charged pixel histories. The event-6625 pad is present with `model_used = 1` in the completed corrected response. The 270 remaining filtered controls have a maximum error of 0.39523 ENC. After the three-frame inverse, the maximum is 0.39136 ENC. Both are below the configured 0.5 ENC tolerance. Solver-refinement differences are below 0.01 ENC. Ten individual default-step waveforms were reacquired at finer steps to satisfy the current-integral check.

## Completed outputs

The filtered study is `results/ion_csa_0p5_enc_production_4mm/0e84b70bcc70c20c`; the corrected study is `results/ion_csa_0p5_enc_corrected_production_4mm/588da52478312c29`. Both HDF5 response files have `metadata/complete = true`, 8,000 events and 715,442 pixel histories. Each stage passed 56,000 response/reconstruction checks and 8,003 summary checks. The filtered study took 726.08 s and the corrected stage took 90.39 s: 816.47 s (13m36s) for the final validated workflow with cached waveform reuse. The earlier interrupted 0.61 V attempt took 1,351.02 s and produced reusable acquisitions.

The corrected event-weighted mean S/N over positive pixel-frame cells is 6.1825. Its mean noiseless summed-charge bias relative to input ions is −0.6212%. The 50% local-compression diagnostic affects 9.0875% of events, with its separate one-ENC minimum deficit condition. These quantities are computed over the saved full frame histories and should be interpreted with the central-square-pad weighting approximation, calibrated stationary noise and no S/N acceptance cut.

The generated filtered and corrected HTML, LaTeX and PDF reports include the figures, tables, parameters, controls and provenance. `make test physics` passed; both PDF reports compiled, with only missing-author LaTeX warnings.
