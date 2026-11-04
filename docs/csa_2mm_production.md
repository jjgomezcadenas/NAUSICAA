# Ion CSA production at 2 mm / 20 ms

The 8,000-event bb0nu study uses `config/specs/runs/ion_csa_0p5_enc_production_2mm.toml`, `config/specs/response/ion_csa/csa_0p5_enc_2mm.toml` and the established three-frame reconstruction card. The configured modelling-error tolerance is 0.5 ENC. The 0.63 V starting guard includes the highest known bb0nu pilot failure at 0.60035 V; 30 controls were selected per voltage stratum.

The full-sample 0.63 V guard passed with 270 remaining controls. The largest filtered error was 0.45995 ENC; after the three-frame inverse it was 0.44839 ENC. Solver refinement passed. Direct SPICE was used for 2,528 of 2,302,198 charged pixel histories. The filtered and corrected responses each contain 8,000 events, have `metadata/complete = true`, and passed 56,000 event response/reconstruction checks and 8,003 summary checks.

The filtered study is `results/ion_csa_0p5_enc_production_2mm/da95679c2bfc42f2`; the corrected study is `results/ion_csa_0p5_enc_corrected_production_2mm/7a75ad6039a234a0`. The filtered workflow took 905.39 s and the corrected stage took 394.94 s, a combined 1,300.33 s (21m40s). The corrected event-weighted mean S/N over positive pixel-frame cells is 1.5970. Its mean noiseless summed-charge bias relative to input ions is −0.01973%. The 50% local-compression diagnostic affects 0.05% of events, with its separate one-ENC minimum deficit condition. These quantities retain all frames on pixels with input charge, use the documented central-square-pad weighting approximation and stationary-noise model, and apply no S/N acceptance cut.

## Storage

The 2 mm plan has 2,302,198 charged pixel histories and 468.2 million pixel-frame cells, or 203.4 frames per pixel on average. The configured three-second recovery tail alone contributes 150 frames at 20 ms. The filtered file saves dense input, linear, noiseless, noise, recovered, voltage and compression arrays on these channels. The corrected file copies the event contents and replaces four count arrays with their three-frame inverse results. Random noise and recovered-count arrays compress poorly.

| Pitch | Charged pixels | Pixel-frame cells | Plan | Filtered response | Corrected response |
|---|---:|---:|---:|---:|---:|
| 2 mm | 2,302,198 | 468.2 million | 6.91 GB | 7.64 GB | 9.05 GB |
| 3 mm | 1,150,123 | 157.4 million | 2.57 GB | 2.83 GB | 3.32 GB |
| 4 mm | 715,442 | 74.1 million | 1.33 GB | 1.49 GB | 1.72 GB |

The corrected response is the downstream reconstruction product. The filtered file remains its auditable parent under the current schema and is required to rerun or change the inverse reconstruction. A slimmer format could retain the parent response in an archive and store only corrected arrays and provenance in the active result; that is a future storage redesign, not part of this production run.

The generated filtered and corrected HTML, LaTeX and PDF reports include figures, tables, parameters and validation results. `make test physics` passed; both PDFs compiled with only missing-author LaTeX warnings.
