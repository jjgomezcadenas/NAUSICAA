# CSA production runtime investigation (2026-10-30)

The first 8,000-event bb0nu, 3 mm CSA acquisition was stopped after its cost became clear. Its plan and completed waveforms remain under `results/ion_csa_bb0nu_3mm_production/`; no full CSA response has been produced.

## Measured cause

The plan contains 1,150,123 pixel histories. The initial 0.35 V direct-SPICE guard selected 437 pixels plus 24 controls (461 histories). Several controls exceeded the specified 0.1 ENC maximum waveform-error tolerance, including one at 0.65 V. `assess` promotes the guard to the upper edge of that control's coarse voltage stratum, 0.8 V. This selects 123,343 histories, about 268 times the initial acquisition.

On the running machine, four workers acquired about 13.8 histories/s during an active stretch (6,750 histories by 660 s, including 176 s before that stretch). At that rate, the 0.8 V acquisition alone would take roughly 2.5 h of active time, before final reduction and reporting. The prepared `plan.h5` is 2.57 GB; 123,343 prepared current waveforms occupy 13.05 GB; the 7,413 saved raw histories average 416 kB each. Extrapolating raw storage gives approximately 51 GB at 0.8 V. These are disk estimates, not measured final response sizes.

The observed wall time was longer because the computer entered clamshell sleep at 20:49 JST and passed through further sleep intervals until 21:36 JST (`pmset -g log`). LTspice itself reports about 0.1–0.15 s per typical transient. The large history count is the dominant active-computation cost; sleep explains much of the additional clock time.

## Bounded error check

Julia reduced 210 previously acquired waveforms, 30 sampled from each voltage interval between 0.35 and 0.8 V. A second, uniformly sampled set of 180 waveforms was acquired across the complete 8,000-event plan, 60 each in 0.725–0.75, 0.75–0.775 and 0.775–0.8 V. Both checks used the saved `maximum_error_ENC` definition: the maximum per-frame difference between direct SPICE and the linear estimate, divided by the calibrated ENC.

| Minimum linear output (V) | Cases | Above 0.1 ENC | Maximum error (ENC) |
|---|---:|---:|---:|
| 0.700–0.725 | 30 | 22 | 0.1922 |
| 0.725–0.750 | 90 | 18 | 0.1515 |
| 0.750–0.775 | 90 | 0 | 0.0955 |
| 0.775–0.800 | 90 | 0 | 0.0700 |

The data support testing a guard near 0.75 V. Zero failures in 180 sampled cases above 0.75 V does not establish that every one of the remaining tens of thousands of cases meets 0.1 ENC. A 0.75 V guard would select 61,378 histories below threshold plus six higher-voltage controls, and about 25 GB of raw output, so further cost reduction is useful.

## Next gate

Before resuming the full CSA run, validate a finer, adaptive guard with controls distributed over the complete event sample, then benchmark the chosen guard and worker count on a bounded batch. Preserve the existing raw waveforms and source-current files so validated cases can be reused. Choose the error tolerance with the reconstructed observables in view.

## Pilot tolerance comparison

The completed 3 mm pilot contains 100 bb0nu and 100 single-electron events and 2,827 direct-SPICE pad histories. For each history, the saved direct-minus-linear response was passed through the same three-frame inverse used by the corrected production reconstruction. A counterfactual guard retains the saved direct response below the boundary and the linear response above it; controls remain direct. The guard below is the highest observed failing voltage across both the filtered and corrected per-bucket errors, followed by the next representable voltage. Full-plan counts include the controls.

| Error tolerance | Empirical joint guard | Direct histories in 8,000-event plan | bb0nu event charge difference: mean absolute / maximum | bb0nu event-to-event SD |
|---|---:|---:|---:|---:|
| 0.1 ENC | 0.760573 V | 70,350 | 0.0373% / 0.0909% | 0.0178% |
| 0.25 ENC | 0.694702 V | 31,197 | 0.1017% / 0.1787% | 0.0311% |
| 0.5 ENC | 0.627936 V | 13,997 | 0.1641% / 0.3294% | 0.0453% |

For the 0.25 ENC counterfactual, the maximum corrected-count difference in positive input buckets containing at most 10 ions was 0.216 ion in bb0nu, and the maximum integrated difference on one pad was 34.37 ions. The sum of absolute omitted cell differences was at most 0.1793% of an event's initial ion count, an upper bound on the charge difference for any fixed subset of pad-frame cells. These are differences from the saved direct-SPICE reference, not errors against physical truth. The pilot has no defined reconstructed blob regions, so a blob-charge comparison is still required. The event charge result includes all frames and does not include electronics noise in the difference; the same saved noise would cancel in this comparison.

The pilot has 1,536 bb0nu and 1,291 single-electron direct-SPICE histories; its 50% local-compression diagnostic found no affected events in either sample. This confirms that the direct-SPICE guard is much stricter than the compression diagnostic.

The initial audit found a literal 0.1 ENC assertion in corrected-response validation and literal 0.1 ENC presentation text. Both were changed for the bounded run described below. The empirical boundaries above are not guarantees for the full 8,000-event sample. They motivate a distributed validation before changing the production guard.

## 0.5 ENC bounded run

A new card sets the per-bucket limit to 0.5 ENC and uses finer voltage-control intervals around 0.6–0.65 V. One 100-event bb0nu, 3 mm run completed the filtered response, summary, corrected inverse, HTML and two compiled PDF reports. The initial filtered pass took 87.82 s on four acquisition workers, including preparation and final reduction; corrected reconstruction and presentation took another 13.85 s. The run acquired 180 direct histories and 12 solver-refinement histories. Its corrected control maximum was 0.4689 ENC and the solver-refinement maximum was 0.0000291 ENC.

That 100-event control set accepted a 0.60 V guard, whereas the 200-event pilot had observed a 0.5 ENC failure near 0.628 V. The small pilot therefore does not establish the production guard. A further 360 histories were sampled across the complete 8,000-event bb0nu plan, 60 each in 0.60–0.61, 0.61–0.62, 0.62–0.63, 0.63–0.64, 0.64–0.65 and 0.65–0.675 V. Julia applied the same three-frame inverse to each saved direct-minus-linear waveform. The highest failure in this sample was at 0.620012 V for the filtered error and 0.612941 V for the corrected error; none of the 180 sampled histories at or above 0.63 V exceeded 0.5 ENC. Together with the earlier pilot, this supports a **working 0.63 V boundary** for the production estimate, subject to production control validation.

The 8,000-event plan contains 14,334 pad histories below 0.63 V, before adding controls. Of these, 1,437 already have a raw waveform in the interrupted 0.1 ENC study and can be reused after matching source-current, model and simulator hashes. The interrupted run sustained about 13.8 acquisitions/s during an active interval, implying 15–18 minutes for the remaining or full direct-SPICE acquisition. Scaling the measured event-processing work and allowing for preparation, reduction, corrected reconstruction and reporting gives a provisional **35–60 minutes of active runtime** for the 8,000-event 0.5 ENC CSA chain on the same machine with four workers. Laptop sleep adds directly to elapsed wall time. Expect roughly 6 GB of raw waveforms, 1.5 GB of prepared PWL currents and a 2.6 GB plan, plus response files. The estimate is less certain for the event-level finish and summary stages because the benchmark has 100 events.

The corrected-response assertion now reads the tolerance saved in the parent response card; the HTML and LaTeX presentation display that configured value. Focused Julia tests, `make test physics`, both LaTeX builds and their logs passed. The only LaTeX warnings were missing-author warnings.

## Full 8,000-event production outcome

The authorized bb0nu 3 mm production run used the 0.5 ENC card and four LTspice workers. It selected 14,388 direct histories after the full-sample controls expanded the guard to 0.63 V; 12,857 new transients were acquired in the expanded iteration, with matching saved histories reused where available. Five individual histories required a finer solver step to pass the current-integral check. The acquisition stage completed 14,388 requests in 584 s, including cache reuse. The final filtered control error was 0.42424 ENC; the corrected three-frame control error was 0.41787 ENC, both below the card limit. Independent solver refinement passed at 0.01 ENC.

The filtered response contains 8,000 events and 1,150,123 pixel histories. Its 56,000 event/charge-ledger checks passed. The corrected response contains the same events and pads and passed its own 56,000 reconstruction checks. Both HDF5 response files report `metadata/complete = true`; the summaries, HTML and LaTeX were generated from those files. Both LaTeX reports compiled to PDF. The filtered command took 1,077.79 s and corrected reconstruction/presentation took 161.21 s: **1,239.00 s (20m39s) total active wall time**. This measured result supersedes the provisional 35–60 minute estimate above. `make test physics` passed after the production changes.

The corrected 3 mm summary gives 3.6106 event-weighted mean S/N over positive pixel-frame cells and a mean noiseless summed-charge bias of −0.1775% relative to input ions. These quantities use the saved full frame histories, the calibrated stationary noise model and the documented central-pad weighting approximation. The 50% local-compression diagnostic affected 1.15% of bb0nu events; it remains distinct from the direct-SPICE selection threshold.
