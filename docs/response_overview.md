# Ion and electron response overview

The overview compares the completed 8,000-event bb0nu productions. Ion pitches
2, 3 and 4 mm use 20, 30 and 40 ms frames, respectively. Electron pure Xe and
Xe + 15% He branches use the fixed 10 mm SiPM pitch and 250 ns buckets.

For each event, compute the mean and maximum over cells with positive input
ions (ion branch) or incident photons (electron branch). Apply the identical
cell selection to the response and retain signed reconstructed amplitudes.
Each event has equal weight in the mean/maximum distributions. The table gives
the mean and median of each distribution. Pulse amplitudes have no S/N cut.
The complete acquired histories remain available for energy sums.

The ion response is the corrected CSA reconstruction with 0.5 ENC modelling
tolerance. Electron stages are incident photons, detected PE, SiPM charge in
PE equivalents, and electronics-reconstructed PE. Comparisons use the same
event IDs, primary-pair checksums and complementary drift distances.

## Commands

- `make response-overview`: Julia reduction followed by Python presentation and PDF compilation.
- `make response-overview-present`: regenerate presentation from saved HDF5.
- `make response-viewer`: serve the existing viewer in production mode at <http://127.0.0.1:8769/>.
- `make response-overview-test`: paired-mask and projection checks.

The input card is `config/specs/runs/response_overview.toml`. It points to the
three corrected ion responses and the complete SiPM shard manifest. Reduction
streams one event at a time and writes `results/response_overview/response_overview.h5`.
The HTML and PNGs are alongside it; the self-contained LaTeX/PDF report is in
`latex_generated/response_overview/`.

## Interactive display

Select ion pitch, ion input or corrected response, electron gas and response
stage. Navigate by one-based position in the file; both arms show the same
source event. Each arm has independent XY, X–time, Y–time and rotatable 3D
views. XY and time projections sum the selected cell amplitudes. Colour is
linear; signed responses use a diverging scale. Positions use the saved grid
reference and phase, and each arm's time starts at its first occupied bucket.

Changing ion pitch keeps the selected electron branch unchanged. Only bb0nu is
offered because its complete electronics productions are the sources in this
card. The original direct-readout mode remains available through `--readout`.

## Validation

The reducer checks event identity, pair checksum and both drift distances
across all branches. Photon keys are matched to SiPM and electronics cells.
Tests verify that quiet input cells are excluded even if they contain noise,
negative output cells remain in the statistics, signed projection sums are
conserved, and negative lattice indices produce the correct coordinates.

## Electron operating-point comparison (4 November 2026)

The comparison uses E/P = 2.5, J = 0.015 and E/P = 3, J = 0.025 for both gases.
These are the user's corrected approximate readings of Oliveira Figure 6.
The existing E/P = 3, J = 0.01 production remains available as “previous”.

`make electron-operating-points` replays photons, runs the SiPM/ASIC response,
reduces the histories in Julia and generates HTML/LaTeX/PDF. Complete response
shards are reused on restart. `make electron-operating-points-present`
regenerates only presentation. Production index TOMLs are in
`results/electron_operating_points/`; the report is
`results/electron_operating_points/electron_operating_points.html`.

The inclusive input-to-electronics distribution uses sensor charge in PE
equivalents on photon-positive buckets. Full event sums use every acquired
cell. Expected dark-only baseline charge is subtracted before the known
attachment correction. Its stationary recovery factor is
`1 / (1 + dark_rate_Hz * recovery_s / microcell_count)`. Illumination can
further suppress dark avalanche charge through recovery. The report retains
the raw sum, corrected sum and corrected response per deposited MeV.
It gives RMS/mean and Gaussian-equivalent FWHM, with the latter defined as
`2 * sqrt(2 * log(2)) * RMS / mean` over the full distribution.

Refresh the viewer to discover completed operating points in the EL selector.
“Produce PNG” downloads both panels and all eight current views, including zoom
and 3D camera settings, with event identity, operating point and drift labels.

## Statistical SiPM pilot

The compact statistical sensor response is selected by
`config/specs/response/sipm/sipm_statistical_v1.toml`. It uses the same saved
PE-per-bucket sensor and ASIC formats. The implementation, commands and
assumptions are described in [sipm_statistical_implementation.md](sipm_statistical_implementation.md).

`make sipm-statistical-validation` compares it with microscopic simulation and
generates [the HTML report](../results/sipm_statistical_validation/validation.html)
and [PDF](../latex_generated/sipm_statistical_validation/validation.pdf).
`make sipm-statistical-pilot` runs two bb0nu events for each gas through the
sensor and ASIC stages. Full production remains stopped; the existing viewer
continues to use its saved production sources.

## Small blue SiPM with TPB

The S14160-1315PS pilot compares 100 matched bb0nu events at E/P = 2.0 and
2.5 for both gases, with a centred 1.3 mm square sensor at zero gap.
The [model and source notes](sipm_tpb_pilot.md) define square-aperture photon
collection, TPB conversion, the 12.5/25 ns recovery comparison and fixed-charge
electronics noise. [HTML results](../results/sipm_tpb_pilot/sipm_tpb_pilot.html)
and [PDF](../latex_generated/sipm_tpb_pilot/sipm_tpb_pilot.pdf) separate optical
statistics, sensor recovery and electronics clipping.

`make sipm-tpb-pilot` runs fresh timing pilots and the complete report workflow.
The runner preserves completed photon and response histories. Fresh timing
rejects an existing response cache; use `make sipm-tpb-report` to regenerate
presentation from saved reductions. `scripts/production/reduce_sipm_tpb_pilot.jl`
refreshes the numerical comparisons independently of simulation.
