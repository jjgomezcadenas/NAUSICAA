# NAUSICA CSA v3 — nominal 20 s discharge

Frozen provisional design for Vicente Herrero and Rogelio Palomo.
The selected circuit is `length_1200um`. Its fitted small-signal discharge time
is 19.54 s; the design label is 20 s. The saved alpha factors use
the fitted value. The discharge comes from Mf in the transistor circuit.

## Open and run

Open an `.asc` in `schematic/` in LTspice and press Run. Keep the complete
v3 folder together: currents are in `inputs/`, with time in seconds and current
in amperes. Plot V(out), V(in) and I(Is1).

- `NAUSICA_20ms_120ions.asc`: illustrative 120-ion packet over 20 ms.
- `NAUSICA_40ms_400ions.asc`: illustrative 400-ion packet over 40 ms.
- `NAUSICA_20ms_calibration.asc` and `NAUSICA_40ms_calibration.asc`: 100 ions.
- `NAUSICA_20ms_baseline.asc` and `NAUSICA_40ms_baseline.asc`: zero current.
- `NAUSICA_20ms_discharge.asc`: 100 ions and a six-second discharge tail.
- `NAUSICA_noise.asc`: output noise from 0.0001 Hz to 1 MHz.

The exact model is `schematic/180nm_professional.inc`, also provided as
`netlists/corrected.inc`. All eight schematics are electrically checked against
the archived netlists. `netlists/` also includes both time resolutions, the
10/100-ion tail checks, all tested packet sizes and readable simulator logs.
The circuit runs used native macOS LTspice 17.2.4.

## Circuit and readout

Supply 1.8 V, temperature 20 C, external input capacitance 1.25 pF plus
intrinsic model capacitances, feedback capacitance 3 fF. Mf: W=0.36 um,
L=1200 um, gate=0 V. M4 gate=0.6621 V. `parameters.json` lists the actual
W, L and multiplicity of every transistor, every bias and the operating point.
The schematic and netlists preserve the exact values used in the saved runs.

| Quantity | 20 ms frame | 40 ms frame |
|---|---:|---:|
| Fitted discharge time (s) | 19.5441 | 19.5440 |
| alpha | 0.998977195 | 0.997955429 |
| Gain (uV / collected ion) | 40.157 | 41.479 |
| Filtered output noise (mV RMS) | 2.015 | 2.424 |
| ENC (electrons RMS) | 50.18 | 58.44 |

Subtract the zero-input baseline, average the output in centred 5 ms windows
at successive frame boundaries, and calculate
`D_n = -(Vbar_n - alpha*Vbar_(n-1))`, with `alpha=exp(-frame/tau)`.
The recovered ion count is `D_n/g`, where g is the calibrated gain above.
The same averaging and frame-difference filter is applied to the one-sided
SPICE output-noise spectrum. These operations are performed offline in Julia.
The external capacitance and averaging width are common to both frame choices.

The physical input is uniform mean collection over a frame, with an on-axis
constant-speed 3 m/s trajectory through the final 5 mm above a 0.5 mm square
pad. The pad weighting potential is one; the surrounding infinite grounded
plane and infinity are zero, with homogeneous dielectric. The integrated
current per ion is e[1-phi_w(g)]. Collection begins at simulation time 10 ms;
induction begins before collection. The waveform keeps the finite-path charge.

## Saved results and checks

`parameters.json`: complete circuit, readout, physical input, calibration and
convergence results, with source hashes. `operating_point.csv`: DC voltages and
currents. `calibration.csv`, `noise_budget.csv` and `compression.csv`: selected
results from the saved Julia analysis. The latter uses 8000 double-beta and
8000 electron events for each geometry, including empty buckets. Its event
fractions count at least one locally compressed bucket.

`report/report.tex` and its figure describe the complete retention scan and
validation. Numerical checks cover charge conservation, transient resolution,
tail fits, noise integration, gain linearity, and nonlinear-history refinement.
`reference/` saves the design and study cards. `manifest.json` records file
checksums and their repository sources. The HDF5 numerical data remain in
the repository's local results directories.

## Reproduction in the repository

    python scripts/electronics/feedback/run_retention_comparison.py --simulator /path/to/LTspice
    python notes/physics/tools/package_csa_v3.py
