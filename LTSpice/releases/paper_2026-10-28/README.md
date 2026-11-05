# NAUSICA circuit files for Vicente Herrero — manuscript snapshot, 28 October 2026

This directory contains the exact LTspice netlists and input waveforms already
run for the current paper. The .cir files are runnable LTspice netlists. Open a
.cir file in LTspice and run it, or execute `LTspice -b path/to/file.cir` with the
appropriate executable name on your system. Keep each folder's files together.
The calculation used macOS LTspice 17.2.4. The source circuit was supplied by
Rogelio Palomo; the selected transistor dimensions and operating point are those
used in the manuscript. `Nausica_Image.jpg` shows the circuit topology.
`nausica.pdf` is the compiled paper, including Vicente Herrero as second author.

## Start here

- `schematic/NAUSICA_20ms.asc`: editable graphical circuit with the paper's
  100-ion calibration current loaded from `single_N100_r2.pwl`.
- `schematic/NAUSICA_baseline.asc`: matching zero-current reference.
- `schematic/NAUSICA_noise.asc`: same circuit with the paper's noise analysis.
- `schematic/NAUSICA_injection.asc`: same circuit with the short injection test.
- `schematic/180nm_professional.inc`: the exact model card used in the paper.
  Keep it beside the .asc and .pwl files. Open an .asc and press Run.

The editable schematics retain Rogelio's v3 drawing layout. Their parameters,
analysis settings and input sources have been synchronized with the archived
paper runs. Their parsed electrical elements and directives are checked against
the corresponding .cir files. The old exponential test current is replaced by
the actual ion-current PWL file. The PWL columns are time in seconds and current
in amperes; the source direction is from ground into the pad input.

The exact batch inputs and recorded run logs are also included:

- `ion_frames/single_N100_r2.cir`: 100 ions uniformly collected over a 20 ms frame,
  with a long tail for measuring discharge. Plot V(out) and I(Is1).
- `ion_frames/single_N0_r2.cir`: matched zero-input baseline. Subtract its output
  from the signal run before measuring the pulse voltage.
- `noise/scale_4_gate_0.5384774789_noise_400.cir`: output noise, 0.1 Hz–1 MHz.
  LTspice V(onoise) is an amplitude spectral density in V/sqrt(Hz).
- `injection/scale_4_gate_0.5384774789_injection_N100_r2.cir`: 100 elementary
  charges injected with a 1 us plateau and 10 ns edges, for the appendix test.

## What is included

`ion_frames/` contains both temporal resolutions of the archived packet runs,
including zero input, calibration and larger packets, the counts bracketing the
paper's illustrative S/N=5 level, and an eight-frame sequence. Each .cir refers
to its local .pwl input current. Collection is from 2 to 22 ms for an isolated
packet; induction starts before collection. r1/r2 use maximum steps of 5/2.5 us.
The input preserves the finite-path induced charge (about 0.9984 e per ion).
Its benchmark is a central trajectory at 3 m/s through a 5 mm gap above a 0.5 mm
square pad in an infinite grounded plane, with uniform gas permittivity, unit
pad weighting potential, zero potential on the rest of the plane and at infinity.
The IFG field is not solved by this benchmark.

`noise/` contains 200 and 400 frequency points per decade for convergence.
`injection/` contains zero, 10- and 100-charge calibration runs at both final
resolutions (r1/r2: 1/0.5 us maximum timestep).
Every folder has the exact `corrected.inc` used in its runs, including the
thermal and flicker noise model. Saved `.run.txt` files record simulator output.
`reference/` contains the design configuration and generated manuscript values.
The historical 10 ms source circuit is included there with a `.txt` suffix for
provenance; the runnable `ion_frames/` decks use the paper's 20 ms inputs.
`manifest.json` records each copied file's source and SHA-256 checksum.

## Circuit parameters

Supply 1.8 V; temperature 20 deg C; external input capacitance 1.25 pF plus
intrinsic MOS capacitances; feedback capacitance 3 fF; feedback gate 10 mV.
M3 is W=200 um, L=2 um; M4 is W=2.88 um, L=1.44 um. Other dimensions are
specified in every netlist. The archived gates are 0.538477478921 V (M3) and
0.556803417206 V (M4). These retain the original runs exactly. The paper displays
0.5385 and 0.5568 V; the separate rounding check changes gain by 0.27% and ENC
by 0.06%. Their long archived decimal strings are numerical operating points,
not a requirement on physical bias-source precision.

## How these runs give the paper's results

The raw transient output is processed outside LTspice. At successive boundaries
20 ms apart, the output is averaged over centred 300 us windows. The readout is
D_n = -(Vbar_n - alpha Vbar_(n-1)), with alpha = exp(-20 ms/tau0).
Fitting the current-free output tail gives tau0 = 2.454158 s. Applying this
readout to the uniform-arrival waveform calibrates the gain to 27.119 uV/ion.
The same averaging and difference transfer functions are applied to the output
noise power spectrum before integrating over 0.1 Hz–1 MHz. The resulting RMS
noise is 3.056 mV, or ENC=112.7 elementary charges. This noise is evaluated at
the stationary DC operating point. It is independent of packet population in
the linear model. Integrating the unprocessed output spectrum gives a different
measurement from this frame readout.

For the injection appendix, a trapezoidal filter with 300 us rise and 300 us flat
top is applied to both the injected-charge response and the noise spectrum.
The saved calculation gives ENC=55.4 elementary charges. The frame readout and
injection test use different filters and input waveforms.

## Reproduce reduction in the repository

The existing Julia-to-HDF5 analysis and Python presentation scripts remain in
the repository. From its root, run:

    python scripts/run_discharge_pipeline.py --julia /path/to/julia --simulator /path/to/LTspice
    julia --project=. scripts/analyze_m3_scan.jl

The manuscript snapshot uses results/ltspice/discharge_calibration/20ms/
(paper_readout.h5 for generated values) and results/ltspice/m3_scan/analysis.h5
for the injection comparison. Noise integration, baseline subtraction and
filtering are implemented by these analysis scripts; pressing Run in LTspice
produces the underlying transients or noise spectra.

Regenerate this package with:

    python notes/physics/tools/package_paper_ltspice.py
