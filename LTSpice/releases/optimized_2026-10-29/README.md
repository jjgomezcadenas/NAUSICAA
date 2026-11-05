# Optimized NAUSICA simulation — 29 October 2026

Snapshot for Vicente Herrero and Rogelio Palomo. The selected circuit is
`length_400um_gate0_rounded` from the feedback/averaging optimization. The paper's
previous circuit remains in the separate `LTSpice/paper_2026-10-28/` snapshot.

## Open and run

- `schematic/NAUSICA_40ms_400ions.asc`: complete editable circuit, 400 ions / 40 ms.
- `schematic/NAUSICA_20ms_120ions.asc`: 120 ions / 20 ms.
- `schematic/NAUSICA_40ms_baseline.asc` and `NAUSICA_20ms_baseline.asc`: matching zero-input references.
- `schematic/NAUSICA_noise.asc`: total output-noise spectrum, 0.1 Hz–1 MHz.
- `schematic/180nm_professional.inc`: exact transistor model card used in these runs.

Keep the entire directory together. Each current source reads a PWL file from
`../inputs/`: first column time in seconds, second column current in amperes.
Open an .asc in LTspice and press Run. Plot V(out), V(in) and I(Is1).
The exact batch netlists, readable run logs and their local model card are in
`netlists/`. Both temporal resolutions and all tested input populations are
included. The graphical circuits are electrically checked against those netlists.

## Selected parameters

Mf W=0.36 um, L=400 um, gate=0 V; M4 gate=0.6164 V. Other transistor dimensions
and biases are those of the selected paper circuit and are explicit in every
netlist. Supply 1.8 V, temperature 20 C, external input capacitance 1.25 pF plus
intrinsic model capacitances, feedback capacitor 3 fF. Model-card thermal and
flicker noise are included. Numerical settings disable waveform compression.

The readout averages the output over centred 5 ms windows at boundaries 20 or
40 ms apart, then takes D_n = -(Vbar_n - alpha*Vbar_(n-1)), with
alpha=exp(-frame/discharge_time). Averaging and differencing are performed in
the Julia analysis, after subtracting the matched zero-input output baseline.
They are applied to the signal and to the noise spectrum consistently.
`parameters.json` and `parameters.h5` contain the fitted discharge times, gains,
noise, ENC, neighbouring-frame responses and source hashes from the saved study.
`report/report.tex` and its scientific figures describe the full scan and checks.

The current uses uniform mean collection times within the specified bucket,
a central trajectory through the final 5 mm at 3 m/s, and the 0.5 mm square-pad
weighting benchmark in an infinite grounded plane with uniform gas permittivity.
The selected pad has unit weighting potential, the surrounding plane and
infinity zero. The finite-path integral e[1-phi_w(g)] is retained. Each isolated
packet begins collection at simulation time 2 ms. Both 2.5 and 5 us maximum
transient steps are included; the noise scan uses 200 and 400 points per decade.

The longer feedback retention will be tested with full sparse pad histories,
including zero-charge buckets and high-charge blob deposits. Those tests and ADC
range/saturation handling are the next study. The netlists here preserve the
already completed isolated-packet calculation, including inputs through 2000 ions.

## Reproduction

From the repository root:

    PYTHONPATH=src python scripts/run_feedback_pipeline.py --julia /path/to/julia --simulator /path/to/LTspice
    python notes/physics/tools/package_optimized_ltspice.py

Raw acquisitions and all scan variants are in `results/ltspice/feedback_optimization/`.
`manifest.json` identifies source files and SHA-256 checksums. The circuit study
used macOS LTspice 17.2.4.
