# NAUSICA concept figure

`nausica_concept.tex` is the editable TikZ source. The paper includes the vector
`nausica_concept.pdf`; `nausica_concept.png` is a review preview.

Rebuild from `notes/physics` with:

```sh
sh figures/build_nausica_concept.sh
```

This requires TeX Live with TikZ and latexmk, and Poppler's `pdftoppm`.
Temporary LaTeX files and the build log stay in `work/figure-build`.

The figure is an original schematic for this paper, made on 23 October 2026.
Geometry and nominal voltages follow `config/detector.toml` at the repository
root. Holes are circular; the drawing shows a section through their axes and
compresses vertical distances. Arrows illustrate ion motion, not a calculated
field solution. Pad width is a design parameter and is not fixed numerically
by this drawing. The spaces between pads are not shown as grounded metal.

The circuit is an equivalent ideal CSA, with input current directed into its
inverting input, parallel feedback capacitance and effective resistance, and
sampling followed by an ADC. The resistor represents the MOS discharge path;
the schematic does not specify a transistor implementation, CMOS process, ADC
location, capacitance value or sampling period. The non-inverting input is
shown at the zero reference used for this conceptual model.

The inherited ITACA image remains in `img/topmetal_diagram.png` with its original
provenance in `img/SOURCES.txt`.

## Transistor CSA manuscript figures (26 October 2026)

`make spice-paper` generates `../generated/spice_circuit.pdf`,
`spice_ion_response.pdf`, `spice_noise.pdf`, `spice_devices_table.tex` and
`spice_values.tex` from the saved 1.25 pF / 300 us transistor study. The circuit
connections, device dimensions and gate biases are checked against the simulated
netlist. The waveform is the uniform one-ion mean arrival template; the noise
uses the same frame filter as the signal. No injection example enters the paper.
Source hashes and calculation conditions are in
`../generated/spice_paper_provenance.json`; PNG previews are scratch artifacts in
`../work/spice-paper/`. The generator does not invoke LTspice or compile LaTeX.
# ENC versus averaging time

The last active manuscript figure is `generated/spice_averaging.pdf`.
`make spice-averaging` runs `scripts/build_averaging_scan.jl` on the saved
1.25 pF SPICE waveforms and spectra, writes
`results/ltspice/averaging_scan.h5`, then renders the figure with
`scripts/plot_averaging_scan.py`. `make spice-averaging-plot` regenerates
the presentation alone. The scan uses 180, 200, 250, 300 and 350 us centered
averages for uniform arrivals over 10 ms, recalibrating the signal at each width.
Source hashes and numerical convergence checks are stored in HDF5.
The self-contained report source is `latex_generated/spice_averaging/report.tex`.
The previous `spice_noise.pdf` remains a diagnostic and is no longer included
in the paper. The manuscript and report have not been compiled for this change.
