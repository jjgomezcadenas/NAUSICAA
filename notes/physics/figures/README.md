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
