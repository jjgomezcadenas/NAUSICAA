# NAUSICAA: ion signal physics

Reproducible calculations and a LaTeX note for positive molecular Xe2+ ions
collected by a focused CMOS pad plane. Each ion carries **+e**.

## First milestone

The implemented model covers detector geometry and the revised nominal fields,
constant-mobility transport, thermal diffusion, and an analytical single-ion
weighting-potential benchmark. It does **not** yet solve the actual IFG field,
simulate realistic event tracks, or evaluate the feedback transistor.

The revised design uses 200 V/cm in the bulk, 5 kV across a 10 mm IFG, and
3 kV across the final 5 mm gap. The pads are at 0 V, IC at +3 kV and gate at
+8 kV. These values supersede the older voltage scheme in the ITACA preprint.

## Reproduce

Use Python >=3.11 (the tested environment is recorded in results/physics_benchmarks.json):

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
make PYTHON=.venv/bin/python test physics readout
make PYTHON=.venv/bin/python note
```

The last command needs a TeX installation with `latexmk`, `pdflatex`, and
`bibtex`. It builds `notes/physics/build/nausica.pdf`. The LaTeX project can also
be opened in Codex at `notes/physics`; run the Python generation step first.
If dependencies are already available, set `PYTHON` to that interpreter.

- `config/detector.toml`: SI inputs, source labels, and unresolved assumptions.
- `src/ion_signal/physics.py`: reusable equations, current convention and models.
- `scripts/build_physics.py`: figures, LaTeX tables/macros and numerical outputs.
- `tests/test_physics.py`: independent surface integration, charge accounting,
  boundary values, drift/diffusion consistency, scaling and convergence.
- `notes/physics/nausica.tex`: the physics note and derivations.
- `notes/physics/references.bib`: primary sources.
- `notes/physics/generated/`: generated figures, tables and numerical macros.
- `results/physics_benchmarks.json`: numerical results, configuration hash and
  software versions.
- `docs/decisions.md`: provenance, limitations and next milestones.

Generated benchmark outputs are checked in so the note is inspectable without
first installing Python. Regenerate them after changing any input or equation.
The note's explanatory prose describes the baseline design; review it alongside
the generated tables whenever changing that design.
Build intermediates and scratch files are ignored. No random sampling is used
in milestone one. Future event simulations must record seeds and data provenance.

## Interpretation

The square-pad benchmark assumes a 0.5 mm square surrounded by an infinite
grounded conductor, a central normal trajectory and uniform local speed.
It is an analytical reference, not a prediction of the full perforated IFG.
The finite final-gap charge is kept as q Delta(phi_w); it is not silently
normalized to e. Full-geometry fields, measured transport and event samples
are the next steps before realistic packet distributions and electronics.

## Analytical front-end model

`make readout` generates a self-contained interactive model at
`notes/physics/readout/analytical_model.html`, numerical results in
`results/readout/benchmarks.json`, and SPICE input-current files in
`results/readout/spice/`. It combines the finite-gap square-pad current with
specified one-pad arrival distributions and an ideal CSA with linear feedback.
The original calculator supplied by Rogelio is preserved under
`notes/physics/readout/`; its transistor and noise estimates are not used.
See `notes/physics/readout/README.md` for assumptions, controls, and SPICE handoff.
