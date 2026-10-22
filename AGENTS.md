# NAUSICAA working conventions

This repository develops a reproducible physical input-current model for a
positive-ion CMOS detector. Establish the physics before implementing electronics.

- Xe2+ means the molecular ion Xe_2^+, with charge +e, not atomic Xe^{2+}.
- Use SI internally; show units in configuration keys, outputs and figure axes.
- The revised design has 5 kV across the 10 mm IFG and 3 kV across the 5 mm
  pad gap, with 200 V/cm in the bulk. Older paper voltage schemes are superseded.
- Physical electric fields and weighting fields are different boundary-value
  problems. Fixed geometry weighting potentials do not depend on applied bias.
- State boundary conditions and approximations alongside every numerical result.
  The central square pad in a grounded plane is a benchmark, not a solved IFG.
- Preserve total signal charge: integral I dt = q Delta(phi_w). Do not append a
  second collection impulse. Do not silently renormalize finite-path integrals.
- A 10 cm straight track along drift is an extreme test, not the event ensemble.
- Keep physical calculations in src/ion_signal; scripts reproduce figures,
  tables and JSON outputs. LaTeX includes generated values rather than copies.
- Use meaningful numerical/physical tests: independent quadrature, boundary
  conditions, charge conservation, scaling and convergence.
- Run `make test physics` with a Python >=3.11 environment containing the pinned
  scientific dependencies. Run `make note` to compile the PDF with latexmk.
- Inspect the compiled PDF after material edits. Keep scratch files in work/.
- Cite primary sources in notes/physics/references.bib. Record user-provided
  revisions as design inputs, not as published measurements.
- Do not change the existing remote or publish/push without a user request.
