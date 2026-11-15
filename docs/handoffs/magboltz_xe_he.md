# Magboltz calculation for Xe–15% He — handoff

Date: 2026-11-15. Branch: `ITACA`. Requested by the user; written by the agent
of the 2026-11-14/15 session.

## Goal

Calculate, with the Magboltz Fortran code itself, the electron drift velocity
and the transverse and longitudinal diffusion in xenon + 15% helium (mole
fraction) at the ITACA operating conditions, after first reproducing the
published values. The result decides the drift field and enters
`config/itaca.toml` and the ITACA detector page
(`docs/handoffs/itaca_detector_page.md`).

## Rules that apply (read `AGENTS.md` first)

- Do not start any calculation or job without the user's explicit OK. Getting
  the source, reading it, compiling it and writing the scripts need no run.
  Before the first run, tell the user exactly what will run (points, collision
  count, expected time, where) and wait.
- No destructive step and no push without the user's explicit OK. Commit only
  when the user says "commit".
- Use the Fortran code directly (user's instruction). Do not substitute
  PyBoltz or Garfield++; PyBoltz may be proposed later as a cross-check.
- Scratch files go in `work/`. Results that are kept go under `results/` with
  their inputs and provenance; the user keeps `results/` untracked on purpose.
- Julia performs calculations and saves to HDF5; Python makes plots. Here the
  calculation is Magboltz: a thin driver may run it and parse its output, then
  save a table (HDF5) with inputs, code version and source checksum.
- State boundary conditions and approximations next to every number. Cite
  formulas and data to primary papers.
- Do not include the user's name, project name or other personal details in
  any request to an outside service.

## Operating conditions

| Quantity | Value |
|---|---|
| Gas | Xe + 15% He by mole (85/15) |
| Temperature | 290 K |
| Pressure | 17.5 bar (proposed); also 15 bar for comparison with the current cards |
| Drift field | 300 V/cm (decided); scan 200–500 V/cm |
| Magnetic field | none |

Reduced drift field at 300 V/cm: 17.1 V/cm/bar at 17.5 bar, 20.0 V/cm/bar at
15 bar. In field per atom, with ideal-gas density at 290 K: 0.0686 Td and
0.0801 Td (1 V/cm/bar = 0.004004 Td at 290 K).

Do not confuse this with the reduced EL field (about 1.8–2 kV/cm/bar inside
the EL holes), which is not part of this task.

## Why it matters

In Xe–15% He the velocity rises steeply with field in this region. From the
measured points of McDonald et al.: about 0.61 mm/us at 20 V/cm/bar and about
0.41 mm/us at 17.1 V/cm/bar (slope about 0.05–0.07 mm/us per V/cm/bar; pure
xenon is flatter and faster there, about 0.70–0.77 mm/us). The cards use
0.6 mm/us. Moving from 15 to 17.5 bar at fixed 300 V/cm therefore costs about
a third of the velocity, and the result is sensitive to how gas density is
treated.

## Density: the main subtlety

Magboltz works in field per atom and converts pressure and temperature to
number density as an ideal gas. Xenon at 17.5 bar is denser than that. The
gas model used in this repository
(`simulation/nexus/source/geometries/ItacaChamber.cc`) is: xenon with
compressibility `z = 1 + (0.083 - 0.422/T_r^1.6) P_r/T_r`, `T_r = T/289.733 K`,
`P_r = P/58.420 bar` (NEXT DocDB 1183; z = 0.8987 at 17.5 bar and 290 K);
helium ideal; mixture by mole fraction at the total pressure. The original
partial-number estimate gave a factor of 1.0958 at 17.5 bar and 290 K. The
actual Geant4 material code combines partial mass densities and nominal mass
fractions, implying a number-density factor of about 1.1121. The user selected
the latter for the corrected Magboltz runs on 2026-11-15. The driver records
both formulas explicitly; the selected run cards use `geant4_material`.

Therefore:

1. Report every result against field per atom (Td) as well as V/cm/bar.
2. For our operating point, run Magboltz at the pressure that gives the real
   number density (the "equivalent ideal-gas pressure"), or equivalently at
   the corresponding E/N; say which. Also report the plain-pressure result so
   the size of the correction is visible.
3. Keep the composition by number at 85/15 and state the selected material
   density formula.
4. McDonald et al. plot against plain E/P at 1–9 bar and about room
   temperature; Simón et al. (NEXT-White) use a density-corrected reduced
   field, `E*/P = (E/P)(T/T0)(Z/Z0)` with T0 = 293.15 K (their Eq. 4.1). Keep
   the two conventions apart when comparing.

Magboltz contains no dense-gas physics beyond this scaling, and no measurement
of Xe–He exists above 9 bar (Fernandes et al., JHEP 04 (2020) 034, say the
pressure scaling of these mixtures is pending). Say so with the results.

## Reference values to reproduce first

Data file: `docs/handoffs/data/xe_he_electron_drift_velocity_literature.csv`
(columns: source, location, gas, kind, pressure label, reduced drift field in
V/cm/bar, velocity in mm/us). Text extractions of the papers are in
`work/literature_transport/` (`mcdonald.txt`, `felkai.txt`, `simon.txt`,
`fernandes.txt`; untracked scratch) together with the digitisation script
`extract.py`.

- A. D. McDonald et al. (NEXT), JINST 14 (2019) P08009, arXiv:1902.05544.
  Pure Xe, 10% and 15% He at 1, 3, 6 and 9 bar (labels rounded), 50–300 V/cm,
  drift length 14.1 cm. Velocity only in Fig. 7 top (pure Xe) and Fig. 8 top
  right (85/15); the points in the CSV were digitised from the PDF vector
  paths (about ±0.01 mm/us, ±0.3 V/cm/bar). They used Magboltz 11.8 and state
  that the velocity scales with E/P at all pressures and agrees with Magboltz
  ("very good" for Xe, "excellent" for mixtures; no percentage). Measured
  longitudinal diffusion exceeds Magboltz in the mixtures.
  Magboltz line as plotted, 85/15: 0.19, 0.35, 0.45, 0.62, 1.15, 1.78, 2.71
  mm/us at 10, 15, 17.14, 20, 25, 30, 40 V/cm/bar. Pure Xe: 0.32, 0.60, 0.70,
  0.74, 0.80, 0.86, 0.89.
- R. Felkai et al. (NEXT), NIM A 905 (2018) 82, arXiv:1710.05600. Magboltz
  11.1 at 15 bar, 300–500 V/cm. Transverse diffusion 2.5 mm/sqrt(m) for 15%
  He (10.5 pure Xe, 3.5 for 10% He); longitudinal "of the order of
  4 mm/sqrt(m)" for 10–15% He, with the 15% mixture peaking near 400 V/cm.
  Velocity read by eye from their Fig. 8 (±0.03): 85/15 about 0.63, 1.38,
  2.10 mm/us at 300, 400, 500 V/cm; pure Xe about 0.75 at 300 V/cm. Text:
  85/15 at 400 V/cm is 1.7 times the pure-xenon velocity.
- A. Simón et al. (NEXT), JINST 13 (2018) P07013, arXiv:1804.01680. Pure Xe:
  0.96799 ± 0.00017 ± 0.00406 mm/us at 7.178 bar, 22.29 °C, density-corrected
  reduced field 53.64 V/cm/bar; data lie 1.2% above Magboltz on average. At
  those conditions D_L = 318.9 and D_T = 1279 um/sqrt(cm).

Acceptance for the validation step: our Magboltz reproduces the plotted
Magboltz lines of McDonald within the digitisation error, and Felkai's values
within the reading error. Report any difference with the Magboltz version
used (theirs: 11.8 and 11.1; cross sections change between versions).

## Steps

1. Source preparation: the CERN page at https://magboltz.web.cern.ch/magboltz/
   labels the downloaded `magboltz-11.21.f` as version 11.21; its internal
   header says 11.20. The source is in `magboltz/source/` and has SHA-256
   `efa3244aa7c033872864997ee005b50d3732a59080e3957d13be9bf6d9569cbc`.
   No redistribution terms were found, so the source stays untracked.
2. The downloaded source's `READIN` routine specifies the input fields.
   `magboltz/MagboltzDriver.jl` prepares the cards and parses the output;
   `magboltz/test/runtests.jl` checks those conversions.
3. The source was compiled with Homebrew gfortran 16.2 and legacy fixed-form
   flags to `magboltz/build/magboltz-11.21`. No Magboltz calculation has run.
   On the Hyperion cluster, run jobs only with the user's OK and never on a
   login node.
   With thermal motion enabled, this source multiplies `NMAX` by 40 million;
   the prepared pilot's `NMAX=2` requests 80 million collisions.
4. Write a driver that writes the input for one point, runs Magboltz, and
   parses drift velocity and the transverse and longitudinal diffusion
   coefficients with their statistical errors. Convert diffusion to the units
   used in the repository cards (mm·sqrt(bar)/sqrt(cm), with
   sigma = D·sqrt(L/P)) and to mm/sqrt(m) at the stated pressure; give the
   conversion in the code and the report.
5. With the user's OK, time a single point; choose the collision count from
   the statistical error; then run the validation set (McDonald at 6 and 9
   bar for both gases; Felkai at 15 bar).
6. With the user's OK, run the scan: 85/15 at 290 K, 200–500 V/cm in steps
   that resolve the steep region (about 20 points), at 17.5 bar (plain and
   density-equivalent) and at 15 bar; pure Xe at 290 K for reference.
7. Save a table (HDF5) with inputs, results, errors, Magboltz version and
   checksum. Make the plots in Python from that table: velocity against
   reduced drift field with the measured points; velocity and maximum drift
   time over 150 cm against field at 17.5 bar; diffusion against field.
   Numerics inside the plots, units and conditions stated.
8. Report to the user: validation table (ours against published), the value
   and uncertainty of the velocity at 300 V/cm and 17.5 bar, the field needed
   for 0.6 mm/us, the diffusion coefficients compared with the card values
   (1.6 transverse, 0.75 longitudinal, in mm·sqrt(bar)/sqrt(cm)), and the
   limits stated above.

## Out of scope

- Ion mobility: Magboltz does not calculate it. The ion velocity is handled
  from the literature in the detector-page task.
- EL yield and the EL region.
- Changing any card or `config/itaca.toml`: report numbers; the user decides.

## What exists

| Item | Path |
|---|---|
| Gas density model | `simulation/nexus/source/geometries/ItacaChamber.cc` |
| Current card values (0.6 mm/us; 1.6 and 0.75 diffusion) | `config/readout/bb0nu_detector_response_3mm_mainframe20.toml` |
| Earlier estimate and its wording | `notes/physics/itaca_poo.tex` (section with Table "User design choices") |
| Bibliography entries (McDonald, Felkai) | `notes/physics/references.bib` |
| Literature points | `docs/handoffs/data/xe_he_electron_drift_velocity_literature.csv` |
| Paper text extractions (scratch) | `work/literature_transport/` |
