# ITACA operating conditions: TOML, HTML page and plots — handoff

Date: 2026-11-15. Branch: `ITACA` (created from `main` at `c86a63108`; work
here). Requested by the user; written by the agent of the 2026-11-14/15
session.

## Final configuration — read this first (2026-11-15, supersedes the rest where they differ)

The user fixed the configuration. `config/itaca.toml` has been updated to it
by the session agent (same schema and keys as the draft, except the drift
velocity and diffusion keys, listed below). Take the values from that file;
do not re-derive design choices.

**Scope change, user's instruction:** no plots on the page for now. All the
numbers are defined, so the page is tables: the operating conditions and the
derived values (masses, drift voltage, the full potential scheme, EL voltage
and reduced EL field, photons per electron, drift times, diffusion). The
plots listed under "Deliverable 3" are not wanted now. Sections describing
17.5 bar, a 320 cm diameter, a 150 cm length, 300 V/cm and "open" drift
velocities are history.

### Inputs (in `config/itaca.toml`)

All values are confirmed by the user, the ion drift velocity included. The
`[status]` table has been removed from the file at the user's request: do not
show a status column or an "open values" section on the page. `[source]`
remains and every value has an entry.

| Quantity | Value | Status |
|---|---|---|
| Gas | Xe + 15% He by mole; xenon 90.2616% Xe-136 | decided |
| Pressure; temperature | 16 bar; 290 K | decided |
| Fiducial radius; length | 180 cm (diameter 360 cm); 130 cm | decided |
| Drift field | 350 V/cm | decided |
| Electron drift velocity | 0.5684 mm/us ± 2.67% (Magboltz) | decided |
| Electron transverse diffusion | 2.0115 mm/sqrt(m) ± 2.46% (Magboltz) | decided |
| Electron longitudinal diffusion | 3.4227 mm/sqrt(m) ± 4.61% (Magboltz) | decided |
| Electron lifetime | 100 ms | decided |
| Ion drift velocity | 16.5 cm/s (literature mobility; Magboltz does not give it) | decided |
| EL depth; voltage | 5 mm; 15 kV | decided |
| ELM | holes 2.5 mm diameter at 4 mm pitch; transfer gap 10 mm at 300 V/cm; light guides 100 mm, index 1.5, 0.5 mm radial gas gap; SiPM 3 mm side; TPB 0.65; forward fraction 0.5 | decided |
| EL yield model; variance ratio J | 139 photons/kV, threshold 0.80 kV/cm/bar; 0.030 | decided |
| Ion collection potentials | 8 kV / 3 kV / 0 V; SiPM plane 0 V | decided |
| W; Fano | 21.9 eV; 0.15 | decided |

Key changes in the file: `drift.electron_drift_velocity_mm_us`,
`drift.electron_transverse_diffusion_mm_sqrtm`,
`drift.electron_longitudinal_diffusion_mm_sqrtm` and their
`_relative_error` keys are new; the two `_mm_sqrtbar_sqrtcm` keys now hold the
Magboltz values converted to card units (0.8046 and 1.3691; the earlier cards
had 1.6 and 0.75); `drift.ion_drift_velocity_cm_s` and `_range` replace the
two `_provisional_range` keys (a single value, no range); `gas.density_model.mixing_rule` is
`additive_partial_volumes`; `el.sipm_plane_V` and the section `[elm]` are new.

The Magboltz runs used the real-gas density of the Geant4 material before its
correction (equivalent pressure 17.624 bar; 17.367 bar after it). With the
corrected density the velocity is expected about 3% higher, within the quoted
error. State this on the page next to the value.

### Derived values (session agent's calculation; use as test references)

Gas, 16 bar and 290 K, with the corrected density formula of the section
"Correction to the gas density" below:

| Quantity | Value |
|---|---|
| Xenon compressibility z | 0.9074 |
| Gas density | 83.50 kg/m^3 (xenon 83.07, helium 0.43) |
| Number density relative to STP, N/N0 | 16.14 |
| Equivalent ideal-gas pressure | 17.367 bar |
| Fiducial volume | 13.23 m^3 |
| Area of each end | 10.18 m^2 |
| Total gas mass | 1105 kg |
| Xenon mass | 1099 kg |
| Xe-136 mass | 994 kg |
| Helium mass | 5.7 kg |
| Track length relative to pure xenon at 15 bar, 290 K | +11% |

Drift region:

| Quantity | Value |
|---|---|
| Drift voltage (350 V/cm over 130 cm) | 45.5 kV |
| Reduced drift field | 21.9 V/cm/bar (plain pressure); 0.0807 Td |
| Maximum electron drift time | 2.29 ms |
| Maximum attachment loss at 100 ms lifetime | 2.3% |
| Electron transverse spread at 130 cm | 2.29 mm |
| Electron longitudinal spread at 130 cm | 3.90 mm (6.9 us) |
| Maximum ion drift time | 7.9 s |
| Ion thermal transverse spread at 130 cm | 1.36 mm |

Potentials (sensors at ground; see "Potentials" below):

| Electrode | Potential |
|---|---|
| Pads | 0 V |
| Pad-facing IFG face | 3 kV |
| TPC-facing IFG face | 8 kV |
| Transfer gap, drift side (8 kV + drift voltage) | 53.5 kV |
| EL entry electrode (+ 0.3 kV across the 10 mm transfer gap) | 53.8 kV |
| EL exit electrode (+ 15 kV across the EL) | 68.8 kV |
| SiPM plane | 0 V |

EL region and ELM:

| Quantity | Value |
|---|---|
| Voltage across the 5 mm EL | 15 kV |
| EL field | 30 kV/cm |
| Reduced EL field | 1.875 kV/cm/bar |
| Photons per electron | 1195 (the earlier cards gave 1251) |
| Light guide length | 100 mm (the simulated value, kept; policy "fixed_length") |
| Average grading field along the guides (68.8 kV over 100 mm) | 6.88 kV/cm (0.43 kV/cm/bar; EL threshold 0.80); stored as `elm.grading_field_kV_cm`, check it against the potentials |

The potential sequence and the grading field follow `docs/elm.tex` (table of
regions and the equation for the grading field). The ELM card records
31.5 and 46.5 kV for its electrodes: those were computed for 200 V/cm and a
1.2 kV cathode and are not the potentials of this configuration.

## Goal

Centralise the ITACA detector description and operating conditions in one
file, `config/itaca.toml`, and generate from it:

1. `notes/html/itaca_detector.html`, linked first in `notes/html/index.html`;
2. the plots on that page;
3. (later, not in this task unless asked) the values used by production jobs
   and by the Geant4 geometry macro.

`AGENTS.md` must then point to the file and the page as the place where the
ITACA conditions live.

## Rules that apply (read `AGENTS.md` first)

- Do not start any job (production run, Slurm job, long computation) without
  the user's explicit OK. Small plotting scripts and the Julia calculation of
  derived values (seconds) are part of this task; anything that reruns the
  production chain is not.
- No destructive step (rm, overwrite of existing data, git push/reset) without
  the user's explicit OK. A question from the user is not an instruction.
- Commit only when the user says "commit". Never push without a request.
- Julia calculates and saves numbers to HDF5; Python reads them and makes the
  page, plots and tables. Do not repeat physical calculations in Python.
- HTML pages follow `notes/html/event_analysis/mc_fast_trigger.html` and
  `sim_detector_response.html`: short items, concrete names, units stated,
  no jargon, no history. Plots carry their numerics inside (values, units,
  selections).
- Never use option menus; ask in plain dialogue. Fix typos in the user's text
  without asking; never rephrase it without asking.
- The user decides every physical design value. Where this document says
  "open", record the value as open; do not choose one.

## Decisions taken by the user (2026-11-15)

| Quantity | Value | Status |
|---|---|---|
| Gas | Xe + 15% He (mole fraction); xenon enriched to 90.2616% Xe-136 (isotopic abundances of NEXUS `GXeEnriched`) | decided |
| Temperature | 290 K | decided |
| Pressure | 17.5 bar | proposed by the user as the configuration to consider; treat as the working value and mark "proposed" until confirmed |
| Fiducial volume | cylinder, diameter 320 cm, length 150 cm (arXiv:2604.02443, Table 1) | decided |
| Drift field | 300 V/cm | decided (keep) |
| EL region | 5 mm depth; 16 kV across it; holes 2.5 mm diameter at 4 mm pitch | proposed with the pressure |
| EL light yield | same as now: about 1251 photons per electron | follows from 16 kV at 17.5 bar |
| EL variance ratio J | 0.030 (same number as the current cards) | decided |
| SiPM | 3 mm side on 4 mm pitch, to limit saturation of the electronics | decided |
| Electron lifetime | 100 ms | decided |
| Electron diffusion coefficients | 1.6 (transverse) and 0.75 (longitudinal) mm·sqrt(bar)/sqrt(cm) | decided (keep the card values) |
| Ion collection potentials | 8 kV (TPC-facing IFG face), 3 kV (pad-facing IFG face), 0 V (pads) | decided (keep) |
| Electron drift velocity | under investigation | open |
| Ion drift velocity | under investigation | open |

Rationale for 17.5 bar and 16 kV, in the user's terms: helium at 15% lowers
the xenon mass and lengthens the tracks; 17.5 bar at 290 K recovers the mass
and track length of pure xenon at 15 bar. The SiPM electronics saturates if
too many photons arrive, so more light is not wanted; 16 kV keeps the present
yield.

### Potentials: the sensors are at ground

Read `manuscripts/itaca3/itaca3.tex`, the paragraph beginning "The cathode is
instrumented with an Ion Focusing Grid (IFG)". The pads (0 V) and the SiPM
plane (0 V) are at ground. The TPC-facing IFG electrode is at 8 kV. The high
voltage is on the anode side: first EL electrode = 8 kV + drift voltage;
second EL electrode = first + EL voltage; the potential is graded back to
ground along the light guides. For 300 V/cm over 150 cm and 16 kV across the
EL: 53 kV and 69 kV. (The manuscript's example is 200 V/cm: 38 kV and 54 kV.)
Do not write "cathode at 45 kV".

### Two different reduced fields

Name them separately everywhere, each with its unit:

- reduced drift field: 300 V/cm / 17.5 bar = 17.1 V/cm/bar;
- reduced EL field: 32 kV/cm / 17.5 bar = 1.83 kV/cm/bar.

## Numbers already calculated in the session (check them, then use them as test references)

Gas density model (as implemented in
`simulation/nexus/source/geometries/ItacaChamber.cc`): xenon as a real gas
with the second-virial compressibility factor of NEXUS
(`z = 1 + (0.083 - 0.422/T_r^1.6) P_r/T_r`, `T_r = T/289.733 K`,
`P_r = P/58.420 bar`; NEXT DocDB 1183); helium as an ideal gas; mixture =
each gas at the total pressure weighted by its mole fraction. Enriched xenon
molar mass 135.685 g/mol; Xe-136 mass fraction 0.9041. Fiducial volume
12.064 m^3.

| Xe–15% He | density (kg/m^3) | xenon (kg) | Xe-136 (kg) | helium (kg) |
|---|---|---|---|---|
| 15 bar, 290 K | 78.9 | 948 | 857 | 4.5 |
| 16 bar, 290 K | 84.7 | 1017 | 920 | 4.8 |
| 17 bar, 290 K | 90.6 | 1088 | 984 | 5.1 |
| 17.5 bar, 290 K | 93.6 | 1124 | 1016 | 5.3 |
| 18 bar, 290 K | 96.6 | 1159 | 1048 | 5.4 |
| 15 bar, 273.15 K | 85.4 | 1025 | 927 | — |

Taking xenon at its partial pressure instead gives about 1.5% less mass
(1105 kg at 17.5 bar); state this as the uncertainty of the mixing rule.

References for comparison: ITACA paper, pure xenon at 15 bar and 300 K,
1050 kg quoted (Table 1; Table 5 gives 85.9 kg/m^3); pure enriched xenon at
15 bar and 290 K, 1115 kg by the model above. Track length scales inversely
with xenon atom density (helium counted as 2/54 of a xenon atom per atom);
relative to pure xenon at 15 bar and 290 K the mixture gives +17% (15 bar),
+9% (16), +2% (17), −1% (17.5), −4% (18 bar). The Monte Carlo files used
until now (`NextPrecdr`, natural xenon, 15 bar, NEXUS default 273.15 K) are
denser than all of these; mention this only as provenance.

EL yield model (current cards; Fernandes et al., JHEP 04 (2020) 034, Table 1,
15% He: 139 ± 1 photons/kV, threshold 0.80 ± 0.02 kV/cm/bar, measured near
1.2 bar): photons per electron = 139 × (V_EL[kV] − 0.80 × P[bar] × d[cm]).
At 17.5 bar and d = 0.5 cm: 17.5 kV → 1460 (E/P 2.00); 16 kV → 1251
(E/P 1.83); 15 kV → 1112 (E/P 1.71). At 15 bar, 15 kV → 1251. Use the
existing Julia function `el_photon_yield` (see
`scripts/production/run_framed_detector_response.jl`).

Drift velocities: literature status, not decisions.

- Electrons. Measured points of McDonald et al. and Magboltz lines are in
  `docs/handoffs/data/xe_he_electron_drift_velocity_literature.csv` (digitised
  from the figures, about ±0.01 mm/us; the Felkai and Simón entries were read
  by eye). In Xe–15% He the velocity is about 0.61 mm/us at 20 V/cm/bar and
  about 0.41 mm/us at 17.1 V/cm/bar; with a real-gas density correction the
  17.5 bar point corresponds to about 15.8 V/cm/bar and 0.33 mm/us. Present
  0.33–0.42 mm/us at 300 V/cm as a range, "provisional". No measurement exists
  above 9 bar. A separate task computes this with Magboltz
  (`docs/handoffs/magboltz_xe_he.md`); leave a place on the plot for its curve.
- Ions (Xe2+ means the dimer Xe_2^+). Reduced mobility in xenon: about 0.65
  cm^2/V/s (Neves et al., J. Chem. Phys. 133 (2010) 124316, as plotted in
  Cortez et al., JINST 12 (2017) P06012, Fig. 3; Cortez measure 0.67 ± 0.02);
  0.70–0.74 quoted in the abstract of Neves et al., NIM A 580 (2007) 66; the
  ITACA papers use 0.74 (table in
  `/Users/jjgomezcadenas/Papers/itaca/srcx/detector_design.tex`). The Neves
  papers were not read in full. Helium raises the mobility by about 17% by
  Blanc's law, using an unmeasured Langevin estimate of 15 cm^2/V/s for the
  ion in helium. With the real gas density (N/N0 = 17.83 at 17.5 bar, 290 K)
  and 300 V/cm: 10.9 cm/s (0.65, no helium effect), 12.8 cm/s (0.65 with
  helium), 14.5 cm/s (0.74 with helium). Present 10.9–14.5 cm/s as a range,
  "provisional". No measurement exists at 10–15 bar.

Diffusion: the card coefficients (1.6 and 0.75 mm·sqrt(bar)/sqrt(cm)) stay.
Felkai et al. (NIM A 905 (2018) 82) give 2.5 mm/sqrt(m) transverse and about
4 mm/sqrt(m) longitudinal for 15% He at 15 bar (Magboltz), i.e. the card is
larger transversely and smaller longitudinally. State the source of each on
the page; do not change the card values.

## Inconsistencies found in the repository (the reason for this task)

| Quantity | ITACA paper | Repository |
|---|---|---|
| Temperature | 300 K | 293.15 K in the detector-response cards; 290 K in `simulation/nexus/macros/geometries/ItacaChamber.mac` |
| Drift field | 200 V/cm | 300 V/cm |
| Electron drift velocity | 1 mm/us | 0.6 in `config/readout/*detector_response*.toml`; 1.0 in `config/specs/detector/electrons/elm_4mm_250ns_ep2_shared.toml` |
| Electron lifetime | — | 30 ms and 15 ms in those two cards |
| SiPM size / pitch | 6 mm / 10 mm | 3 mm / 4 mm |
| Pressure | 15 bar | 15 bar in 39 cards under `config/` |

Existing cards are historical and must not be edited.

## Deliverable 1: `config/itaca.toml`

- Sections: `gas`, `geometry`, `drift`, `el`, `sipm`, `ion_readout`, `energy`
  (W = 21.9 eV, Fano 0.15 from the current cards), plus `status` and `source`
  tables keyed by `"section.key"`.
- Units in key names (`pressure_bar`, `field_V_cm`, `voltage_kV`, ...).
- Only the physical detector and its operating conditions. No run choices
  (event counts, seeds, frame margin, time buckets, cuts, paths) and no
  electronics model parameters (CSA, VH1, SiPM response, ASIC cards stay).
- Status values: `decided`, `proposed`, `provisional`, `open`. Every value has
  a source string (paper with table/figure, card path, or "user design input,
  2026-11-15").
- Whether ion pitch (3 mm), EL pitch (4 mm) and SiPM side (3 mm) belong in
  the file as nominal design values was asked and not answered: ask the user.
- Show the draft to the user before building on it.

## Deliverable 2: derived values (Julia → HDF5)

A short Julia script reads `config/itaca.toml`, computes the derived
quantities and saves them with the file's SHA-256: gas density, xenon /
Xe-136 / helium mass, fiducial volume, reduced drift field, reduced EL field,
EL photons per electron, electrode potentials, electron and ion maximum drift
times (as ranges while the velocities are provisional), transverse and
longitudinal electron spread and thermal ion spread at full drift (use the
existing functions, e.g. `sigma_t_ion_mm`). Add tests against the reference
numbers in this document. Put the code where the repository conventions place
it (trace the existing modules first; do not add physics to
`src/detector_response`).

## Deliverable 3: the page and its plots (Python)

Page sections, following the tables of arXiv:2604.02443 (Tables 1, 2, 3, 5;
text extraction in `work/literature_transport/itaca.txt`, PDF in
`papers/ItacaRevisited.pdf`):

1. Gas. 2. Geometry (fiducial volume; vertical stack of Table 2). 3. Drift
region (field, potentials, velocities, drift times, diffusion, lifetime).
4. EL region. 5. SiPM plane. 6. Ion readout. 7. Energy (W, Fano).
8. Open values, each with what is needed to close it.

Plots, in this priority (the user asked for the session's tables as plots):

1. Xenon mass and relative track length against pressure (15–18 bar), for
   290 K and 273.15 K, with the two pure-xenon references as horizontal lines
   and the operating point marked.
2. Electron drift velocity against reduced drift field: McDonald measured
   points (Xe–15% He and pure Xe) and Magboltz lines from the CSV; vertical
   bands for 15 bar and 17.5 bar at 300 V/cm, from the plain to the
   density-corrected reduced field.
3. Electron drift at 17.5 bar against drift field (250–450 V/cm): velocity
   and maximum drift time as a band; anode-side potentials on a second axis.
4. Ion drift velocity against drift field: band from 10.9 to 14.5 cm/s at
   300 V/cm scaled linearly in field; second axis for maximum drift time and
   the length of a 20 ms frame.
5. EL photons per electron against EL voltage at 15 and 17.5 bar, with the
   three options marked and the reduced EL field on a second axis.
6. Electron and ion diffusion against drift length (card coefficients and the
   Felkai values).
7. Potential along the detector axis.

Plots 1–5 first. One command regenerates derived values, page and plots
(a `make` target, following the existing ones).

## Deliverable 4: `AGENTS.md`

Add a short section near the top: the ITACA operating conditions are in
`config/itaca.toml` and `notes/html/itaca_detector.html`; read them before any
work that uses detector conditions; statements of 15 bar further down describe
historical productions. Show the wording to the user before editing; keep the
edit small.

## Not in this task (agreed direction, to be done later)

New-style detector cards will carry `operating_conditions =
"config/itaca.toml"` and no physical keys; a function `read_detector_card`
reads both files and returns the flat keys the code uses now; a card with both
a reference and a physical key is an error unless listed under `[overrides]`
with a reason; cards without the reference behave as now. The Geant4 geometry
macro will be generated from the same file. Design the TOML so that this is
straightforward, but do not implement it here.

## What exists

| Item | Path |
|---|---|
| Geant4 geometry and gas model | `simulation/nexus/source/geometries/ItacaChamber.cc`, macro in `simulation/nexus/macros/geometries/` (currently 15 bar, 290 K) |
| Current operating values | `config/readout/bb0nu_detector_response_3mm_mainframe20.toml`, `config/specs/detector/electrons/elm_4mm_250ns_ep2_shared.toml`, `config/specs/response/sipm/elm_pitch4mm_3015ps_tau12p5ns_recovered.toml` |
| Principle of operation and earlier velocity estimates | `notes/physics/itaca_poo.tex` |
| Manuscript | `manuscripts/itaca3/itaca3.tex` |
| ITACA paper | `papers/ItacaRevisited.pdf` (arXiv:2604.02443) |
| Literature text extractions (untracked scratch) | `work/literature_transport/` |
| Digitised electron drift velocities | `docs/handoffs/data/xe_he_electron_drift_velocity_literature.csv` |
| Documentation index and style examples | `notes/html/index.html`, `notes/html/event_analysis/` |

## Report back

To the user, with: the draft TOML, the list of values still open, the page and
plots, the test results, and any number in this document that could not be
reproduced.

## Correction to the gas density (2026-11-15, later the same day)

The mixture density in this document and in the first version of
`ItacaChamber.cc` weighted the two gas densities by mole fraction. That is not
consistent with an 85/15 mixture by mole. The geometry source now uses
additive partial volumes consistently:

    rho = P (x_xe M_xe + x_he M_he) / ((x_xe z_xe + x_he) R T)

The values in the table "Xe–15% He" above are 1.4–1.6% too high. Use these
(fiducial volume 12.064 m^3) as the test references instead:

| Xe–15% He | density (kg/m^3) | xenon (kg) | Xe-136 (kg) | helium (kg) | equivalent ideal-gas pressure (bar) |
|---|---|---|---|---|---|
| 15 bar, 290 K | 77.9 | 935 | 845 | 4.9 | 16.195 |
| 16 bar, 290 K | 83.5 | 1002 | 906 | 5.2 | 17.367 |
| 17 bar, 290 K | 89.2 | 1070 | 968 | 5.6 | 18.551 |
| 17.5 bar, 290 K | 92.1 | 1105 | 999 | 5.8 | 19.148 |
| 18 bar, 290 K | 95.0 | 1140 | 1030 | 5.9 | 19.749 |
| 15 bar, 273.15 K | 84.0 | 1008 | 911 | 5.2 | 16.450 |

The number density relative to an ideal gas is 1/(x_xe z_xe + x_he): N/N0 =
16.14 at 16 bar and 290 K. The pure-xenon references are unchanged. The
relative track lengths and ion velocities quoted above shift by the same
1.4–1.6% and must be recomputed from the corrected density.

The configuration is now fixed: see "Final configuration" at the top.
