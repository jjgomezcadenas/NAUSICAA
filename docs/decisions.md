# Model decisions and outstanding inputs

## Scope and provenance (2026-10-22)

This milestone implements the first physics stage requested by J. J.
Gomez-Cadenas in the ion-track discussion. The previous chat estimates were
re-derived and checked rather than adopted as detector predictions.

- Molecular Xe2+ has charge +e. Atomic Xe^{2+} is a different species.
- Revised nominal fields: 200 V/cm bulk; 5 kV/cm in 10 mm holes; 6 kV/cm in
  the 5 mm final gap. Gate +8 kV, IC +3 kV, pads at zero.
- A circular 0.5 mm hole on 2 mm pitch has open fraction 0.049087. The
  entrance flux-capacity check is applied to the bulk-to-hole transition.
- The source preprint is arXiv:2604.02443v1. Its old IFG voltage scheme is
  superseded by the user's correction; its geometrical/transport context is
  useful but is not an instruction and does not validate the revised device.
- 1.5 m is the bulk reference drift, not track arc length. 293.15 K is a
  stated room-temperature assumption. Mobility is inferred from v=0.1 m/s
  at E=200 V/cm, not fitted to a transport measurement.
- The solid-angle single-pad formula assumes grounded surrounding metal.
  Conducting material between actual 0.5 mm pads on 2 mm pitch is not yet
  established. A disk benchmark illustrates shape dependence without bounding
  errors from unspecified surrounding material.

## What is verified

Independent integration of the half-space Poisson kernel reproduces the
square-pad formula. Its derivative, limiting boundary values, integral charge,
finite-start charge, velocity scaling and sample-spacing convergence are tested.
Thermal diffusion is checked against Einstein drift and first-passage moments
for a uniform field. These checks validate the benchmark implementation only.

## Required before a full IFG calculation

1. Pad outline, pitch, guard electrodes and all surrounding conductor potentials.
2. IFG construction, hole wall material/metallization, edge profiles and support.
3. Dielectric thicknesses, permittivities and whether surface charging matters.
4. Mobility and longitudinal/transverse diffusion versus reduced field at the
   operating pressure and temperature; purity and ion species stability.
5. Mechanical alignment and relevant finite-array boundaries.

## Subsequent physics milestones

Solve the physical and weighting fields separately. Validate boundary values,
mesh/domain convergence and charge conservation. Generate trajectory-dependent
single-ion kernels, including neighboring-pad signals and interception losses.
Then propagate actual NEXT/Geant4 energy deposits through drift, diffusion and
focusing to obtain pad-specific arrival distributions and current waveforms.
Use curved synthetic tracks only as labelled diagnostic examples.

The straight 10 cm track aligned with drift is an extreme benchmark. Real
events have isotropic global orientations and scattered paths, nonuniform
ionization, folds and endpoint blobs. A single-electron toy energy budget is
2.2 MeV body +0.3 MeV blob; a two-electron 2.5 MeV toy has 1.9 MeV in the two
bodies +0.3 MeV at each endpoint. These splits require fluctuations and actual
decay kinematics before quantitative event predictions.

Report pad charge, peak current, central-90% charge interval, full width at
half maximum where meaningful, repeated-packet separations and charge per
readout frame. Distinguish global event span from each pad's arrival intervals.

## Electronics question retained for the later milestone

For ideal linear exponential feedback and uniform arrivals lasting T, the
uncorrected amplitude fraction is (tau/T)*(1-exp(-T/tau)). Roughly tau=10T
gives 5% deficit and tau=50T gives 1%. These are analytical reference values,
not an established property of Rogelio's circuit. Packet shape, sampling and
decay correction determine the relevant requirement.

The screenshot has Cf=3 fF, an Mf geometry labelled L=100 um and W=0.36 um,
gate bias 10 mV, input about 16.7 mV and output about 1.0785 V. Its drain
current annotation about 1.11 pA does not establish a feedback decay time.
Actual netlist/process models and pulse tests are needed; numerical leakage,
including LTspice gmin, must be checked at the fA scale. The screenshot and
third-party model files are not redistributed in this repository.

## Analytical front-end extension (2026-10-22)

At the user's request, the paper covers the conceptual design of NAUSICA and
provides physical inputs for Rogelio's circuit simulations. The first analytical
readout model connects the existing finite-gap pad kernel to a linear CSA.
It uses Cf=3 fF as a starting design input from the supplied HTML calculator.
Feedback times are swept as independent parameters; no transistor bias-to-Rf
calibration or noise prediction is implied. One-pad arrivals are specified
single-ion, uniform, Gaussian, and double-Gaussian benchmarks, not an event
ensemble. Gaussian timing uses thermal bulk diffusion at the stipulated drift
speed. The 260-ion example is an illustrative population inspired by the ITACA
occupancy study, not a new reconstruction or full-packet measurement.

The original HTML is preserved in notes/physics/readout/rogelio_reference.html.
The offline adaptation and SPICE current exports are generated by make readout.
The detailed assumptions and numerical checks are documented alongside them.
