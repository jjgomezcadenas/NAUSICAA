# Model decisions and outstanding inputs

## Electron readout extension (2026-10-29)

Design input from the user: the single-electron or double-electron MC track
ionises the gas, producing electron and positive-ion trails. The existing
transport/pixel studies cover positive ions. Extend the simulation to the
electrons drifting to the anode, retaining both readout populations from the
same ionisation deposits.

Petit source audit:

- Loaded definitions: `../Petit/src/itaca_functions.jl:17` and `:24`;
  `../Petit/src/Petit.jl:27` includes that file.
- Electron widths in mm: `sigma_t_mm(Lcm,Pbar;dtmm)` =
  `dtmm * sqrt(Lcm/Pbar)`, and `sigma_l_mm(Lcm,Pbar;dlmm)` =
  `dlmm * sqrt(Lcm/Pbar)`.
- The numeric inputs are drift distance in cm and pressure in bar. Coefficient
  units consistent with this convention are mm sqrt(bar/cm).
- Pure-Xe pair: `(dtmm,dlmm)=(3.5,0.9)`, the function defaults and active
  `scripts/itaca_single_track_analysis.jl` CLI defaults.
- Xe–He pair: `(1.6,0.75)`, explicitly identified by the Xe–He variables in
  `../Petit/pluto/archive/itaca_st.jl:372–379`; the loaded
  `select_events_itaca` helper also defaults to this pair (`:447–449`).
- `get_sigma("electron",...)` (`itaca_functions.jl:556`) calls those two
  functions. Its electron branch uses L, P and the supplied coefficient pair;
  the `edrift` argument enters the ion branch. The MC diffusion routine then
  draws independent x/y Gaussians with sigma_t and a z Gaussian with sigma_l.

User revision on 2026-10-29: operating pressure is 15 bar and the Xe–He case
has a helium fraction of 0.15. These are design inputs, recorded with the Petit
coefficient pairs in `config/electron_transport.toml`. At this pressure the
coefficients multiplying sqrt(L_cm) are the Petit values divided by sqrt(15).
The formula evaluation and Float64 reference checks are saved to
`results/electron_diffusion/electron_diffusion.h5`; the generated report is
`latex_generated/electron_diffusion/electron_diffusion.pdf`.
Electron transport and readout timing need their own configuration and cache
identity. The inspected code uses supplied diffusion coefficients for the
electron branch rather than a field-dependent coefficient table.

### Paired drift and ELCC readout (2026-10-29)

User design inputs: sample one longitudinal event location in a TPC with
L_max = 1.5 m. Ions drift L to the cathode and electrons drift L_max-L to
the anode. Generate and save both diffused clouds before their respective
grid mappings, preserving the common event and ionisation-deposit identity.

The electron readout uses 5 mm diameter cylindrical holes on 10 mm pitch,
with 5 mm depth. Assign arriving electrons to holes using the existing
geometric collection approach. The reduced EL field is 3 kV/cm/bar; at the
design pressure of 15 bar this is 45 kV/cm and a 22.5 kV potential rise
along electron travel through the hole. TOML records entry-grid and exit
anode potentials as 0 and 22,500 V in a relative voltage gauge.

Electroluminescence produces VUV photons along electron passage through the
hole. The requested optical model accepts forward photons reaching the
circular SiPM covering the hole without intersecting a wall. Calculate
emitted and incident photon counts as intermediate quantities; save the
photoelectron response per SiPM and time bucket.

User refinement: compute geometric photon acceptance analytically and sample
that fraction, avoiding photon ray generation in each hole. The readout output
is detected photoelectrons per SiPM per time bucket. Electron
entry counts are intermediate physical inputs. The electron drift conversion
is 1 microsecond per mm; initially use 2 microsecond (2 mm) buckets. SiPM PDE
is configurable in TOML, initially 0.30. The inspected AXEL paper specifies
VUV-sensitive SiPMs but gives no numerical PDE. These inputs are recorded in
`config/electron_transport.toml`.

For isotropic emission uniformly along a centred axial electron path in a
cylindrical hole of radius R and depth h, with a SiPM disk covering the exit
plane and absorbing walls, the direct acceptance averaged over depth is
f = (h + R - sqrt(h*h + R*R))/(2*h). For R=2.5 mm and h=5 mm this is
0.190983 of all emitted photons. This is an axial-path benchmark; lateral
electron paths require their corresponding solid-angle average. The first
implementation holds the fitted EL yield at its mean and samples detected PE
with Binomial(round(N_electrons × yield), f × PDE). PE are assigned to the
parent electron's 2 microsecond arrival bucket.

EL source audit, 2026-10-29:

- `papers/ELCC/NEXT-EL-Xe-He.pdf`, NEXT, arXiv:1906.03984v3,
  Table 1: Y = A * max(E/P - X0,0) * P * h, photons per electron,
  with E/P in kV/cm/bar, P in bar, h in cm, A in photons/electron/kV.
  Pure Xe: A=137 +/- 1, X0=0.73 +/- 0.01. Xe-15% He:
  A=139 +/- 1, X0=0.80 +/- 0.02. At E/P=3, P=15, h=0.5,
  the respective central yields are 2332.425 and 2293.5 photons/electron.
  The mixture yield is 1.67% below pure Xe for these matched fits.
  Measurements were at approximately 1.2 bar with an additional systematic
  uncertainty around 5%; using reduced-field scaling at 15 bar is an
  extrapolation. Their absolute scale was normalised to an external pure-Xe
  result. Table 1 and the conclusion give approximately 3% reduction for
  15% He at E/P=2.5; a discussion paragraph instead says 6%, an internal
  inconsistency. Use the explicit fit parameters when reporting calculations.
- `papers/ELCC/EL-Next-DEMO.pdf`, NEXT, arXiv:2604.15905v1,
  Eq. (2): same yield form, quoting conventional A=140 and X0=0.83
  (2278.5 photons/electron at our operating point). Measurements span
  2.0--9.4 bar and show approximately 5% relative slope increase above
  5 bar. Absolute fitted A values are approximately 211--221, with an
  overall approximately 30% calibration systematic. This does not supply
  a measured 15 bar correction or a Xe-He correction.

The matched Table 1 pure-Xe and Xe-15% He fits are the initial configurable
prescriptions, preserving their provenance and the reduced-field scaling
approximation.

### Three-cloud propagation prototype and concurrency (2026-10-29)

The prototype `scripts/transport/propagate_cloud_stream.jl` writes one HDF5 cloud stream
for each combination of source sample (bb0nu, single electron) and carrier
model (positive ion, pure-Xe electron, Xe-15% He electron). It uses Petit's
per-carrier Gaussian diffusion and a 100-bin axis histogram. The event
position is sampled once from a source/event/base-seed identity and reused
in all three clouds; ion and electron distances sum to 150 cm. Diffusion
seeds also include the carrier model. The event position spans 0--150 cm.
Each process writes its own file. Pure-Xe clouds record helium fraction zero;
the Xe--He and ion clouds record the 15% mixture used by the detector.

`scripts/transport/benchmark_three_clouds.py` ran 100 events in each of six streams,
including Julia startup, HDF5 reading, diffusion and compressed HDF5 writing:
1 worker 39.03 s; 2 workers 19.98 s; 4 workers 13.59 s; 6 workers 7.63 s.
Total cloud output was 64.2 MB at
every concurrency. The six individual process peak resident-memory values
were approximately 0.78-0.83 GB each. Their sum, 4.90 GB, is an upper bound
on simultaneous peak resident memory, not a directly sampled total.

The numerical benchmark is saved in `results/three_cloud_benchmark/benchmark.h5`;
the cloud HDF5 files are in per-worker subdirectories. The generated PDF is
`latex_generated/three_cloud_benchmark/three_cloud_benchmark.pdf`.
`tests/event_analysis/check_three_cloud_benchmark.py` checks that all six
streams conserve counts, paired distances and produced carrier counts, and
that cloud arrays match exactly across worker counts.

This is a propagation and concurrency prototype. EL collection and PE
time-bucketing are subsequent readout stages. A production run still needs
the cache/resume layout and a full-sample resource check.

### Julia parallel launcher, 500-event sample (2026-10-29)

`scripts/transport/propagate_clouds_parallel.jl` now starts up to six independent Julia
processes, one per sample/cloud stream, with `--events`, `--workers` and
`--output` arguments. It writes `propagation_summary.h5` after all streams
finish. Run used:

```
julia --project=. scripts/transport/propagate_clouds_parallel.jl --events 500 --workers 6 --output results/three_cloud_sample_500
```

At 500 events per stream, six workers took 15.33 s including process startup,
reading, diffusion and HDF5 writes. The six cloud files total 310.3 MB.
Per-process peak RSS ranged from 0.819 to 0.864 GB; their 5.05 GB sum is an
upper bound on simultaneous peak RSS. The run summary is
`results/three_cloud_sample_500/propagation_summary.h5` and the generated
report is `latex_generated/three_cloud_sample_500/three_cloud_sample_500.pdf`.

Validation checked all 500 events in each sample for common paired drift
distances, carrier totals and sampled/retained/outside accounting. Every
saved position, count and drift distance in the first 100 events of all six
streams exactly matched the earlier single-worker benchmark.

### Full 8,000-event propagation (2026-10-29)

The six-process Julia launcher processed 8,000 events from each source sample,
producing 48,000 event/cloud records across six HDF5 files under
`results/three_clouds_8000/`. Wall time was 155.73 s and total compressed
cloud storage was 5,016.5 MB. Individual process peak RSS values were
0.95--1.00 GB; the sum of individual peaks is 5.82 GB, a conservative
upper bound on simultaneous peak RSS. The saved numerical summary is
`results/three_clouds_8000/propagation_summary.h5`; the generated report is
`latex_generated/three_clouds_8000/three_clouds_8000.pdf`.
The full event-by-event validator passed for both samples: each event's three
clouds share the same sampled location and produced carrier count, ion and
electron drift distances sum to 150 cm, and sampled counts equal retained
plus out-of-range counts in every cloud.

### Direct six-stream grid readout (2026-10-29)

`scripts/collection/sample_collect_parallel.jl` runs six independent Julia workers.
For each source hit, the worker samples its positive-ion or electron Gaussian
positions and increments the containing pad/EL hole and 2 mm longitudinal
bucket immediately. A sparse signed-index map extends the saved axes around
every retained sample. Individual carrier coordinates and diffusion display
histograms are not saved. The previous binned clouds and two-stage readouts
remain in `results/three_clouds_8000/` and `results/grid_readout_8000/` for
audit; their pad/time distributions are superseded.

Ion buckets span 20 ms at the configured drift speed. Electron buckets are
specified as 2 mm in `config/electron_transport.toml` and span 2 microseconds
at 1 mm per microsecond. Electron arrivals per hole/bucket are converted to
PE using the Table 1 mean EL yield, axial direct-light acceptance, PDE 0.30
and an exact binomial draw. Each PE stays in its parent arrival bucket. The
pure-Xe and Xe–He grids share one event phase. The direct HDF5 histories and
numerical summary are in `results/grid_readout_direct_8000/`; the generated
report is `latex_generated/grid_readout_direct/grid_readout.pdf`.
Run `make grid-readout` for the complete calculation, validation and report,
or `make present-grid-readout` to regenerate presentation from saved HDF5.

The full direct run processed 8,000 events per sample in 55.97 s with six
Julia processes; the six readout files occupy 977.1 MB. The event-wise
validator checked ion conservation, sampled/retained/outside accounting,
bucket widths, shared electron-grid phases and ensemble PE yields. For the
extended electron event 3790, the corrected ion readout has 66 occupied
time slices and 2,181 occupied pad/time cells, compared with 15 and 702
from the binned intermediate representation.

The two-panel viewer in `src/cloud_viewer/` reads the direct histories. Run
`work/event-viewer-venv/bin/python scripts/viewers/serve_cloud_viewer.py --port 8767`
and open `http://127.0.0.1:8767/`. Choose double beta or single electrons,
then pure Xe or Xe–15% He for the PE panel. The ion panel comes from the
same event in the selected sample. Each panel gives complementary drift
lengths, signal count, XY/X–time/Y–time heatmaps and a rotatable 3D view.
Navigation uses one-based positions in the file: 100 selects the hundredth
event, stored under `events/99`. Previous/next buttons and arrow keys step
through the selected sample. The URL retains the sample, gas and position.

The earlier ion-only pixel-charge workflow also has a direct sampler in
`scripts/collection/run_pixel_charge_direct.jl`. It replays each event's Gaussian draws
once to calculate the exact retained-carrier x/y centre and minimum z, then
replays the same draws into its configured pad and frame grid. This preserves
that viewer's grid-reference convention without a carrier-coordinate cache.
The 2, 3 and 4 mm full-sample histories and regenerated reports are in the
matching `results/pixel_charge_1mm*_direct/` and
`latex_generated/pixel_charge_1mm*_direct/` directories. Their old binned
counterparts remain available for audit. Electronics studies previously run
from the binned pixel histories retain their original provenance; they need
new input selection and acquisition before conclusions are updated.

Local source read: J. Hikida et al., “AXEL: High-pressure xenon gas TPC for
neutrinoless double beta decay search”, NIM A 1080 (2025) 170706,
doi:10.1016/j.nima.2025.170706, `papers/ELCC/Axel-ELCC.pdf`.
Figure 2 and section 2 describe field focusing, EL in holes at approximately
3 kV/cm/bar and SiPM detection behind a ground mesh. Section 4 describes
4.5 mm holes covered by 4.7 mm round SiPMs. Our 5/10/5 mm geometry is the
user's design input. This paper supplies no quantitative EL yield equation.

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
