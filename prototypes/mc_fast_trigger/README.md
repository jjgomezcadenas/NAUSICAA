# Fast-trigger prototypes

Open `results/index.html` for both studies. Julia saves numerical results in
HDF5; Python builds HTML and figures from those results.

## Gamma deposits outside the frame

`fast_trigger_gamma_outside_frame.jl` processes the full 1 mm bb0nu and single-electron samples.
Its TOML card controls the event count, placement seed, detector Z interval,
vertex fiducial, frame margin and gamma energy threshold. Uniformly place the
primary vertex in 0–1500 mm and translate all hits together. Only events with
all deposits inside this detector interval enter efficiency denominators.
Apply the 50–1450 mm vertex cut, then reject an event if any gamma deposit
outside the primary box + 10 mm has deposited energy above 3 keV.

The primary box uses direct primary-electron hits. Gamma energy includes all
descendants of its first gamma ancestor, once. Delta rays belong to the
gamma-free track. D is the closest gamma-hit/primary-hit separation; R and Z
are components of that same pair. Outside-frame distances use the maximum
per-hit excursion. Particle ancestry comes from checksum-verified raw NEXUS
files. Numerical tables retain all events; hit arrays are saved for examples.

## Minimum energy

`fast_trigger_emin.jl` evaluates the analytic first and second moments for
1, 3, 5 and 10 keV deposits using the existing EL-yield and ELM-acceptance functions. Its
TOML selects the production cards, gas, energy, drift distances and numbers of
acquired SiPM–bucket cells. Fano, attachment, EL, optical, PDE, dark and ASIC
variance components are saved separately. It assumes linear recovered PE and
checks clipping and a conservative instantaneous microcell occupancy loss.
Acquired-cell count is a sensitivity parameter, not a simulated occupancy.

## Commands

Run from the repository root (use the available Julia executable):

```sh
export JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia"
julia --project=. prototypes/mc_fast_trigger/test_fast_trigger_gamma_outside_frame.jl
julia --project=. prototypes/mc_fast_trigger/fast_trigger_gamma_outside_frame.jl
MPLCONFIGDIR=/tmp/itaca-mpl work/event-viewer-venv/bin/python prototypes/mc_fast_trigger/present_fast_trigger_gamma_outside_frame.py
julia --project=. prototypes/mc_fast_trigger/test_fast_trigger_emin.jl
julia --project=. prototypes/mc_fast_trigger/fast_trigger_emin.jl
MPLCONFIGDIR=/tmp/itaca-mpl work/event-viewer-venv/bin/python prototypes/mc_fast_trigger/present_fast_trigger_emin.py
```

Each study writes `results/<study>.h5`, `results/<study>.html` and figures in
`results/<study>/`. HDF5 outputs stay local. The Python commands regenerate
presentation without rerunning calculations. These prototypes leave production
cards and libraries unchanged.

## Ion minimum energy

`fast_trigger_ion_emin.jl` and its TOML use the matched 4 mm / 40 ms linear
reconstructed-noise covariance. One bucket on each of K independent pads and
K consecutive buckets on one pad are separate cases. A central submatrix of
the saved 160-frame covariance includes all temporal correlations. No SPICE,
nonlinear response or ADC model is run. Ion attachment is zero. The same
presentation script adds these results to `results/fast_trigger_emin.html`.

```sh
julia --project=. prototypes/mc_fast_trigger/test_fast_trigger_ion_emin.jl
julia --project=. prototypes/mc_fast_trigger/fast_trigger_ion_emin.jl
MPLCONFIGDIR=/tmp/itaca-mpl work/event-viewer-venv/bin/python prototypes/mc_fast_trigger/present_fast_trigger_emin.py
```

The event-selection script now calls the MC reader and fast-trigger functions
in `src/event_analysis/`. The earlier MeV EventHits/read_event interface remains
available for existing transport studies. New fast-trigger inputs use
MCEventHits/read_mc_event with keV energies. All event decisions, placements and
frame limits are saved in `all_events`; displayed events also retain their hits.
