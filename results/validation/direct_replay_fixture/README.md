# Legacy diffusion regression fixture

The cloud and readout HDF5 files contain single-electron events 0 and 3790 from
the old histogram-centre collection workflow. They retain the original event
arrays and metadata; event_count is set to two. Source identities are recorded
in manifest.json and the cleanup metadata archives.

Run from the repository root:

```sh
JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia" julia --project=. tests/event_analysis/check_direct_replay.jl
```

The test replays the original RNG and compares binned and direct assignments.
Original source hit tables and direct-readout files remain local.
