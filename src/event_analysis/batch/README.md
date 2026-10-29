# Historical pixel batch helpers

`run_pixel_charge.jl` defines the original cloud-to-pixel batch calculation called
by `scripts/events/analyze_mc_selection.jl --pixel-charge`.
`pixel_diffusion_cache.jl` provides its versioned HDF5 cloud cache. The including
command supplies `ROOT`, the loaded EventAnalysis interface, HDF5, TOML and SHA.
`tests/event_analysis/test_diffusion_cache.jl` checks the cache contract.

The current direct ion-only pipeline starts at
`scripts/collection/run_pixel_charge_direct.jl`; the six-stream ion/PE pipeline
starts at `scripts/collection/sample_collect_parallel.jl`.
