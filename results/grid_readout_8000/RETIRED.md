The large legacy cache files were retired on 2026-11-08. Summary results remain.
See docs/reviews/storage_cleanup_2026-11-08.md and results/validation/direct_replay_fixture.
Historical regeneration uses scripts/transport/propagate_clouds_parallel.jl
followed by scripts/collection/propagate_grids_parallel.jl; both accept --events,
--workers and --output. The grid launcher also accepts --clouds.
