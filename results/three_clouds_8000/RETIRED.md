The large binned cloud files were retired in the approved 2026-11-08 storage
cleanup. Their metadata and SHA-256 identities are in results/storage_cleanup/
2026-11-08. Summary results remain here. Events 0 and 3790 survive in
results/validation/direct_replay_fixture for the legacy replay regression.

Historical regeneration: scripts/transport/propagate_clouds_parallel.jl
with `--events 8000 --workers 6 --output results/three_clouds_8000`. Current production uses direct physical assignment;
see docs/elm_matched_production.md.
