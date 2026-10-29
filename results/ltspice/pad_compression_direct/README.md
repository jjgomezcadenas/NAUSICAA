# Event compression with direct-collected histories

Configuration: `config/pad_history_viewer_direct.toml`.
Completed result identity: `8b0d8b859f56`.

The physical calculation uses the corrected 2 mm / 20 ms and 4 mm / 40 ms
ion histories, with the same optimized CSA and mean-current model as before.
`plan.h5` records the complete conservative candidate set and fixed-seed
controls. `statistics.h5` holds per-event flags and affected charge for each
relative-response-loss/absolute-deficit criterion. `presentation.h5` holds
validated report tables. Individual cases retain their full sampled histories,
reduced curves, currents, exact netlists and model/simulator provenance.

```sh
python scripts/run_pad_compression.py \
  --config config/pad_history_viewer_direct.toml --workers 4
python scripts/cache_standard_pad_histories.py \
  --config config/pad_history_viewer_direct.toml
python scripts/view_pad_histories.py
```

The batch reuses matching completed case files and checks a reproducible subset
at finer current/solver resolution. The standard-menu cache is separate from
the candidate/control statistical selection. The historical binned-input runs
remain under `results/ltspice/pad_compression/`.

Report source and figures:
`latex_generated/pad_compression_direct/8b0d8b859f56/`.
Before/after report:
`latex_generated/pad_history_direct_comparison/report.tex`.
