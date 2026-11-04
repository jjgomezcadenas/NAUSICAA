# Ion CSA v3 pilot results

`latest.toml` identifies the completed 100-event-per-sample, 3 mm / 30 ms pilot.
Its `outputs.json` lists the two versioned response HDF5 files, selected IDs,
checksums and matrix-axis convention. `summary.csv` gives the principal saved
quality, compression and charge diagnostics.

Open `ion_csa.html` inside the completed study for the standalone interactive
display. The display embeds Plotly and is generated locally; it and the large
numerical/case files stay outside Git. The LaTeX source and scientific figures
are under `latex_generated/ion_csa_v3_pilot/<study_id>/`.

Reproduce the calculation and display with:

```sh
work/event-viewer-venv/bin/python scripts/production/run_ion_csa_pilot.py
```

[Calculation, schema and validation](../../docs/ion_csa_response.md).
