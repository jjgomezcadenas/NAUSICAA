# Ten-bucket sequence supplied by the user

The requested counts are [200, 200, 100, 400, 10, 10, 0, 5000, 200, 0].
This example uses 40 ms buckets, uniform collection times within each bucket,
5 ms centred averaging windows and the optimized grounded-feedback-gate CSA.
The current is the charge-preserving finite-gap planar-pad benchmark used in
the paper: square pad side 0.5 mm, gap 5 mm, speed 3 m/s. Each event starts at
the quiescent baseline, and no random noise is added to these mean transients.

The requested history and three controls remove bucket 8, bucket 9, or both.
Differences between paired control runs isolate the marginal contribution of
the 200 ions after the large deposit. Both transient time and current-knot
resolution are refined; all eight charge and eight transient checks pass.

`config/pad_sequence_example.toml` records all settings. Julia functions in
`src/ion_signal/pad_sequence_example.jl` generate the finite-gap PWL currents
(`build`) and reduce the SPICE waveforms to sampled voltages, frame estimates
and paired probe responses (`analyze`). `inputs.h5` stores current arrays;
the four case HDF5 files store acquired LTspice waveforms; `analysis.h5` stores
all quantities presented in the report and provenance. Exact model copies,
netlists, current files and readable run logs are under `decks/`.

From the repository root:

```sh
PYTHONPATH=src python scripts/run_pad_sequence_example.py --simulator /path/to/LTspice
python scripts/present_pad_sequence_example.py
```

The first command creates missing inputs, acquires or reuses matching SPICE
results and runs the Julia reduction. The second command renders the saved
HDF5 results. Report source and plots are in
`latex_generated/pad_sequence_example/`; no report compilation is performed.

The complete sequence remains near linear operation. Signal in neighbouring
buckets reflects induction before collection and the centred averaging windows.
The paired runs show only about a one-percent reduction in the subsequent
200-ion contribution. Exact values and the nominal quiescent ENC reference
are saved in `analysis.h5`.
