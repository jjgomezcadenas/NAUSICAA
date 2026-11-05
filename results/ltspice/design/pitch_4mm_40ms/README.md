# 4 mm pitch / 40 ms circuit exercise

User design inputs: 4 mm pitch and a mean of 400 ions per bucket; comparison
uses 2 mm / 20 ms with the updated mean of 120 ions. Pad geometry, external
capacitance, circuit and operating point stay those of the paper. The local
square-pad finite-gap benchmark and uniform mean arrivals are retained; no
charge renormalization. Noise is stationary at the DC operating point.

Reproduce from the repository root:

```sh
PYTHONPATH=src python scripts/run_pitch_exercise.py --julia /path/to/julia --simulator /path/to/LTspice
```

Julia input generation and reduction: src/ion_signal/pitch_exercise.jl.
HDF5 inputs, raw acquisitions and analysis retain hashes and configuration.
The existing noise spectrum is reused only after its circuit/model hashes match.
Both transient resolutions, current PWLs, netlists and readable logs are saved.
Native LTspice raw/log scratch is ignored; raw HDF5 retains the numerical data.
Report source and figure: latex_generated/pitch_4mm_40ms/.
No manuscript values or previous 20 ms results are changed.
