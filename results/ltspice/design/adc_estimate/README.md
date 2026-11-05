# Preliminary ADC estimate

User-provided populations at 10 ms are scaled to 20 ms at fixed arrival
rate. The saved rounded-circuit gain, calibrated discharge and noise set
the voltage and quantization scales. This calculation is a first-pass
linear estimate; no new 4000-ion transient is run.

```sh
julia --project=. scripts/analyze_adc_estimate.jl
julia --project=. tests/test_adc_estimate.jl
python scripts/present_adc_estimate.py
```

Inputs: `config/adc_estimate.toml`, with a hash of the source HDF5 circuit
analysis. Physical reduction: `src/ion_signal/adc_estimate.jl`. Output:
`analysis.h5`. Report: `latex_generated/adc_estimate/report.tex`.

Quantization assumes independent uniform errors in the two voltage readings,
and assigns no averaging reduction. The ADC signal span is the full voltage
range, including baseline and retained charge. The corrected frame increment
alone does not determine that range. All values and choices remain separate
from the manuscript until reviewed.
