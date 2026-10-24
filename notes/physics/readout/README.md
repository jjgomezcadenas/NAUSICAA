# Analytical readout model

Open `analytical_model.html` in a browser. It is self-contained and works offline.
The interface adapts Rogelio's supplied CSA calculator by retaining its circuit
drawing while replacing the calculation with the Python model in
`src/ion_signal/readout.py`. `rogelio_reference.html` is the unchanged source
received through WhatsApp; that original uses external libraries.

From the repository root:

```sh
make PYTHON=/Users/jjgomezcadenas/miniconda3/bin/python3 test physics readout
```

Any Python >=3.11 environment with `requirements.txt` installed can be used.
`config/detector.toml` controls transport and pad geometry;
`config/readout.toml` controls the specified arrival distributions and CSA.
Outputs are `results/readout/benchmarks.json`, current files under
`results/readout/spice/`, and the offline HTML file.

## Physical assumptions

The single-ion current uses the existing central-square-pad weighting potential
in an infinite grounded plane and constant gap velocity. It starts at the gap
entrance and stops at collection. Its integral is exactly
`e * (1 - phi_w(gap))`, not artificially increased to `e`. Each computational
bin contains the exact endpoint-potential charge difference.

One-pad arrival benchmarks include a single ion, a uniform 10 ms packet,
a thermally diffused compact packet after 1.5 m drift, and two packets after
0.75 m drift with 40 ms between their centres. Gaussian probabilities are
integrated over each arrival bin, then placed at the bin centre for convolution.
Tails outside eight standard deviations are recorded, not renormalized.
These are expected currents for specified populations; they do not include
event morphology, fluctuations, the IFG collection map, final-region diffusion,
or an independently demonstrated distribution of ions per pad. The illustrative
260-ion packet must not be interpreted as a new event-ensemble result.

The CSA is an ideal amplifier with feedback capacitance Cf and linear resistance
Rf, obeying `dV/dt + V/tau = -I/Cf` and `tau = Rf*Cf`. Voltages are relative to
the DC baseline, initially zero. The time evolution is exact for each bin's
average current. The infinite-tau limit is the ideal integrator.
Output is not clipped. No bandwidth, finite gain, input-capacitance effect,
noise, power or stability is predicted by this first model.

Changing Cf at fixed tau in the interface changes Rf accordingly. Charge and
capacitance controls use exact linear scaling of precomputed responses. The
plots use fewer points than the numerical calculation; peak metrics are
computed at full resolution. Sampling is every 10 ms at a stated phase.
The finite plotting window can end before a long discharge has completed;
the JSON records residual output rather than treating it as lost charge.

## Paper outputs

The same generator produces `notes/physics/generated/readout_values.tex`,
`single_ion_response_table.tex` and two vector figures: `single_ion_current.pdf`
(full-gap and near-pad current, with distance and time-to-collection axes) and
`single_ion_response.pdf` (CSA voltage). The current figure shades the final
90% of the full-path charge; the finite-gap signal is not renormalized.
They use the configured single-ion case and the three shortest configured
feedback times, together with the ideal-integrator limit. Numerical metrics
are calculated at full resolution. The paper's `make note` target regenerates
both physics and readout outputs before compiling LaTeX.

## SPICE handoff

Each `.pwl` file contains two columns: time in seconds and current in amperes.
Time is shifted to start at zero. Use a current source directed into the CSA
input, for example `Iion 0 in PWL FILE="diffused_packet.pwl"` in an LTspice
netlist. File syntax and source orientation should be checked in the chosen
SPICE simulator. The original physical time offset, reference ion count and
integrated exported charge are stored in `benchmarks.json`.

The PWL interpolation slightly smooths bin-average currents. The generator
checks its integral against the model charge to a relative tolerance of 1e-4
and reports the error without rescaling. The browser's ion-count control does
not modify these reference files; edit the configuration and regenerate to
change the exported benchmark. A SPICE timestep convergence check is still
required, especially near collection.

## Provenance and validation

Original supplied file:
`/Users/jjgomezcadenas/Library/Containers/net.whatsapp.WhatsApp/Data/tmp/documents/84E8ABA1-46C5-4E11-AF9E-FFFDFD42F103/x_fab_180nm__Nausicaa_csa_designer_simulator.html`.
The starting Cf=3 fF is taken from that calculator. No claim is made that its
transistor constants are foundry-calibrated or that its simplified circuit
has been validated. Its artificial resistance cap, unresponsive controls,
fixed phase margin and inconsistent noise formulas are not used.

Tests check finite-path charge, arrival-tail accounting, packet convolution
charge, ideal integration, the closed-form rectangular-current response and
decay, linear scaling, sampling phase, and convergence to independent continuous
quadrature of the single-ion current with the CSA response kernel.
