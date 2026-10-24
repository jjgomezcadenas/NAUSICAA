# Rogelio's 180 nm CSA

The seven original files in this directory were received from Rogelio and saved
unchanged in commit `7766b7d`. The `.asc` schematic and `.inc` transistor models
are the input sources. The three PDFs and SQLite `.db` contain supplied results.
Do not replace these sources with derived simulations.

## Reproduce the new study

From the repository root, with the Python scientific dependencies installed:

```sh
make PYTHON=/path/to/python3 LTSPICE=/path/to/LTspice.app/Contents/MacOS/LTspice ltspice-study
make PYTHON=/path/to/python3 test physics
```

The runner was verified with **native macOS LTspice 17.2.4**, downloaded from the
[official LTspice page](https://www.analog.com/en/resources/design-tools-and-calculators/ltspice-simulator.html).
The version-26 macOS Wine runner stalled even on a resistor test; it was not used
for the reported calculations. The simulator is not bundled in git. Source,
configuration and simulator hashes are recorded in `results/ltspice/study.json`.

`scripts/build_ltspice_study.py` converts the schematic's wires, labels and known
standard-symbol pins into a text netlist. It checks the MOS/feedback connections
against an independently read node list. It preserves all original device sizes,
biases, models and active directives for the original case. A second original
run tightens the numerical settings and limits the timestep. The native solver
requires a text netlist, not a schematic, for batch execution.

The study then uses 5 pF **external pad capacitance**, in addition to intrinsic
transistor capacitances, at 20 C. The paper's analytical 5 pF was treated as total
input capacitance: this distinction is explicit and should be retained when
fitting circuit parameters. No gain, pole or feedback resistance is imposed on
the transistor circuit. Cf=3 fF and tau=0.3 s are retained in the *measurement
estimator* to compare with the existing paper.

Inputs are uniform mean ion arrivals over 10 ms, using the central square pad in
a grounded plane and the nominal constant final-gap drift speed. The PWL integral
is checked against q Delta(phi_w), without renormalization. Individual stochastic
arrival fluctuations and track distributions are not included in this test.

Outputs: `results/ltspice/README.md`, JSON, gzip CSV waveforms, uncompressed noise
spectra, AC response, diagnostic plot, portable `.cir`/`.pwl` decks, and UTF-8 run
logs. Full simulator `.raw` and UTF-16 `.log` files are reproducible and ignored by
git. The original transient CSV is decimated for display; integral/peak results
are calculated from the full simulation. Other waveform CSVs retain the solver
points. All calculations are in tracked Python code.

## Findings and questions for Rogelio

- The native simulation reproduces approximately 16 fC input charge and the
  supplied operating voltage/recovery behavior. The schematic comment
  `10 ke = 16 fC` is a factor-of-ten charge conversion error.
- M5 uses `Lcs/Wcs/Mp`; its separate `Lp5/Wp5/Mp5` parameters are unused. This
  remains unchanged until the intended device size is established.
- With the adopted 5 pF pad and estimator, the model requires about 1,600 ions
  per 10 ms frame for S/N > 5. This is much less sensitive than the illustrative
  linear model in the paper. See JSON for the adjacent-count check and convergence.
- The supplied model's `1overf` output-noise columns are all zero. This is not a
  prediction of zero physical flicker noise. Establish foundry/model provenance,
  intended noise model and calibrated flicker parameters before circuit claims.
- M4 and M6 account for most of the calculated variance. The model also reports
  DC supply power; see `supply_power_W` in the JSON. These are concrete quantities
  for circuit optimization and a future check against the intended PDK.
- Small-signal stationary noise is used. The nominal estimator discharge time is
  not a fitted MOS discharge time. Further work should examine feedback retention,
  optimal sampling, output loading and transistor parameter variation.

Nothing has been sent to Rogelio, and these results have not changed the paper.
