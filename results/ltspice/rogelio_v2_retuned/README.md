# v2 with restored output baseline

Exploratory bias-only adjustment, separate from Rogelio's supplied proposal.
The M4 width and length and all other transistor sizes remain unchanged.
Julia bisects the M4 gate voltage to match the output baseline of the paper
circuit at the matched input capacitance and temperature.

`tuning_raw.h5` preserves each SPICE operating point; `tuning.h5` stores the
selected gate and complete search. The generated configuration is
`config/rogelio_v2_retuned.toml`. `raw.h5` and `analysis.h5` contain the same
ion-response, noise, convergence and threshold measurements used for the
supplied proposal. Reproduction and reporting commands are in
`../rogelio_v2/README.md`.
