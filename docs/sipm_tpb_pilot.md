# S14160-1315PS and TPB pilot

Compare the same first 100 bb0nu events at E/P = 2.0 and 2.5 kV/cm/bar,
for pure Xe and Xe + 15% He. The sensor is a centred 1.3 x 1.3 mm square
at the exit of the 5 mm diameter, 5 mm deep EL hole. Buckets are 250 ns.

## Optical model

Emission is isotropic and uniform along the hole axis, with absorbing walls.
For square half-side a and depth coordinate u, the accepted fraction is
atan(a^2 / (u sqrt(u^2 + 2a^2))) / pi. Integrate over the hole depth.
The square fits inside the old circular aperture, so existing incident-photon
histories can be independently thinned with probability A_square/A_circle.
Binomial composition preserves the original EL fluctuations and gives the
same marginal distribution as collecting directly onto the square. Save new
VUV photon histories and their parent hashes, geometry, seeds and code.

TPB covers the sensor surface. Adopt complete VUV absorption and a probability
of visible re-emission of 0.62 at 175 nm, an approximate interpolation of the
intrinsic QE curve in Benson, Orebi Gann and Gehman (2018), Figure 19.
The plot gives about 0.65 at 170 nm and 0.50 at 190 nm; its quoted typical
relative QE uncertainty is about 8%. This is distinct from the sample-dependent
external wavelength-shifting efficiency plotted on the same figure.
Use the user's forward fraction 0.5, without reflected-light recovery.
For the pilot, approximate blue PDE by its 0.32 peak value, giving effective
VUV detection probability 0.62 x 0.5 x 0.32 = 0.0992. Sampling this product
as one binomial draw reproduces all three independent Bernoulli stages.
Assign converted photons to their parent bucket (prompt-TPB approximation).

## Sensor and electronics

Use 7,284 microcells, gain 3.6e5, and 120 kHz total dark rate at the datasheet
operating point (25 C, 4 V overvoltage). The proposed provisional recharge
constant is 12.5 ns, with a 25 ns sensitivity run. The 12.5 ns number is a
measured output-pulse fall-time proxy, not a direct microcell recharge fit;
Hamamatsu's technical note demonstrates restored pulse height at 50 ns spacing.
The user approved this provisional constant and sensitivity check on 2026-11-04.

Keep the ASIC limit at 1 mA. Preserve the earlier input charge noise:
0.5 x 2.55e6 / 3.6e5 = 3.5417 PE RMS for this gain. Keep signed noise.
The statistical model retains the established uniform microcell illumination
and constant PDE during recovery approximations. Time substeps resolve recharge.

## Outputs and checks

Julia saves new photon histories, standard sensor/ASIC shards, and reduced
statistics under results/sipm_tpb_pilot. Python generates an HTML comparison
and a self-contained LaTeX/PDF under latex_generated/sipm_tpb_pilot.
Report incident/detected/sensor/ASIC totals, pulse distributions before clipping,
clipped cells/events/charge, mean and maximum positive-input pulses, and
event-total widths after known attachment correction and normalization by
deposited MC energy. Separate optical statistics, recovery and clipping.
Widths are forward-response widths, without a fitted recovery correction.

Check quadrature, binomial moments, geometry matching, dark-rate scaling,
charge/noise ledgers, source event checksums and statistical time convergence.

## Clipping-tail audit

`scripts/production/audit_sipm_tpb_clipping.jl` reads the completed Xe, E/P=2
pilot and writes `results/sipm_tpb_pilot/clipping_audit.h5`. It compares the
sample-integrated charge loss with the distribution of individual event losses,
identifies the largest losses and their drift distances, and measures the
contribution of the worst events to the full response variance. Adding the
saved clipping loss back to the ASIC sum recovers the signed noisy input sum
and isolates the effect of the upper range limit. The unclipped subset width
is reported with its event count; that subset preferentially contains less
concentrated tracks. Full-sample RMS includes the low-response tail.

## Sources

- https://doi.org/10.1140/epjc/s10052-018-5807-z
- https://arxiv.org/abs/1201.2018
- https://www.hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/s14160-1310ps_etc_kapd1070e.pdf
- https://www.hamamatsu.com.cn/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/mppc_kapd9008e.pdf
- https://indico.bnl.gov/event/20473/contributions/84622/attachments/51960/88885/BIC%2520SiPM%2520Choice%2520Intro.pdf
