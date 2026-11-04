# Ion and SiPM response production specification

## First complete case

The first complete response case uses the existing 8,000-event `bb0nu` primary
realisation. Its ion input is the **3 mm pad / 30 ms frame** history. Its
electron inputs are the two **250 ns** incident-photon histories, for pure Xe
and Xe–15% He. They are alternative gas models with the same primary-pair draw
and event placement. The source indices are
`results/simulations/production_pilot/grid_scan_summary_238ddc5a75e0fd9b.h5`
and
`results/simulations/production_pilot/electron_250ns_summary_464425eac32a2c34.h5`.
The 2 and 4 mm ion configurations and `single_electrons` await review of this
first case.

## Response graph and units

| Node | Input | Main output | Time unit |
|---|---|---|---|
| Ion CSA | Sparse positive ions per pad/frame | Signed reconstructed ion charges per selected pad/frame | 30 ms |
| SiPM sensor | Sparse incident VUV photons per SiPM/bucket | Integer detected PE and real-valued sensor PE equivalent per selected SiPM/bucket | 250 ns |
| Electron ASIC | Sensor PE equivalent per SiPM/bucket | Signed, clipped readout PE equivalent per selected SiPM/bucket | 250 ns |

One ion carries `+e`. The CSA converts to SI current and voltage, then
calibrates its output in elementary-charge units. One fully recovered SiPM
avalanche carries `gain × e`; sensor PE equivalent is sensor charge divided
by that charge. The ASIC's bucket-average current is multiplied by the bucket
duration and divided by that same charge to obtain readout PE equivalent.
Recovery, Gaussian noise and clipping make the latter two values real rather
than integer. Negative ASIC estimates remain signed; saturation flags mark
the upper 1 mA limit.

## Saved data contracts

Every node uses the original sample name, source event ID and gas identity
where applicable. The ion response uses sorted signed pad indices and a full
time vector for each pad with input charge. It saves input ions, noiseless
response, sampled noise, signed reconstructed counts and frame covariance.
The HDF5/Python matrices are `(pad, frame)`; Julia reads `(frame, pad)`.
Zero-input frames, leading quiet frames and the recovery interval remain in
the selected histories. Spatial pixels without input charge are absent.

The electron response keeps the source's signed `(pixel_ix, pixel_iy, bucket)`
keys. The sensor file's main columns are `detected_pe::Int64` and
`sensor_pe_equiv::Float64` for each saved SiPM–bucket. The ASIC file's main
column is `readout_pe_equiv::Float64`, with `saturated::Bool`. Both stages
retain the selected SiPMs' global photon interval plus four quiet buckets at
each end, including cells with zero incident photons. A remote photon adds
its receiving SiPM only. The source photon files are referenced by path and
content hash; production files do not duplicate their event groups.

Each output records the exact parsed cards, source content hash, model/code
hashes, Julia version, RNG policy, selected event IDs, event count, bucket
duration, units and a completion flag. Sensor and ASIC files are independent
cache nodes. Changing only the ASIC card reuses the sensor file. A run manifest
lists the ion and both gas-specific electron response nodes, their hashes and
the shared realisation/pair identities.

## Adopted first response models

The ion branch uses the v3 circuit, 5 ms centred voltage averages, finite-gap
square-pad weighting benchmark, uniform mean ion arrival times within each
frame, a three-frame inverse and stationary correlated noise from the
quiescent one-sided SPICE spectrum. The ion card and reconstruction card are
`config/specs/response/ion_csa/csa_v1.toml` and
`config/specs/response/ion_csa/reconstruction_v1.toml`.

The SiPM branch uses binomial PDE 0.30, 7,854 effective microcells, 73.4 ns
recovery, uniform arrival times within each 250 ns bucket, uniform microcell
assignment and area-scaled primary dark pulses. The ASIC adds independent
0.5 PE RMS Gaussian noise per bucket and clips bucket-average current above
1 mA. The approved sensor and ASIC cards are
`config/specs/response/sipm/sipm_v1.toml` and
`config/specs/response/sipm_electronics/readout_v1.toml`. Crosstalk,
afterpulsing, intrinsic gain spread, recovery-dependent PDE and ASIC pulse
shaping are deferred model extensions. The simulated sensor and ASIC values
are assigned to their trigger bucket.

## Validation and execution gates

1. Reject stale CSA calibration if its saved card, source deck, transistor
   model, derived inputs or simulator hash differs from the current resources.
2. Validate source schemas, complete flags, exact event sets, source hashes,
   common primary-pair checksums and complementary drift lengths.
3. Verify every ion cell and event ledger, signed bucket alignment, finite
   output, temporal covariance and reconstruction additivity.
4. Verify `0 ≤ detected_pe ≤ incident_photons` for every SiPM–bucket, sensor
   charge bounds, complete selected-channel windows, ASIC noise and clipping
   behavior. Record clipped-cell fractions with explicit denominators and
   charge/PE-equivalent loss as a separate quantity.
5. Benchmark a bounded subset and an interrupted-run replay before all 8,000
   events. Save versioned HDF5 results atomically or in complete shards with
   hashes and a manifest. Julia computes numerical summaries; Python presents
   the saved results in HTML and a generated LaTeX/PDF report.

The physical source and interface detail is in `notes/physics/csa_design.tex`,
`docs/ion_csa_response.md`, `docs/sipm_electronics_model.tex` and
`docs/sipm_response.md`. This specification fixes the first production data
contract; those notes retain the derivations and validation evidence.
