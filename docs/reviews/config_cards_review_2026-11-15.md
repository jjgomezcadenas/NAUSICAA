# Review of the cards under `config/` against `config/itaca.toml`

Date: 2026-11-15. Branch `ITACA`. Scope: every `.toml` under `config/` (274
files), their physical keys, their consumers, and what is needed so that the
operating conditions are loaded from one place.

## 1. Inventory

| Folder | Cards | Schemas |
|---|---|---|
| `config/specs/runs` | 86 | `sipm_run`, `sipm_production_run`, `detector_run`, `ion_csa_pipeline_run`, `ion_csa_run`, `ion_csa_reconstruction_run` |
| `config/` (top level) | 36 | study and tool cards, mostly without a schema; `config/detector.toml`, `config/itaca.toml` |
| `config/readout` | 35 | framed detector response, electronics, reconstruction, hit output, VH1 and track cards (no schema key) |
| `config/specs/generated/*` | 40 (8 × 5) | frozen copies of detector, electrons, ions, realisation and run cards of past productions |
| `config/specs/detector/configurations` | 12 | `detector_configuration` (pointers to realisation, ions, electrons cards) |
| `config/specs/detector/electrons` | 11 | `electron_detector_spec` |
| `config/specs/detector/ions` | 3 | `ion_detector_spec` |
| `config/specs/response/sipm`, `sipm_electronics`, `ion_csa` | 11, 2, 11 | `sipm_spec`, `sipm_electronics_spec`, `ion_csa_spec`, `vh1_csa_spec` |
| `config/studies/electronics/*` | 23 | LTspice study cards |
| `config/specs/realisations` | 3 | `realisation_spec` |

How the cards are reached: 149 are named in code, scripts, tests or the
Makefile; 50 only from other cards; 15 only from documentation; 60 from
nowhere. "From nowhere" includes the cards passed on the command line to
`run_hit_production.jl` and `run_vh1_pilot.jl` (all the `bb0nu_vh1_*` cards of
the current production), so it does not mean unused; it means the repository
has no record of which card produced which result except the copies saved in
the HDF5 outputs.

## 2. Contradictions with `config/itaca.toml`

No card other than `config/itaca.toml` carries the current conditions. The
cards of the current production chain (the `bb0nu_vh1_*` run cards name them)
hold:

| Quantity | `config/itaca.toml` | `config/readout/bb0nu_detector_response_3mm_mainframe20.toml` | `config/specs/detector/electrons/elm_4mm_250ns_ep2_shared.toml` |
|---|---|---|---|
| Pressure | 16 bar | 15 | 15 |
| Temperature | 290 K | 293.15 | — |
| Drift length | 130 cm | 150 (`chamber_length_mm`) | — |
| Drift field | 350 V/cm | 300 | — |
| Electron drift velocity | 0.5684 mm/µs | 0.6 | 1.0 (`drift_speed_mm_us`) |
| Electron lifetime | 100 ms | 30 | 15 (`lifetime_ms`) |
| Transverse diffusion, card units | 0.8046 | 1.6 | 1.6 (`[xe_he15]`) |
| Longitudinal diffusion, card units | 1.3691 | 0.75 | 0.75 |
| Ion drift velocity | 16.5 cm/s | 15 (`ion_drift_velocity_mm_ms` 0.15) | — |
| EL voltage | 15 kV over 5 mm (1.875 kV/cm/bar) | `EL_reduced_field_kV_cm_bar` 2.0 (15 kV at 15 bar) | 2.0 |
| TPC-facing IFG | 8 kV | 8 kV | `cathode_potential_V` 1200 |
| EL entry and exit | 53.8 and 68.8 kV | — | 31.5 and 46.5 kV (200 V/cm, 1.2 kV cathode) |
| Grading field | 6.88 kV/cm | — | 4.65 |

Across all 274 cards the distributions are: pressure 15 bar in 28 cards and
16 bar in one; temperature 293.15 K in 22 and 290 K in one; drift field
200 V/cm in 11 and 300 in 7 (350 in one); electron velocity 1.0 mm/µs in 19
and 0.6 in 7; lifetime 15 ms in 19 and 30 in 7; chamber length 150 cm in 20;
reduced EL field 3.0 in 12, 2.0 in 6, 2.5 in 2. Also `W_eV` 32.0 and `fano`
0.19 appear once each (`config/pixel_charge*.toml`, the early pixel studies),
and the pure-xenon EL fit (137, 0.73) sits beside the Xe-He fit (139, 0.80)
in every ELM card.

These cards are historical productions and are not to be edited. The
contradiction is structural: every card repeats the operating conditions, so
a new condition cannot be adopted without writing new cards for every stage.

## 3. One quantity, several names

The same physical quantity appears under different keys and units:

| Quantity | Keys found |
|---|---|
| Pressure | `pressure_bar`, `pressure_Pa` |
| Drift field | `drift_field_V_cm`, `field_V_cm`, `bulk_field_V_m`, `field_V_m` |
| Electron drift velocity | `electron_drift_velocity_mm_us`, `drift_speed_mm_us`, `electron_drift_speed_m_s` |
| Ion drift velocity | `ion_drift_velocity_mm_ms`, `ion_drift_velocity_cm_s`, `drift_speed_m_s`, `bulk_speed_m_s` |
| Electron lifetime | `electron_lifetime_ms`, `lifetime_ms` |
| Diffusion | `electron_dT_mm_sqrtbar_sqrtcm`, `d_T_mm_sqrtbar_per_sqrtcm`, `electron_transverse_diffusion_mm_sqrtbar_sqrtcm`, `electron_transverse_diffusion_mm_sqrtm` |
| W, Fano | `W_eV`, `w_eV`; `fano_factor`, `fano` |
| EL yield fit | `EL_gain_photons_e_kV`, `yield_gain_photons_per_electron_kV`, `gain_photons_per_electron_kV` |
| EL variance | `EL_variance_ratio`, `j_variance_ratio`, `variance_ratio_J` |
| Hole and pad pitch | `hole_pitch_mm`, `EL_pitch_mm`, `pitch_mm`, `pitch_m`, `default_pitch_mm`, `pad_pitch_mm`, `ion_pitch_mm` |
| Optics | `tpb_efficiency` / `tpb_conversion_efficiency`; `guide_index` / `guide_refractive_index`; `gap_index` / `gap_refractive_index` |
| SiPM side | `sipm_side_mm`, `side_m` |
| Helium | `helium_fraction`, `helium_mole_fraction` |
| TPC-facing IFG potential | `ion_TPC_face_V`, `gate_voltage_V`, `tpc_face_ifg_V`, and `cathode_potential_V` (a different value, 1.2 kV) |

Each name has its own reader in the scripts (`run_framed_detector_response.jl`,
`run_detector_pilot.jl`, `sipm_pipeline.jl`, `sample_collect_stream.jl`, and
others), so a value can be changed in one card and silently stay in another.

## 4. Proposal

1. **One reader for the operating conditions.** New-style detector and run
   cards carry `operating_conditions = "config/itaca.toml"` and no physical
   keys. A function `read_detector_card(path)` in `src/itaca/` returns the
   flat keys the scripts use now, taking the physical values from
   `config/itaca.toml` and the options (pitches, buckets, SiPM device,
   electronics, cuts) from the card. A card with both the reference and a
   physical key is an error unless the key is listed under `[overrides]` with
   a reason. Cards without the reference load as they do now, so historical
   productions still replay.
2. **One vocabulary.** New cards use the key names of `config/itaca.toml`.
   The reader maps them to the names the code expects, in one table, so the
   scripts need no change at first.
3. **A manifest of the active cards**, `config/active.toml`, listing the cards
   of the current production chain. Everything not listed is historical and
   stays untouched.
4. **Tests**, run by `make test-itaca`:
   - every card listed in `config/active.toml` carries
     `operating_conditions` and no physical key outside `[overrides]`;
   - `read_detector_card` on each active card returns pressure, temperature,
     field, velocities, lifetime, diffusion, EL depth and voltage, potentials,
     W and Fano equal to `config/itaca.toml`;
   - every card under `config/` parses, and every card that declares a
     `schema` has a known one;
   - every path a card names (other cards, SPICE summaries, calibrations)
     exists, and every recorded SHA-256 matches;
   - the derived values (`results/itaca_detector/itaca_conditions.h5`)
     agree with the card within the tests of `tests/itaca/`.
5. **Migration order.** Detector response card first (it feeds both arms),
   then the ELM card, then the run cards of the VH1 chain; the SiPM, ASIC,
   CSA and VH1 response cards keep their physical keys because they describe
   devices, not the detector.
6. **Retirement candidates**, for the user's decision, not for action here:
   the 60 unreferenced cards (list below), the pure-xenon `[pure_xe]` blocks
   in the ELM cards, and `config/detector.toml` (the 2026-10-22 ion benchmark,
   200 V/cm, 15 bar), which still feeds `scripts/benchmarks/`.

## 5. Cards referenced from nowhere

`config/pad_history_viewer.toml`, `config/studies/electronics/histories/pad_history_viewer.toml`,
`config/specs/detector/electrons/elcc_250ns_ep2_square1315.toml`,
`config/specs/detector/electrons/elcc_250ns_ep2p5_square1315.toml`;
under `config/readout/`: `bb0nu_detector_response.toml`,
`bb0nu_detector_response_3mm_frame20.toml`, `bb0nu_hit_output_100.toml`,
`bb0nu_reconstruction_channel_histograms_100.toml`,
`bb0nu_reconstruction_ion4mm_20ms_100.toml`, the eight `bb0nu_tracks_ion4mm_20ms_500*.toml`,
the eleven `bb0nu_vh1_*.toml` and `single_electrons_vh1_hit_output_full_5481_cf20f.toml`
(command-line cards of the current production), `single_electrons_detector_response_3mm_frame20.toml`;
under `config/specs/runs/`: `electron_250ns*.toml` (3), `elm_bb0nu_tau*_pilot100.toml` (2),
`elm_matched_csa_smoke.toml`, `elm_matched_ions_pitch3.toml`, `elm_matched_sipm_smoke.toml`,
`elm_pitch4mm_bb0nu_tau*_pilot100.toml` (2), `ion_csa_0p5_enc_*.toml` (7),
`ion_csa_bb0nu_3mm_production.toml`, `ion_csa_nonlinear_matched_pitch2_full.toml`,
`ion_csa_v3_*.toml` (3), `sipm_bb0nu_ep2p5.toml`, `sipm_bb0nu_ep3.toml`,
`sipm_bb0nu_statistical_*_pilot100.toml` (8).

Method: a scratch script matched each card's file name against the text of
`src/`, `scripts/`, `tests/`, `Makefile`, `config/`, `docs/` and `notes/`.
Cards given only on a command line or in a shell history are not seen.

## 6. Outcome (2026-11-15, branch `new-config`)

The proposal of Section 4 was replaced by a scheme designed from scratch with
the user: a detector set of eight cards under `config/detector/` (gas,
geometry, drift, voltages, el, sipm, ions, energy, with the index
`itaca.toml`), response-model cards unchanged under `config/specs/response/`,
and run cards under `config/runs/` that name the detector set, list their
`[overrides]` explicitly and carry the run parameters. The framed production
chain (`run_framed_detector_response.jl`, `run_hit_production.jl`,
`run_framed_electronics.jl`, `run_framed_reconstruction.jl`) reads run cards
only and writes `resolved.toml` with the hashes of every input beside its
results. `config/itaca.toml` (the single card of the morning) is superseded and
removed. The scheme is documented in `notes/html/detector_configuration.html`.

Check: the run card `config/runs/bb0nu_vh1_legacy_conditions_pilot100.toml`,
the default detector plus nine overrides, reproduced the detector response of
the current bb0nu production on 145 events and the VH1 hits of its first 100
accepted events bit for bit, before and after the scripts were cleaned of the
old card vocabulary. The cards under `config/readout/` and
`config/specs/runs/` are not edited and no longer drive any script; their
productions keep the card copies saved in their HDF5 outputs. Retiring them
remains a separate decision.
