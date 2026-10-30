# Detector and electronics response pipeline plan

## Purpose and boundary

Produce matched, sparse ion-pad and photon-SiPM histories from each MC sample,
then simulate the two electronics branches independently. A top-level run card
selects the detector outputs and response models. Every output records the
configuration and upstream identities that produced it.

The detector stage ends with **ions per pad and time bucket** and **VUV photons
reaching each SiPM per time bucket**. The response stage converts incident
photons to photoelectrons using the SiPM model and applies the electronics to
both branches. The electron branch supplies the energy measurement; the ion
branch has its own electronics and reconstruction observables.

## Shared event realisation

- Read source hits in their stored mm, MeV and ns units. For each source hit,
  sample one ionisation-pair count using its deposited energy, W and the Fano
  prescription. That count creates equal numbers of initial ions and electrons.
- Sample one event location. The ion and electron drift lengths sum to the
  configured TPC length. Store event IDs, lengths, seeds, source hash, exact
  input card and sampler/code identities in a small `realisation.h5`.
- Reproduce the pair count for a hit from a seed keyed by the realisation,
  event ID and source row. Ion and photon branches use the same draw. Record a
  per-event pair-count checksum in each output and verify equality when joining
  branches.
- Source hit positions remain in the reduced hit HDF5. Sample carrier
  diffusion directly into physical grid cells and time buckets. Do not save
  individual ion, electron or photon positions or per-hit pair-count arrays.
- The realisation identity depends on source content, selected event IDs,
  W/F settings, event placement, seed and sampler version. Grid and response
  parameters do not enter that identity. If a future gas model changes W or F,
  it requires its own realisation identity.

## Stage 1: detector transport and collection

### Ion branch

The ion detector card supplies the ion transport and pad-grid settings,
including pitch and time-bucket width. Replay the shared pair counts, diffuse
ions, and save a versioned `ions.h5` with sparse `(event, pad_x, pad_y,
time_bucket, ion_count)` histories, axes, charge accounting and provenance.
Electron-grid, EL and SiPM-response parameters cannot invalidate this file.

### Electron and photon branches

One electron-detector TOML contains shared EL-grid geometry, E/P, electron
time-bucket width and optical boundaries, plus separately labelled pure-Xe and
Xe--15%He diffusion and EL-yield prescriptions. One invocation produces two
versioned files, `photons_pure_xe.h5` and `photons_xe_he15.h5`. Each saves
sparse `(event, sipm_x, sipm_y, time_bucket, photon_count)` histories and its
gas-model identity. Both use the same source events and initial pair draws;
they are alternative gas-model outputs.

Electron attachment uses the configured lifetime and each electron's drift
time. It changes the electron branch only. EL production uses the
field-dependent mean yield and J; direct optical collection uses the analytic
acceptance and its sampling fluctuation. Generate counts without materialising
individual photons. No PDE or electronics response enters these files. SiPM
active area and optical alignment are detector-stage geometry; SiPM PDE, gain
and sensor noise belong to stage 2.

A detector-configuration TOML selects one ion detector card and the electron
detector card that produces the two photon files. It also identifies the shared
realisation. The detector stage may run the three branches concurrently.

## Stage 2: sensor and electronics response

A final-run TOML selects a detector configuration and independent response
cards:

- `sipm.toml`: photon-to-photoelectron conversion and SiPM sensor behaviour;
- `ion_csa.toml`: ion-pad electronics;
- `sipm_electronics.toml`: electronics following the SiPM.

The ion response consumes `ions.h5`. The SiPM and its electronics consume each
selected photon file and produce separate pure-Xe and Xe--He response outputs.
Changing the CSA reruns only the ion response. Changing PDE, SiPM gain or
SiPM electronics reruns only the affected electron response. A combined
reconstruction or comparison reads the two response branches and preserves
their shared primary-ionisation correlation.

## Specifications and generated directories

```text
config/specs/
  realisations/
    bb0nu_8000_seed2026.toml
    single_electrons_8000_seed2026.toml
  detector/
    ions/pitch3mm_bucket30ms.toml
    electrons/elcc_pitch10mm_bucket2us.toml  # pure Xe and Xe--15%He sections
    configurations/bb0nu_ion3mm_elcc.toml
  response/
    sipm/sipm_v1.toml
    ion_csa/csa_v1.toml
    sipm_electronics/readout_v1.toml
  runs/
    bb0nu_ion3mm_elcc_response_v1.toml

results/simulations/
  bb0nu/<realisation_id>/
    realisation.h5
    detector/
      ions/<ion_spec_name>__<hash>/ions.h5
      photons/<electron_spec_name>__<hash>/
        photons_pure_xe.h5
        photons_xe_he15.h5
    response/
      ions/<ion_output_id>/<csa_spec_id>/ion_response.h5
      electrons/<photon_output_id>/<sipm_response_id>/electron_response.h5
    runs/<final_graph_id>/
      manifest.toml
      reconstruction.h5
  single_electrons/<realisation_id>/
    ...
```

These are proposed paths and schemas; the implementation must not mistake them
for existing files. Directory labels support browsing. The resolver derives
identities from canonical parsed TOML values, source content and model version,
locates existing HDF5 files, and checks their saved provenance. TOML comments
and key order do not affect identities. A generated list command presents the
named specs and available outputs, so selection never requires opening result
manifests. Generated manifests snapshot the exact inputs and paths used.

## Dependency and reuse rules

| Change | Detector files to regenerate | Response files to regenerate |
| --- | --- | --- |
| Ion pitch or bucket | Ion file | Ion response |
| Electron grid or E/P | Photon files affected by that setting | Their SiPM responses |
| Ion CSA | None | Ion response |
| SiPM PDE, gain or electronics | None | SiPM responses |
| Shared W/F, source events or placement | Ion and both photon files | Both response branches |

Every node's cache identity includes only its actual physical inputs, upstream
identity and implementation version. Shared conditions such as pressure enter
both branch identities when they affect both branches. Outputs are written
atomically and marked complete only after all required events and ledgers are
saved. The final graph references branch files instead of copying them.

## Implementation and validation sequence

1. Define the TOML contracts, canonical identity calculation, result resolver
   and new HDF5 schemas. Preserve existing PE/readout files for audit.
2. Implement deterministic shared Fano draws and event placement. Verify that
   ion and both photon branches obtain identical per-event pair-count checksums
   regardless of worker count or grid choices.
3. Regenerate one ion baseline under the new realisation; existing ion files
   used deterministic counts and cannot be paired with Fano-fluctuated photons.
4. Implement electron attachment, gas-specific EL fluctuations, analytic
   optical collection and sparse incident-photon histories. Validate means,
   variances, charge ledgers, seeds and time-bucket assignment with independent
   numerical references.
5. Implement independent SiPM, ion-CSA and SiPM-electronics consumers. Verify
   that changing one card invalidates only the nodes shown in the dependency
   table and that mismatched realisations are rejected.
6. Benchmark a small event subset before the full samples. Save numerical
   validation in HDF5; use Python to present it in a generated LaTeX report
   that compiles to PDF.

## Physical inputs still to settle

Oliveira et al. give J versus E/P for pure Xe in a 10 bar planar EL region.
The sampling prescription and J(E/P) for Xe--15%He at our 15 bar hole geometry
need an explicit, sourced model or a recorded approximation before production
runs. The Fano-compatible integer pair sampler must also be selected and
validated for low-energy source hits.
