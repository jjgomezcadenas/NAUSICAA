# Production provenance fixes — 2026-11-06

Follow-up to the independent [repository review](repository_review_2026-11-05.html).
Scope: W-01, W-10 and R-01. The original review remains the dated baseline.

| Finding | Implemented change | Verification |
|---|---|---|
| W-01 | `CSAResources.verify` checks the actual 1200 µm production deck, model and linear waveform against acquisition hashes. New model/run identities and response metadata include these hashes. | Mutating each resource independently rejects reuse; a fresh three-event production completes acquisition, correction and publication. |
| W-10 | Matched and SiPM reduction read saved cards and shard provenance. Presentation uses the saved realization card and electron branch. Missing workflow timing is reported as unavailable. | Streaming reduction succeeds with the original card path unavailable and reproduces the saved numerical observables. |
| R-01 | `primary_sampling.jl` owns pair and drift draws; `primary_identity.jl` owns seed and identity derivation. Both arms check `primary_code_hash`. Schema `nausicaa.event_realisation.v2` includes that implementation in its directory identity. | Collection-code edits leave primary identity stable; sampler edits change it. Both arms reject mismatches. All 8,000 matched events replay with identical pair counts/checksums and complementary drift distances. |

## Existing productions and new runs

Historical ledgers and response arrays stay in their original directories. The
new common-ledger writer replays source hits and verifies any same-card legacy
ledger before publishing the v2 file. RNG stream derivation, Float32 sampling,
Fano prescription and placement are unchanged. Runtime and distribution-package
versions are part of the sampler identity. No individual carrier positions are saved.

Create or reuse the new shared ledger before either collection arm:

```sh
julia --project=. scripts/production/run_common_realisation.jl config/specs/runs/elm_matched_ions_pitch2.toml
```

`run_matched_ion_pitch.jl` performs this step automatically and checks the new ion
histories against the completed 4 mm reference. The completed 3 mm production keeps its original
provenance and requires no numerical rerun.

The CSA analysis file pins the waveform hash (`raw_sha256`). The waveform records
the acquired production deck hash and transistor-model hash. Validation follows
that chain rather than recording current resource hashes without an independent
reference. Historical saved-study registration remains supported after resource
validation; new productions store the explicit resource manifest.

## Checks and reproduction

```sh
make production-provenance-test
make ion-csa-test
make test physics PYTHON=work/event-viewer-venv/bin/python
julia --project=. tests/event_analysis/runtests.jl
```

The streaming integration test uses the saved local ELM pilot dataset. Generated
simulation HDF5 files and scratch logs remain outside Git. The matched report is
regenerated from existing physical histories: all 162 numerical datasets match
their pre-change hashes exactly. The HTML/TeX/PDF report also remains identical.
Machine-readable evidence is in `results/provenance_validation/`.

R-02's low-mean integer-sampling law and the separate detector-window/selection
questions are outside this implementation. The sampler docstring now describes
its actual integer-moment law.

The separate `make note` manuscript target still reports its existing LaTeX
error (`Missing $ inserted` in `notes/physics/nausica.tex`). The generated matched
analysis report builds successfully; manuscript text was left unchanged.
