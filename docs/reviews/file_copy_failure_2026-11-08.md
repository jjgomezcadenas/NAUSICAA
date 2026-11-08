# Publication-copy failure and guards — 2026-11-08

## Finding

The 2 mm physics calculation and reconstruction completed and passed validation.
Publication then called `Base.cp(parent, temporary; force=true)` on an
8,793,451,765-byte corrected HDF5. The copy grew to 341,698,407,669 bytes and
failed with `ENOSPC`. The temporary file has been deleted after verification of
the corrected source against its saved SHA-256. No physical histories were lost.

The installed Julia 1.12.6 implementation uses a 32-bit signed return value for
`jl_fs_sendfile`; its copy loop uses that value as the bytes transferred. The
[C wrapper](https://github.com/JuliaLang/julia/blob/v1.12.6/src/jl_uv.c#L545)
and [Julia loop](https://github.com/JuliaLang/julia/blob/v1.12.6/base/filesystem.jl#L218)
show this path. For transfers above the counter's range, the reported count wraps.
The loop also has no zero-progress exit.

A controlled test on this machine made only two individual calls, with an OS
file-size ceiling of 110% of the source. The scratch files were removed afterward:

| Call | Bytes actually written | Bytes reported |
|---|---:|---:|
| First: 4 GiB + 64 KiB | 4,295,032,832 | 65,536 |
| Second: 4 GiB | 4,294,967,296 | 0 |

The recorded results are in
`results/copy_guard_validation/counter_probe.toml`; the bounded probe source and
its Python launcher invocation are retained in the task's `work/copy_guard/`
scratch record/tool history. The unbounded `Base.cp` loop was never run in this
experiment.

For the actual 2 mm file, the first reported count is 203,517,173 bytes. The
apparent remainder is exactly 8,589,934,592 bytes (8 GiB), whose 32-bit count is
zero. The next loop iteration therefore copies the same remainder again without
advancing its source offset. This accounts for the runaway copy. The earlier
3 mm error is also consistent with a wrapped negative count after a successful
3,236,869,732-byte transfer.

Our workflow compounded the runtime bug: it had neither a copy-size bound nor
failure cleanup. The earlier 3 mm copy error should have prompted this check
before the 2 mm production.

## Implemented protections

`src/ion_signal/guarded_files.jl` provides shared copy/publication operations:

- Copy with reads and writes of at most 8 MiB, independent of `sendfile`.
- Require every write to make the expected progress; stop exactly at the
  source's recorded length. Reject a source that grows, shrinks or changes.
- Verify exact byte length and SHA-256 before any metadata edit.
- Allow at most 10% total size growth during publication. Check after every
  flushed metadata update and again after validation, before publication.
- Check available disk space before copying and during each block, keeping a
  1 GiB reserve. The initial budget includes the metadata-growth allowance.
- Write into a uniquely owned temporary beside the target and publish by rename.
  Delete that temporary on exceptions or interruption; preserve the source and
  any pre-existing destination. A hard process kill can leave a temporary, but
  cannot revive the unbounded copy loop.

The CSA publisher, photon archival/pilot workflows and production-index copies
use the helper. CSA source snapshots use it as well. Historical recovery now
checks a temporary's length before hashing and uses the shared guarded publisher.

An exact copy permits **zero** extra bytes. The 10% allowance applies only to
intentional publication metadata. For genuinely generated datasets, an output
budget must come from the dataset dimensions and storage representation: their
source files can legitimately be much smaller.

## Verification

`tests/test_guarded_files.jl` checks ordinary and empty copies, checksum failure,
insufficient disk budget, self-copy rejection, oversized writes, validation
failure, cancellation cleanup, preservation of an existing target, and the
inclusive 10% boundary. These tests are included in `make ion-csa-test`.

The complete CSA test suite (including 20 guard checks) and all 73 Python tests
plus physics checks pass. Guarded publication of the actual 2 mm response
completed successfully: 8,793,451,765 source bytes became 8,793,456,974 published
bytes, an increase of 5,209 bytes for metadata (0.00005924%). All 8,000 events and
897,306,359 ions passed validation. Publication/registration took 200.21 s,
recorded separately from the original 1,063.56 s numerical calculation.

Machine-readable evidence: `results/copy_guard_validation/publication.toml`.
The final [2 mm HTML report](../../latex_generated/ion_readout/ion_csa_nonlinear_pitch2/index.html)
and its PDF were generated successfully. No simulation rerun was required.

The separate manuscript `make note` target still reports its pre-existing LaTeX
build error. The generated 2 mm report compiled successfully; manuscript text was
not changed.
