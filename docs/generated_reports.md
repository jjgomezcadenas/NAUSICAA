# Generated reports

Open [the report catalogue](../latex_generated/index.html). It groups generated
notes by scientific subject, includes a text filter, and links TeX, PDF, HTML,
source notes, run cards, saved results and presentation scripts when available.

## Layout

| Subject | Collections |
|---|---:|
| `archive/` | 1 |
| `detector/` | 1 |
| `electron_readout/` | 9 |
| `ion_readout/` | 52 |
| `literature/` | 1 |
| `production/` | 8 |
| `reconstruction/` | 6 |
| `transport/` | 4 |

Authored sources remain in `docs/` and `notes/physics/`. Numerical results remain
in `results/`. Each old top-level report directory is a relative symlink to its
new location; existing bookmarks, source inputs and recorded run-card paths
continue to work. The production-formulas PDF in `docs/` is also a compatibility
link. Current report names include `csa_design`, `sipm_pilot` and
`sipm_validation`, with aliases for their previous filenames.

The catalogue distinguishes standalone documents by their TeX document class.
A PDF plot is an asset of its report, not a separate note. It includes externally
authored notes whose PDF is delivered under `latex_generated/`.

## Regeneration

```sh
make report-catalogue
make report-clean
make elm-note-present
make elm-report
make production-formulas
```

`report-catalogue` reads `docs/generated_reports.json`, checks every generated
link, and writes `latex_generated/index.html` and `catalogue.json`. Add an entry
to the registry when adding a report family. Its optional configuration summary
is curated; full run settings and event counts remain in the report and original
cards. Presentation and complete-workflow commands are shown where verified.

`report-clean` moves compilation intermediates into `work/latex/`. The active
presentation commands use `scripts/documentation/build_report.py` to place all
new LaTeX intermediates there and copy only the finished PDF to the report
folder. The helper checks unresolved references/citations and reports layout
warnings. For a generated report:

```sh
python3 scripts/documentation/build_report.py latex_generated/ion_readout/csa_design/csa_design.tex
```

For an authored source with repository-relative inputs:

```sh
python3 scripts/documentation/build_report.py docs/elm.tex --cwd . --destination latex_generated/detector/elm/elm.pdf
```

Scientific calculations are separate from these presentation commands. The ELM
worked example has `make elm-note` for calculations plus presentation, and
`make elm-note-present` to reuse its saved HDF5 result. The ELM event pilot has
`make elm-pilot` and `make elm-report`, respectively.

## Migration record

The reorganization moved 81 existing generated directories after verifying
byte-for-byte contents, plus the production-formulas PDF. It preserves scientific
manifests, hashes, run identities and numerical files. The ion-track TeX needed
its relative figure paths updated for the extra directory level; its presenter
now resolves the report directory before calculating relative paths. The
production-formulas bibliography now uses the same repository-root convention
as the other authored model notes.

`archive/grid_readout/` is identified as historical by the documented transition
from histogram-centre collection to direct carrier sampling in
`docs/diffusion_sampling_fix.md`. Other study outputs remain under their
scientific subject; their date alone does not make them historical.

Active ELM, SiPM, combined-response, pixel-charge and Makefile report paths use
the new directories. Older workflow scripts and saved TOML report paths continue
through the compatibility links. Their scientific cards and source manifests
were not rewritten. Remaining older presenters may still compile in place;
`make report-clean` collects those intermediates. Migrating those builders to
the shared helper can be done when that workflow is next edited.

Validation used presentation-only builds for ELM, the ELM pilot, production
formulas, the SiPM model/statistical validation, MC selection, CSA design and
ion tracks. CSA design retains its existing overfull paragraph near its
reproduction/provenance section. Other checked builds had no reported layout or
reference warnings. No physical simulation was repeated for this reorganization.
