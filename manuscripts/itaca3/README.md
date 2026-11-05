# ITACA III manuscript

Edit `itaca3.tex`; build with `make -C manuscripts/itaca3` from the repository root.
The compiled manuscript is `manuscripts/itaca3/itaca3.pdf`. Build intermediates
are written to `work/itaca3_build/` at the repository root.

This directory contains the class, bibliography style, bibliography, active
includes, generated numerical macros and figures required to compile the
manuscript. It also builds independently after copying this directory:

```sh
make BUILD_DIR=work/build
```

## Source and template

The initial manuscript is copied from `notes/physics/nausica.tex` on
2026-11-05. The original manuscript and its assets remain unchanged. Source
paths and SHA-256 hashes of the original manuscript, supplied template ZIP and
copied dependencies are recorded in `provenance.json`. Generated values are
copied snapshots of the existing calculations.

The supplied `Download+the+journal+article+template+package+(December+2024+version).zip`
provides `sn-jnl.cls` and `sn-mathphys-num.bst`. The current
[EPJ C submission guidelines](https://link.springer.com/journal/10052/submission-guidelines)
recommend this Springer Nature template with the `iicol` option. The manuscript
uses `pdflatex,iicol,sn-mathphys-num`: two columns and numbered references.
The template files are copied unchanged.

## Adaptations in this copy

- Replaced the article layout, custom margins and running headers with the
  journal class; converted the title, authors, affiliations and abstract to
  its front-matter commands. The title and author list retain their original wording.
- Retained the active manuscript, appendices, equations, data and citations.
  Omitted the inactive historical material after `\end{document}`.
- Corrected the unmatched math delimiter in the introduction:
  `$\bbonu\ lifetimes in excess of $10^{28}$ years` becomes
  `\bbonu\ lifetimes in excess of $10^{28}$ years`.
- Wrapped long displayed equations without changing their expressions, sized
  the weighting-field drawing to a column, and used full-width floats for the
  paired diffusion plots, NAUSICA concept, single-ion-current plots and Topmetal
  comparison table.
- Kept the editable includes during drafting. The template's single-TeX-file
  submission recommendation can be addressed by expanding these inputs when
  assembling the final submission package.

The corresponding-author designation/email and final submission declarations
remain for the authors. The scientific text is the existing working draft;
the ELM material has not yet been incorporated.
