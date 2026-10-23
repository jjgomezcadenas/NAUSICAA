#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
mkdir -p work/figure-build
latexmk -pdf -interaction=nonstopmode -halt-on-error \
    -outdir=work/figure-build figures/nausica_concept.tex \
    > work/figure-build/build.log 2>&1
cp work/figure-build/nausica_concept.pdf figures/nausica_concept.pdf
pdftoppm -scale-to 2000 -png -singlefile \
    figures/nausica_concept.pdf figures/nausica_concept
