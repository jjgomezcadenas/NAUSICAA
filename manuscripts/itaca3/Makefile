# Run from manuscripts/itaca3. Override BUILD_DIR for a standalone checkout.
BUILD_DIR ?= ../../work/itaca3_build

.PHONY: all clean
all:
	latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=$(BUILD_DIR) itaca3.tex
	cp $(BUILD_DIR)/itaca3.pdf itaca3.pdf

clean:
	latexmk -c -outdir=$(BUILD_DIR) itaca3.tex
