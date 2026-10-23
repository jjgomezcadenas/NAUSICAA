PYTHON ?= python3
export PYTHONPATH := $(CURDIR)/src
export MPLCONFIGDIR := $(CURDIR)/work/matplotlib-cache
export XDG_CACHE_HOME := $(CURDIR)/work/cache

.PHONY: all physics readout test note check
all: test readout note

physics:
	$(PYTHON) scripts/build_physics.py

readout:
	$(PYTHON) scripts/build_readout.py

test:
	$(PYTHON) -m unittest discover -s tests -v

note: physics
	cd notes/physics && latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build nausica.tex

check: test readout note
