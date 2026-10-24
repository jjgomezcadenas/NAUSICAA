PYTHON ?= python3
export PYTHONPATH := $(CURDIR)/src
export MPLCONFIGDIR := $(CURDIR)/work/matplotlib-cache
export XDG_CACHE_HOME := $(CURDIR)/work/cache

.PHONY: all physics readout realistic-readout gain-figure arrival-cases noise-study capacitance-check test note check
all: test readout note

physics:
	$(PYTHON) scripts/build_physics.py

readout:
	$(PYTHON) scripts/build_readout.py

realistic-readout:
	$(PYTHON) scripts/build_realistic_csa.py

gain-figure:
	$(PYTHON) scripts/build_gain_figure.py

arrival-cases:
	$(PYTHON) scripts/build_arrival_cases.py

noise-study: arrival-cases
	$(PYTHON) scripts/build_noise_study.py

capacitance-check:
	$(PYTHON) scripts/build_capacitance_check.py

test:
	$(PYTHON) -m unittest discover -s tests -v

note: physics readout gain-figure noise-study
	cd notes/physics && latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build nausica.tex

check: test readout note
