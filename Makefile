PYTHON ?= python3
JULIA ?= julia
EVENT_JULIA_DEPOT ?= $(CURDIR)/work/julia-depot:$(HOME)/.julia
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

note: physics readout spice-paper
	cd notes/physics && latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build nausica.tex

check: test readout note

.PHONY: ltspice-study
ltspice-study:
	$(PYTHON) scripts/build_ltspice_study.py --simulator "$(LTSPICE)"
	$(PYTHON) scripts/analyze_ltspice_study.py

.PHONY: shaping-study
shaping-study:
	$(PYTHON) scripts/build_shaping_study.py --simulator "$(LTSPICE)"

.PHONY: ion-shaping-study
ion-shaping-study:
	$(PYTHON) scripts/build_shaping_study.py --include-ion-frames --simulator "$(LTSPICE)"

.PHONY: spice-paper
spice-paper: spice-averaging-plot
	$(PYTHON) scripts/build_spice_paper.py

.PHONY: spice-averaging spice-averaging-plot
spice-averaging:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --project=. scripts/build_averaging_scan.jl
	$(PYTHON) scripts/plot_averaging_scan.py

spice-averaging-plot:
	$(PYTHON) scripts/plot_averaging_scan.py
