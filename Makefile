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

.PHONY: julia-env test-events extract-hits
julia-env:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --project=. -e 'using Pkg; Pkg.instantiate()'

test-events:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --project=. tests/test_event_io.jl

.PHONY: test-event-core event-review
test-event-core:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/event_analysis/runtests.jl

event-review:
	$(PYTHON) scripts/build_event_review.py

.PHONY: validate-event-step1 report-event-step1
# Reader/reconstruction checkpoint; input reduced files must already be v3.
validate-event-step1:
	mkdir -p work/event_analysis
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_event_io.jl > work/event_analysis/step1_io.log 2>&1 || { cat work/event_analysis/step1_io.log; exit 1; }
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/event_analysis/runtests.jl > work/event_analysis/step1_core.log 2>&1 || { cat work/event_analysis/step1_core.log; exit 1; }
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_event_analysis.jl > work/event_analysis/step1_legacy.log 2>&1 || { cat work/event_analysis/step1_legacy.log; exit 1; }
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/validate_event_step1.jl
	$(MAKE) report-event-step1
	$(MAKE) event-review

report-event-step1:
	$(PYTHON) scripts/report_event_step1.py
	cd latex_generated/event_analysis_step1 && latexmk -pdf -interaction=nonstopmode -halt-on-error report.tex

.PHONY: test-event-diffusion validate-event-step2 report-event-step2
test-event-diffusion:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/event_analysis/runtests.jl --step2

validate-event-step2:
	mkdir -p work/event_analysis
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/event_analysis/runtests.jl --step2 > work/event_analysis/step2_tests.log 2>&1 || { cat work/event_analysis/step2_tests.log; exit 1; }
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/validate_event_step2.jl
	$(MAKE) report-event-step2
	$(MAKE) event-review

report-event-step2:
	$(PYTHON) scripts/report_event_step2.py
	cd latex_generated/event_analysis_step2 && latexmk -pdf -interaction=nonstopmode -halt-on-error report.tex

.PHONY: test-event-collection validate-event-step3 report-event-step3
test-event-collection:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/event_analysis/runtests.jl --step3

validate-event-step3:
	mkdir -p work/event_analysis
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/event_analysis/runtests.jl --step3 > work/event_analysis/step3_tests.log 2>&1 || { cat work/event_analysis/step3_tests.log; exit 1; }
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/validate_event_step3.jl
	$(MAKE) report-event-step3
	$(MAKE) event-review

report-event-step3:
	$(PYTHON) scripts/report_event_step3.py
	cd latex_generated/event_analysis_step3 && latexmk -pdf -interaction=nonstopmode -halt-on-error report.tex

.PHONY: test-event-ensemble validate-event-step4 report-event-step4

VIEWER_PYTHON ?= $(if $(wildcard work/event-viewer-venv/bin/python),work/event-viewer-venv/bin/python,$(PYTHON))
.PHONY: view-events check-event-viewer
view-events:
	$(VIEWER_PYTHON) scripts/view_events.py

check-event-viewer:
	$(VIEWER_PYTHON) tests/event_viewer/check_session.py

test-event-ensemble:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/event_analysis/runtests.jl --step4

validate-event-step4:
	mkdir -p work/event_analysis
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/event_analysis/runtests.jl --step4 > work/event_analysis/step4_tests.log 2>&1 || { cat work/event_analysis/step4_tests.log; exit 1; }
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/scan_event_analysis.jl
	$(MAKE) report-event-step4
	$(MAKE) event-review

report-event-step4:
	$(PYTHON) scripts/report_event_step4.py
	cd latex_generated/event_analysis_step4 && latexmk -pdf -interaction=nonstopmode -halt-on-error report.tex

extract-hits:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --project=. scripts/extract_hits.jl events/0nubb.next.h5 events/electron_0nubb_energy.next.h5

.PHONY: spice-averaging spice-averaging-plot
spice-averaging:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --project=. scripts/build_averaging_scan.jl
	$(PYTHON) scripts/plot_averaging_scan.py

spice-averaging-plot:
	$(PYTHON) scripts/plot_averaging_scan.py

.PHONY: test-ion-analysis ion-analysis ion-report ion-tracks
test-ion-analysis:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --project=. tests/test_event_analysis.jl

ion-analysis:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --project=. scripts/analyze_hits.jl

ion-report:
	$(PYTHON) scripts/plot_hit_analysis.py
	mkdir -p work/ion_tracks_latex
	cd latex_generated/ion_tracks && latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=../../work/ion_tracks_latex report.tex
	cp work/ion_tracks_latex/report.pdf latex_generated/ion_tracks/report.pdf

ion-tracks: ion-analysis
	$(MAKE) ion-report

.PHONY: spice2-reproduction spice2-plots
spice2-reproduction:
	$(PYTHON) scripts/reproduce_spice2.py --simulator "$(LTSPICE)"
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --project=. scripts/analyze_spice2.jl
	$(PYTHON) scripts/plot_spice2.py

spice2-plots:
	$(PYTHON) scripts/plot_spice2.py

.PHONY: spice2-enc spice2-enc-plots
spice2-enc:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(PYTHON) scripts/run_spice2_enc.py --simulator "$(LTSPICE)" --julia "$(JULIA)"
	$(PYTHON) scripts/plot_spice2_enc.py

spice2-enc-plots:
	$(PYTHON) scripts/plot_spice2_enc.py

.PHONY: ramo-note ramo-note-render
ramo-note:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --project=. scripts/check_ramo_note.jl
	$(PYTHON) scripts/write_ramo_note.py

ramo-note-render:
	$(PYTHON) scripts/write_ramo_note.py
