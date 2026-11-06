PYTHON ?= python3
JULIA ?= julia
# Resolved lazily by simulator targets; an explicit LTSPICE assignment takes precedence.
LTSPICE ?= $(shell PYTHONPATH="$(CURDIR)/src" $(PYTHON) -m ion_signal.spice_runtime --path)
EVENT_JULIA_DEPOT ?= $(CURDIR)/work/julia-depot:$(HOME)/.julia
export PYTHONPATH := $(CURDIR)/src
export MPLCONFIGDIR := $(CURDIR)/work/matplotlib-cache
export XDG_CACHE_HOME := $(CURDIR)/work/cache

.PHONY: all physics readout realistic-readout gain-figure arrival-cases noise-study capacitance-check test note check
all: test readout note

physics:
	$(PYTHON) scripts/benchmarks/build_physics.py

readout:
	$(PYTHON) scripts/benchmarks/build_readout.py

realistic-readout:
	$(PYTHON) scripts/benchmarks/build_realistic_csa.py

gain-figure:
	$(PYTHON) scripts/benchmarks/build_gain_figure.py

arrival-cases:
	$(PYTHON) scripts/benchmarks/build_arrival_cases.py

noise-study: arrival-cases
	$(PYTHON) scripts/benchmarks/build_noise_study.py

capacitance-check:
	$(PYTHON) scripts/benchmarks/build_capacitance_check.py

test:
	$(PYTHON) -m unittest discover -s tests -v

note: physics readout spice-paper
	cd notes/physics && latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build nausica.tex

check: test readout note

.PHONY: sipm-test sipm-pilot sipm-report
PHOTON_PILOT_EVENTS ?= 500
.PHONY: photon-production-pilot photon-production-report
photon-production-pilot:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/benchmark_photon_production.jl $(PHOTON_PILOT_EVENTS)
	$(MAKE) photon-production-report

photon-production-report:
	$(PYTHON) scripts/production/present_photon_production_pilot.py

sipm-test:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_sipm_response.jl
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_sipm_statistical.jl
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_sipm_statistical_pipeline.jl
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_sipm_tpb.jl

.PHONY: sipm-tpb-pilot sipm-tpb-report
sipm-tpb-pilot:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/run_sipm_tpb_pilot.jl 12.5 512
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/run_sipm_tpb_pilot.jl 25.0 512
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/run_sipm_tpb_pilot.jl 12.5 1024
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/reduce_sipm_tpb_pilot.jl
	$(MAKE) sipm-tpb-report

sipm-tpb-report:
	$(PYTHON) scripts/production/present_sipm_tpb_pilot.py

.PHONY: sipm-statistical-validation sipm-statistical-report sipm-statistical-pilot
sipm-statistical-validation: sipm-test
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/benchmark_sipm_statistical.jl
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/benchmark_sipm_statistical_followup.jl
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/reduce_sipm_statistical_comparison.jl
	$(MAKE) sipm-statistical-report

sipm-statistical-report:
	$(PYTHON) scripts/production/present_sipm_statistical.py
	$(PYTHON) scripts/documentation/build_report.py docs/sipm_electronics_model.tex --cwd . --destination latex_generated/electron_readout/sipm_response/sipm_electronics_model.pdf

sipm-statistical-pilot:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/run_sipm_production.jl config/specs/runs/sipm_bb0nu_statistical_pilot.toml

.PHONY: sipm-production-pilot sipm-production-pilot-report
sipm-production-pilot:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/benchmark_sipm_production.jl
	$(MAKE) sipm-production-pilot-report

sipm-production-pilot-report:
	$(PYTHON) scripts/production/present_sipm_production_pilot.py

.PHONY: electron-ep2-pilot sipm-clipping-pilot-analysis sipm-clipping-pilot-report
electron-ep2-pilot:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/run_electron_ep2_pilot.jl
	$(MAKE) sipm-clipping-pilot-analysis
	$(MAKE) sipm-clipping-pilot-report

sipm-clipping-pilot-analysis:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/reduce_sipm_clipping_pilots.jl

sipm-clipping-pilot-report:
	$(PYTHON) scripts/production/present_sipm_clipping_pilots.py
	$(PYTHON) scripts/production/present_sipm_production_pilot.py

sipm-pilot: sipm-test
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/run_sipm_pilot.jl
	$(MAKE) sipm-report

sipm-report:
	$(PYTHON) scripts/production/present_sipm_pilot.py
	$(PYTHON) scripts/documentation/build_report.py latex_generated/electron_readout/sipm_response/sipm_pilot.tex
	$(PYTHON) scripts/documentation/build_report.py docs/sipm_electronics_model.tex --cwd . --destination latex_generated/electron_readout/sipm_response/sipm_electronics_model.pdf

.PHONY: ltspice-study
ltspice-study:
	$(PYTHON) scripts/electronics/reproduction/initial/build_ltspice_study.py --simulator "$(LTSPICE)"
	$(PYTHON) scripts/electronics/reproduction/initial/analyze_ltspice_study.py

.PHONY: shaping-study
shaping-study:
	$(PYTHON) scripts/electronics/reproduction/initial/build_shaping_study.py --simulator "$(LTSPICE)"

.PHONY: ion-shaping-study
ion-shaping-study:
	$(PYTHON) scripts/electronics/reproduction/initial/build_shaping_study.py --include-ion-frames --simulator "$(LTSPICE)"

.PHONY: spice-paper
spice-paper: spice-averaging-plot
	$(PYTHON) scripts/paper/build_spice_paper.py

.PHONY: julia-env test-events extract-hits
julia-env:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --project=. -e 'using Pkg; Pkg.instantiate()'

test-events:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --project=. tests/test_event_io.jl

.PHONY: test-event-core event-review
test-event-core:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/event_analysis/runtests.jl

event-review:
	$(PYTHON) scripts/documentation/build_event_review.py

.PHONY: validate-event-step1 report-event-step1
# Reader/reconstruction checkpoint; input reduced files must already be v3.
validate-event-step1:
	mkdir -p work/event_analysis
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_event_io.jl > work/event_analysis/step1_io.log 2>&1 || { cat work/event_analysis/step1_io.log; exit 1; }
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/event_analysis/runtests.jl > work/event_analysis/step1_core.log 2>&1 || { cat work/event_analysis/step1_core.log; exit 1; }
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_event_analysis.jl > work/event_analysis/step1_legacy.log 2>&1 || { cat work/event_analysis/step1_legacy.log; exit 1; }
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/validation/validate_event_step1.jl
	$(MAKE) report-event-step1
	$(MAKE) event-review

report-event-step1:
	$(PYTHON) scripts/validation/report_event_step1.py
	$(PYTHON) scripts/documentation/build_report.py latex_generated/reconstruction/event_analysis_step1/report.tex

.PHONY: test-event-diffusion validate-event-step2 report-event-step2
test-event-diffusion:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/event_analysis/runtests.jl --step2

validate-event-step2:
	mkdir -p work/event_analysis
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/event_analysis/runtests.jl --step2 > work/event_analysis/step2_tests.log 2>&1 || { cat work/event_analysis/step2_tests.log; exit 1; }
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/validation/validate_event_step2.jl
	$(MAKE) report-event-step2
	$(MAKE) event-review

report-event-step2:
	$(PYTHON) scripts/validation/report_event_step2.py
	$(PYTHON) scripts/documentation/build_report.py latex_generated/reconstruction/event_analysis_step2/report.tex

.PHONY: test-event-collection validate-event-step3 report-event-step3
test-event-collection:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/event_analysis/runtests.jl --step3

validate-event-step3:
	mkdir -p work/event_analysis
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/event_analysis/runtests.jl --step3 > work/event_analysis/step3_tests.log 2>&1 || { cat work/event_analysis/step3_tests.log; exit 1; }
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/validation/validate_event_step3.jl
	$(MAKE) report-event-step3
	$(MAKE) event-review

report-event-step3:
	$(PYTHON) scripts/validation/report_event_step3.py
	$(PYTHON) scripts/documentation/build_report.py latex_generated/reconstruction/event_analysis_step3/report.tex

.PHONY: test-event-ensemble validate-event-step4 report-event-step4

VIEWER_PYTHON ?= $(if $(wildcard work/event-viewer-venv/bin/python),work/event-viewer-venv/bin/python,$(PYTHON))
.PHONY: view-events check-event-viewer
view-events:
	$(VIEWER_PYTHON) scripts/viewers/view_events.py

check-event-viewer:
	$(VIEWER_PYTHON) tests/event_viewer/check_session.py

test-event-ensemble:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/event_analysis/runtests.jl --step4

validate-event-step4:
	mkdir -p work/event_analysis
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/event_analysis/runtests.jl --step4 > work/event_analysis/step4_tests.log 2>&1 || { cat work/event_analysis/step4_tests.log; exit 1; }
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/validation/scan_event_analysis.jl
	$(MAKE) report-event-step4
	$(MAKE) event-review

report-event-step4:
	$(PYTHON) scripts/validation/report_event_step4.py
	$(PYTHON) scripts/documentation/build_report.py latex_generated/reconstruction/event_analysis_step4/report.tex

RAW_EVENTS_DIR ?= $(HOME)/NAUSICAA
extract-hits:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --project=. scripts/events/extract_hits.jl "$(RAW_EVENTS_DIR)/0nubb.next.h5" "$(RAW_EVENTS_DIR)/electron_0nubb_energy.next.h5"

.PHONY: spice-averaging spice-averaging-plot
spice-averaging:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --project=. scripts/electronics/readout/build_averaging_scan.jl
	$(PYTHON) scripts/electronics/readout/plot_averaging_scan.py

spice-averaging-plot:
	$(PYTHON) scripts/electronics/readout/plot_averaging_scan.py

.PHONY: test-ion-analysis ion-analysis ion-report ion-tracks
test-ion-analysis:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --project=. tests/test_event_analysis.jl

ion-analysis:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --project=. scripts/events/analyze_hits.jl

ion-report:
	$(PYTHON) scripts/events/plot_hit_analysis.py
	mkdir -p work/ion_tracks_latex
	$(PYTHON) scripts/documentation/build_report.py latex_generated/ion_readout/ion_tracks/report.tex

ion-tracks: ion-analysis
	$(MAKE) ion-report

.PHONY: spice2-reproduction spice2-plots
spice2-reproduction:
	$(PYTHON) scripts/electronics/reproduction/spice2/reproduce_spice2.py --simulator "$(LTSPICE)"
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --project=. scripts/electronics/reproduction/spice2/analyze_spice2.jl
	$(PYTHON) scripts/electronics/reproduction/spice2/plot_spice2.py

spice2-plots:
	$(PYTHON) scripts/electronics/reproduction/spice2/plot_spice2.py

.PHONY: spice2-enc spice2-enc-plots
spice2-enc:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(PYTHON) scripts/electronics/reproduction/spice2/run_spice2_enc.py --simulator "$(LTSPICE)" --julia "$(JULIA)"
	$(PYTHON) scripts/electronics/reproduction/spice2/plot_spice2_enc.py

spice2-enc-plots:
	$(PYTHON) scripts/electronics/reproduction/spice2/plot_spice2_enc.py

.PHONY: ramo-note ramo-note-render
ramo-note:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --project=. scripts/paper/check_ramo_note.jl
	$(PYTHON) scripts/paper/write_ramo_note.py

ramo-note-render:
	$(PYTHON) scripts/paper/write_ramo_note.py

.PHONY: mc-selection present-mc-selection
mc-selection:
	mkdir -p work/mc_selection
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/events/analyze_mc_selection.jl
	$(MAKE) present-mc-selection

present-mc-selection:
	MPLCONFIGDIR="$(CURDIR)/work/matplotlib" $(PYTHON) scripts/events/present_mc_selection.py

.PHONY: pixel-charge present-pixel-charge
pixel-charge:
	mkdir -p work/pixel_charge
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/collection/run_pixel_charge_direct.jl config/pixel_charge.toml
	$(MAKE) present-pixel-charge

present-pixel-charge:
	MPLCONFIGDIR="$(CURDIR)/work/matplotlib" $(PYTHON) scripts/collection/present_pixel_charge.py --results-dir results/pixel_charge_direct --latex-dir latex_generated/ion_readout/pixel_charge_direct

.PHONY: pixel-charge-1mm
pixel-charge-1mm:
	mkdir -p work/pixel_charge_1mm
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/collection/run_pixel_charge_direct.jl config/pixel_charge_1mm.toml
	MPLCONFIGDIR="$(CURDIR)/work/matplotlib" $(PYTHON) scripts/collection/present_pixel_charge.py --results-dir results/pixel_charge_1mm_direct --latex-dir latex_generated/ion_readout/pixel_charge_1mm_direct

.PHONY: pixel-charge-1mm-3mm
pixel-charge-1mm-3mm:
	mkdir -p work/pixel_charge_1mm_3mm
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/collection/run_pixel_charge_direct.jl config/pixel_charge_1mm_3mm.toml
	MPLCONFIGDIR="$(CURDIR)/work/matplotlib" $(PYTHON) scripts/collection/present_pixel_charge.py --results-dir results/pixel_charge_1mm_3mm_direct --latex-dir latex_generated/ion_readout/pixel_charge_1mm_3mm_direct

.PHONY: pixel-charge-1mm-4mm
pixel-charge-1mm-4mm:
	mkdir -p work/pixel_charge_1mm_4mm
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/collection/run_pixel_charge_direct.jl config/pixel_charge_1mm_4mm.toml
	MPLCONFIGDIR="$(CURDIR)/work/matplotlib" $(PYTHON) scripts/collection/present_pixel_charge.py --results-dir results/pixel_charge_1mm_4mm_direct --latex-dir latex_generated/ion_readout/pixel_charge_1mm_4mm_direct

.PHONY: grid-readout present-grid-readout
grid-readout:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/event_analysis/test_grid_readout.jl
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/collection/sample_collect_parallel.jl --events 8000 --workers 6 --output results/grid_readout_direct_8000
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/collection/summarize_grid_readout.jl results/grid_readout_direct_8000
	$(PYTHON) tests/event_analysis/check_grid_readout.py results/grid_readout_direct_8000 --events 8000
	$(MAKE) present-grid-readout

present-grid-readout:
	MPLCONFIGDIR="$(CURDIR)/work/matplotlib" $(PYTHON) scripts/collection/present_grid_readout.py

.PHONY: grid-readout-3mm present-grid-readout-3mm
grid-readout-3mm:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/event_analysis/test_grid_readout.jl
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/collection/sample_collect_parallel.jl --events 8000 --workers 6 --ion-card config/pixel_charge_1mm_3mm.toml --output results/grid_readout_direct_3mm_8000
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/collection/summarize_grid_readout.jl results/grid_readout_direct_3mm_8000
	$(PYTHON) tests/event_analysis/check_grid_readout.py results/grid_readout_direct_3mm_8000 --events 8000 --ion-pitch-mm 3 --ion-bucket-ms 30
	$(MAKE) present-grid-readout-3mm

present-grid-readout-3mm:
	MPLCONFIGDIR="$(CURDIR)/work/matplotlib" $(PYTHON) scripts/collection/present_grid_readout.py --input results/grid_readout_direct_3mm_8000/grid_readout_summary.h5 --output latex_generated/production/grid_readout_direct_3mm
	$(MAKE) pdf-grid-readout-3mm

.PHONY: pdf-grid-readout pdf-grid-readout-3mm pdf-electron-diffusion pdf-mc-selection
pdf-grid-readout:
	$(PYTHON) scripts/documentation/build_report.py latex_generated/production/grid_readout_direct/grid_readout.tex

pdf-grid-readout-3mm:
	$(PYTHON) scripts/documentation/build_report.py latex_generated/production/grid_readout_direct_3mm/grid_readout.tex

pdf-electron-diffusion:
	$(PYTHON) scripts/documentation/build_report.py latex_generated/transport/electron_diffusion/electron_diffusion.tex

pdf-mc-selection:
	$(PYTHON) scripts/documentation/build_report.py latex_generated/reconstruction/mc_selection/mc_selection.tex

.PHONY: test-script-layout
test-script-layout:
	$(VIEWER_PYTHON) tests/script_layout/check_migration.py

.PHONY: ion-csa-test ion-csa-pilot csa-design-note
ion-csa-test:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_ion_csa_nonlinear.jl
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_ion_csa_response.jl
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_ion_csa_correction.jl
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_ion_csa_pipeline.jl

ion-csa-pilot:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/run_ion_csa.jl config/specs/runs/ion_csa_v3_pilot_all.toml

.PHONY: ion-csa-nonlinear-calibrate ion-csa-nonlinear-pilot ion-csa-nonlinear-report
ion-csa-nonlinear-calibrate:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/calibrate_ion_csa_nonlinear.jl

ion-csa-nonlinear-pilot:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/run_ion_csa.jl config/specs/runs/ion_csa_nonlinear_pilot.toml
	$(MAKE) ion-csa-nonlinear-report

ion-csa-nonlinear-report:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/reduce_ion_csa_nonlinear_report.jl
	$(VIEWER_PYTHON) scripts/production/present_ion_csa_nonlinear.py

.PHONY: ion-csa-nonlinear-production-100
ion-csa-nonlinear-production-100:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/run_ion_csa.jl config/specs/runs/ion_csa_nonlinear_production_100.toml
	$(VIEWER_PYTHON) scripts/production/present_ion_csa_production_timing.py

csa-design-note:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/reduce_csa_design.jl
	$(VIEWER_PYTHON) scripts/production/present_csa_design.py
	$(PYTHON) scripts/documentation/build_report.py latex_generated/ion_readout/csa_design/csa_design.tex

.PHONY: response-overview response-overview-present response-viewer response-overview-test
response-overview:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/reduce_response_overview.jl
	$(MAKE) response-overview-present

response-overview-present:
	$(VIEWER_PYTHON) scripts/production/present_response_overview.py

RESPONSE_VIEWER_CARD ?= config/specs/runs/elm_matched_response_viewer.toml
response-viewer:
	$(VIEWER_PYTHON) scripts/viewers/serve_cloud_viewer.py --production $(RESPONSE_VIEWER_CARD) --port 8769

response-overview-test:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_response_overview.jl
	$(VIEWER_PYTHON) tests/test_production_viewer.py

.PHONY: electron-operating-points electron-point-2p5 electron-point-3 electron-operating-points-present
electron-operating-points:
	$(MAKE) -j2 electron-point-2p5 electron-point-3
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/reduce_electron_operating_points.jl
	$(MAKE) electron-operating-points-present

electron-point-2p5:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/run_electron_operating_point.jl ep2p5

electron-point-3:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/run_electron_operating_point.jl ep3

electron-operating-points-present:
	MPLCONFIGDIR=work/matplotlib $(VIEWER_PYTHON) scripts/production/present_electron_operating_points.py

.PHONY: elm-test elm-pilot elm-report
elm-test:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_elm_pilot.jl

elm-pilot:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/run_elm_pilot.jl
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/reduce_elm_pilot.jl
	$(MAKE) elm-report

elm-report:
	$(PYTHON) scripts/production/present_elm_pilot.py

# Generated-document catalogue and presentation-only note builds.
.PHONY: report-catalogue report-clean elm-note elm-note-present production-formulas
report-catalogue:
	$(PYTHON) scripts/documentation/catalogue_reports.py

report-clean:
	$(PYTHON) scripts/documentation/catalogue_reports.py --clean-intermediates

elm-note:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/documentation/calculate_elm_example.jl
	$(MAKE) elm-note-present

elm-note-present:
	$(VIEWER_PYTHON) scripts/documentation/present_elm_example.py
	$(PYTHON) scripts/documentation/build_report.py docs/elm.tex --cwd . --destination latex_generated/detector/elm/elm.pdf
	$(MAKE) report-catalogue

production-formulas:
	$(PYTHON) scripts/documentation/build_report.py docs/production_formulas.tex --cwd . --destination latex_generated/production/formulas/production_formulas.pdf
	$(MAKE) report-catalogue

.PHONY: elm-sipm-test elm-sipm-pilot elm-sipm-report
elm-sipm-test:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_elm_sipm_pilot.jl
elm-sipm-pilot:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/run_elm_sipm_pilot.jl
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/reduce_elm_sipm_pilot.jl
	$(MAKE) elm-sipm-report
elm-sipm-report:
	$(PYTHON) scripts/production/present_elm_sipm_pilot.py

.PHONY: elm-pitch4-test elm-pitch4-pilot elm-pitch4-report
elm-pitch4-test:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_elm_pitch4mm_pilot.jl
elm-pitch4-pilot:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/run_elm_pitch4mm_pilot.jl
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/reduce_elm_pitch4mm_pilot.jl
	$(MAKE) elm-pitch4-report
elm-pitch4-report:
	$(PYTHON) scripts/production/present_elm_pitch4mm_pilot.py

.PHONY: elm-pitch4-production elm-pitch4-production-report elm-streaming-test
elm-pitch4-production:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/run_elm_pitch4mm_production.jl
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/reduce_elm_pitch4mm_production.jl
	$(MAKE) elm-pitch4-production-report
elm-pitch4-production-report:
	$(PYTHON) scripts/production/present_elm_pitch4mm_production.py
elm-streaming-test:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_elm_streaming_reduction.jl

.PHONY: elm-matched-production elm-matched-analysis elm-matched-report elm-matched-test
elm-matched-production:
	$(PYTHON) scripts/production/run_elm_matched_production.py --kind full --julia "$(JULIA)"
	$(MAKE) elm-matched-analysis
elm-matched-analysis:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. scripts/production/reduce_matched_analysis.jl
	$(MAKE) elm-matched-report
elm-matched-report:
	$(PYTHON) scripts/production/present_matched_analysis.py
elm-matched-test:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_matched_realisation.jl

.PHONY: rogelio-notes
rogelio-notes:
	$(PYTHON) scripts/documentation/build_report.py LTSpice/notes/calculos_para_rogelio.tex --destination latex_generated/ion_readout/rogelio_notes/calculos_para_rogelio.pdf
	$(PYTHON) scripts/documentation/build_report.py LTSpice/notes/ramo_para_rogelio.tex --destination latex_generated/ion_readout/rogelio_notes/ramo_para_rogelio.pdf
	$(MAKE) report-catalogue

.PHONY: production-provenance-test
production-provenance-test:
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_matched_realisation.jl
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_primary_ledger.jl
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_ion_csa_pipeline.jl
	JULIA_DEPOT_PATH="$(EVENT_JULIA_DEPOT)" $(JULIA) --startup-file=no --project=. tests/test_elm_streaming_reduction.jl
	$(VIEWER_PYTHON) tests/test_matched_report_provenance.py
