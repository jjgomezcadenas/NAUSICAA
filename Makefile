PYTHON ?= python3
export PYTHONPATH := $(CURDIR)/src
export MPLCONFIGDIR := $(CURDIR)/work/matplotlib-cache
export XDG_CACHE_HOME := $(CURDIR)/work/cache

.PHONY: all physics test note check
all: test note

physics:
	$(PYTHON) scripts/build_physics.py

test:
	$(PYTHON) -m unittest discover -s tests -v

note: physics
	cd notes/physics && latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build main.tex

check: test note
