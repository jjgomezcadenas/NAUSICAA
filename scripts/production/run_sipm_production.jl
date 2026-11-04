include("sipm_production.jl")
card=isempty(ARGS) ? "config/specs/runs/sipm_bb0nu_production.toml" : only(ARGS)
SiPMProduction.run(card)
