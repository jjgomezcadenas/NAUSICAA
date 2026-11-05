#!/usr/bin/env julia
module MatchedAnalysis
using HDF5,TOML,Statistics,SHA
include("reduce_matched_electrons.jl")
include("../../src/ion_signal/ion_csa_pipeline.jl")
const ROOT=MatchedElectronAnalysis.ROOT
const OUT=joinpath(ROOT,"results/elm_matched_production")
sha(path)=bytes2hex(sha256(read(path)))

"""Reduce both matched arms with the same positive-input-cell pulse convention.

Electron energy diagnostics use the existing validated streaming reduction.
Ion pulse means/maxima preserve signed corrected estimates on the input-positive
cells; total response ratios include all stored quiet/recovery frames. Charge
compression and signed noise remain separate from the electron ASIC clipping.
"""
function main()
    graph=TOML.parsefile(joinpath(OUT,"manifest.toml"))
    times=TOML.parsefile(joinpath(OUT,"run_timings.toml"))
    MatchedElectronAnalysis.main("full")
    target=joinpath(OUT,"matched_analysis.h5")
    h5open(target*".tmp","w") do out
      attributes(out)["schema"]="nausicaa.matched_response_analysis.v1"
      meta=create_group(out,"metadata")
      for (k,v) in pairs((;complete=false,manifest_path=joinpath(OUT,"manifest.toml"),
          manifest_sha256=sha(joinpath(OUT,"manifest.toml")),ion_path=graph["ion_response"],
          electron_analysis_path=joinpath(OUT,"electron_analysis.h5"),
          reduction_sha256=sha(@__FILE__),run_timings_toml=read(joinpath(OUT,"run_timings.toml"),String)))
        meta[string(k)]=v
      end
      columns=Dict{String,Vector{Float64}}()
      h5open(graph["ion_response"],"r") do source
        ids=sort!(parse.(Int64,collect(keys(source["events"]))))
        h5open(joinpath(OUT,"electron_analysis.h5"),"r") do electrons
            @assert ids==Int64.(read(electrons["tau_12p5/xe_he15/event_id"]))
        end
        for (j,id) in enumerate(ids)
            event=IonCSAProduction.read_event(source,id)
            n=event.input_ions;q=Float64.(event.noiseless_ions);r=Float64.(event.reconstructed_ions)
            pos=n.>0;initial=sum(n);totalq=sum(q);totalr=sum(r)
            flags=read(source["events/$id/compression_level"])
            snr=[q[i]/event.noise_sigma_ions[i[1]] for i in findall(pos)]
            values=(;event_id=id,ion_drift_cm=event.ion_drift_cm,electron_drift_cm=event.electron_drift_cm,
                initial_ions=initial,positive_cells=count(pos),pixels=length(event.pixel_ix),
                event_mean_input=mean(n[pos]),event_max_input=maximum(n[pos]),
                event_mean_noiseless=mean(q[pos]),event_max_noiseless=maximum(q[pos]),
                event_mean_reconstructed=mean(r[pos]),event_max_reconstructed=maximum(r[pos]),
                event_mean_snr=mean(snr),total_noiseless=totalq,total_reconstructed=totalr,
                noiseless_ratio=totalq/initial,reconstructed_ratio=totalr/initial,
                compressed50=any(flags.>=3),spice_pixels=count(>(0),event.model_used),
                negative_input_cells=count(<(0),r[pos]))
            for (k,v) in pairs(values);push!(get!(columns,string(k),Float64[]),Float64(v));end
            if j%1000==0;println("Matched ion analysis $j/$(length(ids))");flush(stdout);end
        end
        g=create_group(out,"ions")
        for (k,v) in columns;g[k]=v;end
        for key in ("event_mean_input","event_max_input","event_mean_noiseless","event_max_noiseless",
                    "event_mean_reconstructed","event_max_reconstructed","event_mean_snr",
                    "initial_ions","noiseless_ratio","reconstructed_ratio")
            g["average_"*key]=mean(columns[key]);g["median_"*key]=median(columns[key])
        end
        g["event_count"]=length(ids);g["compression50_event_percent"]=100mean(columns["compressed50"])
        g["noiseless_total_bias_percent"]=100(sum(columns["total_noiseless"])/sum(columns["initial_ions"])-1)
        g["reconstructed_total_bias_percent"]=100(sum(columns["total_reconstructed"])/sum(columns["initial_ions"])-1)
        g["noiseless_ratio_rms_percent"]=100std(columns["noiseless_ratio"];corrected=false)
        g["reconstructed_ratio_rms_percent"]=100std(columns["reconstructed_ratio"];corrected=false)
        g["response_card"]=read(joinpath(ROOT,"config/specs/response/ion_csa/csa_0p5_enc_4mm.toml"),String)
        g["reconstruction_card"]=read(joinpath(ROOT,"config/specs/response/ion_csa/reconstruction_v1.toml"),String)
      end
      corrected=joinpath(dirname(graph["ion_manifest"]),"execution/pitch4/corrected")
      summary=only(filter(isfile,[joinpath(d,"summary.h5") for d in readdir(corrected;join=true)]))
      h5open(summary,"r") do source
          v=create_group(out,"validation");a=attributes(source["validation"])
          for key in ("passed","guard_V","maximum_remaining_control_error_ENC","tolerance_ENC")
              v[key]=read(a[key])
          end
          v["maximum_solver_refinement_error_ENC"]=maximum(read(source["validation/refinement_error_ENC"]))
          v["source_summary"]=summary;v["source_sha256"]=sha(summary)
          @assert Bool(read(v["passed"]))
      end
      write(meta["complete"],true)
    end
    mv(target*".tmp",target;force=true);println(target)
end
end
if abspath(PROGRAM_FILE)==@__FILE__;MatchedAnalysis.main();end
