#!/usr/bin/env julia
# Statistics of all accepted events of a run kept as per-batch files.
#   julia --project=. scripts/production/prepare_batch_statistics.jl results/readout/<run_label>
# For each batch under results/readout/<run_label>/batches/ it reads
# <output_subdir>/features.h5 (prepare_batch_features.jl) or, when the table is
# absent, computes the same features from hits.h5 and the selection controls.
# Writes results/readout/<run_label>/<output_subdir>/statistics/events.h5 with
# one row per accepted event, the energy fits and the selection counts. No
# detector or electronics calculation is repeated. present_batch_statistics.py
# makes the page.
include("prepare_batch_features.jl")

"""Gaussian peak fits of a quantity in slices of another; NaN where a slice has too few events."""
function sliced_fits(values,by,edges;minimum_events=200)
    n=length(edges)-1;out=Dict("edges"=>collect(Float64.(edges)),"entries"=>zeros(Int64,n),"mu"=>fill(NaN,n),"sigma"=>fill(NaN,n),"mu_error"=>fill(NaN,n),"sigma_error"=>fill(NaN,n))
    for k in 1:n
        x=[v for (v,b) in zip(values,by) if isfinite(v)&&edges[k]<=b<edges[k+1]];out["entries"][k]=length(x)
        length(x)<minimum_events&&continue
        fit=fit_peak(x);out["mu"][k]=fit.mu;out["sigma"][k]=fit.sigma;out["mu_error"][k]=fit.mu_error;out["sigma_error"][k]=fit.sigma_error
    end
    out
end

function prepare_batch_statistics(run_dir::AbstractString)
    run_dir=abspath(run_dir);batch_root=joinpath(run_dir,"batches");isdir(batch_root)||error("No batches under $run_dir")
    tags=sort(filter(t->occursin(r"^b\d{5}$",t),readdir(batch_root)));isempty(tags)&&error("No batch folders under $batch_root")
    subdirs=unique(batch_output_subdir(joinpath(batch_root,t)) for t in tags);length(subdirs)==1||error("Expected one hit output folder per batch, found $subdirs")
    subdir=only(subdirs);out=joinpath(run_dir,subdir,"statistics");mkpath(out);target=joinpath(out,"events.h5");isfile(target)&&error("Preserve existing statistics: $target")
    started=time();tables=Dict{String,Any}[]
    for t in tags
        saved=joinpath(batch_root,t,subdir,"features.h5")
        push!(tables,isfile(saved) ? read_batch_features(saved) : compute_batch_features(joinpath(batch_root,t)))
        last(tables)["batch"]==t||error("Batch label differs in $t")
        println("batch $t: $(sum(length(x["events"]["event_number"]) for x in tables)) events so far, $(round(time()-started;digits=1)) s");flush(stdout)
    end
    first_table=first(tables)
    for x in tables, k in ("track_analysis","detector","ion_hit_cut_keV","electron_hit_cut_keV","species","run_label","settings")
        x[k]==first_table[k]||error("Setting $k differs in batch $(x["batch"])")
    end
    columns=Dict{String,Any}(k=>reduce(vcat,[x["events"][k] for x in tables]) for k in keys(first_table["events"]))
    allunique(columns["event_number"])||error("Duplicate event numbers across batches")
    dc=TOML.parse(first_table["detector"])
    budget=16*1024^2+2Base.summarysize(columns);GuardedFiles.space_check(out,budget,1024^3);tmp=target*".tmp"
    try
        h5open(tmp,"w") do f
            f["schema"]="itaca.batch_statistics.v1";f["species"]=first_table["species"];f["run_label"]=first_table["run_label"];f["batches"]=tags
            f["hit_files"]=[relpath(joinpath(batch_root,t,subdir,"hits.h5"),out) for t in tags];f["hit_file_sha256"]=[x["hit_file_sha256"] for x in tables]
            f["track_analysis"]=first_table["track_analysis"];f["detector"]=first_table["detector"]
            f["ion_hit_cut_keV"]=first_table["ion_hit_cut_keV"];f["electron_hit_cut_keV"]=first_table["electron_hit_cut_keV"]
            for (k,v) in columns;f["events/$k"]=v;end
            for k in keys(first_table["selection"])
                v=Int64[x["selection"][k] for x in tables];f["selection/batch/$k"]=v;f["selection/total/$k"]=sum(v)
            end
            for (k,v) in first_table["settings"];f["selection/settings/$k"]=v;end
            drift_edges=range(0,dc["chamber_length_mm"]/10;length=7)
            for (key,drift) in (("ion_energy_keV","ion_drift_cm"),("electron_energy_keV","electron_drift_cm"))
                values=Float64.(columns[key]);fit=fit_peak(filter(isfinite,values))
                for (k,v) in Base.pairs(fit);f["energy_fits/$key/$k"]=v;end
                for (k,v) in sliced_fits(values,columns[drift],drift_edges);f["energy_fits_by_drift/$key/$k"]=v;end
            end
            f["time_s"]=time()-started;flush(f);GuardedFiles.check_size(tmp,budget)
        end
        Base.rename(tmp,target)
    catch
        isfile(tmp)&&rm(tmp);rethrow()
    end
    println("STATISTICS $target events=$(length(columns["event_number"])) batches=$(length(tags)) time_s=$(round(time()-started;digits=1))")
    target
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1||error("Usage: prepare_batch_statistics.jl results/readout/<run_label>")
    prepare_batch_statistics(ARGS[1])
end
