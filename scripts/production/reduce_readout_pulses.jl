#!/usr/bin/env julia
# Reduce sparse ion/pixel and photon/SiPM histories to event pulse statistics.
using HDF5, SHA, Statistics

const ROOT=normpath(joinpath(@__DIR__,"../.."))

"""Read a scalar HDF5 dataset as a string."""
saved_string(group,name)=String(read(group[name]))

"""Return one event's count, mean and maximum over positive readout cells.

One cell is a pixel–frame pulse for ions or a SiPM–bucket pulse for photons.
All occupied cells enter. An empty event has zero metrics and is excluded from
the mean/maximum distributions, while remaining in the event-count ledger.
"""
function pulse_metrics(counts::Vector{Int64})
    all(>(0),counts) || error("Sparse readout contains a nonpositive cell")
    isempty(counts) && return (Int64(0),Int64(0),0.0,Int64(0))
    total=foldl(Base.Checked.checked_add,counts;init=Int64(0))
    return (Int64(length(counts)),total,Float64(total)/length(counts),maximum(counts))
end

"""Read one complete sparse history file and return aligned event summaries.

The saved `count` vector is already indexed by transverse pixel or SiPM and
longitudinal time bucket. Verify its event-level charge/light ledger while
reducing; retain event IDs and exact integer totals.
"""
function reduce_file(path::String,kind::String,expected_events::Int,expected_total::Int64)
    h5open(path,"r") do file
        schema=kind=="ions" ? "nausicaa.ion_pad_file.pilot.v1" :
            "nausicaa.photon_sipm_file.pilot.v1"
        String(read(attributes(file)["schema"]))==schema || error("Readout schema mismatch: $path")
        m=file["metadata"]
        Bool(read(m["complete"])) && Int(read(m["event_count"]))==expected_events ||
            error("Incomplete readout or event-count mismatch: $path")
        ids=sort!(parse.(Int64,collect(keys(file["events"]))))
        length(ids)==expected_events || error("Missing event groups: $path")
        occupied=zeros(Int64,expected_events)
        totals=zeros(Int64,expected_events)
        event_mean=zeros(Float64,expected_events)
        event_max=zeros(Int64,expected_events)
        for (i,id) in enumerate(ids)
            group=file["events/$(id)"]
            Int64(read(group["event_id"]))==id || error("Event ID mismatch: $path:$id")
            occupied[i],totals[i],event_mean[i],event_max[i]=
                pulse_metrics(Vector{Int64}(read(group["count"])))
            if kind=="ions"
                totals[i]==Int64(read(group["initial_pairs"])) ||
                    error("Ion ledger mismatch: $path:$id")
            else
                totals[i]<=Int64(read(group["emitted_photons"])) ||
                    error("Photon ledger mismatch: $path:$id")
            end
        end
        sum(totals)==expected_total || error("Aggregate signal mismatch: $path")
        return (;ids,occupied,totals,event_mean,event_max,
                source_sha=saved_string(m,"source_sha256"),
                code_hash=saved_string(m,"code_hash"),
                card_hash=saved_string(m,"card_hash"),
                signal_unit=saved_string(m,"signal_unit"))
    end
end

"""Use common logarithmic edges for one carrier and statistic across samples.

The endpoints extend beyond the observed extrema, so all positive event
values enter exactly one bin. Eighty bins preserve comparable log-scale plots.
"""
function log_edges(values::Vector{Float64};bins::Int=80)
    all(>(0),values) || error("Histogram values must be positive")
    low=minimum(values)*0.95
    high=maximum(values)*1.05
    return 10.0 .^ range(log10(low),log10(high);length=bins+1)
end

"""Count event values into supplied edges with explicit under/overflow checks."""
function histogram(values::Vector{Float64},edges::Vector{Float64})
    counts=zeros(Int64,length(edges)-1)
    for value in values
        index=searchsortedlast(edges,value)
        1<=index<length(edges) || error("Histogram value outside saved edges")
        counts[index]+=1
    end
    sum(counts)==length(values) || error("Histogram entry ledger failed")
    return counts
end

"""Reduce six ion branches and four 250 ns photon branches.

Write per-event mean and maximum pulse counts, exact signal totals, shared
histogram edges, summary statistics and source/code identities to one HDF5
file. The presentation stage only reads these numerical results.
"""
function reduce_scan(index_path::String,electron_index_path::String,output_path::String)
    cases=Dict{Tuple{String,String,String},Any}()
    cards=Dict{String,String}()
    index_sha=bytes2hex(sha256(read(index_path)))
    electron_index_sha=bytes2hex(sha256(read(electron_index_path)))
    electron_bucket_us=Ref(0.0)
    h5open(index_path,"r") do index
        String(read(attributes(index)["schema"]))=="nausicaa.detector_grid_scan.v1" ||
            error("Grid-scan index schema mismatch")
        Bool(read(index["metadata/complete"])) || error("Incomplete grid-scan index")
        for name in ("realisation_card","ion_card_2mm",
                     "ion_card_3mm","ion_card_4mm","run_id")
            cards[name]=saved_string(index["metadata"],name)
        end
        h5open(electron_index_path,"r") do electron_index
            String(read(attributes(electron_index)["schema"]))==
                "nausicaa.electron_rebin_summary.v1" || error("Electron index schema mismatch")
            Bool(read(electron_index["metadata/complete"])) || error("Incomplete electron index")
            cards["electron_card"]=saved_string(electron_index["metadata"],"electron_card")
            cards["electron_run_id"]=saved_string(electron_index["metadata"],"run_id")
            saved_string(electron_index["metadata"],"realisation_card")==
                cards["realisation_card"] || error("Primary realisation cards differ")
            electron_bucket_us[]=Float64(read(electron_index["metadata/bucket_width_us"]))
            electron_bucket_us[]==0.25 || error("Expected 0.25 µs photon buckets")
            for sample in ("bb0nu","single_electrons")
                g=index[sample];eg=electron_index[sample]
                events=Int(read(g["events"]))
                events==Int(read(eg["events"])) || error("Event counts differ: $sample")
                for field in ("realisation_id","source_sha256")
                    saved_string(g,field)==saved_string(eg,field) ||
                        error("$field differs: $sample")
                end
                Int64(read(g["pairs"]))==Int64(read(eg["initial_pairs"])) ||
                    error("Primary pair totals differ: $sample")
                for pitch in ("2","3","4")
                    branch=g["ions/$pitch"]
                    path=saved_string(branch,"path")
                    cases[("ions",sample,pitch)]=(
                        data=reduce_file(path,"ions",events,Int64(read(branch["ions"]))),
                        path=path,label="$(pitch) mm",frame=Float64(read(branch["bucket_width_ms"])))
                end
                for gas in ("pure_xe","xe_he15")
                    path=saved_string(eg,"$(gas)_path")
                    cases[("photons",sample,gas)]=(
                        data=reduce_file(path,"photons",events,
                            Int64(read(eg["$(gas)_incident"]))),
                        path=path,label=gas,frame=electron_bucket_us[])
                end
            end
        end
    end
    mkpath(dirname(output_path))
    temporary=output_path*".tmp"
    isfile(temporary) && rm(temporary)
    h5open(temporary,"w") do out
        attributes(out)["schema"]="nausicaa.readout_pulses.v2"
        m=create_group(out,"metadata")
        m["complete"]=false
        m["scan_index"]=abspath(index_path)
        m["scan_index_sha256"]=index_sha
        m["electron_index"]=abspath(electron_index_path)
        m["electron_index_sha256"]=electron_index_sha
        m["reducer_sha256"]=bytes2hex(sha256(read(@__FILE__)))
        for (name,body) in cards
            m[name]=body
        end
        m["event_mean_definition"]="total positive cell count / number of positive cells in each event"
        m["event_max_definition"]="largest positive cell count in each event"
        m["ion_cell"]="one pixel and one time frame"
        m["photon_cell"]="one SiPM and one time bucket"
        m["selection"]="all occupied cells; no count threshold; all events with pulses"
        groups=create_group(out,"cases")
        edges_group=create_group(out,"histogram_edges")
        for kind in ("ions","photons")
            kg=create_group(groups,kind)
            eg=create_group(edges_group,kind)
            for field in ("mean","max")
                all_values=Float64[]
                for ((case_kind,_,_),case) in cases
                    case_kind==kind || continue
                    data=case.data
                    values=field=="mean" ? data.event_mean : Float64.(data.event_max)
                    append!(all_values,values[data.occupied.>0])
                end
                eg[field]=log_edges(all_values)
            end
            for sample in ("bb0nu","single_electrons")
                sg=create_group(kg,sample)
                labels=kind=="ions" ? ("2","3","4") : ("pure_xe","xe_he15")
                for label in labels
                    case=cases[(kind,sample,label)]
                    data=case.data
                    g=create_group(sg,label)
                    g["source_file"]=case.path
                    g["source_sha256"]=data.source_sha
                    g["code_hash"]=data.code_hash
                    g["card_hash"]=data.card_hash
                    g["signal_unit"]=data.signal_unit
                    g["label"]=case.label
                    g["bucket_duration"]=case.frame
                    g["bucket_unit"]=kind=="ions" ? "ms" : "us"
                    g["event_id"]=data.ids
                    g["occupied_cells"]=data.occupied
                    g["signal_total"]=data.totals
                    g["event_mean"]=data.event_mean
                    g["event_max"]=data.event_max
                    selected=data.occupied.>0
                    values_mean=data.event_mean[selected]
                    values_max=Float64.(data.event_max[selected])
                    g["events"]=length(data.ids)
                    g["events_with_charge"]=count(selected)
                    g["total_occupied_cells"]=sum(data.occupied)
                    g["total_signal"]=sum(data.totals)
                    g["mean_of_event_means"]=mean(values_mean)
                    g["median_of_event_means"]=median(values_mean)
                    g["mean_of_event_maxima"]=mean(values_max)
                    g["median_of_event_maxima"]=median(values_max)
                    g["hist_mean_counts"]=histogram(values_mean,read(eg["mean"]))
                    g["hist_max_counts"]=histogram(values_max,read(eg["max"]))
                end
            end
        end
        write(m["complete"],true)
    end
    mv(temporary,output_path;force=true)
    return output_path
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==3 || error("Usage: reduce_readout_pulses.jl GRID_SCAN_INDEX.h5 ELECTRON_INDEX.h5 OUTPUT.h5")
    println(reduce_scan(ARGS[1],ARGS[2],ARGS[3]))
end
