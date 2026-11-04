#!/usr/bin/env julia
using HDF5,TOML,Statistics,SHA
const ROOT=normpath(joinpath(@__DIR__,"../.."))
resolve(p)=isabspath(p) ? p : joinpath(ROOT,p)

"""Correct the full acquired sensor sum for expected dark baseline and attachment.

The stationary dark-only mean charge per avalanche is 1/(1+rate*tau/Ncells).
Subtract its expected total over all acquired cells, then multiply by exp(t/lifetime).
Optical illumination can reduce dark avalanche charge further through recovery;
this operational baseline correction uses the independently measurable dark-only mean.
Drift distance is cm, velocity mm/us, lifetime ms, and returned charge is PE-equivalent.
"""
function corrected_total(total,cells,dt,sensor,electron,L_cm)
    area=get(sensor,"shape","circle")=="square" ? sensor["side_m"]^2 : pi*sensor["diameter_m"]^2/4
    rate=sensor["dark_reference_rate_Hz"]*area/sensor["dark_reference_area_m2"]
    dark=cells*rate*dt/(1+rate*sensor["recovery_s"]/sensor["microcell_count"])
    drift_ms=L_cm*10/electron["drift_speed_mm_us"]/1000
    ((total-dark)*exp(drift_ms/electron["lifetime_ms"]),dark)
end

"""Save values, full-range histogram and moment-based relative RMS/FWHM.

FWHM is 2sqrt(2log(2))*RMS, a Gaussian-equivalent width. It is not a peak fit.
"""
function distribution(parent,key,values)
    g=create_group(parent,key);g["values"]=values
    lo,hi=extrema(values);lo==hi && (lo-=.5;hi+=.5)
    edges=collect(range(lo,hi;length=101));counts=zeros(Int64,100)
    for x in values;counts[clamp(searchsortedlast(edges,x),1,100)]+=1;end
    mu=mean(values);sigma=std(values;corrected=false)
    g["edges"]=edges;g["counts"]=counts;g["mean"]=mu;g["median"]=median(values)
    g["rms"]=sigma;g["relative_rms_percent"]=100sigma/mu
    g["gaussian_equivalent_fwhm_percent"]=100sqrt(8log(2))*sigma/mu
end

"""Stream pulse and energy statistics for old and corrected E/P productions.

Inclusive pulses use photon-positive SiPM–bucket cells. Event totals use every
acquired cell, including dark counts. Match source energy by event ID and retain
both corrected PE totals and response per deposited MeV as a diagnostic.
"""
function main(;limit=typemax(Int),tags=("legacy","ep2p5","ep3"),output_dir="results/electron_operating_points")
    points=[("legacy","results/sipm_response_production/0246ff441c441632/manifest.toml")]
    for tag in filter(!=("legacy"),tags)
        index=TOML.parsefile(resolve("results/electron_operating_points/$tag.toml"))
        index["complete"] || error("Incomplete operating point")
        push!(points,(tag,index["sipm_manifest"]))
    end
    energies=h5open(resolve("events/0nubb_hitdf.h5"),"r") do f
        Dict(zip(Int.(read(f["events/event_id"])),Float64.(read(f["events/energy_MeV"]))))
    end
    outdir=resolve(output_dir);mkpath(outdir)
    h5open(joinpath(outdir,"electron_operating_points.h5.tmp"),"w") do out
        attributes(out)["schema"]="nausicaa.electron_operating_points.v1"
        out["reducer_sha256"]=bytes2hex(sha256(read(@__FILE__)))
        for (tag,path) in points
            manifest=TOML.parsefile(resolve(path));manifest["complete"] || error("Incomplete response")
            for gas in ("pure_xe","xe_he15")
                println("Reducing $tag/$gas");flush(stdout)
                g=create_group(out,"$tag/$gas");g["manifest_toml"]=read(resolve(path),String)
                ids=Int[];totals=Float64[];corrected=Float64[];responses=Float64[];dark=Float64[]
                clipped=Int64[];pulses=Int64[];output=Float64[];energy=Float64[]
                histogram=Dict{Int,Int64}();pulse_sum=0.;pulse_count=0;binwidth=25.
                for sh in manifest["shards"]
                    length(ids)>=limit && break
                    sh["gas"]==gas && sh["sample"]=="bb0nu" || continue
                    h5open(resolve(sh["sensor"]),"r") do sf
                        h5open(resolve(sh["asic"]),"r") do af
                            sensor=TOML.parse(String(read(sf["metadata/sensor_card_toml"])))
                            dt=read(sf["metadata/bucket_width_s"])
                            h5open(String(read(sf["metadata/source_path"])),"r") do ph
                                electron=TOML.parse(String(read(ph["metadata/branch_card_toml"])))
                                if !haskey(g,"electron_card_toml")
                                    g["electron_card_toml"]=String(read(ph["metadata/branch_card_toml"]))
                                    g["sensor_card_toml"]=String(read(sf["metadata/sensor_card_toml"]))
                                    g["electronics_card_toml"]=String(read(af["metadata/electronics_card_toml"]))
                                    g["reduced_field_kV_cm_bar"]=electron["reduced_field_kV_cm_bar"]
                                    g["j_variance_ratio"]=electron["j_variance_ratio"]
                                end
                                for id in sh["event_ids"]
                                    length(ids)>=limit && break
                                    s=sf["events/$id"];a=af["events/$id"];p=ph["events/$id"]
                                    read(s["pair_checksum"])==read(p["pair_checksum"])==read(a["pair_checksum"]) || error("Pair mismatch")
                                    lookup=Dict(k=>i for (i,k) in enumerate(zip(read(s["pixel_ix"]),read(s["pixel_iy"]),read(s["bucket"]))))
                                    positions=[lookup[k] for k in zip(read(p["pixel_ix"]),read(p["pixel_iy"]),read(p["bucket"]))]
                                    charge=read(s["sensor_pe_equiv"]);v=charge[positions]
                                    for x in v
                                        x>=0 && isfinite(x) || error("Invalid sensor charge")
                                        b=floor(Int,x/binwidth);histogram[b]=get(histogram,b,0)+1
                                    end
                                    pulse_sum+=sum(v);pulse_count+=length(v)
                                    total=sum(charge);c,d=corrected_total(total,length(charge),dt,sensor,electron,read(s["electron_drift_cm"]))
                                    push!(ids,id);push!(totals,total);push!(corrected,c);push!(responses,c/energies[id]);push!(dark,d)
                                    push!(pulses,length(v));push!(clipped,count(identity,Bool.(read(a["saturated"])[positions])))
                                    push!(output,Float64(read(a["readout_pe_total"])));push!(energy,energies[id])
                                end
                            end
                        end
                    end
                end
                g["event_id"]=ids;g["expected_dark_pe"]=dark;g["occupied_pulses"]=pulses;g["clipped_pulses"]=clipped
                g["clipped_pulse_fraction"]=sum(clipped)/sum(pulses)
                g["readout_charge_retention"]=sum(output)/sum(totals)
                h=create_group(g,"inclusive_sensor_pe");lastbin=maximum(keys(histogram))
                h["edges"]=collect(0:lastbin+1).*binwidth
                h["counts"]=[get(histogram,b,Int64(0)) for b in 0:lastbin]
                h["mean"]=pulse_sum/pulse_count;h["pulse_count"]=pulse_count
                distribution(g,"total_sensor_pe",totals)
                distribution(g,"attachment_corrected_pe",corrected)
                distribution(g,"corrected_pe_per_MeV",responses)
                distribution(g,"source_energy_MeV",energy)
                distribution(g,"total_readout_pe",output)
            end
        end
        # Common edges make overlaid event distributions directly comparable.
        for gas in ("pure_xe","xe_he15"),key in ("total_sensor_pe","attachment_corrected_pe","corrected_pe_per_MeV","total_readout_pe")
            groups=[out["$tag/$gas/$key"] for (tag,_) in points]
            vectors=[read(g["values"]) for g in groups]
            lo,hi=extrema(vcat(vectors...));lo==hi && (lo-=.5;hi+=.5)
            edges=collect(range(lo,hi;length=101))
            for (g,v) in zip(groups,vectors)
                counts=zeros(Int64,100)
                for x in v;counts[clamp(searchsortedlast(edges,x),1,100)]+=1;end
                write(g["edges"],edges);write(g["counts"],counts)
            end
        end
        out["complete"]=true
    end
    mv(joinpath(outdir,"electron_operating_points.h5.tmp"),joinpath(outdir,"electron_operating_points.h5");force=true)
end
if abspath(PROGRAM_FILE)==@__FILE__
    main()
end
