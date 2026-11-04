#!/usr/bin/env julia
"""Matched square-sensor/TPB pilots using nested-aperture photon thinning."""
module TPBPilot
using HDF5,TOML,Random,Distributions,Statistics
include("benchmark_sipm_production.jl")
const EA=Prod.EA
const PIPE=Prod.Pilot
const S=Prod.S
const OUT_TPB=joinpath(ROOT,"results/sipm_tpb_pilot")
toml_string(d)=sprint(io->TOML.print(io,d;sorted=true))
function save_card(path,d)
    full=joinpath(ROOT,path);mkpath(dirname(full));write(full,toml_string(d));path
end

"""Thin every photon cell into a contained aperture, preserving event and EL ledgers.

Probability is A_new/A_parent for the uniform axial-emission prescription.
Each photon is an independent trial, so binomial composition gives exactly
the direct smaller-aperture marginal law. Retain occupied cells and owned
coordinate vectors. The RNG belongs to the caller.
"""
function thin_history(rng,h,p)
    isfinite(p)&&0<=p<=1 || throw(ArgumentError("Aperture ratio must be in [0,1]"))
    counts=Int64[rand(rng,Binomial(n,p)) for n in h.cells.count]
    keep=findall(>(0),counts);c=h.cells
    cells=EA.SparseCells(c.pixel_ix[keep],c.pixel_iy[keep],c.bucket[keep],counts[keep])
    EA.PhotonSiPMHistory(h.placement,h.gas_model,h.grid,cells,h.initial_pairs,
        h.surviving_electrons,h.emitted_photons)
end

"""Write 100-event VUV histories for a square sensor and a standard input index.

Use archived circular-aperture photon sources for both gases. Identity includes
source content, geometry, event selection, thinning seed and implementation.
The source deposited charge, attachment, EL fluctuations and bucket assignments
are retained. Metadata totals are recomputed for the selected events.
"""
function square_photons(tag;side_mm=1.3,n=100)
    sourceidx=tag=="ep2" ? "results/simulations/production_pilot/electron_250ns_ep2_pilot100.h5" :
        "results/simulations/production_pilot/electron_250ns_summary_ea65f602bf417abb.h5"
    code=S.content_hash(read(@__FILE__,String)*read(joinpath(ROOT,"src/event_analysis/grid_readout.jl"),String))
    paths=Dict{String,String}()
    h5open(joinpath(ROOT,sourceidx),"r") do index
        for gas in ("pure_xe","xe_he15")
            source=String(read(index["bb0nu/$(gas)_path"]))
            hash=Prod.sha(source)
            h5open(source,"r") do input
                read(input["metadata/complete"]) || error("Incomplete parent photons")
                old=TOML.parse(String(read(input["metadata/branch_card_toml"])))
                old["sipm_diameter_mm"]==old["hole_diameter_mm"] || error("Parent aperture")
                ac=EA.axial_light_acceptance(Float32(old["hole_diameter_mm"]/2),Float32(old["hole_depth_mm"]))
                aq=EA.axial_square_acceptance(side_mm,old["hole_depth_mm"],old["hole_diameter_mm"]/2)
                ratio=Float64(aq)/Float64(ac)
                ids=sort!(parse.(Int64,collect(keys(input["events"]))))[1:n]
                card=deepcopy(old);delete!(card,"sipm_diameter_mm")
                card["sipm_shape"]="square";card["sipm_side_mm"]=side_mm;card["sipm_plane_gap_mm"]=0.0
                card["light_acceptance"]=Float64(aq)
                card["optical_model"]="uniform_axial_emission_absorbing_walls"
                save_card("config/specs/detector/electrons/elcc_250ns_$(tag)_square1315.toml",card)
                identity=S.content_hash(hash*toml_string(card)*string(ids)*code*"seed2026")
                parent=dirname(dirname(dirname(dirname(source))))
                target=joinpath(parent,"detector/photons",identity,"photons_$(gas).h5")
                paths[gas]=target
                if isfile(target)
                    h5open(target,"r") do f
                        read(f["metadata/complete"]) && read(f["metadata/card_hash"])==identity || error("Stale square cache")
                    end
                    return
                end
                mkpath(dirname(target))
                h5open(target*".tmp","w") do output
                    attributes(output)["schema"]="nausicaa.photon_sipm_file.pilot.v1"
                    m=create_group(output,"metadata")
                    replaced=Set(["complete","card_hash","code_hash","branch_card_toml","event_count",
                        "initial_pairs","surviving_electrons","emitted_photons","incident_photons","elapsed_s"])
                    for k in keys(input["metadata"])
                        k in replaced || (m[k]=read(input["metadata"][k]))
                    end
                    Prod.put(m,(;complete=false,card_hash=identity,code_hash=code,branch_card_toml=toml_string(card),
                        event_count=n,parent_photon_path=source,parent_photon_sha256=hash,
                        thinning_probability=ratio,square_acceptance=aq,circle_acceptance=ac,
                        thinning_seed=2026,thinning_source=read(@__FILE__,String)))
                    totals=zeros(Int64,4);started=time()
                    for id in ids
                        original=input["events/$id"];h=EA.read_photon_history(original)
                        seed=S.event_seed(identity,id,"square_aperture:2026")
                        result=thin_history(Xoshiro(seed),h,ratio)
                        g=create_group(output,"events/$id");EA.write_photon_history(g,result)
                        g["pair_checksum"]=read(original["pair_checksum"]);g["thinning_seed"]=seed
                        totals .+= [result.initial_pairs,result.surviving_electrons,
                            result.emitted_photons,sum(result.cells.count)]
                    end
                    for (k,v) in zip(("initial_pairs","surviving_electrons","emitted_photons","incident_photons"),totals);m[k]=v;end
                    m["elapsed_s"]=time()-started;write(m["complete"],true)
                end
                mv(target*".tmp",target)
            end
        end
    end
    mkpath(OUT_TPB);write(joinpath(OUT_TPB,".gitignore"),"*.h5\n")
    idx=joinpath(OUT_TPB,"photons_$(tag).h5")
    h5open(idx,"w") do f
        f["metadata/complete"]=true;f["metadata/parent_index"]=sourceidx
        for (gas,p) in paths;f["bb0nu/$(gas)_path"]=p;end
    end
    relpath(idx,ROOT)
end

"""Create explicit sensor, ASIC and run cards for a selected recovery sensitivity.

Blue PDE=0.32, TPB conversion=0.62, forward fraction=0.5; the library samples
the product. ASIC input charge noise is preserved from 0.5 PE at gain 2.55e6.
The card records that recharge is provisional and the blue PDE is a peak-value
approximation. All stage outputs use the existing immutable production layout.
"""
function response_cards(tag,index;tau_ns=12.5,substeps=512)
    label="tpb_$(tag)_tau$(replace(string(tau_ns),'.'=>'p'))_n$(substeps)"
    sensor=TOML.parsefile(joinpath(ROOT,"config/specs/response/sipm/sipm_statistical_v1.toml"))
    delete!(sensor,"diameter_m")
    merge!(sensor,Dict("shape"=>"square","side_m"=>.0013,"pde"=>.32,"microcell_pitch_m"=>15e-6,
        "microcell_count"=>7284,"gain"=>3.6e5,"recovery_s"=>tau_ns*1e-9,
        "dark_reference_rate_Hz"=>120000.,"dark_reference_area_m2"=>.0013^2,
        "statistical_substeps"=>substeps,"reference_overvoltage_V"=>4.0,
        "tpb"=>Dict("conversion_efficiency"=>.62,"forward_fraction"=>.5)))
    sensor["provenance"]=Dict("device"=>"Hamamatsu S14160-1315PS, 25 C, +4 V",
        "datasheet"=>"https://www.hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/s14160-1310ps_etc_kapd1070e.pdf",
        "pde"=>"0.32 peak blue PDE used as pilot approximation to TPB-spectrum average",
        "tpb"=>"Benson et al. 2018, Fig.19: approximate interpolation of intrinsic QE to 175 nm; complete absorption; user forward factor 0.5",
        "tpb_source"=>"https://doi.org/10.1140/epjc/s10052-018-5807-z",
        "recovery"=>"Provisional recharge proxy: 12.5 ns measured pulse decay; 25 ns sensitivity; not a direct recharge measurement",
        "illumination"=>"Uniform microcells, constant PDE during recovery; prompt TPB; absorbing surroundings")
    cp=save_card("config/specs/response/sipm/s14160_1315ps_tpb_tau$(replace(string(tau_ns),'.'=>'p'))_n$(substeps).toml",sensor)
    asic=TOML.parsefile(joinpath(ROOT,"config/specs/response/sipm_electronics/readout_v1.toml"))
    asic["noise_rms_pe"]=.5*2.55e6/3.6e5
    asic["provenance"]["noise"]="Same input charge RMS as 0.5 PE at gain 2.55e6; rescaled to gain 3.6e5"
    ep=save_card("config/specs/response/sipm_electronics/readout_s14160_1315ps.toml",asic)
    pilot=TOML.parsefile(joinpath(ROOT,"config/specs/runs/sipm_response_statistical_ep2_pilot.toml"))
    merge!(pilot,Dict("photon_index"=>index,"sipm_card"=>cp,"electronics_card"=>ep,
        "results_directory"=>"results/sipm_tpb_pilot","report_directory"=>"latex_generated/sipm_tpb_pilot"))
    pp=save_card("config/specs/runs/sipm_response_$(label).toml",pilot)
    run=TOML.parsefile(joinpath(ROOT,"config/specs/runs/sipm_bb0nu_statistical_ep2_pilot100.toml"))
    merge!(run,Dict("photon_index"=>index,"pilot_card"=>pp,"results_directory"=>"results/sipm_tpb_pilot/$(label)",
        "report_directory"=>"latex_generated/sipm_tpb_pilot"))
    save_card("config/specs/runs/sipm_bb0nu_statistical_$(label)_pilot100.toml",run)
    label
end

"""Run both operating points and gas models, with optional recharge/time-step scans."""
function main(;tau_ns=12.5,substeps=512)
    labels=[response_cards(tag,square_photons(tag);tau_ns,substeps) for tag in ("ep2","ep2p5")]
    benchmark_response(tags=labels,output_dir=joinpath(OUT_TPB,"tau$(replace(string(tau_ns),'.'=>'p'))_n$(substeps)"))
end
end
if abspath(PROGRAM_FILE)==@__FILE__
    TPBPilot.main(tau_ns=length(ARGS)>=1 ? parse(Float64,ARGS[1]) : 12.5,
        substeps=length(ARGS)>=2 ? parse(Int,ARGS[2]) : 512)
end
