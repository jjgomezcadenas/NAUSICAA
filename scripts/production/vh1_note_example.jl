#!/usr/bin/env julia
# One simulated VH1 pad for the note: readings, codes and reconstruction.
#
# Writes results/readout/vh1_note/example_pad.csv (per reading: time, noiseless
# and noisy voltage), example_counts.csv (per bucket: input, reconstructed,
# predicted sigma) and readings_scan.csv (per-bucket noise against readings per
# frame). No event data; counts are a fixed profile.
#
#     JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia" julia --startup-file=no --project=. scripts/production/vh1_note_example.jl
using Random,LinearAlgebra,Statistics,DelimitedFiles
const ROOT=normpath(joinpath(@__DIR__,"../.."))
include(joinpath(ROOT,"src/event_analysis/EventAnalysis.jl"))
include(joinpath(ROOT,"src/detector_response/DetectorResponse.jl"))
include(joinpath(ROOT,"src/electronics/ElectronicsResponse.jl"))
include(joinpath(ROOT,"src/reconstruction/SignalReconstruction.jl"))
include(joinpath(ROOT,"src/vh1/VH1Readout.jl"))
using .EventAnalysis,.DetectorResponse,.ElectronicsResponse,.SignalReconstruction,.VH1Readout
const ER=ElectronicsResponse;const V=VH1Readout
out=joinpath(ROOT,"results/readout/vh1_note");mkpath(out)
s=load_vh1_settings(ROOT,joinpath(ROOT,"config/specs/response/ion_csa/vh1_device_changes.toml");minimum_samples=4096)
withk(s,k,w)=VH1Settings(s.gain_V_per_ion,s.charge_fraction,s.feedback_capacitance_F,s.v_ref_V,s.window_V,s.frame_ms,k,Float32(w),s.reference_frames,s.adc,s.noise_sample_us,s.noise,s.current_pulse,s.provenance)
counts=Int64[0,0,80,350,900,1500,1200,600,200,40,0,0];nb=length(counts)
t=IonReadoutTime(1000f0,20f0,nb)
channels=ReadoutChannels(Int32[1],Int32[1],Float32[0],Float32[0])
d=IonDetectorResponse(1,EventPlacement(1,50f0,100f0),EventFrame(mc_vertex_frame,(0f0,0f0,0f0),(1f0,1f0,1f0)),
    DetectorVectors{IonReadoutTime}(channels,t,reshape(counts,nb,1)),sum(counts),0,ReadoutLosses(0,0,0,0,0,0),SparseCells(Int64[],Int64[],Int64[],Int64[]))
ft=vh1_frame_time(t,s);grid=V.internal_grid(ft,s)
packet=ER.prepare_bucket_current(s.current_pulse,t.bucket_ms)
current=ER.prepared_ion_current(counts,t,packet,grid)
clean=V.frame_readings!(zeros(V.readings(ft,s)),current,ft,s)
r=sim_vh1_electronics(Xoshiro(7),Xoshiro(8),d,s)
noisy=vec(frame_voltages(r,s))
plan=prepare_vh1_reconstruction(ft,t,s);ri=reco_vh1_ions(r,plan)
times=[ (last(V.window_range(i,s))*Float64(s.noise_sample_us)*1e-3) for i in 1:V.readings(ft,s)]  # ms from acquisition start
open(joinpath(out,"example_pad.csv"),"w") do io
    println(io,"reading_end_ms,noiseless_V,noisy_V")
    for i in eachindex(times);println(io,times[i],",",clean[i],",",noisy[i]);end
end
sigma=sqrt.(diag(plan.ion_covariance))
open(joinpath(out,"example_counts.csv"),"w") do io
    println(io,"bucket_start_ms,input,reconstructed,predicted_sigma")
    for k in 1:nb;println(io,(t.start_ms-ft.start_ms)+(k-1)*t.bucket_ms,",",counts[k],",",ri.ions[k,1],",",sigma[k]);end
end
t30=IonReadoutTime(1000f0,20f0,30)
open(joinpath(out,"readings_scan.csv"),"w") do io
    println(io,"conversions_per_frame,window_ms,median_sigma_ions,max_sigma_ions")
    for (k,w) in ((1,20.0),(1,5.0),(2,10.0),(4,5.0),(10,2.0),(20,1.0))
        sw=withk(s,k,w);pw=prepare_vh1_reconstruction(vh1_frame_time(t30,sw),t30,sw);sg=sqrt.(diag(pw.ion_covariance))
        println(io,k,",",w,",",median(sg),",",maximum(sg))
    end
end
println("frames ",ft.frame_count," readings ",V.readings(ft,s)," reference end at ",ft.start_ms+ft.frame_ms," ms; bucket 1 starts at ",t.start_ms)
println("input ",counts);println("reconstructed ",round.(vec(ri.ions);digits=1));println("sigma ",round.(sigma;digits=1))
println("saturated ",r.saturated[1]," valid ",all(ri.valid))
