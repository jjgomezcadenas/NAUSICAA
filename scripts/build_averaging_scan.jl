# Reuse saved transistor simulations; no new LTspice run is needed.
using HDF5, SHA
include(joinpath(@__DIR__, "../src/ion_signal/averaging_scan.jl"))
using .AveragingScan
root=normpath(joinpath(@__DIR__, ".."))
data=joinpath(root,"results/ltspice/capacitance/pad_1.25pF")
paths=[joinpath(data,name) for name in ("ions_N1_r2.csv.gz","ions_N1_r1.csv.gz","noise_400.csv","noise_200.csv")]
t,v=read_columns(paths[1],["physical_time_s","baseline_subtracted_output_V"])
tc,vc=read_columns(paths[2],["physical_time_s","baseline_subtracted_output_V"])
f,a=read_columns(paths[3],["frequency","V(onoise)"])
fc,ac=read_columns(paths[4],["frequency","V(onoise)"])
widths=[180.,200.,250.,300.,350.]*1e-6
frame=.01; tau=.3; alpha=exp(-frame/tau)
gains=[waveform_gain(t,v,w,frame,alpha) for w in widths]
sigmas=[noise_rms(f,a,w,frame,alpha) for w in widths]
qerr=[abs(noise_rms(f,a,w,frame,alpha;order=16)/s-1) for (w,s) in zip(widths,sigmas)]
serr=[abs(noise_rms(fc,ac,w,frame,alpha)/s-1) for (w,s) in zip(widths,sigmas)]
terr=[abs(waveform_gain(tc,vc,w,frame,alpha)/g-1) for (w,g) in zip(widths,gains)]
@assert maximum(qerr)<1e-6
@assert maximum(serr)<.005
@assert maximum(terr)<.005
# Independent integral identity: white PSD with non-overlapping box averages.
white=noise_rms(f,ones(length(f)),widths[4],frame,alpha)
@assert abs(white/sqrt((1+alpha^2)/widths[4]/2)-1)<.002
out=joinpath(root,"results/ltspice/averaging_scan.h5")
h5open(out,"w") do h
    h["averaging_s"]=widths; h["gain_V_per_ion"]=gains
    h["noise_rms_V"]=sigmas; h["enc_electrons_rms"]=sigmas./gains
    h["quadrature_relative_error"]=qerr; h["spectral_relative_error"]=serr
    h["transient_relative_error"]=terr
    attrs=attributes(h)
    attrs["frame_s"]=frame; attrs["correction_tau_s"]=tau
    attrs["external_capacitance_F"]=1.25e-12; attrs["feedback_capacitance_F"]=3e-15
    attrs["temperature_K"]=293.15; attrs["frequency_low_Hz"]=first(f); attrs["frequency_high_Hz"]=last(f)
    attrs["assumptions"]="Uniform mean ion arrivals over 10 ms; centered continuous voltage averages; stationary one-sided SPICE noise at the operating point; supplied-model flicker terms zero; intrinsic MOS capacitances included. Input current: central square pad in grounded plane, finite 5 mm path, no charge renormalization."
    sources=vcat(paths,[joinpath(data,"study.json"),@__FILE__,joinpath(root,"src/ion_signal/averaging_scan.jl")])
    h["source_paths"]=[relpath(p,root) for p in sources]
    h["source_sha256"]=[bytes2hex(sha256(read(p))) for p in sources]
end
for (w,g,s) in zip(widths,gains,sigmas)
    println("$(w*1e6) us: ENC=$(s/g) electrons RMS, gain=$(g*1e6) uV/ion")
end
println("Maximum relative checks: quadrature=$(maximum(qerr)), spectrum=$(maximum(serr)), transient=$(maximum(terr))")
