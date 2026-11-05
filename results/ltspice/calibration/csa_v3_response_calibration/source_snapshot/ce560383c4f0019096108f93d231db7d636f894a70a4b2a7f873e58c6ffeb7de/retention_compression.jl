module RetentionCompression
using HDF5
include("pad_compression.jl")
const PC=PadCompression

"""Use the same calibrated kernel with an explicitly bounded, independently measured baseline tolerance."""
function PC.kernel(c::PC.Context,name)
 get!(c.kernels,name) do
  s=PC.sample(c,name);k=PC.P.kernel(c.root,c.readout,s["frame_s"])
  noise=h5open(joinpath(c.root,c.readout["circuit_analysis"]),"r") do h
   group=h[c.readout["circuit_name"]]["$(round(Int,1000s["frame_s"]))ms"]
   j=read(attributes(group)["best_index"]);read(group["noise_rms_V"])[j]
  end
  baseline_error=h5open(joinpath(c.root,c.readout["circuit_results"]),"r") do h
   maximum(abs.(read(h["$(round(Int,1000s["frame_s"]))ms/single_N0_r2/V(out)"]).-k.baseline))
  end
  limit=c.card["baseline_validation_max_error_V"]
  @assert 0<limit<=.005*k.gain
  @assert baseline_error<limit
  (;k...,noise,enc=noise/k.gain,baseline_error)
 end
end
end
