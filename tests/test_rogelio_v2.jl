using Test,HDF5,TOML
root=normpath(joinpath(@__DIR__,".."))
@testset "Rogelio v2 physical and numerical checks" begin
 h5open(joinpath(root,"results/ltspice/rogelio_v2/analysis.h5"),"r") do h
  @test read(attributes(h["reproduction/original"])["max_relative_difference_supplied"])<1e-6
  for name in ("v2_native_corrected","v2_matched","v2_m3_large")
   g=h[name];a=attributes(g)
   @test read(a["gain_V_per_ion"])>0
   @test read(a["noise_rms_V"])>0
   @test maximum(abs.(read(g["charge_relative_error"])))<1e-8
   @test maximum(read(g["signal_refinement_error"]))<1e-3
   @test read(a["spectral_relative_error"])<1e-4
   @test read(a["quadrature_relative_error"])<1e-5
   @test abs(sum(read(g["device_variance_fraction"]))-1)<1e-3
   @test read(a["tail_fit_max_relative_error"])<1e-4
   @test abs(read(a["tau_amplitude_relative_change"]))<.01
   @test read(a["threshold_below_snr"])<5<read(a["threshold_above_snr"])
   counts=read(g["counts"]);n=read(a["verified_N_snr5"])
   @test (n-1 in counts) && (n in counts)
   gain=read(g["gain_by_count_V_per_ion"])
   @test abs(gain[findfirst(==(10),counts)]/first(gain)-1)<1e-3
  end
 end
 h5open(joinpath(root,"results/ltspice/rogelio_v2/raw.h5"),"r") do h
  for name in ("v2_native_corrected","v2_matched","v2_m3_large")
   # Verify that every transistor remains in the translated netlist.
   for device in ("M1","M2","M3","M4","M5","M6","Mf")
    @test haskey(h[name*"/op"],"Id($device)")
   end
  end
 end
 h5open(joinpath(root,"results/ltspice/rogelio_v2_retuned/analysis.h5"),"r") do h
  g=h["v2_retuned"];a=attributes(g)
  @test abs(read(a["baseline_V"])-0.878931543)<1e-6
  @test maximum(abs.(read(g["charge_relative_error"])))<1e-8
  @test maximum(read(g["signal_refinement_error"]))<1e-3
  @test read(a["spectral_relative_error"])<1e-4
  @test read(a["quadrature_relative_error"])<1e-5
  @test read(a["threshold_below_snr"])<5<read(a["threshold_above_snr"])
  @test abs(sum(read(g["device_variance_fraction"]))-1)<1e-3
 end
end
