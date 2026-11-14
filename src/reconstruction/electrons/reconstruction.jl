"""
    reco_electron_signal(response, settings) -> ReconstructedPE

Correct the measured ASIC charge using the stationary mean SiPM inverse at the
original bucket duration. The illuminated microcell count sets the recovery
correction. Subtract the full device's stationary mean dark charge first when
requested; fluctuations and negative estimates are retained. This mean inverse
is validated against changing illumination for the pilot's charge range.
Clipped samples and nonpositive inverse denominators return NaN/valid=false.
The function uses no Monte Carlo counts or simulated microcell states.
"""
function reco_electron_signal(r::SiPMASICResponse,s::ElectronReconstructionSettings)
 0<s.illuminated_microcells<=s.sipm.microcells||throw(ArgumentError("Illuminated cells"))
 dt=Float64(r.time.bucket_us)*1e-6
 dt>0||throw(ArgumentError("Bucket duration"))
 dark=s.subtract_mean_dark ? s.sipm.dark_rate_Hz*dt/(1+s.sipm.dark_rate_Hz*s.sipm.recovery_s/s.sipm.microcells) : 0.
 alpha=s.sipm.recovery_s/(s.illuminated_microcells*dt)
 out=fill(NaN32,size(r.response_pe));valid=falses(size(out))
 for j in eachindex(out)
  m=Float64(r.response_pe[j])-dark;den=1-alpha*m
  if !r.clipped[j] && isfinite(m) && den>0
   out[j]=m/den;valid[j]=true
  end
 end
 stop=Float32(Float64(r.time.start_us)+Float64(r.time.bucket_us)*r.time.bucket_count)
 ReconstructedPE(r.event_id,r.channels,r.time,stop,out,valid,copy(r.clipped))
end
