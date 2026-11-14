"""Positive reconstructed hits in vertex-relative mm, with charge and energy in keV."""
struct ReconstructedHits
 x_mm::Vector{Float32}
 y_mm::Vector{Float32}
 z_mm::Vector{Float32}
 charge::Vector{Float32}
 energy_keV::Vector{Float32}
end

"""Sum signed vectors before hit selection. Correct PE for mean attachment.
Return NaN for a complete-event estimate when any input bucket is invalid;
retain the valid-bucket sum separately. Energy uses W for ions and
W/(PDE × EL yield × optical acceptance) for attachment-corrected PE.
"""
function event_charge_estimates(ions::ReconstructedIons,pe::ReconstructedPE,placement,detector,pde)
 correction=exp(10Float64(placement.electron_drift_cm)/detector.electron_transport.electron_drift_velocity_mm_us/1000/detector.electron_lifetime_ms)
 ion_scale=Float64(detector.W_eV)/1000
 pe_scale=ion_scale/(pde*detector.el_photons_per_electron*detector.optical_acceptance)
 iv=ions.valid .& isfinite.(ions.ions);ev=pe.valid .& isfinite.(pe.pe)
 qi=sum(Float64,ions.ions[iv];init=0.);qe=correction*sum(Float64,pe.pe[ev];init=0.)
 (;attachment_correction=correction,ion_keV_per_count=ion_scale,pe_keV_per_count=pe_scale,
 ions_valid_sum=qi,pe_valid_sum=qe,ions_total=all(iv) ? qi : NaN,pe_total=all(ev) ? qe : NaN,
 ion_energy_keV=all(iv) ? qi*ion_scale : NaN,electron_energy_keV=all(ev) ? qe*pe_scale : NaN,
 invalid_ion_buckets=count(.!iv),invalid_pe_buckets=count(.!ev))
end

"""Convert valid Q > 0 buckets to hits. Both tracks use vertex-relative Z.
For electrons, use the already rebinned signal and the actual final bucket end.
The PE correction is a constant for the event, matching the transport model.
"""
function reconstruct_hits(signal,placement,detector;charge_correction=1.,keV_per_count,minimum_energy_keV=0.)
 ion=signal isa ReconstructedIons
 values=ion ? signal.ions : signal.pe
 selected=findall(signal.valid .& isfinite.(values) .& (values.>0))
 x=Float32[];y=Float32[];z=Float32[];q=Float32[];energy=Float32[]
 vertex_z=10Float64(placement.ion_drift_cm)
 for index in selected
  k,j=Tuple(index)
  if ion
   time=Float64(signal.time.start_ms)+(k-.5)*signal.time.bucket_ms
   zd=time*detector.ion_transport.ion_drift_velocity_mm_ms
  else
   lo=Float64(signal.time.start_us)+(k-1)*signal.time.bucket_us
   hi=min(lo+signal.time.bucket_us,Float64(signal.acquired_end_us))
   zd=detector.chamber_length_mm-(lo+hi)/2*detector.electron_transport.electron_drift_velocity_mm_us
  end
  charge=Float64(values[index])*charge_correction
  charge*keV_per_count>minimum_energy_keV||continue
  push!(x,signal.channels.x_mm[j]);push!(y,signal.channels.y_mm[j]);push!(z,zd-vertex_z)
  push!(q,charge);push!(energy,charge*keV_per_count)
 end
 ReconstructedHits(x,y,z,q,energy)
end

"""Z centres of all acquired buckets, including empty ones, in vertex-relative mm."""
function hit_z_centres(signal,placement,detector)
 z0=10Float64(placement.ion_drift_cm)
 if signal isa ReconstructedIons
  return Float32[(Float64(signal.time.start_ms)+(k-.5)*signal.time.bucket_ms)*detector.ion_transport.ion_drift_velocity_mm_ms-z0 for k in 1:signal.time.bucket_count]
 end
 out=Float32[]
 for k in 1:signal.time.bucket_count
  lo=Float64(signal.time.start_us)+(k-1)*signal.time.bucket_us
  hi=min(lo+signal.time.bucket_us,Float64(signal.acquired_end_us))
  push!(out,detector.chamber_length_mm-(lo+hi)/2*detector.electron_transport.electron_drift_velocity_mm_us-z0)
 end
 out
end
