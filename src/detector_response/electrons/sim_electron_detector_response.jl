"""
    sim_electron_detector_response(rng,hits,e_ion_pairs,placement,frame,settings;
                                   max_vector_bytes=256*1024^2)

Apply the configured attachment at the event drift distance, diffuse surviving
electrons, and calculate EL and incident photons for each physical hole/bucket.
Reuse sample_attached_electrons, sample_el_photons and sample_incident_photons.
Compute photons in excluded cells as well to measure optical signal losses.
Return blue or VUV photon vectors before PDE. Photon time is the parent arrival
bucket. Simulation does not correct attachment; later signal reconstruction
uses exp(t/lifetime), leaving the attachment fluctuations in the reconstructed PE.
"""
function sim_electron_detector_response(rng::AbstractRNG,h::MCEventHits,pairs::Vector{Int64},
        p::EventPlacement,f::EventFrame,s::DetectorResponseSettings;max_vector_bytes::Int=256*1024^2)::ElectronDetectorResponse
    _check(h,pairs,p,f,s)
    tr=s.electron_transport;g=s.electron_grid
    t=electron_readout_time(f,p,s.chamber_length_mm,tr.electron_drift_velocity_mm_us,s.electron_bucket_us)
    ch=readout_channels(f,g)
    drift_ms=10f0*p.electron_drift_cm/tr.electron_drift_velocity_mm_us/1000f0
    surviving=sample_attached_electrons(rng,pairs,drift_ms,s.electron_lifetime_ms)
    cells=_collect(rng,h,surviving,p,g,t,tr.electron_drift_velocity_mm_us,tr.sigma_t_mm,tr.sigma_l_mm,s.chamber_length_mm,false)
    photons=Dict{NTuple{3,Int64},Int64}();emitted=Int64(0);incident=Int64(0);eout=Int64(0)
    selected=Set(zip(Int64.(ch.ix),Int64.(ch.iy)))
    for key in sort!(collect(keys(cells)))
        n=sample_el_photons(rng,cells[key],s.el_photons_per_electron,s.el_variance_ratio)
        m=sample_incident_photons(rng,n,s.optical_acceptance)
        emitted+=n;incident+=m
        m>0&&(photons[key]=m)
        (!(key[1:2] in selected)||!(1<=key[3]<=t.bucket_count)) && (eout+=cells[key])
    end
    v,loss,excluded=_vectors(photons,ch,t;max_vector_bytes)
    ElectronDetectorResponse(h.event_id,p,f,v,sum(pairs),sum(pairs)-sum(surviving),eout,
        emitted,emitted-incident,loss,excluded)
end
