"""
    sim_ion_detector_response(rng,hits,e_ion_pairs,placement,frame,settings;
                              max_vector_bytes=256*1024^2)

Diffuse the supplied ions and assign them to physical pads and t0-based buckets.
Read every channel overlapping the frame in full. Record all diffusion outside
those channels or the complete bucket interval, then exclude it from vectors.
The same pair vector is also used by sim_electron_detector_response. No
attachment, electronics, reconstruction or zero suppression is applied here.
"""
function sim_ion_detector_response(rng::AbstractRNG,h::MCEventHits,pairs::Vector{Int64},
        p::EventPlacement,f::EventFrame,s::DetectorResponseSettings;max_vector_bytes::Int=256*1024^2)::IonDetectorResponse
    _check(h,pairs,p,f,s)
    tr=s.ion_transport;g=s.ion_grid
    t=ion_readout_time(f,p,s.chamber_length_mm,tr.ion_drift_velocity_mm_ms,s.ion_bucket_ms)
    ch=readout_channels(f,g)
    cells=_collect(rng,h,pairs,p,g,t,tr.ion_drift_velocity_mm_ms,tr.sigma_t_mm,tr.sigma_l_mm,s.chamber_length_mm,true)
    v,loss,excluded=_vectors(cells,ch,t;max_vector_bytes)
    IonDetectorResponse(h.event_id,p,f,v,sum(pairs),loss.excluded_total,loss,excluded)
end
