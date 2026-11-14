"""
    sim_sipm_readout(sipm_rng, asic_rng, detector, settings)

Apply the existing statistical illuminated-microcell model and ASIC response
to every SiPM in the frame. Photons retain their original time buckets; dark
activity covers the whole device. Return PE-equivalent measured charge and
simulation totals. Reconstruction, attachment correction and rebinning follow.
"""
function sim_sipm_readout(sipm_rng::AbstractRNG,asic_rng::AbstractRNG,
                         detector::ElectronDetectorResponse,s::SiPMReadoutSettings)
    v=detector.photons;dt=Float64(v.time.bucket_us)*1e-6
    isapprox(dt,s.asic.bucket_width_s;rtol=1e-6)||throw(ArgumentError("ASIC/photon bucket mismatch"))
    0<s.illuminated_microcells<=s.sipm.microcells||throw(ArgumentError("Illuminated microcells"))
    sizeof(v.counts)<=256*1024^2||throw(ArgumentError("Response vector exceeds 256 MiB budget"))
    response=Matrix{Float32}(undef,size(v.counts));clipped=falses(size(v.counts))
    detected=Int64(0);dark=Int64(0);charge=0.;lost=0.;nclip=Int64(0)
    q0=SiPMResponse.charge_per_pe_C(s.sipm)
    for j in axes(v.counts,2)
        a=SiPMResponse.simulate_channel_illuminated(sipm_rng,view(v.counts,:,j),dt,s.sipm,s.illuminated_microcells;options=s.statistical)
        b=SiPMElectronics.apply_electronics(asic_rng,a.charge_C,q0,s.asic)
        response[:,j]=Float32.(b.current_A.*(dt/q0));clipped[:,j]=b.saturated
        detected+=sum(a.detected_pe);dark+=sum(a.dark_pulses);charge+=sum(a.charge_C)/q0
        lost+=sum(b.current_unclipped_A.-b.current_A)*(dt/q0);nclip+=count(b.saturated)
    end
    (SiPMASICResponse(detector.event_id,v.channels,v.time,response,clipped),
     SiPMReadoutControls(detected,dark,charge,nclip,lost))
end
