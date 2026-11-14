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
    n=size(v.counts,1);M=s.illuminated_microcells;total=s.sipm.microcells
    lit=SiPMResponse.Sensor(s.sipm.pde,M,s.sipm.gain,s.sipm.recovery_s,s.sipm.dark_rate_Hz*M/total)
    unlit=M<total ? SiPMResponse.Sensor(0.,total-M,s.sipm.gain,s.sipm.recovery_s,s.sipm.dark_rate_Hz*(total-M)/total) : nothing
    lit_buffers=SiPMResponse.StatisticalBuffers(n);dark_buffers=SiPMResponse.StatisticalBuffers(n)
    empty_counts=zeros(Int64,n);charge_buffer=zeros(n);noise_buffer=zeros(n)
    for j in axes(v.counts,2)
        a=SiPMResponse.simulate_channel_statistical(Xoshiro(rand(sipm_rng,UInt64)),view(v.counts,:,j),dt,lit;options=s.statistical,buffers=lit_buffers)
        detected+=sum(a.detected_pe);dark+=sum(a.dark_pulses)
        if unlit===nothing
            copyto!(charge_buffer,a.charge_C)
        else
            b=SiPMResponse.simulate_channel_statistical(Xoshiro(rand(sipm_rng,UInt64)),empty_counts,dt,unlit;options=s.statistical,buffers=dark_buffers)
            charge_buffer .= a.charge_C .+ b.charge_C
            dark+=sum(b.dark_pulses)
        end
        charge+=sum(charge_buffer)/q0
        randn!(asic_rng,noise_buffer)
        channel_loss=0.
        for i in 1:n
            noise=s.asic.noise_rms_pe*noise_buffer[i]
            current=(charge_buffer[i]+q0*noise)/s.asic.bucket_width_s
            limited=min(current,s.asic.current_limit_A)
            response[i,j]=Float32(limited*(dt/q0))
            clipped[i,j]=current>s.asic.current_limit_A
            nclip+=clipped[i,j];channel_loss+=current-limited
        end
        lost+=channel_loss*(dt/q0)
    end
    (SiPMASICResponse(detector.event_id,v.channels,v.time,response,clipped),
     SiPMReadoutControls(detected,dark,charge,nclip,lost))
end
