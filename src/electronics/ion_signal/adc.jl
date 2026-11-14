function check_adc(s)
    1<=s.bits<=16 && isfinite(s.minimum_V) && isfinite(s.maximum_V) &&
    s.maximum_V>s.minimum_V && isfinite(s.sample_us) && s.sample_us>0 &&
    isfinite(s.input_noise_rms_V) && s.input_noise_rms_V>=0 || throw(ArgumentError("Invalid ADC settings"))
end
"""
    digitise_adc(rng, voltage_V, settings)

Add independent ADC input noise, then quantise to UInt16 with the reviewed
floor convention. Return counts and below/above-range sample totals. CSA noise
must already be present in voltage_V. Zero ADC noise consumes no random numbers.
"""
function digitise_adc(rng::AbstractRNG,voltage_V::AbstractVector,s::ADCSettings)
    out=Vector{UInt16}(undef,length(voltage_V))
    digitise_adc!(out,rng,voltage_V,s)
end

"""Write ADC counts into a supplied channel column; return range counters."""
function digitise_adc!(out,rng,voltage_V,s)
    length(out)==length(voltage_V)||throw(DimensionMismatch("ADC output"))
    check_adc(s);all(isfinite,voltage_V)||throw(ArgumentError("Nonfinite voltage"))
    width=(Float64(s.maximum_V)-s.minimum_V)/2^s.bits
    lo=0;hi=0
    for j in eachindex(voltage_V)
        v=Float64(voltage_V[j])+(s.input_noise_rms_V==0 ? 0. : s.input_noise_rms_V*randn(rng))
        lo+=v<s.minimum_V;hi+=v>s.maximum_V
        out[j]=UInt16(floor(Int,clamp((v-s.minimum_V)/width,0.,Float64(2^s.bits-1))))
    end
    (;counts=out,below_range=Int64(lo),above_range=Int64(hi))
end
"""
    ion_adc_time(time, pulse, acquisition, adc)

Cover all uniform within-bucket arrivals, the complete single-ion current pulse,
and requested baseline/tail. The returned start uses the common event t0.
"""
function ion_adc_time(t::IonReadoutTime,pulse::IonCurrentPulse,a::IonAcquisitionSettings,adc::ADCSettings)
    check_adc(adc)
    all(x->isfinite(x)&&x>=0,(a.baseline_ms,a.signal_tail_ms))||throw(ArgumentError("Invalid baseline/tail"))
    first(pulse.time_s)<=last(pulse.time_s)<=0||throw(ArgumentError("Single-ion pulse ends at collection"))
    start=Float32(Float64(t.start_ms)+1000first(pulse.time_s)-a.baseline_ms)
    stop=Float64(t.start_ms)+Float64(t.bucket_ms)*t.bucket_count+a.signal_tail_ms
    n=ceil(Int64,(stop-Float64(start))*1000/adc.sample_us)+1
    IonADCTime(start,adc.sample_us,n)
end
