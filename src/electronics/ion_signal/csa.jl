const E_CHARGE=1.602176634e-19
function interpolate(t,y,x)
    x<=first(t)&&return first(y);x>=last(t)&&return last(y)
    j=searchsortedlast(t,x);a=(x-t[j])/(t[j+1]-t[j]);(1-a)*y[j]+a*y[j+1]
end
function primitive(t,y)
    sumq=vcat(0.,cumsum(diff(t).*(y[1:end-1].+y[2:end])./2))
    function integral(x)
        x<=first(t)&&return 0.;x>=last(t)&&return last(sumq)
        j=searchsortedlast(t,x);d=x-t[j]
        sumq[j]+d*y[j]+d*d*(y[j+1]-y[j])/(2(t[j+1]-t[j]))
    end
end
"""Single-ion current averaged over one uniform-arrival bucket."""
struct PreparedBucketCurrent
    time_s::Vector{Float64}
    current_A::Vector{Float64}
end

"""Integrate the single-ion pulse once for a given collection-bucket duration."""
function prepare_bucket_current(p::IonCurrentPulse,bucket_ms)
    all(diff(p.time_s).>0)&&all(p.current_A.>=0)||throw(ArgumentError("Ion pulse"))
    dt=Float64(bucket_ms)*.001;charge=primitive(p.time_s,p.current_A)
    ut=sort!(unique(vcat(p.time_s,dt.+p.time_s)))
    PreparedBucketCurrent(ut,[(charge(x)-charge(x-dt))/dt for x in ut])
end

"""Build a pad's current using a prepared uniform-arrival packet. Arrays are owned."""
function prepared_ion_current(counts,t,p::PreparedBucketCurrent,a)
    length(counts)==t.bucket_count && all(counts.>=0)||throw(ArgumentError("Ion bucket counts"))
    dt=Float64(t.bucket_ms)*.001;offset=(Float64(t.start_ms)-a.start_ms)*.001
    stop=(a.sample_count-1)*Float64(a.sample_us)*1e-6
    active=findall(>(0),counts);isempty(active)&&return IonCurrentPulse([0.,stop],[0.,0.])
    ut=p.time_s;ui=p.current_A
    times=sort!(unique(vcat(0.,stop,[offset+(j-1)*dt+x for j in active for x in ut])))
    first(times)>=-1e-10 && last(times)<=stop+1e-10||throw(ArgumentError("ADC interval misses current"))
    # Only overlapping bucket pulses can contribute. Keep their original order.
    vals=zeros(length(times))
    for (k,x) in enumerate(times)
        low=max(1,floor(Int,(x-offset-last(ut))/dt)+1)
        high=min(length(counts),ceil(Int,(x-offset-first(ut))/dt)+1)
        value=0.
        for j in low:high
            counts[j]==0&&continue
            u=x-offset-(j-1)*dt
            first(ut)<=u<=last(ut)&&(value+=counts[j]*interpolate(ut,ui,u))
        end
        vals[k]=value
    end
    IonCurrentPulse(times,vals)
end
"""Common current knots and interpolation weights for every pad in one event.
Each row stores only overlapping bucket contributions, in bucket order.
"""
struct EventCurrentPlan
    time_s::Vector{Float64}
    bucket_rows::Vector{Vector{Int}}
    selected::BitVector
    pulse::IonCurrentPulse
    row_start::Vector{Int}
    bucket::Vector{Int}
    weight::Vector{Float64}
end

"""Prepare all bucket shifts once. The returned current buffer is reused per pad."""
function prepare_event_current(t,p::PreparedBucketCurrent,a)
    dt=Float64(t.bucket_ms)*.001;offset=(Float64(t.start_ms)-a.start_ms)*.001
    max_rows=Int128(2)+Int128(t.bucket_count)*length(p.time_s)
    overlaps=ceil(Int,(last(p.time_s)-first(p.time_s))/dt)+3
    max_rows*(160+16overlaps)<=64*1024^2||throw(ArgumentError("Current preparation exceeds 64 MiB budget"))
    stop=(a.sample_count-1)*Float64(a.sample_us)*1e-6
    times=sort!(unique(vcat(0.,stop,[offset+(j-1)*dt+x for j in 1:t.bucket_count for x in p.time_s])))
    first(times)>=-1e-10 && last(times)<=stop+1e-10||throw(ArgumentError("ADC interval misses current"))
    starts=Int[1];buckets=Int[];weights=Float64[]
    for x in times
        low=max(1,floor(Int,(x-offset-last(p.time_s))/dt)+1)
        high=min(t.bucket_count,ceil(Int,(x-offset-first(p.time_s))/dt)+1)
        for j in low:high
            u=x-offset-(j-1)*dt
            if first(p.time_s)<=u<=last(p.time_s)
                push!(buckets,j);push!(weights,interpolate(p.time_s,p.current_A,u))
            end
        end
        push!(starts,length(buckets)+1)
    end
    row_lookup=Dict(x=>k for (k,x) in enumerate(times))
    bucket_rows=[[row_lookup[offset+(j-1)*dt+x] for x in p.time_s] for j in 1:t.bucket_count]
    EventCurrentPlan(times,bucket_rows,falses(length(times)),IonCurrentPulse(Float64[],Float64[]),starts,buckets,weights)
end

"""Fill the shared current buffer from one pad's ion counts; consume before reuse."""
function event_current!(p::EventCurrentPlan,counts)
    fill!(p.selected,false);p.selected[1]=true;p.selected[end]=true
    for j in eachindex(counts)
        if counts[j]>0
            for k in p.bucket_rows[j];p.selected[k]=true;end
        end
    end
    n=count(p.selected);resize!(p.pulse.time_s,n);resize!(p.pulse.current_A,n);i=0
    for k in eachindex(p.time_s)
        p.selected[k]||continue
        value=0.
        for q in p.row_start[k]:p.row_start[k+1]-1
            value+=counts[p.bucket[q]]*p.weight[q]
        end
        i+=1;p.pulse.time_s[i]=p.time_s[k];p.pulse.current_A[i]=value
    end
    p.pulse
end

"""
    ion_current(counts, time, pulse, adc_time)

Average the single-ion current over uniform collection times in each input
bucket, as in the existing CSA prescription. Sum shifted packet currents and
return local seconds from the first ADC sample. Current knots include all pulse
edges shifted by bucket starts and ends; each interval is linear in the result.
"""
function ion_current(counts::AbstractVector{Int64},t::IonReadoutTime,p::IonCurrentPulse,a::IonADCTime)
    prepared_ion_current(counts,t,prepare_bucket_current(p,t.bucket_ms),a)
end

"""Constant gain and decay time recovered from an exponential response table."""
struct ExponentialCSA
    gain_V_C::Float64
    tau_s::Float64
    baseline_V::Float64
end
function exponential_parameters(s::LinearCSAResponse)
    t=s.time_s;y=s.voltage_per_charge_V_C
    length(t)==length(y)>=2 && first(t)==0 && all(diff(t).>0)||throw(ArgumentError("Linear response table"))
    all(<(0),y)||return nothing
    g=-first(y);rate=-log(last(y)/first(y))/last(t)
    isfinite(rate)&&rate>0||return nothing
    maximum(abs.(y./(-g.*exp.(-rate.*t)).-1))<=1e-10||return nothing
    ExponentialCSA(g,1/rate,s.baseline_V)
end

"""
    exponential_voltage!(voltage, current, sample_times, model)

Exact exponential response to each straight segment of the supplied current.
Split intervals at current knots and requested sample times. Small h/tau uses
series for the two integrated-current weights to avoid cancellation. Output is
voltage in V; all integration times are local seconds. No intermediate sampling
grid, current renormalisation or approximation of narrow current pulses is used.
"""
function exponential_voltage!(out,c,ts,s::ExponentialCSA)
    length(out)==length(ts)||throw(DimensionMismatch("Voltage samples"))
    first(ts)>=first(c.time_s)&&last(ts)<=last(c.time_s)+1e-12||throw(ArgumentError("Current interval"))
    now=first(c.time_s);x=0.;j=1
    for (k,target) in enumerate(ts)
        target>=now||throw(ArgumentError("Sorted sample times required"))
        while now<target
            while j<length(c.time_s)-1 && c.time_s[j+1]<=now;j+=1;end
            stop=min(target,c.time_s[j+1]);h=stop-now
            h>0||error("Current integration did not advance")
            slope=(c.current_A[j+1]-c.current_A[j])/(c.time_s[j+1]-c.time_s[j])
            i0=c.current_A[j]+slope*(now-c.time_s[j]);di=slope*h;u=h/s.tau_s
            a=exp(-u)
            if u<1e-3
                A=1+u*(-.5+u*(1/6+u*(-1/24+u/120)))
                B=.5+u*(-1/6+u*(1/24+u*(-1/120+u/720)))
            else
                A=-expm1(-u)/u;B=(1-A)/u
            end
            x=a*x+s.gain_V_C*h*(i0*A+di*B);now=stop
        end
        out[k]=s.baseline_V-x
    end
    out
end

"""
    linear_csa_voltage(current, sample_time_s, model)

Return the linear CSA voltage at the requested times. An exponential response
uses exact integration over straight current segments. Other response tables
use the convolution reference. Production prepares the exponential parameters
once and reuses output arrays through the private workspace.
"""
function linear_csa_voltage(c::IonCurrentPulse,ts::Vector{Float64},s::LinearCSAResponse)
    model=exponential_parameters(s)
    model===nothing && return convolution_voltage_reference(c,ts,s)
    last(s.time_s)>=last(ts)-first(c.time_s)||throw(ArgumentError("Linear calibration duration"))
    exponential_voltage!(Vector{Float64}(undef,length(ts)),c,ts,model)
end

"""Prepared convolution for a uniform ADC grid and a supplied response table.

The kernel transform and FFT arrays are reusable. Refinement divides each ADC
interval for a convergence check against the exact exponential calculation.
"""
struct PreparedCSAConvolution{F,I}
    refinement::Int
    h::Float64
    n::Int
    charge::Vector{Float64}
    spectrum::Vector{ComplexF64}
    kernel::Vector{ComplexF64}
    voltage::Vector{Float64}
    forward::F
    inverse::I
end
function prepare_csa_convolution(ts,s;refinement=4)
    length(ts)>=2 && all(diff(ts).>0)||throw(ArgumentError("Sample times"))
    refinement>=1||throw(ArgumentError("Convolution refinement"))
    step=ts[2]-ts[1]
    maximum(abs.(diff(ts).-step))<1e-9*step||throw(ArgumentError("Uniform samples required"))
    h=step/refinement;n=refinement*(length(ts)-1);m=nextpow(2,2n)
    charge=zeros(m);kernel=zeros(m)
    for j in 1:n;kernel[j]=interpolate(s.time_s,s.voltage_per_charge_V_C,(j-.5)*h);end
    spectrum=zeros(ComplexF64,m÷2+1)
    PreparedCSAConvolution(refinement,h,n,charge,spectrum,rfft(kernel),zeros(m),
       plan_rfft(charge;flags=FFTW.ESTIMATE),plan_irfft(spectrum,m;flags=FFTW.ESTIMATE))
end
function convolution_voltage!(out,c,ts,s,p::PreparedCSAConvolution)
    first(ts)>=first(c.time_s)&&last(ts)<=last(c.time_s)+1e-12||throw(ArgumentError("Current interval"))
    last(s.time_s)>=last(ts)-first(c.time_s)||throw(ArgumentError("Linear calibration duration"))
    q=primitive(c.time_s,c.current_A);fill!(p.charge,0.)
    for j in 1:p.n;p.charge[j]=q(first(ts)+j*p.h)-q(first(ts)+(j-1)*p.h);end
    mul!(p.spectrum,p.forward,p.charge)
    p.spectrum .*= p.kernel
    mul!(p.voltage,p.inverse,p.spectrum)
    fill!(out,s.baseline_V)
    for j in 2:length(ts);out[j]+=p.voltage[p.refinement*(j-1)];end
    out
end
function convolution_voltage_reference(c::IonCurrentPulse,ts::Vector{Float64},s::LinearCSAResponse;refinement=4)
    p=prepare_csa_convolution(ts,s;refinement)
    convolution_voltage!(similar(ts),c,ts,s,p)
end
"""Prepare the selected pad's noiseless voltage; domain failures return to the caller."""
function pad_voltage(counts,t,s::IonElectronicsSettings,at)
    current=ion_current(counts,t,s.csa.current_pulse,at)
    ts=collect(0:at.sample_count-1).*(Float64(at.sample_us)*1e-6)
    v=all(iszero,counts) ? fill(s.csa.linear_response.baseline_V,length(ts)) : linear_csa_voltage(current,ts,s.csa.linear_response);mode=linear_csa
    if s.csa.mode==hybrid_nonlinear_csa && minimum(v)<s.csa.nonlinear_below_V
        out=nonlinear_waveform(s.csa.nonlinear_response,current.time_s,current.current_A;sample_time_s=ts)
        out.valid||throw(DomainError(minimum(v),"Nonlinear voltage calibration exceeded; stop and review this pad"))
        v=out.voltage_V;mode=hybrid_nonlinear_csa
    end
    (;current,time_s=ts,voltage_V=v,mode)
end
"""
    sim_ion_electronics_response(noise_rng, adc_rng, detector, settings)

Simulate every pad in the frame using the reviewed current, CSA, voltage-noise
and ADC interfaces. Return codes and simulation counters. Limit the output to
256 MiB before allocation. A nonlinear domain failure stops this calculation
and reports the pad for calibration review; no invalid codes are saved. Production
never launches SPICE. Peak-current exceedances are recorded by the run script.
"""
function sim_ion_electronics_response(noise_rng::AbstractRNG,adc_rng::AbstractRNG,
                                     d::IonDetectorResponse,s::IonElectronicsSettings)
    sim_ion_electronics_prepared(noise_rng,adc_rng,d,s,CSAWorkspace())
end

"""Run-owned arrays and tables. Nothing is shared between concurrent workers.

Noise entries are keyed by covariance array and required FFT period. Each entry
covers every lag requested by its event. The noise cache is limited to 64 MiB;
working voltage arrays are limited by the 256 MiB ADC output budget.
"""
mutable struct CSAWorkspace
    noise::Dict{Tuple{UInt,Int},Any}
    pulse::Any
    bucket_ms::Float32
    packet::Any
    linear::Any
    exponential::Any
    convolution::Any
    times::Vector{Float64}
    voltage::Vector{Float64}
    measured::Vector{Float64}
    peak_current_A::Vector{Float64}
end
CSAWorkspace()=CSAWorkspace(Dict{Tuple{UInt,Int},Any}(),nothing,0f0,nothing,nothing,nothing,nothing,Float64[],Float64[],Float64[],Float64[])

"""Simulate with run-owned preparation; public responses own their output arrays."""
function sim_ion_electronics_prepared(noise_rng,adc_rng,d,s,w::CSAWorkspace)
    v=d.ions;at=ion_adc_time(v.time,s.csa.current_pulse,s.acquisition,s.adc)
    at.sample_us==s.csa.noise.sample_us||throw(ArgumentError("CSA noise sampling differs"))
    n=at.sample_count;nc=size(v.counts,2)
    Int128(n)*nc*sizeof(UInt16)<=256*1024^2||throw(ArgumentError("ADC vectors exceed 256 MiB budget"))
    key=(objectid(s.csa.noise.covariance_V2),nextpow(2,max(2,2(n-1))))
    if !haskey(w.noise,key)
        pn=prepare_csa_noise(n,s.csa.noise)
        entry=(pn,s.csa.noise)
        bytes=Base.summarysize(entry)
        bytes<=64*1024^2||throw(ArgumentError("Prepared noise exceeds 64 MiB"))
        Base.summarysize(w.noise)+bytes>64*1024^2&&empty!(w.noise)
        # Retain the covariance array too, so its object identity cannot be recycled.
        w.noise[key]=entry
    end
    pn=w.noise[key][1];resize!(pn.output,n)
    if w.pulse!==s.csa.current_pulse || w.bucket_ms!=v.time.bucket_ms
        w.packet=prepare_bucket_current(s.csa.current_pulse,v.time.bucket_ms)
        w.pulse=s.csa.current_pulse;w.bucket_ms=v.time.bucket_ms
    end
    if w.linear!==s.csa.linear_response
        w.exponential=exponential_parameters(s.csa.linear_response);w.linear=s.csa.linear_response;w.convolution=nothing
    end
    if length(w.times)!=n || (n>1 && length(w.times)>1 && w.times[2]!=Float64(at.sample_us)*1e-6);w.convolution=nothing;end
    resize!(w.times,n);resize!(w.voltage,n);resize!(w.measured,n);resize!(w.peak_current_A,nc)
    for i in 1:n;w.times[i]=(i-1)*Float64(at.sample_us)*1e-6;end
    last(s.csa.linear_response.time_s)>=last(w.times)||throw(ArgumentError("Linear calibration duration"))
    fill!(w.peak_current_A,0.)
    current_plan=prepare_event_current(v.time,w.packet,at)
    out=Matrix{UInt16}(undef,n,nc);lin=0;nonlin=0;lo=0;hi=0
    for j in axes(v.counts,2)
        counts=view(v.counts,:,j);voltage=w.voltage
        if all(iszero,counts)
            fill!(voltage,s.csa.linear_response.baseline_V);lin+=1
        else
            current=event_current!(current_plan,counts)
            w.peak_current_A[j]=maximum(current.current_A)
            if w.exponential===nothing
                w.convolution===nothing && (w.convolution=prepare_csa_convolution(w.times,s.csa.linear_response))
                convolution_voltage!(voltage,current,w.times,s.csa.linear_response,w.convolution)
            else
                exponential_voltage!(voltage,current,w.times,w.exponential)
            end
            if s.csa.mode==hybrid_nonlinear_csa && minimum(voltage)<s.csa.nonlinear_below_V
                wave=nonlinear_waveform(s.csa.nonlinear_response,current.time_s,current.current_A;sample_time_s=w.times)
                wave.valid||throw(DomainError(minimum(voltage),"Nonlinear voltage calibration exceeded; stop and review this pad"))
                copyto!(voltage,wave.voltage_V);nonlin+=1
            else
                lin+=1
            end
        end
        noise=sample_csa_noise!(noise_rng,pn)
        for i in 1:n;w.measured[i]=voltage[i]+noise[i];end
        adc=digitise_adc!(view(out,:,j),adc_rng,w.measured,s.adc)
        lo+=adc.below_range;hi+=adc.above_range
    end
    (IonADCResponse(d.event_id,v.channels,at,out),IonElectronicsControls(lin,nonlin,0,lo,hi))
end
