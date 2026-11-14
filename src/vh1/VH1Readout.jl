"""VH1 ion readout: reset two-capacitor CSA with one integrating conversion per frame.

An independent realisation beside ElectronicsResponse and SignalReconstruction,
reusing their current, noise, ADC and least-squares pieces without changing them.
Physics and numbers follow notes/physics/nausicaa_VH.tex, Sections 7 and 13:
the frame reading is V_REF minus the gain per ion times the accumulated charge,
with no discharge and no nonlinear branch. Load DetectorResponse,
ElectronicsResponse and SignalReconstruction first.
"""
module VH1Readout
using Random, LinearAlgebra, Statistics, TOML, SHA, DelimitedFiles
using ..DetectorResponse, ..ElectronicsResponse, ..SignalReconstruction
const ER=ElectronicsResponse
const SR=SignalReconstruction

"""Frames of equal length; start_ms is measured from the common event t0. The
first `reference_frames` frames end before any induced current."""
struct VH1FrameTime
    start_ms::Float32
    frame_ms::Float32
    frame_count::Int32
    reference_frames::Int32
end

"""Circuit, converter and noise description of one VH1 readout.

gain_V_per_ion is the output change per collected ion, the charge fraction
times e/Cf. Each frame is read conversions_per_frame times, at the ends of
equal sub-intervals; each reading is the mean of the output over the last
window_ms of its sub-interval (a sigma-delta integrating that long). One
conversion with window_ms equal to frame_ms averages the whole frame; a window
equal to noise_sample_us is a point sample. Noise is generated at noise_sample_us from the
SPICE spectrum and averaged like the signal. The
window gives the output voltages inside which the amplifier is linear; frames
outside it are flagged, not corrected. provenance records the SPICE sources.
"""
struct VH1Settings
    gain_V_per_ion::Float64
    charge_fraction::Float64
    feedback_capacitance_F::Float64
    v_ref_V::Float64
    window_V::Tuple{Float64,Float64}
    frame_ms::Float32
    conversions_per_frame::Int64
    window_ms::Float32
    reference_frames::Int64
    adc::ADCSettings
    noise_sample_us::Float32
    noise::CSAVoltageNoise
    current_pulse::IonCurrentPulse
    provenance::Dict{String,String}
end

"""Converter codes for one event: reading rows (conversions_per_frame per
frame, in time order) and pad columns, with the noiseless window check per pad. saturated marks pads whose noiseless output left the
linear window in at least one frame."""
struct VH1FrameResponse
    event_id::Int64
    channels::ReadoutChannels
    time::VH1FrameTime
    codes::Matrix{UInt16}
    saturated::BitVector
    adc_below_range::Int64
    adc_above_range::Int64
end

"""Prepared response matrix, covariance and whitened inverse for one event timing."""
struct VH1ReconstructionPlan
    frames::VH1FrameTime
    buckets::IonReadoutTime
    settings::VH1Settings
    response::Matrix{Float64}
    covariance::Matrix{Float64}
    lower::Matrix{Float64}
    inverse::Matrix{Float64}
    ion_covariance::Matrix{Float64}
end

sha(path)=open(io->bytes2hex(sha256(io)),path)

"""Single-ion current of the existing prescription: central 0.5 mm square pad in
a grounded plane, constant speed through the final gap, finite-path charge
retained. Same formula as scripts/production/framed_electronics_settings.jl."""
function single_ion_pulse(pad_side_m,gap_m,speed_m_s,step_s)
    a=pad_side_m/2;g=gap_m;v=speed_m_s
    t=collect(range(-g/v,0.;length=ceil(Int,g/v/step_s)+1))
    current=[ER.E_CHARGE*v*4a^2/(pi*(z^2+a^2)*sqrt(z^2+2a^2)) for z in -v.*t]
    IonCurrentPulse(t,current)
end

"""
    load_vh1_settings(root, cardpath; minimum_samples=4096)

Read a VH1 card, verify the recorded hashes of its SPICE inputs, build the
noise covariance at the internal sampling and return VH1Settings. The gain per
ion must agree with charge_fraction*e/Cf to 1 %.
"""
function load_vh1_settings(root,cardpath;minimum_samples=4096)
    c=TOML.parsefile(cardpath)
    get(c,"schema","")=="nausicaa.vh1_csa_spec.v1"||throw(ArgumentError("VH1 card schema"))
    csv=joinpath(root,c["noise_spectrum_csv"]);summary=joinpath(root,c["spice_summary"])
    hashes=Dict("noise_spectrum_sha256"=>sha(csv),"spice_summary_sha256"=>sha(summary))
    for (k,v) in hashes
        c[k]==v||throw(ArgumentError("$k differs from the card: rerun the VH1 SPICE study or update the card"))
    end
    data=readdlm(csv,',';skipstart=1)
    f=Float64.(data[:,1]);asd=Float64.(data[:,2])
    f[1]<=c["noise_low_Hz"] && f[end]>=c["noise_high_Hz"]||throw(ArgumentError("Spectrum does not cover the card band"))
    gain=Float64(c["gain_V_per_ion"]);fraction=Float64(c["charge_fraction"]);cf=Float64(c["feedback_capacitance_F"])
    abs(gain/(fraction*ER.E_CHARGE/cf)-1)<.01||throw(ArgumentError("gain_V_per_ion inconsistent with charge_fraction*e/Cf"))
    frame_ms=Float32(c["frame_ms"]);sample_us=Float32(c["noise_sample_us"]);window_ms=Float32(c["conversion_window_ms"])
    abs(frame_ms*1000/sample_us-round(frame_ms*1000/sample_us))<1e-6||throw(ArgumentError("frame_ms must be a multiple of noise_sample_us"))
    abs(window_ms*1000/sample_us-round(window_ms*1000/sample_us))<1e-6||throw(ArgumentError("conversion_window_ms must be a multiple of noise_sample_us"))
    k=Int64(c["conversions_per_frame"]);k>=1||throw(ArgumentError("conversions_per_frame must be at least 1"))
    m=round(Int,frame_ms*1000/sample_us);m%k==0||throw(ArgumentError("conversions_per_frame must divide the samples per frame"))
    0<window_ms<=frame_ms/k||throw(ArgumentError("conversion_window_ms must lie in (0, frame_ms/conversions_per_frame]"))
    adc=ADCSettings(c["adc_bits"],Float32(c["adc_minimum_V"]),Float32(c["adc_maximum_V"]),Float32(frame_ms*1000),Float32(c["adc_input_noise_rms_V"]))
    noise=ER.voltage_noise_from_spectrum(f,asd,sample_us,minimum_samples)
    pulse=single_ion_pulse(c["pad_side_m"],c["gap_m"],c["gap_speed_m_s"],c["current_step_s"])
    window=(Float64(c["output_window_V"][1]),Float64(c["output_window_V"][2]))
    window[1]<c["v_ref_V"]<=window[2]||throw(ArgumentError("v_ref_V must lie inside the output window"))
    provenance=merge(hashes,Dict("card"=>relpath(cardpath,root),"card_sha256"=>sha(cardpath),"spice_variant"=>String(c["spice_variant"]),"label"=>String(c["label"])))
    VH1Settings(gain,fraction,cf,Float64(c["v_ref_V"]),window,frame_ms,k,window_ms,Int64(c["reference_frames"]),adc,sample_us,noise,pulse,provenance)
end

"""Frames covering the reference period, every collection bucket and the end of
the last single-ion pulse. The reference frames end where the earliest induced
current begins."""
function vh1_frame_time(t::IonReadoutTime,s::VH1Settings)
    onset=Float64(t.start_ms)+1000first(s.current_pulse.time_s)
    start=onset-s.reference_frames*Float64(s.frame_ms)
    stop=Float64(t.start_ms)+Float64(t.bucket_ms)*t.bucket_count
    n=s.reference_frames+ceil(Int,(stop-onset)/Float64(s.frame_ms)-1e-9)
    VH1FrameTime(Float32(start),s.frame_ms,Int32(n),Int32(s.reference_frames))
end

samples_per_frame(s::VH1Settings)=round(Int,Float64(s.frame_ms)*1000/s.noise_sample_us)
window_samples(s::VH1Settings)=round(Int,Float64(s.window_ms)*1000/s.noise_sample_us)
readings(ft::VH1FrameTime,s::VH1Settings)=Int(ft.frame_count)*s.conversions_per_frame
"""Internal sample indices averaged by reading i: the last window of its sub-interval."""
function window_range(i,s::VH1Settings)
    m=samples_per_frame(s);k=s.conversions_per_frame;sub=m÷k
    n=(i-1)÷k;j=(i-1)%k+1;stop=n*m+j*sub
    (stop-window_samples(s)+1):stop
end

"""Internal sampling grid of the acquisition as an IonADCTime, for the current builder."""
function internal_grid(ft::VH1FrameTime,s::VH1Settings)
    n=samples_per_frame(s)*Int(ft.frame_count)
    IonADCTime(ft.start_ms,s.noise_sample_us,n+1)
end

"""Noiseless frame readings of one pad, in V, from its piecewise-linear current:
the mean of the output at the internal cell midpoints of the conversion window
at the end of each frame. The same operation is used for the unit-ion response,
so the reconstruction inverts exactly what is simulated."""
function frame_readings!(out,current::IonCurrentPulse,ft::VH1FrameTime,s::VH1Settings)
    dt=Float64(s.noise_sample_us)*1e-6
    q=ER.primitive(current.time_s,current.current_A);k=s.gain_V_per_ion/ER.E_CHARGE
    for i in 1:readings(ft,s)
        r=window_range(i,s);acc=0.
        for j in r;acc+=q((j-.5)*dt);end
        out[i]=s.v_ref_V-k*acc/length(r)
    end
    out
end

"""Noise reading per frame from an internal-sampling noise vector: the same window mean."""
function frame_noise!(out,noise,s::VH1Settings)
    for i in eachindex(out)
        r=window_range(i,s);acc=0.
        for j in r;acc+=noise[j];end
        out[i]=acc/length(r)
    end
    out
end

"""
    sim_vh1_electronics(noise_rng, adc_rng, detector, settings) -> VH1FrameResponse

Simulate every pad of a framed ion detector response: build the pad current,
form the noiseless frame readings, add VH1 noise drawn at the internal sampling
and converted like the signal, and quantise once per frame. Pads are independent.
Output size is limited to 256 MiB.
"""
function sim_vh1_electronics(noise_rng::AbstractRNG,adc_rng::AbstractRNG,d::IonDetectorResponse,s::VH1Settings)
    v=d.ions;ft=vh1_frame_time(v.time,s);grid=internal_grid(ft,s)
    nf=Int(ft.frame_count);nc=size(v.counts,2);m=samples_per_frame(s);nr=readings(ft,s)
    Int128(nr)*nc*sizeof(UInt16)<=256*1024^2||throw(ArgumentError("VH1 frame vectors exceed 256 MiB budget"))
    packet=ER.prepare_bucket_current(s.current_pulse,v.time.bucket_ms)
    plan=ER.prepare_event_current(v.time,packet,grid)
    pn=ER.prepare_csa_noise(nf*m,s.noise)
    codes=Matrix{UInt16}(undef,nr,nc);saturated=falses(nc);lo=0;hi=0
    clean=zeros(nr);noisy=zeros(nr);nz=zeros(nr)
    for j in 1:nc
        counts=view(v.counts,:,j)
        if all(iszero,counts)
            fill!(clean,s.v_ref_V)
        else
            frame_readings!(clean,ER.event_current!(plan,counts),ft,s)
        end
        saturated[j]=any(x->x<s.window_V[1]||x>s.window_V[2],clean)
        frame_noise!(nz,ER.sample_csa_noise!(noise_rng,pn),s)
        for i in 1:nr;noisy[i]=clean[i]+nz[i];end
        adc=ER.digitise_adc!(view(codes,:,j),adc_rng,noisy,s.adc)
        lo+=adc.below_range;hi+=adc.above_range
    end
    VH1FrameResponse(d.event_id,v.channels,ft,codes,saturated,lo,hi)
end

"""Converter codes to voltages at the centre of each code interval."""
function frame_voltages(r::VH1FrameResponse,s::VH1Settings)
    lsb=(Float64(s.adc.maximum_V)-s.adc.minimum_V)/2^s.adc.bits
    s.adc.minimum_V .+(Float64.(r.codes).+.5).*lsb
end

"""
    prepare_vh1_reconstruction(frames, buckets, settings) -> VH1ReconstructionPlan

Response matrix: drop of each reading after the reference frame below the
reference reading (the last reading of the reference frame) for one ion in each
bucket, formed with the same window operation as the simulation.
Covariance: VH1 noise at the internal sampling, averaged over each conversion
window, plus the converter quantisation, for reference minus each signal frame. Whitened
least squares as in the existing reconstruction.
"""
function prepare_vh1_reconstruction(ft::VH1FrameTime,t::IonReadoutTime,s::VH1Settings)
    grid=internal_grid(ft,s);nf=Int(ft.frame_count);r=Int(ft.reference_frames);kc=s.conversions_per_frame
    r==1 && nf>r||throw(ArgumentError("One reference frame and at least one signal frame are required"))
    nr=readings(ft,s);i0=r*kc   # the reference reading: last reading of the reference frame
    packet=ER.prepare_bucket_current(s.current_pulse,t.bucket_ms)
    A=zeros(nr-i0,Int(t.bucket_count));reading=zeros(nr)
    for k in 1:Int(t.bucket_count)
        counts=zeros(Int64,Int(t.bucket_count));counts[k]=1
        current=ER.prepared_ion_current(counts,t,packet,grid)
        frame_readings!(reading,current,ft,s)
        for i in i0+1:nr;A[i-i0,k]=reading[i0]-reading[i];end
    end
    groups=[window_range(i,s) for i in i0:nr]
    lsb=(Float64(s.adc.maximum_V)-s.adc.minimum_V)/2^s.adc.bits
    C=SR.averaged_covariance(s.noise.covariance_V2,groups,lsb^2/12)
    L=Matrix(cholesky(Symmetric(C)).L)
    W=L\A;H=qr(W)\(L\Matrix{Float64}(I,size(L)...))
    VH1ReconstructionPlan(ft,t,s,A,C,L,H,Matrix(Symmetric(H*C*H')))
end

"""
    reco_vh1_ions(response, plan) -> ReconstructedIons

Signed ion counts per bucket and pad from the frame codes. A pad is invalid
when its noiseless output left the linear window or a code sits at either end
of the converter range. Checks reuse the existing structure: baseline is the
reference reading, nonlinear is always false, residual is the RMS of the
frame drops after the fit.
"""
function reco_vh1_ions(r::VH1FrameResponse,p::VH1ReconstructionPlan)
    s=p.settings;nr=readings(p.frames,s);i0=Int(p.frames.reference_frames)*s.conversions_per_frame
    size(r.codes,1)==nr||throw(DimensionMismatch("Reading count"))
    nc=size(r.codes,2);nb=Int(p.buckets.bucket_count)
    V=frame_voltages(r,s);Y=zeros(nr-i0,nc);base=zeros(Float32,nc)
    for j in 1:nc
        ref=V[i0,j];base[j]=ref
        for i in i0+1:nr;Y[i-i0,j]=ref-V[i,j];end
    end
    Q=p.inverse*Y
    limit=BitVector([any(c->c==0||c==2^s.adc.bits-1,view(r.codes,:,j)) for j in 1:nc])
    converged=.!(limit.|r.saturated)
    rms=Float32[sqrt(mean(abs2,view(Y,:,j)-p.response*view(Q,:,j))) for j in 1:nc]
    for j in 1:nc;converged[j]||(Q[:,j].=NaN);end
    valid=repeat(reshape(converged,1,:),nb,1)
    time=p.buckets
    SR.ReconstructedIons(r.event_id,r.channels,time,Float32.(Q),BitMatrix(valid),
        SR.IonReconstructionChecks(base,falses(nc),converged,limit,rms))
end

export VH1FrameTime,VH1Settings,VH1FrameResponse,VH1ReconstructionPlan,
       load_vh1_settings,vh1_frame_time,sim_vh1_electronics,frame_voltages,
       prepare_vh1_reconstruction,reco_vh1_ions
end
