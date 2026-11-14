"""
    voltage_noise_from_spectrum(frequency_Hz, asd, sample_us, minimum_samples; oversampling=16)

Reuse the saved one-sided CSA spectrum and its log-linear PSD interpolation.
Integrate each spectral interval analytically, assign its power to the nearest
FFT frequency, and fold frequencies above Nyquist by sampling aliasing. Low
frequencies contribute to the DC bin rather than being discarded. The FFT
period is at least oversampling times the requested vector duration. Compare
with direct covariance integration and a doubled period to check resolution.
No additional analogue filter or high-frequency cutoff is introduced.
"""
function voltage_noise_from_spectrum(f,asd,sample_us,minimum_samples;oversampling=16)
    length(f)==length(asd)>=2 && all(diff(f).>0) && all(f.>0) && all(asd.>0)||throw(ArgumentError("Noise spectrum"))
    m=nextpow(2,max(16,oversampling*minimum_samples));dt=sample_us*1e-6;df=1/(m*dt)
    power=zeros(m÷2+1);j=1;x=Float64(first(f));stop=Float64(last(f))
    k=max(0,floor(Int,x/df+.5))
    while x<stop
        while j<length(f)-1 && x>=f[j+1];j+=1;end
        edge=min(stop,Float64(f[j+1]),(k+.5)*df)
        if edge<=x;k+=1;continue;end
        slope=2log(asd[j+1]/asd[j])/log(f[j+1]/f[j]);p=asd[j]^2
        exponent=slope+1
        area=abs(exponent)<1e-10 ? p*f[j]*log(edge/x) :
             p*f[j]*(x/f[j])^exponent*expm1(exponent*log(edge/x))/exponent
        alias=mod(k,m);alias=min(alias,m-alias);power[alias+1]+=area
        x=edge
        x>=(k+.5)*df*(1-1e-14)&&(k+=1)
    end
    eigenvalues=power.*(m/2);eigenvalues[1]*=2;eigenvalues[end]*=2
    cov=irfft(eigenvalues,m)
    CSAVoltageNoise(Float32(sample_us),cov)
end
"""Reusable FFT buffers and factors for independent pad-noise vectors."""
struct PreparedCSANoise{P}
    sample_count::Int64
    period::Int64
    factors::Vector{Float64}
    spectrum::Vector{ComplexF64}
    voltage::Vector{Float64}
    output::Vector{Float32}
    inverse::P
end

"""
    prepare_csa_noise(sample_count, noise)

Construct the smallest power-of-two circulant covariance that contains every
requested lag and has a nonnegative spectrum. Reflect the supplied positive
lags to close the shorter period. This preserves the requested covariance
exactly; it does not recompute or coarsen the physical noise spectrum. If a
short period is not positive semidefinite, double it up to the original period.
Prepare spectral factors and the inverse FFT once for reuse across pads.
"""
function prepare_csa_noise(n::Int64,s::CSAVoltageNoise)
    c=s.covariance_V2;full=length(c)
    0<n<=full÷2||throw(ArgumentError("Noise covariance period too short"))
    m=nextpow(2,max(2,2(n-1)))
    eig=Float64[]
    while true
        embedding=vcat(c[1:m÷2+1],reverse(c[2:m÷2]))
        eig=real.(rfft(embedding))
        tol=1e-12*max(maximum(abs,eig),eps())
        minimum(eig)>=-tol && break
        m>=full && throw(ArgumentError("Noise covariance is not positive semidefinite"))
        m=min(2m,full)
    end
    # Fourier coefficients of real Gaussian white noise have variance m/2
    # per real/imaginary component; DC and Nyquist are real with variance m.
    factors=sqrt.(max.(eig,0.).*(m/2));factors[1]*=sqrt(2.);factors[end]*=sqrt(2.)
    spectrum=zeros(ComplexF64,m÷2+1);voltage=zeros(m);out=zeros(Float32,n)
    inverse=plan_irfft(spectrum,m;flags=FFTW.ESTIMATE)
    PreparedCSANoise(n,m,factors,spectrum,voltage,out,inverse)
end

"""
    sample_csa_noise!(rng, prepared)

Generate one independent pad vector using prepared spectral factors and buffers.
Return a borrowed Float32 buffer, overwritten by the next call. Callers consume
it before advancing to the next pad. Only the inverse FFT is required.
"""
function sample_csa_noise!(rng::AbstractRNG,p::PreparedCSANoise)
    z=p.spectrum;f=p.factors
    z[1]=f[1]*randn(rng);z[end]=f[end]*randn(rng)
    for j in 2:length(z)-1
        z[j]=f[j]*complex(randn(rng),randn(rng))
    end
    mul!(p.voltage,p.inverse,z)
    for j in eachindex(p.output);p.output[j]=p.voltage[j];end
    p.output
end

"""
    sample_csa_noise(rng, sample_count, noise)

Return an owned Float32 voltage-noise vector with the supplied covariance.
Production prepares the same sampler once per event and reuses its buffers
across pads. Different pads remain independent; temporal covariance is retained.
"""
function sample_csa_noise(rng::AbstractRNG,n::Int64,s::CSAVoltageNoise)
    p=prepare_csa_noise(n,s)
    sample_csa_noise!(rng,p)
end
