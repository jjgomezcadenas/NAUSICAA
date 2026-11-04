"""Stochastic avalanche charge from incident photon buckets and primary dark pulses."""
module SiPMResponse
using Random, Distributions
export Sensor, elementary_charge_C, charge_per_pe_C, dark_rate_Hz, recovered_fraction,
       simulate_channel, axel_mean, axel_inverse

const elementary_charge_C = 1.602176634e-19

"""Effective SiPM: detection probability, microcells, gain, recovery time and total dark rate.

Times are seconds, dark rate Hz, gain electrons per fully recovered avalanche.
The independent primary dark process is stationary at the acquisition boundary.
"""
struct Sensor
    pde::Float64
    microcells::Int
    gain::Float64
    recovery_s::Float64
    dark_rate_Hz::Float64
    function Sensor(pde, microcells::Integer, gain, recovery_s, dark_rate_Hz)
        isfinite(pde) && 0 <= pde <= 1 || throw(ArgumentError("PDE must lie in [0,1]"))
        microcells > 0 || throw(ArgumentError("Microcell count must be positive"))
        isfinite(gain) && gain > 0 || throw(ArgumentError("Gain must be positive"))
        all(x -> isfinite(x) && x >= 0, (recovery_s, dark_rate_Hz)) ||
            throw(ArgumentError("Recovery and dark rate must be finite and nonnegative"))
        new(pde, microcells, gain, recovery_s, dark_rate_Hz)
    end
end
charge_per_pe_C(s::Sensor) = elementary_charge_C * s.gain
dark_rate_Hz(rate, reference_area, diameter) = rate * (pi * diameter^2 / 4) / reference_area

"""Recovered charge fraction after a nonnegative interval; zero recovery time is ideal."""
function recovered_fraction(elapsed::Real, tau::Real)
    elapsed >= 0 && isfinite(tau) && tau >= 0 || throw(ArgumentError("Invalid recovery interval"))
    tau == 0 ? 1.0 : -expm1(-elapsed / tau)
end

"""Simulate one channel in successive buckets, preserving input counts.

PDE uses a binomial draw per bucket. Detected optical triggers and Poisson dark
pulses receive uniform times and microcell indices. Every avalanche resets its
cell and contributes its recovered charge to the trigger bucket. Float64 local
times avoid loss of precision for large signed source-bucket indices. The return
value owns count and charge arrays; the caller owns the RNG. No pulse shaping
or recovery-dependent PDE is applied in this first model.
"""
function simulate_channel(rng::AbstractRNG, incident::AbstractVector{<:Integer},
                          dt::Real, sensor::Sensor)
    isfinite(dt) && dt > 0 || throw(ArgumentError("Bucket duration must be positive"))
    all(>=(0), incident) || throw(ArgumentError("Incident counts must be nonnegative"))
    n = length(incident); optical = zeros(Int64,n); dark = zeros(Int64,n)
    optical_charge = zeros(n); dark_charge = zeros(n)
    # Independent stationary Poisson ages per microcell before the first bucket.
    last = sensor.dark_rate_Hz == 0 ? fill(-Inf, sensor.microcells) :
           -randexp(rng,sensor.microcells) .* (sensor.microcells / sensor.dark_rate_Hz)
    arrivals = Tuple{Float64,Bool}[]
    for b in eachindex(incident)
        optical[b] = rand(rng, Binomial(incident[b],sensor.pde))
        dark[b] = rand(rng, Poisson(sensor.dark_rate_Hz * dt))
        empty!(arrivals)
        for _ in 1:optical[b]; push!(arrivals,(rand(rng)*dt,false)); end
        for _ in 1:dark[b]; push!(arrivals,(rand(rng)*dt,true)); end
        sort!(arrivals;by=first)
        for (offset,is_dark) in arrivals
            time = (b-1)*dt + offset
            cell = rand(rng,1:sensor.microcells)
            q = charge_per_pe_C(sensor) * recovered_fraction(time-last[cell],sensor.recovery_s)
            last[cell] = time
            if is_dark; dark_charge[b] += q; else; optical_charge[b] += q; end
        end
    end
    (; incident_photons=Int64.(incident), detected_pe=optical, dark_pulses=dark,
       optical_charge_C=optical_charge, dark_charge_C=dark_charge,
       charge_C=optical_charge+dark_charge)
end

"""AXEL stationary mean charge in fully recovered avalanche equivalents."""
axel_mean(n::Real, s::Sensor, dt::Real) = n / (1 + n*s.recovery_s/(s.microcells*dt))

"""Stationary mean inverse with an explicit validity flag; clipped data are invalid.

Signed finite noise estimates are retained. The inverse is a diagnostic mean
correction and does not recover the realized incident photon count.
"""
function axel_inverse(m::Real,s::Sensor,dt::Real;clipped::Bool=false)
    isfinite(dt) && dt > 0 || throw(ArgumentError("Bucket duration must be positive"))
    denominator = 1-m*s.recovery_s/(s.microcells*dt)
    valid = !clipped && isfinite(m) && denominator > 0
    (; pe=valid ? m/denominator : NaN, valid)
end
include("sipm_statistical.jl")
export StatisticalOptions, simulate_channel_statistical, occupancy_moments, recover_totals

end
