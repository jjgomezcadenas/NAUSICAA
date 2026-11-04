"""Bucket-integrating ASIC response in current units."""
module SiPMElectronics
using Random
export Electronics, apply_electronics

struct Electronics
    bucket_width_s::Float64
    noise_rms_pe::Float64
    current_limit_A::Float64
    function Electronics(dt, noise, limit)
        all(x->isfinite(x)&&x>0,(dt,limit)) && isfinite(noise) && noise>=0 ||
            throw(ArgumentError("Require positive duration/range and nonnegative noise"))
        new(dt,noise,limit)
    end
end

"""Convert owned sensor charges to average current, add Gaussian PE noise and clip above range.

The sensor charge input is preserved. Negative noise estimates are retained.
Returns the noise, unclipped/clipped current and saturation flags in owned arrays.
"""
function apply_electronics(rng::AbstractRNG,charge_C::AbstractVector{<:Real},q0::Real,e::Electronics)
    all(x->isfinite(x)&&x>=0,charge_C) && isfinite(q0) && q0>0 ||
        throw(ArgumentError("Require nonnegative finite charge and positive PE calibration"))
    noise_pe = e.noise_rms_pe .* randn(rng,length(charge_C))
    current_unclipped_A = (charge_C .+ q0 .* noise_pe) ./ e.bucket_width_s
    saturated = current_unclipped_A .> e.current_limit_A
    current_A = min.(current_unclipped_A,e.current_limit_A)
    (;noise_pe,current_unclipped_A,current_A,saturated)
end
end
