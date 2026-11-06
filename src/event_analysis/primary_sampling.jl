# Shared primary-pair and event-placement sampling. Readout-independent.
"""
    sample_integer_moments(rng, mean::Float64, variance_ratio::Float64)::Int64

Sample a nonnegative integer with target mean `mean` and target variance
`variance_ratio*mean`, for `0 ≤ variance_ratio < 1`. The adjacent-integer
law has the minimum possible variance `f*(1-f)`, where `f` is the fractional
part of `mean`. Mix that law with `Binomial(ceil(mean/(1-ratio)), mean/n)`
to reach the target variance when feasible. Below the integer variance floor,
use the adjacent law. `mean=0` returns zero. Inputs are preserved; the caller
owns the RNG. Invalid or unrepresentable parameters raise an error.
"""
function sample_integer_moments(rng::AbstractRNG,mean::Float64,variance_ratio::Float64)::Int64
    isfinite(mean)&&0<=mean<typemax(Int64) &&
        isfinite(variance_ratio)&&0<=variance_ratio<1 ||
        throw(ArgumentError("Require representable mean ≥ 0 and variance ratio in [0,1)"))
    mean==0 && return Int64(0)
    lo=floor(Int64,mean);f=mean-Float64(lo)
    vfloor=f*(1-f);target=variance_ratio*mean
    adjacent()=lo+Int64(rand(rng)<f)
    target<=vfloor+eps(mean) && return adjacent()
    nfloat=ceil(mean/(1-variance_ratio))
    nfloat<typemax(Int64) || throw(OverflowError("Binomial trial count overflows Int64"))
    n=Int64(nfloat)
    vbin=mean*(1-mean/Float64(n))
    mix=clamp((target-vfloor)/(vbin-vfloor),0.0,1.0)
    return rand(rng)<mix ? Int64(rand(rng,Binomial(n,mean/Float64(n)))) : adjacent()
end

"""
    sample_fano_pairs(rng, hits::EventHits, W_eV::Float32, F::Float32)::Vector{Int64}

Sample one initial electron–ion pair count for each ordered source hit. The
per-hit mean is `1e6*energy_MeV/W_eV`. Use `sample_integer_moments` with
target variance `F*mean`, bounded below by the integer variance floor. The same resulting vector supplies both carrier branches.
`F` is a sub-Poisson variance ratio in [0,1). Preserve hits and consume only
the supplied RNG. Invalid inputs and sum overflow raise errors.
"""
function sample_fano_pairs(rng::AbstractRNG,hits::EventHits,
        W_eV::Float32,F::Float32)::Vector{Int64}
    hits.frame==source_frame || throw(ArgumentError("Pair sampling requires source-frame hits"))
    _validate_hit_vectors(hits.x_mm,hits.y_mm,hits.z_mm,hits.energy_MeV,
                          hits.time_ns,hits.particle_id,hits.hit_id)
    isfinite(W_eV) && W_eV>0 || throw(ArgumentError("W_eV must be finite and positive"))
    isfinite(F)&&0f0<=F<1f0 || throw(ArgumentError("F must be in [0,1)"))
    counts=Vector{Int64}(undef,length(hits.energy_MeV))
    for i in eachindex(counts)
        mean=Float64(hits.energy_MeV[i])*1e6/Float64(W_eV)
        counts[i]=sample_integer_moments(rng,mean,Float64(F))
    end
    foldl(Base.Checked.checked_add,counts;init=Int64(0))
    return counts
end

"""
    sample_drift_cm(rng::AbstractRNG, min_cm::Float32, max_cm::Float32)::Float32

Sample the common event drift uniformly in [min_cm,max_cm). Bounds are finite,
nonnegative and ordered; equal bounds return that fixed distance without
advancing rng. A nonzero interval consumes one Float32 uniform draw. Roundoff
at the upper bound is assigned to prevfloat(max_cm), preserving the half-open
interval. This pure boundary makes random state explicit to the caller.

# Example
```julia
using Random
rng = Xoshiro(42)
@assert 10f0 <= sample_drift_cm(rng,10f0,150f0) < 150f0
@assert sample_drift_cm(rng,20f0,20f0) === 20f0
```
"""
function sample_drift_cm(rng::AbstractRNG,min_cm::Float32,max_cm::Float32)::Float32
    all(isfinite,(min_cm,max_cm)) && 0f0 <= min_cm <= max_cm ||
        throw(ArgumentError("Drift bounds must be finite, nonnegative and ordered"))
    min_cm == max_cm && return min_cm
    return min(min_cm+(max_cm-min_cm)*rand(rng,Float32),prevfloat(max_cm))
end

