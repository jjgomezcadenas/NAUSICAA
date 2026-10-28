"""
    ion_counts(hits::EventHits, ionisation_energy_eV::Float32)

Convert deposited MeV into expected positive-ion counts using the common mean
pair-creation energy W in eV: expected[j] = energy_MeV[j] × 10⁶ / W. Apply
Petit's per-hit `round(Int64, expected[j])` rule (nearest integer, ties to even)
to obtain the counts used for Monte Carlo diffusion. Each ion has charge +e.

Return an owned NamedTuple with expected::Vector{Float32}, sampled::Vector{Int64},
event_energy_MeV::Float32, total_expected::Float32, mean_per_hit::Float32 and
rounding_difference::Float32 = sum(sampled) − total_expected. Reductions use
Float32; the sum of integer counts is checked for Int64 overflow. Rounding is
deterministic. Statistical ionisation fluctuations are a separate model.

W must be finite and positive. Validate current hit values before calculating;
unrepresentable expectations/counts or totals raise ArgumentError/OverflowError.
Preserve the input arrays and row order. Source: Petit transform_hits_df in
src/itaca_functions.jl; explicit W and charge accounting are the adaptations.

# Example
```julia
h = EventHits(1,source_frame,zeros(Float32,2),zeros(Float32,2),zeros(Float32,2),
              Float32[0.001,0],zeros(Float32,2),Int32[7,7],Int32[0,1])
p = ion_counts(h,20f0)
@assert p.sampled == Int64[50,0] && p.mean_per_hit ≈ 25f0
```
"""
function ion_counts(hits::EventHits, ionisation_energy_eV::Float32)
    _validate_hit_vectors(hits.x_mm,hits.y_mm,hits.z_mm,hits.energy_MeV,
                         hits.time_ns,hits.particle_id,hits.hit_id)
    _check_w(ionisation_energy_eV)
    expected = hits.energy_MeV .* 1f6 ./ ionisation_energy_eV
    all(x -> isfinite(x) && x < Float32(typemax(Int64)), expected) ||
        throw(ArgumentError("Ion expectations must be finite and representable as Int64"))
    sampled = round.(Int64,expected)
    integer_total = foldl(Base.Checked.checked_add,sampled;init=Int64(0))
    energy = sum(hits.energy_MeV)
    total = sum(expected)
    isfinite(energy) && isfinite(total) || throw(ArgumentError("Event totals overflow Float32"))
    return (;expected,sampled,event_energy_MeV=energy,total_expected=total,
             mean_per_hit=total/Float32(length(expected)),
             rounding_difference=Float32(integer_total)-total)
end

"""Validate the common positive finite mean pair-creation energy in eV."""
function _check_w(w::Float32)
    isfinite(w) && w > 0f0 || throw(ArgumentError("ionisation_energy_eV must be finite and positive"))
    return nothing
end

"""
    sigma_t_ion_mm(T_K::Float32, L_cm::Float32, E_V_cm::Float32)::Float32

Petit's thermal positive-ion width: sqrt(2(kB/e) T L/E) × 10, in mm, with
kB/e = 8.6173f-5 V/K as in the source. Inputs are temperature in K, common event
drift distance in cm and field magnitude in V/cm. The factor 10 converts cm
to mm. The model assumes a uniform field and thermal, singly charged ions.
The baseline uses this same width longitudinally and transversely.

Require finite positive T and E and nonnegative L. Zero L gives zero width.
Return Float32, preserving inputs; nonfinite results raise ArgumentError.
Source: Petit src/itaca_functions.jl:sigma_t_ion_mm. Detector SI fields are
converted explicitly by the caller: E_V_cm = E_V_m / 100.

# Example
```julia
s = sigma_t_ion_mm(293.15f0,100f0,200f0)
@assert 1.58f0 < s < 1.60f0
@assert sigma_t_ion_mm(293.15f0,0f0,200f0) === 0f0
```
"""
function sigma_t_ion_mm(T_K::Float32,L_cm::Float32,E_V_cm::Float32)::Float32
    all(isfinite,(T_K,L_cm,E_V_cm)) && T_K > 0f0 && E_V_cm > 0f0 && L_cm >= 0f0 ||
        throw(ArgumentError("Require finite T > 0 K, L ≥ 0 cm and E > 0 V/cm"))
    sigma = sqrt(2f0*8.6173f-5*T_K*(L_cm/E_V_cm))*10f0
    isfinite(sigma) || throw(ArgumentError("Thermal width overflows Float32"))
    return sigma
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

"""
    sample_grid_offset_mm(rng::AbstractRNG, pitch_mm::Float32)::NTuple{2,Float32}

Return independent uniform x/y offsets in [0,pitch_mm), in mm. Consume two
Float32 uniform draws from the supplied RNG. The same pair is retained for all
z slices of an event by the later pad stage. Pitch must be finite and positive;
roundoff at the upper bound is assigned to prevfloat(pitch_mm).

# Example
```julia
using Random
dx,dy = sample_grid_offset_mm(Xoshiro(42),2f0)
@assert 0f0 <= dx < 2f0 && 0f0 <= dy < 2f0
```
"""
function sample_grid_offset_mm(rng::AbstractRNG,pitch_mm::Float32)::NTuple{2,Float32}
    isfinite(pitch_mm) && pitch_mm > 0f0 || throw(ArgumentError("pitch_mm must be finite and positive"))
    return ntuple(_ -> min(pitch_mm*rand(rng,Float32),prevfloat(pitch_mm)),2)
end
