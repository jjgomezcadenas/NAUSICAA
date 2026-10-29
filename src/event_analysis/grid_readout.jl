"""
    el_photon_yield(gain::Float32, threshold::Float32, reduced_field::Float32,
                    pressure_bar::Float32, depth_mm::Float32)::Float32

Evaluate the fitted mean EL photons per electron. Gain is photons/electron/kV;
threshold and reduced field are kV/cm/bar, pressure is bar and depth is mm.
The result is zero below threshold. The user selects the pure-Xe or Xe-He
Table 1 fit in TOML. Require finite nonnegative gain and threshold, positive
pressure and depth, and nonnegative reduced field.
"""
function el_photon_yield(gain::Float32, threshold::Float32, reduced_field::Float32,
                         pressure_bar::Float32, depth_mm::Float32)::Float32
    all(isfinite,(gain,threshold,reduced_field,pressure_bar,depth_mm)) &&
        gain >= 0f0 && threshold >= 0f0 && reduced_field >= 0f0 &&
        pressure_bar > 0f0 && depth_mm > 0f0 || throw(ArgumentError("Invalid EL yield inputs"))
    yield = gain * max(reduced_field-threshold,0f0) * pressure_bar * (depth_mm/10f0)
    isfinite(yield) || throw(ArgumentError("EL yield overflows Float32"))
    return yield
end

"""
    axial_light_acceptance(radius_mm::Float32, depth_mm::Float32)::Float32

Return the direct fraction of isotropic VUV photons reaching an exit-plane
circular SiPM through a cylindrical hole. Assume emission uniformly along the
hole axis, a SiPM covering the hole, absorbing walls and no mesh loss. The
depth average is (depth+radius-hypot(depth,radius))/(2*depth). Both lengths
are positive finite mm. The returned probability lies in (0,0.5).
"""
function axial_light_acceptance(radius_mm::Float32, depth_mm::Float32)::Float32
    all(isfinite,(radius_mm,depth_mm)) && min(radius_mm,depth_mm)>0f0 ||
        throw(ArgumentError("Hole radius and depth must be positive finite mm"))
    # Rationalized form avoids cancellation for apertures much larger than h.
    value=radius_mm/(depth_mm+radius_mm+hypot(depth_mm,radius_mm))
    0f0 < value < 0.5f0 || throw(ArgumentError("Invalid optical acceptance"))
    return value
end

"""
    sample_detected_pe(rng, electrons, photons_per_electron, acceptance, pde)::Int64

Draw an integer PE count from a binomial distribution. The trial count is
round(Int64, electrons × fitted mean photons/electron), ties to even; the
success probability is the product of direct geometrical acceptance and SiPM
PDE. This treats the fitted photon yield as a fixed mean per input cell and
samples optical detection fluctuations. Inputs are not mutated. Require
nonnegative electrons/yield and probabilities in [0,1]. Detect overflow.
"""
function sample_detected_pe(rng::AbstractRNG,electrons::Int64,
                            photons_per_electron::Float32,acceptance::Float32,
                            pde::Float32)::Int64
    electrons >= 0 && isfinite(photons_per_electron) && photons_per_electron >= 0f0 &&
        isfinite(acceptance) && 0f0 <= acceptance <= 1f0 &&
        isfinite(pde) && 0f0 <= pde <= 1f0 || throw(ArgumentError("Invalid PE inputs"))
    trials_f = Float64(electrons)*Float64(photons_per_electron)
    isfinite(trials_f) && trials_f < typemax(Int64) || throw(OverflowError("Photon trial count overflows Int64"))
    trials = round(Int64,trials_f)
    return Int64(rand(rng,Binomial(trials,Float64(acceptance)*Float64(pde))))
end

"""
    collect_grid_history(diffused::DataFrame; reference_mm, pitch_mm,
                         z_bin_width_mm, grid_offset_mm, drift_speed_mm_ms,
                         drift_direction, minimum_xy_span_mm=160f0)

Map the complete retained cloud to a square lattice of holes and time slices.
The lattice reference comes from the common undiffused event; one random phase
per physical grid is supplied by the caller. Every occupied diffusion cell is
assigned to its containing x/y grid cell and relative arrival slice using
Petit's adapted pad-axis and occupancy functions. Input counts are Int64 and
positions Float32 mm. Return the owned axes and sparse count history.
"""
function collect_grid_history(diffused::DataFrame;reference_mm::NTuple{3,Float32},
        pitch_mm::Float32,z_bin_width_mm::Float32,
        grid_offset_mm::NTuple{2,Float32},drift_speed_mm_ms::Float32,
        drift_direction::Int,minimum_xy_span_mm::Float32=160f0)
    axes=pad_axes(diffused;reference_mm,pitch_mm,z_bin_width_mm,
                  grid_offset_mm,minimum_xy_span_mm)
    history=compute_occupancy_analysis(diffused;x_edges_mm=axes.x_edges_mm,
        y_edges_mm=axes.y_edges_mm,z_edges_mm=axes.z_edges_mm,
        drift_speed_mm_ms,drift_direction)
    return (;axes,history)
end

"""
    detect_pe_history(rng, electron_history; photons_per_electron,
                      acceptance, pde)

Sample one PE count for each occupied hole/time cell using `sample_detected_pe`.
Return a sparse DataFrame (pad_ix, pad_iy, slice_index, pe), all Int64. Zero-PE
cells are omitted. The per-cell EL yield is evaluated from the total electrons
entering that cell; no photon arrays are retained. Preserve the input history.
"""
function detect_pe_history(rng::AbstractRNG,electron_history;
        photons_per_electron::Float32,acceptance::Float32,pde::Float32)
    counts=DataFrame(pad_ix=Int64[],pad_iy=Int64[],slice_index=Int64[],pe=Int64[])
    for row in eachrow(electron_history.counts)
        pe=sample_detected_pe(rng,row.ions,photons_per_electron,acceptance,pde)
        pe>0 && push!(counts,(row.pad_ix,row.pad_iy,row.slice_index,pe))
    end
    return counts
end

"""Return the signed cell index on a Float32 edge lattice.

The Float64 quotient supplies an initial integer estimate. Neighbour checks
use the same Float32 edge arithmetic as the saved axes; a carrier on an
internal edge belongs to the cell on its right. This works before the finite
axis extent is known and permits a single sampling/collection pass.
"""
@inline function _direct_lattice_index(x::Float32,base::Float32,width::Float32)::Int64
    estimate=(Float64(x)-Float64(base))/Float64(width)
    isfinite(estimate) && typemin(Int64)+2 < estimate < typemax(Int64)-2 ||
        throw(ArgumentError("Grid index is outside Int64 range"))
    index=floor(Int64,estimate)
    while x < base+Float32(index)*width
        index-=1
    end
    while x >= base+Float32(index+1)*width
        index+=1
    end
    return index
end

"""
    diffuse_to_grid(rng, hits::EventHits, carriers_per_hit::Vector{Int64};
        sigma_t_mm::Float32, sigma_l_mm::Float32, nsigma::Float32,
        reference_mm::NTuple{3,Float32}, grid_offset_mm::NTuple{2,Float32},
        pitch_mm::Float32, bucket_width_mm::Float32,
        drift_speed_mm_ms::Float32, drift_direction::Int,
        minimum_xy_span_mm::Float32=160f0)

Sample and collect one source-frame event without storing individual carrier
positions. Follow Petit's hit-major, carrier-minor Gaussian draw order: x, y,
then z, each with `randn(rng,Float32)`. Reject samples outside the same finite
hit-extent ± `nsigma` × width window used by `diffuse_xyz_image_mc`. Each
retained sample increments its containing transverse lattice cell and 2D
position-derived longitudinal bucket by one. The lattice reference and phase
are caller supplied; the transverse reference plus phase is an edge. The z
reference is an edge. The output extends each axis to cover all retained
samples and at least `minimum_xy_span_mm` transversely.

Return owned sparse `counts::DataFrame` with `pad_ix`, `pad_iy`,
`slice_index`, `ions` as Int64, explicit Float32 x/y/z edges and relative
`time_edges_ms`, plus `sampled_ions`, `retained_ions`, `outside_ions`,
`total_ions`, reference and phase. `slice_index` increases with relative
arrival time; `drift_direction=+1` makes increasing z later, while -1 makes
decreasing z later. Missing cells and slices have zero count. The ledger
satisfies sampled = retained + outside = sum(counts.ions) + outside.

Validate source-frame hits, nonnegative counts, finite widths/coordinates,
positive pitch/bucket/speed/nsigma, and direction ±1 before drawing. Inputs
retain their values. The result has no dependence on a display histogram's
bin count. The event's physics grid uses Float32 edge arithmetic; Float64 is
used only for an initial integer-index estimate.
"""
function diffuse_to_grid(rng::AbstractRNG,hits::EventHits,carriers_per_hit::Vector{Int64};
        sigma_t_mm::Float32,sigma_l_mm::Float32,nsigma::Float32,
        reference_mm::NTuple{3,Float32},grid_offset_mm::NTuple{2,Float32},
        pitch_mm::Float32,bucket_width_mm::Float32,
        drift_speed_mm_ms::Float32,drift_direction::Int,
        minimum_xy_span_mm::Float32=160f0)
    hits.frame==source_frame || throw(ArgumentError("Direct collection requires source-frame hits"))
    _validate_hit_vectors(hits.x_mm,hits.y_mm,hits.z_mm,hits.energy_MeV,
                          hits.time_ns,hits.particle_id,hits.hit_id)
    length(carriers_per_hit)==length(hits.x_mm) ||
        throw(DimensionMismatch("One carrier count is required per hit"))
    all(>=(0),carriers_per_hit) || throw(ArgumentError("Carrier counts must be nonnegative"))
    all(v->isfinite(v)&&v>=0f0,(sigma_t_mm,sigma_l_mm)) ||
        throw(ArgumentError("Diffusion widths must be finite and nonnegative"))
    all(v->isfinite(v)&&v>0f0,(nsigma,pitch_mm,bucket_width_mm,drift_speed_mm_ms)) ||
        throw(ArgumentError("Window, grid widths and speed must be finite and positive"))
    all(isfinite,reference_mm) && all(isfinite,grid_offset_mm) &&
        isfinite(minimum_xy_span_mm) && minimum_xy_span_mm>=0f0 ||
        throw(ArgumentError("Reference, phase and span must be finite"))
    drift_direction in (-1,1) || throw(ArgumentError("drift_direction must be ±1"))
    sampled=foldl(Base.Checked.checked_add,carriers_per_hit;init=Int64(0))
    window=(_diffusion_axis(hits.x_mm,sigma_t_mm,1,nsigma),
            _diffusion_axis(hits.y_mm,sigma_t_mm,1,nsigma),
            _diffusion_axis(hits.z_mm,sigma_l_mm,1,nsigma))
    bases=(reference_mm[1]+grid_offset_mm[1],
           reference_mm[2]+grid_offset_mm[2],reference_mm[3])
    all(isfinite,bases) || throw(ArgumentError("Grid base overflows Float32"))
    widths=(pitch_mm,pitch_mm,bucket_width_mm)
    sums=Dict{NTuple{3,Int64},Int64}()
    retained=Int64(0)
    minidx=[typemax(Int64),typemax(Int64),typemax(Int64)]
    maxidx=[typemin(Int64),typemin(Int64),typemin(Int64)]
    for hit in eachindex(carriers_per_hit)
        hx,hy,hz=hits.x_mm[hit],hits.y_mm[hit],hits.z_mm[hit]
        for _ in 1:carriers_per_hit[hit]
            x=hx+randn(rng,Float32)*sigma_t_mm
            y=hy+randn(rng,Float32)*sigma_t_mm
            z=hz+randn(rng,Float32)*sigma_l_mm
            if !(window[1][1]<=x<=window[1][2] &&
                 window[2][1]<=y<=window[2][2] &&
                 window[3][1]<=z<=window[3][2])
                continue
            end
            key=(_direct_lattice_index(x,bases[1],widths[1]),
                 _direct_lattice_index(y,bases[2],widths[2]),
                 _direct_lattice_index(z,bases[3],widths[3]))
            sums[key]=Base.Checked.checked_add(get(sums,key,Int64(0)),Int64(1))
            retained=Base.Checked.checked_add(retained,Int64(1))
            for axis in 1:3
                minidx[axis]=min(minidx[axis],key[axis])
                maxidx[axis]=max(maxidx[axis],key[axis])
            end
        end
    end
    if retained==0
        minidx .= 0;maxidx .= 0
    end
    half=ceil(Int64,(minimum_xy_span_mm/2f0)/pitch_mm)
    for axis in 1:2
        minidx[axis]=min(minidx[axis],-half)
        maxidx[axis]=max(maxidx[axis],half-1)
    end
    axes=ntuple(axis->begin
        edges=Float32[bases[axis]+Float32(i)*widths[axis] for i in minidx[axis]:maxidx[axis]+1]
        _check_edges(edges)
        edges
    end,3)
    x_edges_mm,y_edges_mm,z_edges_mm=axes
    nz=length(z_edges_mm)-1
    time_origin_z_mm=drift_direction==1 ? first(z_edges_mm) : last(z_edges_mm)
    time_edges_ms=drift_direction==1 ?
        (z_edges_mm.-time_origin_z_mm)./drift_speed_mm_ms :
        (time_origin_z_mm.-reverse(z_edges_mm))./drift_speed_mm_ms
    _check_edges(time_edges_ms)
    counts=DataFrame(pad_ix=Int64[],pad_iy=Int64[],slice_index=Int64[],ions=Int64[])
    for (key,n) in sums
        ix=key[1]-minidx[1]+1;iy=key[2]-minidx[2]+1
        iz=drift_direction==1 ? key[3]-minidx[3]+1 : maxidx[3]-key[3]+1
        push!(counts,(ix,iy,iz,n))
    end
    sort!(counts,[:pad_ix,:pad_iy,:slice_index])
    outside=sampled-retained
    return (;counts,x_edges_mm,y_edges_mm,z_edges_mm,time_edges_ms,
            sampled_ions=sampled,retained_ions=retained,outside_ions=outside,
            total_ions=retained,reference_mm,grid_offset_mm,drift_direction,
            time_origin_z_mm)
end

"""
    diffused_carrier_reference(rng, hits, carriers_per_hit;
        sigma_t_mm::Float32, sigma_l_mm::Float32, nsigma::Float32)

Replay the same per-carrier Gaussian draws and finite diffusion window as
`diffuse_to_grid`, reducing the retained samples to their ion-weighted x/y
centre and minimum z. This supplies the earlier pixel viewer's grid reference
before a second direct collection replay, without retaining coordinates or
depending on a display histogram. Return `reference_mm::NTuple{3,Float32}`
and the exact sampled/retained/outside ledger. Float64 sums reduce the many
Float32 coordinates; the returned reference remains Float32 mm. Require at
least one retained carrier. Inputs preserve their values.
"""
function diffused_carrier_reference(rng::AbstractRNG,hits::EventHits,
        carriers_per_hit::Vector{Int64};sigma_t_mm::Float32,
        sigma_l_mm::Float32,nsigma::Float32)
    hits.frame==source_frame || throw(ArgumentError("Reference requires source-frame hits"))
    _validate_hit_vectors(hits.x_mm,hits.y_mm,hits.z_mm,hits.energy_MeV,
                          hits.time_ns,hits.particle_id,hits.hit_id)
    length(carriers_per_hit)==length(hits.x_mm) ||
        throw(DimensionMismatch("One carrier count is required per hit"))
    all(>=(0),carriers_per_hit) || throw(ArgumentError("Carrier counts must be nonnegative"))
    all(v->isfinite(v)&&v>=0f0,(sigma_t_mm,sigma_l_mm)) &&
        isfinite(nsigma) && nsigma>0f0 ||
        throw(ArgumentError("Diffusion widths and nsigma are invalid"))
    sampled=foldl(Base.Checked.checked_add,carriers_per_hit;init=Int64(0))
    window=(_diffusion_axis(hits.x_mm,sigma_t_mm,1,nsigma),
            _diffusion_axis(hits.y_mm,sigma_t_mm,1,nsigma),
            _diffusion_axis(hits.z_mm,sigma_l_mm,1,nsigma))
    sum_x=0.0;sum_y=0.0;minimum_z=Inf32;retained=Int64(0)
    for hit in eachindex(carriers_per_hit)
        hx,hy,hz=hits.x_mm[hit],hits.y_mm[hit],hits.z_mm[hit]
        for _ in 1:carriers_per_hit[hit]
            x=hx+randn(rng,Float32)*sigma_t_mm
            y=hy+randn(rng,Float32)*sigma_t_mm
            z=hz+randn(rng,Float32)*sigma_l_mm
            if window[1][1]<=x<=window[1][2] &&
               window[2][1]<=y<=window[2][2] &&
               window[3][1]<=z<=window[3][2]
                sum_x+=Float64(x);sum_y+=Float64(y)
                minimum_z=min(minimum_z,z)
                retained=Base.Checked.checked_add(retained,Int64(1))
            end
        end
    end
    retained>0 || throw(ArgumentError("Grid reference requires retained carriers"))
    reference_mm=(Float32(sum_x/retained),Float32(sum_y/retained),minimum_z)
    return (;reference_mm,sampled_ions=sampled,retained_ions=retained,
            outside_ions=sampled-retained)
end
