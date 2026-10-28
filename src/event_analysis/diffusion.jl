# Adapted from Petit src/itaca_functions.jl (MIT; petit/LICENSE).
# Source hashes and the list of targeted changes are in PROVENANCE.toml.

"""
Construct owned Float32 histogram edges from the hit extent plus nsigma × sigma.
Positive-span axes have nbins cells; require strictly increasing representable
edges. This checks the precision of the saved axis before any random sampling.
When all coordinates coincide and sigma=0, use one cell whose edges are the
adjacent representable Float32 values. Its centre is the original coordinate.
"""
function _diffusion_axis(values::Vector{Float32},sigma::Float32,nbins::Int,nsigma::Float32)
    lo,hi = extrema(values)
    if lo == hi && sigma == 0f0
        edges = Float32[prevfloat(lo),nextfloat(hi)]
        all(isfinite,edges) || throw(ArgumentError("A finite neighbouring Float32 edge is required"))
        return edges
    end
    lower,upper = lo-nsigma*sigma,hi+nsigma*sigma
    all(isfinite,(lower,upper)) || throw(ArgumentError("Diffusion bounds overflow Float32"))
    lower < upper || throw(ArgumentError("Diffusion bounds collapse at Float32 precision"))
    edges = collect(range(lower,upper;length=nbins+1))
    all(diff(edges) .> 0f0) || throw(ArgumentError("Histogram bins are too fine for Float32 coordinates"))
    return edges
end

"""
Return the half-open histogram cell containing x, including the final upper edge
in the final cell; return 0 outside the saved Float32 edge vector. Callers supply
strictly increasing finite edges. Binary search follows the stored boundaries
even when Float32 rounding makes adjacent widths slightly unequal.
"""
@inline function _histogram_index(x::Float32,edges::Vector{Float32})::Int
    first(edges) <= x <= last(edges) || return 0
    return min(searchsortedlast(edges,x),length(edges)-1)
end

"""
    diffuse_xyz_image_mc(rng::AbstractRNG, hits::EventHits, ions::Vector{Int64};
        sigma_t_mm::Float32, sigma_l_mm::Float32, nbins::Int,
        nsigma::Float32, ionisation_energy_eV::Float32)

Diffuse each supplied ion using Petit's per-carrier Monte Carlo procedure.
For each hit, draw x, y and z independent standard normals in that order;
multiply x/y by sigma_t_mm and z by sigma_l_mm and add the original coordinates.
The baseline supplies the common thermal event width on all three axes. Setting
sigma_l_mm=0f0 gives the approved longitudinal comparison. Each ion consumes
three `randn(rng,Float32)` draws, including axes with zero width.

Accumulate into a dense Int64 histogram spanning each hit-coordinate extent plus
nsigma times its width on each side. Use nbins cells on each positive-span axis.
An axis with identical coordinates and zero width uses one bin centred on that
coordinate, with adjacent representable Float32 values as edges.
Save explicit Float32 edges; internal cells are half open and the final upper
edge belongs to the last cell. Samples outside any axis contribute to
outside_ions. Inputs use the original hit frame, mm and MeV. No per-ion array is
retained; memory scales with the histogram and occupied cells.

Return an owned NamedTuple:
- bins::DataFrame: event_id::Int64, x/y/z/energy::Float32, ions::Int64. Positions
  are geometric bin centres in mm; energy = ions × W / 10⁶ in MeV. Only occupied
  cells are stored, in x-outer/y-middle/z-inner order, following Petit.
- sampled_ions, retained_ions, outside_ions::Int64, with sampled = retained +
  outside and sum(bins.ions) = retained. The full retained table is available
  for subsequent pad counting, including charge below reconstruction thresholds.
- x_edges_mm, y_edges_mm, z_edges_mm::Vector{Float32}: independently owned axes.

Validate hit values, count length/nonnegativity, finite nonnegative widths,
positive nbins/nsigma/W, representable edges and Int64 total counts before RNG
consumption. Invalid arguments raise ArgumentError or DimensionMismatch; count
overflow raises OverflowError. An all-zero count vector returns a typed empty
table and consumes no draws. Input arrays retain their values. Identical inputs,
RNG state and Julia version reproduce the same histogram.

Source: Petit src/itaca_functions.jl:diffuse_xyz_image_mc. Adaptations are Float32,
explicit counts/W/RNG, saved axes and charge ledger, and edge-based assignment
with final-edge inclusion. The same W must be passed to ion_counts and this
function. Ion production and histogram truncation remain separately accounted.

# Example
```julia
using Random
h = EventHits(1,source_frame,Float32[0,1],Float32[0,1],Float32[0,1],
              Float32[0.001,0.001],zeros(Float32,2),Int32[7,7],Int32[0,1])
n = ion_counts(h,20f0)
d = diffuse_xyz_image_mc(Xoshiro(42),h,n.sampled;sigma_t_mm=0.5f0,
    sigma_l_mm=0.5f0,nbins=20,nsigma=3f0,ionisation_energy_eV=20f0)
@assert d.sampled_ions == 100
@assert sum(d.bins.ions) + d.outside_ions == d.sampled_ions
```
"""
function diffuse_xyz_image_mc(rng::AbstractRNG,hits::EventHits,ions::Vector{Int64};
        sigma_t_mm::Float32,sigma_l_mm::Float32,nbins::Int,nsigma::Float32,
        ionisation_energy_eV::Float32)
    _validate_hit_vectors(hits.x_mm,hits.y_mm,hits.z_mm,hits.energy_MeV,
                         hits.time_ns,hits.particle_id,hits.hit_id)
    length(ions) == length(hits.x_mm) || throw(DimensionMismatch("One ion count is required per hit"))
    all(>=(0),ions) || throw(ArgumentError("Ion counts must be nonnegative"))
    all(x -> isfinite(x) && x >= 0f0,(sigma_t_mm,sigma_l_mm)) ||
        throw(ArgumentError("Diffusion widths must be finite and nonnegative"))
    nbins > 0 && isfinite(nsigma) && nsigma > 0f0 ||
        throw(ArgumentError("nbins and nsigma must be positive; nsigma must be finite"))
    _check_w(ionisation_energy_eV)
    sampled_ions = foldl(Base.Checked.checked_add,ions;init=Int64(0))
    x_edges = _diffusion_axis(hits.x_mm,sigma_t_mm,nbins,nsigma)
    y_edges = _diffusion_axis(hits.y_mm,sigma_t_mm,nbins,nsigma)
    z_edges = _diffusion_axis(hits.z_mm,sigma_l_mm,nbins,nsigma)
    nx,ny,nz = length(x_edges)-1,length(y_edges)-1,length(z_edges)-1
    hist_counts = zeros(Int64,nx,ny,nz)

    # Petit's hit-major, carrier-minor loop with three scalar Gaussian draws.
    for i in eachindex(ions)
        x_center,y_center,z_center = hits.x_mm[i],hits.y_mm[i],hits.z_mm[i]
        for _ in 1:ions[i]
            x_diffused = x_center + randn(rng,Float32)*sigma_t_mm
            y_diffused = y_center + randn(rng,Float32)*sigma_t_mm
            z_diffused = z_center + randn(rng,Float32)*sigma_l_mm
            ix = _histogram_index(x_diffused,x_edges)
            iy = _histogram_index(y_diffused,y_edges)
            iz = _histogram_index(z_diffused,z_edges)
            if ix > 0 && iy > 0 && iz > 0
                hist_counts[ix,iy,iz] += 1
            end
        end
    end

    x_centers = x_edges[1:end-1]./2f0 .+ x_edges[2:end]./2f0
    y_centers = y_edges[1:end-1]./2f0 .+ y_edges[2:end]./2f0
    z_centers = z_edges[1:end-1]./2f0 .+ z_edges[2:end]./2f0
    result_x,result_y,result_z = Float32[],Float32[],Float32[]
    result_ions = Int64[]
    for ix in 1:nx, iy in 1:ny, iz in 1:nz
        count = hist_counts[ix,iy,iz]
        if count > 0
            push!(result_x,x_centers[ix]);push!(result_y,y_centers[iy]);push!(result_z,z_centers[iz])
            push!(result_ions,count)
        end
    end
    retained_ions = sum(result_ions)
    bins = DataFrame(event_id=fill(hits.event_id,length(result_ions)),x=result_x,y=result_y,
        z=result_z,energy=Float32.(result_ions).*ionisation_energy_eV./1f6,ions=result_ions)
    all(isfinite,bins.energy) || throw(ArgumentError("Equivalent bin energies overflow Float32"))
    return (;bins,sampled_ions,retained_ions,outside_ions=sampled_ions-retained_ions,
             x_edges_mm=x_edges,y_edges_mm=y_edges,z_edges_mm=z_edges)
end
