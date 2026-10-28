"""
    MCEnergyMaps

Owned energy maps of one MC event, in its supplied coordinate frame. Explicit
Float32 x/y/z edge vectors define every display cell, in mm. Energy matrices
are Float32 MeV sums with Julia shapes (nx,ny), (nx,nz), (ny,nz); the matching
Int64 hit-count matrices count contributing hits. `occupied` is a sparse
DataFrame with Int64 ix/iy/iz (one-based), Float32 energy_MeV and Int64 hit_count.
It includes cells containing zero-energy hits. Empty cells have zero energy/count.
`total_energy_MeV` sums original hit energies; `hit_count` counts original hits.
See [`bin_mc_energy`](@ref) for binning and precision conventions.
"""
struct MCEnergyMaps
    x_edges_mm::Vector{Float32}
    y_edges_mm::Vector{Float32}
    z_edges_mm::Vector{Float32}
    xy_energy_MeV::Matrix{Float32}
    xz_energy_MeV::Matrix{Float32}
    yz_energy_MeV::Matrix{Float32}
    xy_hit_count::Matrix{Int64}
    xz_hit_count::Matrix{Int64}
    yz_hit_count::Matrix{Int64}
    occupied::DataFrame
    total_energy_MeV::Float32
    hit_count::Int64
end

"""Select Petit/StatsBase display edges; a coincident axis gets one 0.5 mm cell.
Convert edges to Float32, extend endpoints by one representable step if needed
to retain all hits, and validate increasing finite edges. The number of bins
is the realised StatsBase count rather than necessarily the requested target.
"""
function _display_edges(values::Vector{Float32}, nbins::Int)
    lo,hi=extrema(values)
    edges = lo==hi ? Float32[lo-0.25f0,hi+0.25f0] :
        Float32.(collect(only(StatsBase.fit(StatsBase.Histogram,values;nbins).edges)))
    if first(edges)>lo;edges[1]=prevfloat(lo);end
    if last(edges)<hi;edges[end]=nextfloat(hi);end
    _check_edges(edges)
    return edges
end

"""
    bin_mc_energy(hits::EventHits; nbins::Int=100) -> MCEnergyMaps

Sum deposited MC energy and contributing hit counts on a common display grid.
Use Petit's StatsBase histogram edge-selection target `nbins` independently on
each coordinate. The realised edges are shared by the three projections and
the occupied 3D voxels. This adapts Petit src/plots.jl to the requested energy
heatmaps. Graph voxelisation is independently controlled by its physical width.

Every hit contributes once to a sparse 3D cell. Each projection sums these cell
energies/counts along the omitted coordinate. Internal bins are half open; the
final upper edge belongs to the final bin. Coincident coordinates use one cell
of width 0.5 mm. Arithmetic and owned physical arrays are Float32; counts are
Int64. Input coordinates, identifiers and energy arrays retain their values.

Require nbins > 0 and valid finite EventHits data. Invalid input or edges that
collapse at Float32 precision raise ArgumentError. Nonfinite summed energies
raise ArgumentError. Projection sums agree with event energy within Float32
reduction rounding; validation compares with independent Float64 references.

# Example
```julia
h = EventHits(1, source_frame, Float32[0,1], zeros(Float32,2), zeros(Float32,2),
    Float32[0.1,0.2], zeros(Float32,2), Int32[1,1], Int32[0,1])
m = bin_mc_energy(h;nbins=10)
@assert m.hit_count == 2 && sum(m.xy_hit_count) == 2
@assert isapprox(sum(m.xy_energy_MeV), 0.3f0)
```
"""
function bin_mc_energy(hits::EventHits;nbins::Int=100)::MCEnergyMaps
    nbins>0 || throw(ArgumentError("nbins must be positive"))
    _validate_hit_vectors(hits.x_mm,hits.y_mm,hits.z_mm,hits.energy_MeV,hits.time_ns,hits.particle_id,hits.hit_id)
    xe,ye,ze = (_display_edges(v,nbins) for v in (hits.x_mm,hits.y_mm,hits.z_mm))
    cells=Dict{NTuple{3,Int64},Tuple{Float32,Int64}}()
    for r in eachindex(hits.energy_MeV)
        key=(_histogram_index(hits.x_mm[r],xe),_histogram_index(hits.y_mm[r],ye),_histogram_index(hits.z_mm[r],ze))
        e,n=get(cells,key,(0f0,Int64(0)))
        cells[key]=(e+hits.energy_MeV[r],n+1)
    end
    occupied=DataFrame(ix=Int64[],iy=Int64[],iz=Int64[],energy_MeV=Float32[],hit_count=Int64[])
    nx,ny,nz=length(xe)-1,length(ye)-1,length(ze)-1
    xy,xz,yz=zeros(Float32,nx,ny),zeros(Float32,nx,nz),zeros(Float32,ny,nz)
    nxy,nxz,nyz=zeros(Int64,nx,ny),zeros(Int64,nx,nz),zeros(Int64,ny,nz)
    for key in sort!(collect(keys(cells)))
        i,j,k=key;e,n=cells[key]
        push!(occupied,(i,j,k,e,n))
        xy[i,j]+=e;xz[i,k]+=e;yz[j,k]+=e
        nxy[i,j]+=n;nxz[i,k]+=n;nyz[j,k]+=n
    end
    total=sum(hits.energy_MeV)
    isfinite(total) && all(isfinite,xy) && all(isfinite,xz) && all(isfinite,yz) ||
        throw(ArgumentError("Display energy sums overflow Float32"))
    return MCEnergyMaps(xe,ye,ze,xy,xz,yz,nxy,nxz,nyz,occupied,total,Int64(length(hits.energy_MeV)))
end
