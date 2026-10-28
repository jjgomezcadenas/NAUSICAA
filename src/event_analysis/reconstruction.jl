"""
    DiffusionParams()
    DiffusionParams(drift_length_cm, sigma_t_mm, sigma_l_mm, bin_width_mm,
                    max_distance_mm, energy_threshold_keV, histogram_nbins, histogram_nsigma)

Per-event diffusion and reconstruction settings attached to a `Tracks` result,
following Petit's `DiffusionParams`. All physical fields are Float32. Drift
length is in cm; widths and connection distance in mm; threshold in keV.
`histogram_nbins::Int` and `histogram_nsigma::Float32` describe the intermediate
diffusion histogram. The zero-argument constructor supplies Petit's neutral
metadata defaults for a standalone graph calculation. The caller supplies
evaluated settings for an event; computation uses its explicit arguments.

This immutable record owns scalar values. The event record is the source of
settings when writing a future run output. See `PROVENANCE.toml` for Petit origin.
"""
struct DiffusionParams
    drift_length_cm::Float32
    sigma_t_mm::Float32
    sigma_l_mm::Float32
    bin_width_mm::Float32
    max_distance_mm::Float32
    energy_threshold_keV::Float32
    histogram_nbins::Int
    histogram_nsigma::Float32
end
DiffusionParams() = DiffusionParams(0f0, 0f0, 0f0, 1f0, 1.5f0, 0f0, 100, 3f0)

"""
    Tracks(voxels::DataFrame, graph::SimpleGraph{Int},
           components::Vector{Vector{Int}}, diffusion::DiffusionParams)

One connected component returned by Petit's graph construction. `voxels` has
Float32 x/y/z (mm) and energy (MeV), with an optional repeated event_id::Int64.
Local graph vertex i identifies row i in this table. `components` contains the
component's row indices in the energy-filtered parent table, before induced
subgraph renumbering. Both conventions follow Petit.

`make_tracks` creates independent component tables, graphs and index arrays,
sorted by decreasing component energy. The direct constructor stores supplied
objects; callers constructing fixtures ensure graph/table alignment. Stored
arrays and graphs are mutable. Backbone calculations preserve their values.
"""
struct Tracks
    voxels::DataFrame
    graph::SimpleGraph{Int}
    components::Vector{Vector{Int}}
    diffusion::DiffusionParams
end

"""Validate a single-event Float32 x/y/z/energy table before graph or binning work."""
function validate_bin_table(df::DataFrame)
    for name in (:x, :y, :z, :energy)
        hasproperty(df, name) || throw(ArgumentError("Missing column $name"))
        eltype(df[!, name]) == Float32 || throw(ArgumentError("$name must have Float32 values"))
        all(isfinite, df[!, name]) || throw(ArgumentError("$name must be finite"))
    end
    all(>=(0f0), df.energy) || throw(ArgumentError("energy must be nonnegative"))
    if hasproperty(df, :event_id)
        eltype(df.event_id) == Int64 || throw(ArgumentError("event_id must be Int64"))
        length(unique(df.event_id)) <= 1 || throw(ArgumentError("Expected one event"))
    end
    return nothing
end

"""
    to_petit_dataframe(hits::EventHits) -> DataFrame

Adapt one event to Petit's column names. Return event_id::Int64, x/y/z/energy/
time::Float32 and particle_id/hit_id::Int32 columns, with units mm/MeV/ns.
Every row corresponds to the original hit at the same index. The frame is
supplied by `hits.frame`; numerical coordinates retain that frame. The result
owns new columns and preserves the input arrays. Validate current hit values
using the same contract as `EventHits` (including unique particle/hit pairs).

# Example
```julia
h = EventHits(1, source_frame, Float32[1], Float32[2], Float32[3],
              Float32[0.01], Float32[0], Int32[7], Int32[0])
d = to_petit_dataframe(h)
@assert d.x == h.x_mm && d.x !== h.x_mm
```
"""
function to_petit_dataframe(hits::EventHits)::DataFrame
    _validate_hit_vectors(hits.x_mm, hits.y_mm, hits.z_mm, hits.energy_MeV,
                         hits.time_ns, hits.particle_id, hits.hit_id)
    return DataFrame(event_id=fill(hits.event_id, length(hits.x_mm)),
        x=hits.x_mm, y=hits.y_mm, z=hits.z_mm, energy=hits.energy_MeV,
        time=hits.time_ns, particle_id=hits.particle_id, hit_id=hits.hit_id)
end

"""
    voxelize_event(event::DataFrame, bin_width_mm::Float32) -> DataFrame

Bin a single-event x/y/z/energy table using Petit's fixed convention:
`floor(Int, coordinate / bin_width_mm)`, with centre `(index + 0.5) * width`.
Coordinates are in mm; energy in MeV. Sum Float32 energy in input-row order for
each occupied cell. Return event_id::Int64 and x/y/z/energy::Float32 columns;
row order follows dictionary iteration. Absent event_id is represented by 0.
When an ions::Vector{Int64} column is present, sum its nonnegative counts into
an additional ions column using checked integer addition. These counts give the
exact reconstruction threshold and charge ledger alongside Petit's energy sums.

Adapted from Petit `src/voxels.jl:voxelize_event`; ion production belongs to its
own stage. This operation preserves input tables and allocates new output
columns. Width must be finite and positive; validate finite Float32 coordinates
and nonnegative energies. Empty inputs yield a typed empty result. Unrepresentable
integer indices raise a conversion error; nonfinite output sums/centres raise
`ArgumentError`. Negative coordinates use mathematical floor.

# Example
```julia
using DataFrames
d = DataFrame(x=Float32[-0.2, 0.2, 0.8], y=zeros(Float32,3),
              z=zeros(Float32,3), energy=Float32[1,2,3])
v = sort(voxelize_event(d, 1f0), :x)
@assert v.x == Float32[-0.5, 0.5]
@assert v.energy == Float32[1,5]
```
"""
function voxelize_event(event::DataFrame, bin_width_mm::Float32)::DataFrame
    validate_bin_table(event)
    isfinite(bin_width_mm) && bin_width_mm > 0f0 || throw(ArgumentError("bin_width_mm must be finite and positive"))
    id = hasproperty(event, :event_id) && nrow(event) > 0 ? event.event_id[1] : Int64(0)
    ix = floor.(Int, event.x ./ bin_width_mm)
    iy = floor.(Int, event.y ./ bin_width_mm)
    iz = floor.(Int, event.z ./ bin_width_mm)
    sums = Dict{Tuple{Int,Int,Int},Float32}()
    hasions = hasproperty(event,:ions)
    ion_sums = Dict{Tuple{Int,Int,Int},Int64}()
    if hasions
        eltype(event.ions) == Int64 && all(>=(0),event.ions) ||
            throw(ArgumentError("ions must be nonnegative Int64 counts"))
        for (key,n) in zip(zip(ix,iy,iz),event.ions)
            ion_sums[key] = Base.Checked.checked_add(get(ion_sums,key,Int64(0)),n)
        end
    end
    for (key, energy) in zip(zip(ix, iy, iz), event.energy)
        sums[key] = get(sums, key, 0f0) + energy
    end
    rows = NamedTuple{(:event_id,:x,:y,:z,:energy),Tuple{Int64,Float32,Float32,Float32,Float32}}[]
    counts = Int64[]
    for ((i,j,k), energy) in sums
        push!(rows, (event_id=id, x=(Float32(i)+0.5f0)*bin_width_mm,
                    y=(Float32(j)+0.5f0)*bin_width_mm, z=(Float32(k)+0.5f0)*bin_width_mm,
                    energy=energy))
        hasions && push!(counts,ion_sums[(i,j,k)])
    end
    result = DataFrame(rows)
    hasions && (result.ions = counts)
    validate_bin_table(result)
    return result
end

"""Euclidean distance in mm between two Float32 Cartesian points (Petit helper)."""
euclidean_distance(x1, y1, z1, x2, y2, z2) = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)

"""
    get_voxel_size_and_distance(sigma_t_mm::Float32;
        bin_sigma_scale=3f0, edge_sigma_scale=6f0, truth_bin_width_mm=0.5f0)

Return `(bin_width_mm, truth_bin_width_mm, max_distance_mm)::NTuple{3,Float32}`.
Use b = bin_sigma_scale × sigma_t_mm and d = edge_sigma_scale × sigma_t_mm;
defaults implement the agreed b = 3σt and d = 6σt for all drift distances.
The return order follows Petit. Arguments are finite and positive and the
products must remain finite. Settings are explicit; this deterministic helper
performs no I/O or sampling. The truth width is independently configurable.
"""
function get_voxel_size_and_distance(sigma_t_mm::Float32;
        bin_sigma_scale::Float32=3f0, edge_sigma_scale::Float32=6f0,
        truth_bin_width_mm::Float32=0.5f0)::NTuple{3,Float32}
    all(x -> isfinite(x) && x > 0f0, (sigma_t_mm,bin_sigma_scale,edge_sigma_scale,truth_bin_width_mm)) ||
        throw(ArgumentError("Widths and scale factors must be finite and positive"))
    b, d = bin_sigma_scale*sigma_t_mm, edge_sigma_scale*sigma_t_mm
    isfinite(b) && isfinite(d) || throw(ArgumentError("Scaled widths overflow"))
    return (b, truth_bin_width_mm, d)
end

"""
    is_single_track(tracks::AbstractVector{Tracks}) -> Bool

Return true exactly when graph construction returned one connected component.
An empty or multi-component event fails. Invoke on reconstructed, diffused
input when integrating the event pipeline; the predicate itself is deterministic
and preserves its input. This is Petit's RECO event-level selection criterion.
"""
is_single_track(tracks::AbstractVector{Tracks}) = length(tracks) == 1
