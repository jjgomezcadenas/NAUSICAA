"""
    align_pad_history(history, reference_mm::NTuple{3,Float32}) -> DataFrame

Align a compute_occupancy_analysis history by the pad containing the selected
track's energy-weighted x/y reference. That pad has relative integer indices
(pad_dx,pad_dy)=(0,0). Keep every ion in its original pad; this is a relabelling
of physical collection cells. The random position of the track centre within
its central pad is retained statistically across events.

Return owned Int64 columns pad_dx, pad_dy, arrival_bin, ions. arrival_bin is
slice_index-1: its zero edge is the first saved relative arrival edge, and bins
follow the history's explicit drift direction. Occupied entries remain sparse;
omitted cells/slices represent zero. The supplied reference must be finite and
covered transversely by the saved axes. Input arrays retain their values.
The third reference coordinate is recorded upstream; longitudinal alignment
uses the saved first arrival edge. This convention defines the ensemble maps.

# Example
```julia
using DataFrames
d = DataFrame(x=Float32[0],y=Float32[0],z=Float32[2],energy=Float32[0.001],ions=Int64[10])
h = compute_occupancy_analysis(d;x_edges_mm=Float32[-3,-1,1,3],y_edges_mm=Float32[-1,1],
    z_edges_mm=Float32[0,2,4],drift_speed_mm_ms=0.1f0,drift_direction=1)
a = align_pad_history(h,(0f0,0f0,0f0))
@assert a.pad_dx == Int64[0] && a.arrival_bin == Int64[1] && a.ions == Int64[10]
```
"""
function align_pad_history(history,reference_mm::NTuple{3,Float32})::DataFrame
    all(isfinite,reference_mm) || throw(ArgumentError("The track reference must be finite"))
    i = _histogram_index(reference_mm[1],history.x_edges_mm)
    j = _histogram_index(reference_mm[2],history.y_edges_mm)
    min(i,j)>0 || throw(ArgumentError("The transverse reference must lie in the saved grid"))
    return DataFrame(pad_dx=history.counts.pad_ix.-Int64(i),pad_dy=history.counts.pad_iy.-Int64(j),
        arrival_bin=history.counts.slice_index.-Int64(1),ions=copy(history.counts.ions))
end

"""
    ensemble_pad_maps(aligned::AbstractVector{DataFrame}; processed_events::Int,
        pitch_mm::Float32, z_bin_width_mm::Float32, drift_speed_mm_ms::Float32,
        minimum_half_cells::Int=40)

Sum aligned accepted-event histories into common transverse and arrival axes.
Each input table has Int64 pad_dx/pad_dy/arrival_bin/ions columns as returned by
align_pad_history; arrival bins and counts are nonnegative. processed_events
includes rejected events and must be positive and ≥ length(aligned).

Return owned arrays in a NamedTuple: x/y_edges_mm::Vector{Float32} (central pad
centred at zero), arrival_edges_mm::Vector{Float32}, time_edges_ms::Vector{Float32},
xy_sum_ions::Matrix{Int64} (x first), arrival_sum_ions::Vector{Int64}; Float32
xy_mean_accepted, xy_mean_processed, arrival_mean_accepted, arrival_mean_processed;
accepted_events::Int64, processed_events::Int64 and total_ions::Int64.

The numerator is all retained charge from accepted events. The accepted mean
divides by accepted_events. The processed mean divides by processed_events, with
rejected events contributing zero selected charge. If none pass, accepted means
are NaN (undefined), processed means and sums are zero, and the counts record this
explicitly. Missing accepted-event cells contribute zero in both denominators.

Axes extend to every occupied cell, with at least [-minimum_half_cells,+minimum_half_cells]
pad centres in x/y and one arrival bin. Widths and speed are finite and positive.
Input tables must already share pitch, slice width and time-order convention;
the run orchestrator supplies these from one resolved detector/card. Exact Int64
addition checks overflow. Mean arithmetic is Float32. Preserve all inputs.

# Example
```julia
using DataFrames
a = DataFrame(pad_dx=Int64[0],pad_dy=Int64[0],arrival_bin=Int64[0],ions=Int64[12])
e = ensemble_pad_maps([a];processed_events=3,pitch_mm=2f0,z_bin_width_mm=2f0,
    drift_speed_mm_ms=0.1f0,minimum_half_cells=0)
@assert e.xy_mean_accepted[1,1] === 12f0 && e.xy_mean_processed[1,1] === 4f0
@assert e.time_edges_ms == Float32[0,20]
```
"""
function ensemble_pad_maps(aligned::AbstractVector{DataFrame};processed_events::Int,
        pitch_mm::Float32,z_bin_width_mm::Float32,drift_speed_mm_ms::Float32,minimum_half_cells::Int=40)
    processed_events > 0 && processed_events >= length(aligned) || throw(ArgumentError("A positive processed count must include every accepted event"))
    all(x->isfinite(x)&&x>0f0,(pitch_mm,z_bin_width_mm,drift_speed_mm_ms)) || throw(ArgumentError("Widths and speed must be finite and positive"))
    minimum_half_cells>=0 || throw(ArgumentError("minimum_half_cells must be nonnegative"))
    xy = Dict{Tuple{Int64,Int64},Int64}(); arrival = Dict{Int64,Int64}()
    total_ions = Int64(0)
    for table in aligned
        all(s->hasproperty(table,s)&&eltype(table[!,s])==Int64,(:pad_dx,:pad_dy,:arrival_bin,:ions)) ||
            throw(ArgumentError("Aligned columns must be Int64"))
        all(>=(0),table.arrival_bin) && all(>=(0),table.ions) || throw(ArgumentError("Arrival bins and counts must be nonnegative"))
        for r in eachrow(table)
            key = (r.pad_dx,r.pad_dy)
            xy[key] = Base.Checked.checked_add(get(xy,key,Int64(0)),r.ions)
            arrival[r.arrival_bin] = Base.Checked.checked_add(get(arrival,r.arrival_bin,Int64(0)),r.ions)
            total_ions = Base.Checked.checked_add(total_ions,r.ions)
        end
    end
    xmin,xmax,ymin,ymax = -minimum_half_cells,minimum_half_cells,-minimum_half_cells,minimum_half_cells
    for (x,y) in keys(xy);xmin=min(xmin,x);xmax=max(xmax,x);ymin=min(ymin,y);ymax=max(ymax,y);end
    zmax = isempty(arrival) ? 0 : maximum(keys(arrival))
    x_edges_mm = (Float32.(xmin:xmax+1).-0.5f0).*pitch_mm
    y_edges_mm = (Float32.(ymin:ymax+1).-0.5f0).*pitch_mm
    arrival_edges_mm = Float32.(0:zmax+1).*z_bin_width_mm
    time_edges_ms = arrival_edges_mm./drift_speed_mm_ms
    foreach(_check_edges,(x_edges_mm,y_edges_mm,arrival_edges_mm,time_edges_ms))
    xy_sum_ions = zeros(Int64,xmax-xmin+1,ymax-ymin+1)
    arrival_sum_ions = zeros(Int64,zmax+1)
    for ((x,y),n) in xy;xy_sum_ions[x-xmin+1,y-ymin+1]=n;end
    for (z,n) in arrival;arrival_sum_ions[z+1]=n;end
    accepted_events = Int64(length(aligned))
    xy_mean_accepted = accepted_events>0 ? Float32.(xy_sum_ions)./Float32(accepted_events) : fill(NaN32,size(xy_sum_ions))
    arrival_mean_accepted = accepted_events>0 ? Float32.(arrival_sum_ions)./Float32(accepted_events) : fill(NaN32,size(arrival_sum_ions))
    return (;x_edges_mm,y_edges_mm,arrival_edges_mm,time_edges_ms,xy_sum_ions,arrival_sum_ions,
        xy_mean_accepted,arrival_mean_accepted,xy_mean_processed=Float32.(xy_sum_ions)./Float32(processed_events),
        arrival_mean_processed=Float32.(arrival_sum_ions)./Float32(processed_events),
        accepted_events,processed_events=Int64(processed_events),total_ions)
end

"""
    write_ensemble(group::Union{HDF5.File,HDF5.Group}, result)

Write an ensemble_pad_maps result to a new/empty HDF5 group using schema
nausicaa.event.ensemble.v1. Save every NamedTuple field with its Julia datatype,
explicit axes, two denominators and exact count numerator. Store the alignment
and matrix axis order as attributes. Preserve the caller's open handle and data.
Existing field names raise the HDF5 creation error. Return nothing.
"""
function write_ensemble(group::Union{HDF5.File,HDF5.Group},result)
    attributes(group)["schema"] = "nausicaa.event.ensemble.v1"
    attributes(group)["alignment"] = "central pad contains selected-track x/y reference; arrival zero is first saved edge"
    attributes(group)["matrix_axes"] = "Julia x,y; HDF5 physical y,x"
    for (name,value) in pairs(result);group[string(name)] = value;end
    return nothing
end

"""
    read_ensemble(group::Union{HDF5.File,HDF5.Group}) -> NamedTuple

Read the explicit axes, exact numerator, means and denominators of an ensemble
v1 group. Return owned arrays, preserving Float32/Int64 storage and the open
file. Validate schema, axis lengths, matrix/vector shapes, positive processed
count and the stored integer charge ledger. An inconsistent group raises error.
"""
function read_ensemble(group::Union{HDF5.File,HDF5.Group})
    read(attributes(group)["schema"]) == "nausicaa.event.ensemble.v1" || error("Unsupported ensemble schema")
    fields = (:x_edges_mm,:y_edges_mm,:arrival_edges_mm,:time_edges_ms,:xy_sum_ions,:arrival_sum_ions,
        :xy_mean_accepted,:arrival_mean_accepted,:xy_mean_processed,:arrival_mean_processed,
        :accepted_events,:processed_events,:total_ions)
    r = NamedTuple{fields}(Tuple(read(group[string(s)]) for s in fields))
    for s in (:x_edges_mm,:y_edges_mm,:arrival_edges_mm,:time_edges_ms)
        eltype(getproperty(r,s))==Float32 || error("Axis type must be Float32")
        _check_edges(getproperty(r,s))
    end
    nx,ny,nz = length(r.x_edges_mm)-1,length(r.y_edges_mm)-1,length(r.arrival_edges_mm)-1
    length(r.time_edges_ms)==nz+1 || error("Time and arrival axes differ")
    for s in (:xy_sum_ions,:xy_mean_accepted,:xy_mean_processed)
        size(getproperty(r,s))==(nx,ny) || error("Invalid transverse shape")
    end
    for s in (:arrival_sum_ions,:arrival_mean_accepted,:arrival_mean_processed)
        length(getproperty(r,s))==nz || error("Invalid arrival shape")
    end
    eltype(r.xy_sum_ions)==eltype(r.arrival_sum_ions)==Int64 || error("Count type must be Int64")
    for s in (:xy_mean_accepted,:xy_mean_processed,:arrival_mean_accepted,:arrival_mean_processed)
        eltype(getproperty(r,s))==Float32 || error("Mean type must be Float32")
    end
    r.processed_events>0 && 0<=r.accepted_events<=r.processed_events || error("Invalid denominators")
    all(s->getproperty(r,s) isa Int64,(:accepted_events,:processed_events,:total_ions)) || error("Denominators and total must be Int64")
    all(>=(0),r.xy_sum_ions) && all(>=(0),r.arrival_sum_ions) || error("Negative ensemble charge")
    sum(r.xy_sum_ions)==sum(r.arrival_sum_ions)==r.total_ions || error("Ensemble charge ledger differs")
    for (counts,accepted,processed) in ((r.xy_sum_ions,r.xy_mean_accepted,r.xy_mean_processed),
                                      (r.arrival_sum_ions,r.arrival_mean_accepted,r.arrival_mean_processed))
        isequal(processed,Float32.(counts)./Float32(r.processed_events)) || error("Processed mean differs from its numerator/denominator")
        expected = r.accepted_events>0 ? Float32.(counts)./Float32(r.accepted_events) : fill(NaN32,size(counts))
        isequal(accepted,expected) || error("Accepted mean differs from its numerator/denominator")
    end
    return r
end
