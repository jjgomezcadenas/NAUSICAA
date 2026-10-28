"""
    collect_ion_cloud(diffused; pitch_mm, z_bin_width_mm, grid_offset_mm,
        drift_speed_mm_ms, drift_direction, minimum_xy_span_mm=160f0)

Collect the full retained ion table before reconstruction. Use its ion-weighted
x/y centre and minimum z as the reference, the supplied transverse grid phase,
and the existing ideal containing-cell hole/pad assignment. Return owned axes
and sparse pad-by-frame history. Input uses Float32 mm/MeV and Int64 ion counts;
all charge reductions in the history are Int64. Inputs are preserved. Require
positive total charge and inherit pad_axes/compute_occupancy_analysis validation.
The speed is mm/ms, frame width is set by z_bin_width_mm/speed, and the sign
sets relative arrival order. No reconstruction selection gates collection.
"""
function collect_ion_cloud(d::DataFrame;pitch_mm::Float32,z_bin_width_mm::Float32,
        grid_offset_mm::NTuple{2,Float32},drift_speed_mm_ms::Float32,drift_direction::Int,
        minimum_xy_span_mm::Float32=160f0)
    total=_diffused_total(d);total>0 || throw(ArgumentError("Collection requires positive total ions"))
    weights=Float32.(d.ions)
    reference_mm=(sum(d.x.*weights)/Float32(total),sum(d.y.*weights)/Float32(total),minimum(d.z))
    axes=pad_axes(d;reference_mm,pitch_mm,z_bin_width_mm,grid_offset_mm,minimum_xy_span_mm)
    history=compute_occupancy_analysis(d;x_edges_mm=axes.x_edges_mm,y_edges_mm=axes.y_edges_mm,
        z_edges_mm=axes.z_edges_mm,drift_speed_mm_ms,drift_direction)
    return (;axes,history)
end

"""
    pixel_frame_statistics(history; n_min::Int64)

Use every positive pixel-frame cell for inclusive counts and event extrema.
Return owned Int64 counts, maximum/minimum (nothing for an empty event), and
minimum_histogram_accepted, true only when the event minimum > n_min. The cut
selects events for the minimum histogram alone. Counts are never integrated
over frames. Preserve the history; negative n_min raises ArgumentError.
"""
function pixel_frame_statistics(history;n_min::Int64)
    n_min>=0 || throw(ArgumentError("n_min must be nonnegative"))
    counts=copy(history.counts.ions[history.counts.ions.>0])
    return (;counts,maximum_ions=isempty(counts) ? nothing : maximum(counts),
        minimum_ions=isempty(counts) ? nothing : minimum(counts),
        minimum_histogram_accepted=!isempty(counts) && minimum(counts)>n_min)
end
