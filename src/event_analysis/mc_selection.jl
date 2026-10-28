"""
    analyze_mc_event(hits::EventHits, settings::MCGraphSettings) -> NamedTuple

Reduce one MC event in source coordinates. Voxelise with the configured width,
apply the inclusive voxel-energy threshold and construct every connected
component using the existing graph and backbone algorithms. Event track length
is the maximum backbone length among components (zero for no components).
Accept exactly one component. All physical outputs are Float32 mm/MeV.

Return `event` (scalar summaries), `voxel_energy_MeV` (all occupied voxels),
`retained_voxel_energy_MeV` and `component_length_mm`. Inputs are preserved;
returned arrays are independently owned. Float64 reference sums and absolute
Float32 reduction errors are included in `event` for validation against the
12.5 keV energy-resolution scale. Errors from voxelisation/graph validation
propagate to the caller. No diffusion or fiducial selection is applied.
"""
function analyze_mc_event(hits::EventHits, settings::MCGraphSettings)
    voxels=voxelize_event(to_petit_dataframe(hits),settings.voxel_width_mm)
    tracks=make_tracks(voxels;method=settings.method,max_distance_mm=resolved_mc_distance(settings),
                       energy_threshold_kev=settings.energy_threshold_keV,k=settings.k)
    lengths=Float32[walk_track_from_extremes(t).total_length for t in tracks]
    keep=voxels.energy .>= settings.energy_threshold_keV*1f-3
    hit_energy=sum(hits.energy_MeV);voxel_energy=sum(voxels.energy)
    hit_ref=sum(Float64,hits.energy_MeV);voxel_ref=sum(Float64,voxels.energy)
    record=(event_id=hits.event_id,n_voxels=Int64(nrow(voxels)),n_retained_voxels=Int64(count(keep)),
        hit_energy_MeV=hit_energy,voxel_energy_MeV=voxel_energy,
        retained_energy_MeV=sum(voxels.energy[keep]),n_components=Int64(length(tracks)),
        track_length_mm=maximum(lengths;init=0f0),accepted=length(tracks)==1,
        hit_reference_MeV=hit_ref,voxel_reference_MeV=voxel_ref,
        hit_reduction_error_MeV=abs(Float64(hit_energy)-hit_ref),
        voxel_reduction_error_MeV=abs(Float64(voxel_energy)-hit_ref),
        energy_closure_MeV=Float64(voxel_energy)-Float64(hit_energy))
    return (event=record,voxel_energy_MeV=copy(voxels.energy),
        retained_voxel_energy_MeV=voxels.energy[keep],component_length_mm=lengths)
end
