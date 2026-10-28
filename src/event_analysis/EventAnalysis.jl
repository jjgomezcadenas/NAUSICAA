"""
Event-level data and geometry for the NAUSICAA analysis.

Load with `include("src/event_analysis/EventAnalysis.jl")`, followed by
`using .EventAnalysis`. Physical hit fields use Float32 with positions in mm,
energy in MeV and MC deposition time in ns.
"""
module EventAnalysis

using DataFrames, Graphs, NearestNeighbors, LinearAlgebra, SparseArrays, Statistics, HDF5, Random
import StatsBase

export CoordinateFrame, source_frame, vertex_frame, EventInfo, EventHits, center_hits,
       read_event, to_petit_dataframe, voxelize_event, DiffusionParams, Tracks, TrackCoords,
       make_tracks, is_single_track, get_voxel_size_and_distance, compute_mc_path,
       find_track_extremes, walk_track_from_extremes, get_raw_path, find_extremes_mst_diameter,
       ion_counts, sigma_t_ion_mm, sample_drift_cm, sample_grid_offset_mm, diffuse_xyz_image_mc,
       select_diffused_event, pad_axes, compute_occupancy_analysis, pad_projections,
       align_pad_history, ensemble_pad_maps, write_ensemble, read_ensemble,
       MCEnergyMaps, bin_mc_energy

include("types.jl")
include("coordinates.jl")
include("../ion_signal/event_io.jl")
include("io.jl")
include("reconstruction.jl")
include("petit/graphs.jl")
include("petit/backbone.jl")
include("petit/path.jl")
include("petit/truth.jl")
include("paths.jl")
include("production.jl")
include("diffusion.jl")
include("collection.jl")
include("ensemble.jl")
include("display.jl")

end
