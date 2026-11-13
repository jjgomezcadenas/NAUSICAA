"""Square channel lattice: pitch and lower XY edge of cell (0,0), in vertex-relative mm."""
struct ReadoutGrid
    pitch_mm::Float32
    edge_mm::NTuple{2,Float32}
    function ReadoutGrid(p::Float32,e::NTuple{2,Float32})
        isfinite(p)&&p>0&&all(isfinite,e) || throw(ArgumentError("Invalid grid"))
        new(p,e)
    end
end
"""All acquired channels, ordered by X index then Y index, with centres in mm."""
struct ReadoutChannels
    ix::Vector{Int32}
    iy::Vector{Int32}
    x_mm::Vector{Float32}
    y_mm::Vector{Float32}
end
"""Complete ion buckets; start_ms is measured from the common event t0."""
struct IonReadoutTime
    start_ms::Float32
    bucket_ms::Float32
    bucket_count::Int32
end
"""Complete electron buckets; start_us is measured from the same event t0 as ions."""
struct ElectronReadoutTime
    start_us::Float32
    bucket_us::Float32
    bucket_count::Int32
end
"""Ion drift velocity in mm/ms and event diffusion widths in mm."""
struct IonTransportSettings
    ion_drift_velocity_mm_ms::Float32
    sigma_t_mm::Float32
    sigma_l_mm::Float32
end
"""Electron drift velocity in mm/μs and event diffusion widths in mm."""
struct ElectronTransportSettings
    electron_drift_velocity_mm_us::Float32
    sigma_t_mm::Float32
    sigma_l_mm::Float32
end
@enum GasMixture Xe_100 Xe_85_He_15
"""Resolved event settings. Times are ms for ions, μs for electrons; lifetime is ms.

Diffusion widths, EL yield and optical acceptance are calculated from the run
card before collection. photon_wavelength identifies :blue_tpb or :vuv. The
existing per-event drift-distance approximation is used for widths and attachment.
"""
struct DetectorResponseSettings
    chamber_length_mm::Float32
    W_eV::Float32
    fano_factor::Float32
    gas_mixture::GasMixture
    ion_grid::ReadoutGrid
    electron_grid::ReadoutGrid
    ion_bucket_ms::Float32
    electron_bucket_us::Float32
    ion_transport::IonTransportSettings
    electron_transport::ElectronTransportSettings
    electron_lifetime_ms::Float32
    el_photons_per_electron::Float32
    el_variance_ratio::Float32
    optical_acceptance::Float32
    photon_wavelength::Symbol
end
"""Counts arranged as time buckets × channels, including zero-signal channels."""
struct DetectorVectors{T<:Union{IonReadoutTime,ElectronReadoutTime}}
    channels::ReadoutChannels
    time::T # IonReadoutTime (ms) or ElectronReadoutTime (μs)
    counts::Matrix{Int64}
end
"""Signal before acquisition cuts and losses. XY and time counts may overlap;
excluded_total counts their union once. Channel counts refer to positive signal.
"""
struct ReadoutLosses
    total::Int64
    excluded_xy::Int64
    excluded_time::Int64
    excluded_total::Int64
    channels_with_signal::Int64
    outside_channels_with_signal::Int64
end
"""Ion vectors and excluded cells for control plots. No individual ion positions are saved."""
struct IonDetectorResponse
    event_id::Int64
    placement::EventPlacement
    frame::EventFrame
    ions::DetectorVectors{IonReadoutTime}
    e_ion_pairs::Int64
    ions_outside_frame::Int64
    losses::ReadoutLosses
    excluded::SparseCells
end
"""Photon vectors at the SiPM, before PDE, with electron and optical loss counts.

uncollected_photons counts optical loss. losses counts photons reaching SiPMs
outside acquired channels or times. excluded stores those aggregate cells.
"""
struct ElectronDetectorResponse
    event_id::Int64
    placement::EventPlacement
    frame::EventFrame
    photons::DetectorVectors{ElectronReadoutTime}
    e_ion_pairs::Int64
    e_lost_attachment::Int64
    electrons_outside_frame::Int64
    EL_photons::Int64
    uncollected_photons::Int64
    losses::ReadoutLosses
    excluded::SparseCells
end
