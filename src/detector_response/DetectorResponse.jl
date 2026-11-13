"""Framed ion and photon response. Load EventAnalysis before this module."""
module DetectorResponse
using Random, Statistics
using ..EventAnalysis
import ..EventAnalysis: sample_fano_pairs
export ReadoutGrid, ReadoutChannels, IonReadoutTime, ElectronReadoutTime,
       IonTransportSettings, ElectronTransportSettings, GasMixture, Xe_100, Xe_85_He_15,
       DetectorResponseSettings, DetectorVectors, ReadoutLosses, IonDetectorResponse,
       ElectronDetectorResponse, acquisition_frame, hits_in_frame, readout_channels,
       ion_readout_time, electron_readout_time, sample_fano_pairs,
       sim_ion_detector_response, sim_electron_detector_response, vector_statistics
include("types.jl")
include("frame.jl")
include("collection.jl")
include("ions/sim_ion_detector_response.jl")
include("electrons/sim_electron_detector_response.jl")
end
