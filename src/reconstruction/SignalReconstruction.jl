"""Reconstruct ion and PE vectors from measured electronics responses."""
module SignalReconstruction
using LinearAlgebra, Statistics
using ..DetectorResponse, ..ElectronicsResponse
const ER=ElectronicsResponse
include("types.jl")
include("ions/reconstruction.jl")
include("ions/direct_inverse.jl")
include("electrons/reconstruction.jl")
include("electrons/rebin.jl")
include("hits.jl")
include("hit_output.jl")
export SavedHits,EventEnergy,HitEventInfo,HitEvent,HitOutputSettings,HitFileMetadata,HitWriter,HitReader,
 saved_hits,open_hit_output,write_hit_event,finish_hit_output,abort_hit_output,open_hit_input,read_hit_event,close_hit_input
export IonReconstructionCalibration,IonReconstructionSettings,IonReconstructionChecks,
 ReconstructedIons,ElectronReconstructionSettings,ReconstructedPE,
 ReconstructedHits,event_charge_estimates,reconstruct_hits,
 adc_voltages,reco_ion_signal,reco_electron_signal,rebin_electron_signal
end
