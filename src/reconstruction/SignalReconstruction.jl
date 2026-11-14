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
export IonReconstructionCalibration,IonReconstructionSettings,IonReconstructionChecks,
 ReconstructedIons,ElectronReconstructionSettings,ReconstructedPE,
 adc_voltages,reco_ion_signal,reco_electron_signal,rebin_electron_signal
end
