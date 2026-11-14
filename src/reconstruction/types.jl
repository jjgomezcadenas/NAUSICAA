"""ADC calibration, single-ion current, CSA voltage models and measured noise covariance."""
struct IonReconstructionCalibration
 adc::ADCSettings
 current_pulse::IonCurrentPulse
 linear_response::LinearCSAResponse
 nonlinear_response::NonlinearCSA
 voltage_noise::CSAVoltageNoise
end
"""Baseline and voltage-average durations in ms, nonlinear voltage threshold,
reference-fit iteration limit and convergence threshold in ions per bucket.
The direct nonlinear inverse does not use the two reference-fit settings."""
struct IonReconstructionSettings
 baseline_ms::Float32
 average_ms::Float32
 nonlinear_below_V::Float32
 max_iterations::Int64
 convergence_ions::Float32
end
"""One value per pad: measured baseline, selected model, fit convergence,
endpoint ADC-code flag and RMS residual of the fitted voltage averages."""
struct IonReconstructionChecks
 baseline_V::Vector{Float32}
 nonlinear::BitVector
 converged::BitVector
 adc_at_limit::BitVector
 voltage_residual_rms_V::Vector{Float32}
end
"""Signed reconstructed ions; rows are collection buckets and columns are pads.
Invalid estimates are NaN and have valid=false. Times retain the common event t0."""
struct ReconstructedIons
 event_id::Int64
 channels::ReadoutChannels
 time::IonReadoutTime
 ions::Matrix{Float32}
 valid::BitMatrix
 checks::IonReconstructionChecks
end
"""SiPM calibration, optically illuminated microcells and optional mean dark subtraction."""
struct ElectronReconstructionSettings
 sipm::SiPM
 illuminated_microcells::Int64
 subtract_mean_dark::Bool
end
"""Signed PE estimates before sensor charge loss. Rows are time buckets and
columns are SiPMs. acquired_end_us retains the actual acquisition end when the
last rebinned group is incomplete. Invalid/clipped estimates remain marked;
attachment correction is a later energy calculation."""
struct ReconstructedPE
 event_id::Int64
 channels::ReadoutChannels
 time::ElectronReadoutTime
 acquired_end_us::Float32
 pe::Matrix{Float32}
 valid::BitMatrix
 asic_clipped::BitMatrix
end
