"""Additional milliseconds before the first induced current and after collection."""
struct IonAcquisitionSettings
    baseline_ms::Float32
    signal_tail_ms::Float32
end
"""ADC start after event t0, spacing in μs, and number of samples. Use local times for arithmetic."""
struct IonADCTime
    start_ms::Float32
    sample_us::Float32
    sample_count::Int64
end
"""Ideal ADC limits and resolution, with separate Gaussian input noise in V RMS.

Codes use floor((V-minimum)/LSB), bounded to 0:2^bits-1. The exact upper
voltage is an endpoint code; values beyond the range are counted separately.
"""
struct ADCSettings
    bits::Int64
    minimum_V::Float32
    maximum_V::Float32
    sample_us::Float32
    input_noise_rms_V::Float32
end
"""Piecewise-linear induced current in A versus seconds.

For a single-ion pulse, collection is at time zero. For a complete pad vector,
times are relative to the first ADC sample. Integrals preserve finite-path charge.
"""
struct IonCurrentPulse
    time_s::Vector{Float64}
    current_A::Vector{Float64}
end
"""Signed impulse response in V/C, tabulated versus elapsed seconds, and baseline V."""
struct LinearCSAResponse
    time_s::Vector{Float64}
    voltage_per_charge_V_C::Vector{Float64}
    baseline_V::Float64
end
"""Stationary voltage covariance at successive ADC lags; element one is variance.

The vector contains a complete even periodic covariance representation for FFT
sampling. Its period must be at least twice the requested event sample count.
The spectral preparation includes the full saved noise band and sampling aliasing.
"""
struct CSAVoltageNoise
    sample_us::Float32
    covariance_V2::Vector{Float64}
end
@enum CSAResponseMode linear_csa hybrid_nonlinear_csa
"""CSA pulse, calibrated linear/nonlinear models, voltage selection threshold and noise."""
struct CSASettings
    mode::CSAResponseMode
    current_pulse::IonCurrentPulse
    linear_response::LinearCSAResponse
    nonlinear_response::NonlinearCSA
    nonlinear_below_V::Float32
    noise::CSAVoltageNoise
end
"""Ion acquisition, CSA and ADC settings supplied independently of event data."""
struct IonElectronicsSettings
    acquisition::IonAcquisitionSettings
    csa::CSASettings
    adc::ADCSettings
end
# Reuse validated numerical implementations under the reviewed public names.
const SiPM=SiPMResponse.Sensor
const SiPMStatisticalOptions=SiPMResponse.StatisticalOptions
const SiPmAsic=SiPMElectronics.Electronics
"""SiPM parameters, illuminated microcells, calculation accuracy and ASIC settings."""
struct SiPMReadoutSettings
    sipm::SiPM
    illuminated_microcells::Int64
    statistical::SiPMStatisticalOptions
    asic::SiPmAsic
end
"""ADC codes for one event: sample rows and pad columns, with channel positions and timing."""
struct IonADCResponse
    event_id::Int64
    channels::ReadoutChannels
    time::IonADCTime
    counts::Matrix{UInt16}
end
"""ASIC charge in fully recovered avalanche units, including signed noise and clipping flags."""
struct SiPMASICResponse
    event_id::Int64
    channels::ReadoutChannels
    time::ElectronReadoutTime
    response_pe::Matrix{Float32}
    clipped::BitMatrix
end
"""Counts of pad models used and samples outside ADC voltage limits; simulation controls."""
struct IonElectronicsControls
    linear_channels::Int64
    nonlinear_channels::Int64
    spice_channels::Int64
    adc_below_range::Int64
    adc_above_range::Int64
end
"""Photodetection, dark, sensor-charge and clipping totals for one event."""
struct SiPMReadoutControls
    detected_pe::Int64
    dark_pulses::Int64
    sipm_charge_pe::Float64
    asic_clipped_buckets::Int64
    asic_charge_lost_pe::Float64
end
