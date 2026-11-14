"""Electronics response of framed detector vectors. Load DetectorResponse first."""
module ElectronicsResponse
using Random, LinearAlgebra, FFTW
using ..DetectorResponse
include("../ion_signal/ion_csa_nonlinear_io.jl")
include("../ion_signal/sipm_response.jl")
include("../ion_signal/sipm_electronics.jl")
include("../ion_signal/bucket_noise_covariance.jl")
const NonlinearCSA=IonCSANonlinearIO.N.NonlinearCSA
const nonlinear_waveform=IonCSANonlinearIO.N.nonlinear_waveform
include("types.jl")
include("ion_signal/adc.jl")
include("ion_signal/noise.jl")
include("ion_signal/csa.jl")
include("electron_signal/sipm.jl")
export IonAcquisitionSettings,IonADCTime,ADCSettings,IonCurrentPulse,LinearCSAResponse,
       NonlinearCSA,CSAVoltageNoise,CSAResponseMode,linear_csa,hybrid_nonlinear_csa,
       CSASettings,IonElectronicsSettings,SiPM,SiPMStatisticalOptions,SiPmAsic,
       SiPMReadoutSettings,IonADCResponse,SiPMASICResponse,IonElectronicsControls,
       SiPMReadoutControls,ion_adc_time,ion_current,linear_csa_voltage,nonlinear_waveform,
       sample_csa_noise,digitise_adc,sim_ion_electronics_response,sim_sipm_readout
end
