# TOML/HDF5 adapters for the reviewed electronics interfaces.
using LinearAlgebra
include("../../src/ion_signal/csa_resources.jl")
const ER=ElectronicsResponse
const B=ER.BucketNoiseCovariance.AveragingScan
"""Load the existing circuit and SiPM cards; prepare numerical tables for this event length."""
function prepare_electronics_run(card)
 setup_started=time()
 c=TOML.parsefile(joinpath(ROOT,card["csa_card"]));resources=CSAResources.verify(ROOT,c)
 k=ER.IonCSANonlinearIO.load_calibration(joinpath(ROOT,card["nonlinear_calibration"]))
 resources["model_sha256"]==k.circuit_sha256||error("CSA circuit/calibration mismatch")
 # Same 0.5 mm pad, 5 mm finite gap and 3 m/s prescription as the saved calibration.
 a=c["pad_side_m"]/2;g=c["gap_m"];v=c["gap_speed_m_s"]
 (a,g,v)==(.00025,.005,3.)||error("Revalidate changed collection geometry")
 t=collect(range(-g/v,0.;length=ceil(Int,g/v/c["current_step_s"])+1))
 current=[ER.E_CHARGE*v*4a^2/(pi*(z^2+a^2)*sqrt(z^2+2a^2)) for z in -v.*t]
 pulse=IonCurrentPulse(t,current)
 adc=ADCSettings(card["adc_bits"],Float32(card["adc_minimum_V"]),Float32(card["adc_maximum_V"]),Float32(card["adc_sample_us"]),Float32(card["adc_input_noise_rms_V"]))
 acquisition=IonAcquisitionSettings(Float32(card["baseline_ms"]),Float32(card["signal_tail_ms"]))
 calpath=joinpath(ROOT,c["calibration_directory"],"analysis.h5")
 rawpath=joinpath(ROOT,c["calibration_directory"],"variants",c["variant"]*".h5")
 f,asd,tau=h5open(calpath,"r") do h
  group=h[c["variant"]];(read(group["frequency_Hz"]),read(group["output_asd_V_sqrtHz"]),read(attributes(group["20ms"])["tau_s"]))
 end
 # Derive the amplitude from the saved linear 100-ion waveform at the end of collection.
 tref=.02;pre=.01
 vr=h5open(rawpath,"r") do h
  group=h["20ms"];signal=group["single_N100_r2"];zero=group["single_N0_r2"]
  ER.interpolate(read(signal["time"]),read(signal["V(out)"]),pre+tref)-
  ER.interpolate(read(zero["time"]),read(zero["V(out)"]),pre+tref)
 end
 # Same uniform finite-gap packet used in existing CSA calibration; analytical weighting potential.
 phi(z)=z<=0 ? 1. : 2/pi*atan(a^2/(z*sqrt(z^2+2a^2)))
 cumulative(u)=ER.E_CHARGE*(phi(clamp(-v*u,0,g))-phi(g))
 ut=sort!(unique(vcat(t,tref.+t)))
 ui=[(cumulative(x)-cumulative(x-tref))/tref for x in ut]
 # Gauss integrate the exponential response to the piecewise-linear input.
 nodes,weights=let eig=eigen(SymTridiagonal(zeros(8),[j/sqrt(4j^2-1) for j in 1:7]));(eig.values,2eig.vectors[1,:].^2);end
 integral=0.
 for j in 1:length(ut)-1
  mid=(ut[j]+ut[j+1])/2;half=(ut[j+1]-ut[j])/2
  integral+=half*sum(w*ER.interpolate(ut,ui,mid+half*x)*exp(-(tref-mid-half*x)/tau) for (x,w) in zip(nodes,weights))
 end
 gain=-vr/(100integral)
 ion_setup_s=time()-setup_started;setup_started=time()
 sc=TOML.parsefile(joinpath(ROOT,card["sipm_card"]));ac=TOML.parsefile(joinpath(ROOT,card["asic_card"]))
 sipm=SiPMReadoutSettings(SiPM(sc["pde"],sc["microcell_count"],sc["gain"],sc["recovery_s"],sc["dark_reference_rate_Hz"]),
 Int64(sc["illuminated_microcell_count"]),SiPMStatisticalOptions(sc["statistical_substeps"],sc["sparse_trigger_limit"];initial_state=Symbol(get(sc,"initial_state","stationary_dark"))),
 SiPmAsic(ac["bucket_width_s"],ac["noise_rms_pe"],ac["current_limit_A"]))
 sipm_setup_s=time()-setup_started
 (;ion_setup_s,sipm_setup_s,sipm,resources,nonlinear=k,frequency=f,asd,tau,gain,rawpath,pulse,adc,acquisition,cache=Dict{Int,Any}())
end

"""Reuse calibration and noise tables for events with the same covariance period.

The cache belongs to one run and its fixed cards. Clear it before retained tables
would exceed 64 MiB. A changed card requires a new prepared run.
"""
function electronics_settings(prepared::NamedTuple,ion_time;max_samples=nothing)
 at=ion_adc_time(ion_time,prepared.pulse,prepared.acquisition,prepared.adc)
 n=isnothing(max_samples) ? at.sample_count : max_samples
 capacity=nextpow(2,n)
 if !haskey(prepared.cache,capacity)
  duration=(capacity-1)*Float64(prepared.adc.sample_us)*1e-6
  lt=collect(range(0.,max(duration,.2);step=prepared.adc.sample_us*1e-6/4))
  last(lt)<duration&&push!(lt,duration)
  linear=LinearCSAResponse(lt,-prepared.gain.*exp.(-lt./prepared.tau),prepared.nonlinear.model.baseline_V)
  noise=ER.voltage_noise_from_spectrum(prepared.frequency,prepared.asd,prepared.adc.sample_us,capacity)
  bytes=sizeof(lt)+sizeof(linear.voltage_per_charge_V_C)+sizeof(noise.covariance_V2)
  bytes<=64*1024^2||error("Prepared CSA tables exceed 64 MiB")
  Base.summarysize(prepared.cache)+bytes>64*1024^2&&empty!(prepared.cache)
  prepared.cache[capacity]=IonElectronicsSettings(prepared.acquisition,
   CSASettings(hybrid_nonlinear_csa,prepared.pulse,linear,prepared.nonlinear.model,Float32(prepared.nonlinear.guard_V),noise),prepared.adc)
 end
 (;prepared...,ion=prepared.cache[capacity])
end
"""Standalone settings adapter; production retains prepare_electronics_run(card)."""
electronics_settings(card,ion_time;max_samples=nothing)=electronics_settings(prepare_electronics_run(card),ion_time;max_samples)
