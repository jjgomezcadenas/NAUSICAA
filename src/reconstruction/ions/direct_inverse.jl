"""Prepared integral of reciprocal gain and the linear collected-charge inverse.

Gain is V per elementary charge. F(x)=integral(1/g(x),x) therefore has units
of elementary charges. The response matrix describes cumulative induced charge
from each collection bucket, averaged at the same times as measured samples.
"""
struct DirectIonInverse
 model::NonlinearCSA
 linear::ER.ExponentialCSA
 integral_gain::Vector{Float64}
 current_response::Matrix{Float64}
 inverse::Matrix{Float64}
end

"""Integral of 1/g across a segment with linearly varying positive gain.
Use log1p and its small-argument limit for nearly constant gain."""
function reciprocal_gain_integral(width,gain,slope)
 u=slope*width/gain
 abs(u)<1e-8 ? width/gain*(1-u/2+u^2/3) : log1p(u)/slope
end

"""Prepare F at the gain-table knots and the cumulative-current response once.
The finite single-ion pulse integral is preserved; it is not renormalised to e.
This solve contains no fitted waveform or iterated gain estimate."""
function prepare_direct_inverse!(p::IonReconstructionPlan)
 p.direct_inverse!==nothing&&return p.direct_inverse
 model=p.calibration.nonlinear_response
 linear=ER.exponential_parameters(p.calibration.linear_response)
 linear===nothing&&throw(ArgumentError("Exponential small-signal calibration required"))
 F=zeros(length(model.knots_V))
 for j in 1:length(F)-1
  width=model.knots_V[j+1]-model.knots_V[j]
  slope=(model.gain_V_per_e[j+1]-model.gain_V_per_e[j])/width
  F[j+1]=F[j]+reciprocal_gain_integral(width,model.gain_V_per_e[j],slope)
 end
 B=zeros(size(p.response));packet=ER.prepare_bucket_current(p.calibration.current_pulse,p.time.bucket_ms)
 for k in 1:p.time.bucket_count
  counts=zeros(Int64,p.time.bucket_count);counts[k]=1
  current=ER.prepared_ion_current(counts,p.time,packet,p.adc_time)
  integral=ER.primitive(current.time_s,current.current_A)
  cumulative=[integral(t)/ER.E_CHARGE for t in p.times]
  baseline=mean(view(cumulative,p.groups[1]))
  for j in 2:length(p.groups);B[j-1,k]=mean(view(cumulative,p.groups[j]))-baseline;end
 end
 # Fixed averaging and an ordinary linear solve; nonlinear noise weighting is deferred.
 inverse=qr(B)\Matrix{Float64}(I,size(B,1),size(B,1))
 p.direct_inverse=DirectIonInverse(model,linear,F,B,inverse)
end

"""Return F(x) and recovery(x)/gain(x). Negative excursions use the small-signal
linear extension, preserving signed estimates. Positive excursions must remain
inside the measured nonlinear table; no extrapolation or clipping is applied."""
function inverse_terms(d::DirectIonInverse,x)
 if x<0
  g=d.linear.gain_V_C*ER.E_CHARGE
  return x/g,x/(d.linear.tau_s*g)
 end
 model=d.model
 j=clamp(searchsortedlast(model.knots_V,x),1,length(model.knots_V)-1)
 width=model.knots_V[j+1]-model.knots_V[j];delta=x-model.knots_V[j]
 slope=(model.gain_V_per_e[j+1]-model.gain_V_per_e[j])/width
 a=delta/width;g=model.gain_V_per_e[j]+slope*delta
 recovery=(1-a)*model.recovery_V_s[j]+a*model.recovery_V_s[j+1]
 d.integral_gain[j]+reciprocal_gain_integral(delta,model.gain_V_per_e[j],slope),recovery/g
end

"""
    direct_ion_inverse(plan, voltage_drop) -> (ions, valid, residual, iterations)

Invert dx/dt=g(x)I/e-r(x) without waveform fitting:
Q(t)/e = F(x(t))-F(x(0)) + integral(r(x)/g(x),dt).
Evaluate F exactly within each linear gain-table segment and use trapezoidal
integration of recovery at the ADC sample times. Then average cumulative charge
and apply the prepared linear current-pulse inverse to obtain collected ions.

voltage_drop is measured baseline minus voltage at every ADC sample. It shares
the plan's sampling times. No random numbers or Monte Carlo charge are used.
Noise-induced nonlinear bias is not calibrated by this operation. A domain
failure returns NaN counts and valid=false. Voltage residuals are not calculated
by this direct inverse; residual=NaN and iterations=0 keep the existing caller
result shape. Forward-voltage checks belong to validation.
"""
function direct_ion_inverse(p::IonReconstructionPlan,x::AbstractVector)
 length(x)==length(p.times)||throw(DimensionMismatch("ADC voltage samples"))
 d=prepare_direct_inverse!(p)
 if !all(isfinite,x)||maximum(x)>last(d.model.knots_V)
  return (;ions=fill(NaN,Int(p.time.bucket_count)),valid=false,residual=NaN,iterations=0)
 end
 cumulative=zeros(length(x));F0,previous=inverse_terms(d,Float64(first(x)));integral=0.
 for j in 2:length(x)
  F,rate=inverse_terms(d,Float64(x[j]))
  integral+=(p.times[j]-p.times[j-1])*(previous+rate)/2
  cumulative[j]=F-F0+integral;previous=rate
 end
 baseline=mean(view(cumulative,p.groups[1]))
 averaged=[mean(view(cumulative,g))-baseline for g in p.groups[2:end]]
 q=d.inverse*averaged
 (;ions=q,valid=all(isfinite,q),residual=NaN,iterations=0)
end
