"""Convert endpoint-coded ADC samples to voltage at the centre of each interval."""
function adc_voltages(r::IonADCResponse,s::ADCSettings)
 ER.check_adc(s)
 s.sample_us==r.time.sample_us||throw(ArgumentError("ADC sampling mismatch"))
 all(<=(2^s.bits-1),r.counts)||throw(ArgumentError("ADC code exceeds bit range"))
 s.minimum_V .+(Float64.(r.counts).+.5).*((Float64(s.maximum_V)-s.minimum_V)/2^s.bits)
end

"""Exact covariance of rectangular averages of sampled stationary voltage noise.

Add quantisation variance to zero lag. The first interval is the baseline;
return covariance of baseline minus each subsequent average. Double cumulative
sums evaluate interval pairs without constructing a sample-by-sample matrix.
"""
function averaged_covariance(covariance,groups,quantisation_variance)
 n=maximum(last,groups);length(covariance)>=n||throw(ArgumentError("Covariance length"))
 c=copy(covariance[1:n]);c[1]+=quantisation_variance
 D=zeros(n+1);D[2]=c[1]/2
 slope=D[2]
 for k in 2:n;slope+=c[k];D[k+1]=D[k]+slope;end
 d(k)=D[abs(k)+1]
 m=length(groups);C=zeros(m,m)
 for i in 1:m,j in 1:i
  a,b=first(groups[i]),last(groups[i]);u,v=first(groups[j]),last(groups[j])
  value=(d(b-u+1)-d(a-u)-d(b-v)+d(a-v-1))/(length(groups[i])*length(groups[j]))
  C[i,j]=C[j,i]=value
 end
 C[2:end,2:end].-C[2:end,1].-transpose(C[1,2:end]).+C[1,1]
end

"""Prepared voltage averages, linear inverse and covariance for one event's timing.
The direct current inverse is prepared when first needed. The nonlinear
sensitivity grid is used only by the reference fit in validation."""
mutable struct IonReconstructionPlan
 time::IonReadoutTime
 adc_time::IonADCTime
 calibration::IonReconstructionCalibration
 settings::IonReconstructionSettings
 groups::Vector{UnitRange{Int}}
 times::Vector{Float64}
 response::Matrix{Float64}
 covariance::Matrix{Float64}
 lower::Matrix{Float64}
 inverse::Matrix{Float64}
 ion_covariance::Matrix{Float64}
 nonlinear_grid::Any
 direct_inverse::Any
end

"""Prepare voltage averages and the unit-ion response for the recorded bucket times.
The first group contains baseline samples before any induced current. Subsequent
nonoverlapping averages include the acquisition tail and the last short group."""
function prepare_ion_reconstruction(at,time,c,s)
 c.adc.sample_us==at.sample_us||throw(ArgumentError("ADC sample mismatch"))
 s.baseline_ms>0&&s.average_ms>0&&s.max_iterations>0&&s.convergence_ions>0||throw(ArgumentError("Reconstruction settings"))
 nb=round(Int,Float64(s.baseline_ms)*1000/at.sample_us)
 width=round(Int,Float64(s.average_ms)*1000/at.sample_us)
 1<=nb<at.sample_count && width>=1||throw(ArgumentError("Averaging intervals"))
 ts=collect(0:at.sample_count-1).*(Float64(at.sample_us)*1e-6)
 onset=(Float64(time.start_ms)-at.start_ms)*.001+first(c.current_pulse.time_s)
 ts[nb]<onset+1e-10||throw(ArgumentError("Baseline contains induced current"))
 groups=vcat([1:nb],[i:min(i+width-1,at.sample_count) for i in nb+1:width:at.sample_count])
 length(groups)-1>=time.bucket_count||throw(ArgumentError("Insufficient voltage averages"))
 A=zeros(length(groups)-1,time.bucket_count)
 packet=ER.prepare_bucket_current(c.current_pulse,time.bucket_ms)
 for k in 1:time.bucket_count
  counts=zeros(Int64,time.bucket_count);counts[k]=1
  current=ER.prepared_ion_current(counts,time,packet,at)
  # Include the final ADC time despite last-bit arithmetic differences.
  current=IonCurrentPulse(vcat(current.time_s, max(last(ts),last(current.time_s))+1e-12),vcat(current.current_A,0.))
  v=linear_csa_voltage(current,ts,c.linear_response)
  for j in 2:length(groups);A[j-1,k]=mean(view(v,groups[1]))-mean(view(v,groups[j]));end
 end
 lsb=(Float64(c.adc.maximum_V)-c.adc.minimum_V)/2^c.adc.bits
 C=averaged_covariance(c.voltage_noise.covariance_V2,groups,lsb^2/12)
 L=Matrix(cholesky(Symmetric(C)).L)
 W=L\A;H=qr(W)\(L\Matrix{Float64}(I,size(L)...))
 IonReconstructionPlan(time,at,c,s,groups,ts,A,C,L,H,Matrix(Symmetric(H*C*H')),nothing,nothing)
end

"""Build piecewise-linear unit currents at current knots and ADC times.
Every integration interval is at most the calibrated nonlinear step duration.
Only sampled voltage averages and their derivatives are retained during a fit."""
function nonlinear_grid!(p)
 p.nonlinear_grid!==nothing&&return p.nonlinear_grid
 c=p.calibration;t=p.time;at=p.adc_time
 packet=ER.prepare_bucket_current(c.current_pulse,t.bucket_ms)
 offset=(Float64(t.start_ms)-at.start_ms)*.001;dt=Float64(t.bucket_ms)*.001
 extra=collect(0:c.nonlinear_response.step_s:last(p.times))
 knots=[offset+(j-1)*dt+x for j in 1:t.bucket_count for x in packet.time_s]
 grid=sort!(unique(vcat(p.times,extra,filter(x->0<=x<=last(p.times),knots))))
 U=zeros(t.bucket_count,length(grid))
 for k in eachindex(grid),j in 1:t.bucket_count
  x=grid[k]-offset-(j-1)*dt
  first(packet.time_s)<=x<=last(packet.time_s)&&(U[j,k]=ER.interpolate(packet.time_s,packet.current_A,x))
 end
 output=zeros(Int,length(grid));weights=zeros(length(grid))
 for g in 2:length(p.groups),sample in p.groups[g]
  k=searchsortedfirst(grid,p.times[sample]);output[k]=g-1;weights[k]=1/length(p.groups[g])
 end
 p.nonlinear_grid=(;grid,U,output,weights)
end

"""Nonlinear rate and its voltage/current derivatives. Positive voltage drops use
the calibrated tables. Negative noise-induced drops use the calibrated linear
small-signal law, allowing signed reconstructed ion estimates."""
function rate_derivatives(x,current,model,linear)
 if x<0
  return linear.gain_V_C*current-x/linear.tau_s,-1/linear.tau_s,linear.gain_V_C
 end
 j=clamp(searchsortedlast(model.knots_V,x),1,length(model.knots_V)-1)
 dx=model.knots_V[j+1]-model.knots_V[j];a=(x-model.knots_V[j])/dx
 gain=(1-a)*model.gain_V_per_e[j]+a*model.gain_V_per_e[j+1]
 recovery=(1-a)*model.recovery_V_s[j]+a*model.recovery_V_s[j+1]
 derivative=(model.gain_V_per_e[j+1]-model.gain_V_per_e[j])/dx*current/ER.E_CHARGE-(model.recovery_V_s[j+1]-model.recovery_V_s[j])/dx
 gain*current/ER.E_CHARGE-recovery,derivative,gain/ER.E_CHARGE
end

"""Predict averaged voltage drops and their ion-count derivatives using RK4.
The physical nonlinear tables and timestep are the same as the CSA simulator;
the derivative equations allow all bucket counts to be fitted simultaneously."""
function nonlinear_prediction(p,q)
 data=nonlinear_grid!(p);t=data.grid;U=data.U;n=length(q);m=size(p.response,1)
 model=p.calibration.nonlinear_response
 linear=ER.exponential_parameters(p.calibration.linear_response)
 linear===nothing&&throw(ArgumentError("Signed inverse requires exponential small-signal calibration"))
 x=0.;J=zeros(n);v=zeros(m);D=zeros(m,n);k1=zeros(n);k2=zeros(n);k3=zeros(n);k4=zeros(n)
 current=transpose(U)*q
 for k in 1:length(t)-1
  h=t[k+1]-t[k];u=current[k];ue=current[k+1];um=(u+ue)/2
  r1,a1,b1=rate_derivatives(x,u,model,linear)
  x2=x+h*r1/2;r2,a2,b2=rate_derivatives(x2,um,model,linear)
  x3=x+h*r2/2;r3,a3,b3=rate_derivatives(x3,um,model,linear)
  x4=x+h*r3;r4,a4,b4=rate_derivatives(x4,ue,model,linear)
  xn=x+h*(r1+2r2+2r3+r4)/6
  if !isfinite(xn)||max(x,x2,x3,x4,xn)>last(model.knots_V)
   return (;valid=false,voltage=v,jacobian=D)
  end
  @inbounds for j in 1:n
   ui=U[j,k];uj=U[j,k+1];mid=(ui+uj)/2
   k1[j]=a1*J[j]+b1*ui
   k2[j]=a2*(J[j]+h*k1[j]/2)+b2*mid
   k3[j]=a3*(J[j]+h*k2[j]/2)+b3*mid
   k4[j]=a4*(J[j]+h*k3[j])+b4*uj
   J[j]+=h*(k1[j]+2k2[j]+2k3[j]+k4[j])/6
  end
  x=xn;row=data.output[k+1]
  if row>0
   w=data.weights[k+1];v[row]+=w*x
   @inbounds for j in 1:n;D[row,j]+=w*J[j];end
  end
 end
 (;valid=true,voltage=v,jacobian=D)
end

"""Fit signed bucket counts with covariance-weighted nonlinear voltage residuals.
Use a decreasing update length when a full step worsens the fit. Failed domain
or convergence checks return valid=false, never a silently accepted partial fit."""
function fit_nonlinear(p,y,q0)
 q=copy(q0);prediction=nonlinear_prediction(p,q)
 prediction.valid||return (;ions=q,valid=false,residual=NaN,iterations=0)
 for iteration in 1:p.settings.max_iterations
  residual=p.lower\(y-prediction.voltage);W=p.lower\prediction.jacobian
  delta=qr(W)\residual
  if maximum(abs,delta)<=p.settings.convergence_ions
   return (;ions=q,valid=true,residual=sqrt(mean(abs2,y-prediction.voltage)),iterations=iteration)
  end
  objective=sum(abs2,residual);accepted=false
  for step in (1.,.5,.25,.125,.0625,.03125)
   candidate=q.+step.*delta;trial=nonlinear_prediction(p,candidate)
   if trial.valid && sum(abs2,p.lower\(y-trial.voltage))<objective
    q=candidate;prediction=trial;accepted=true;break
   end
  end
  accepted||return (;ions=q,valid=false,residual=sqrt(mean(abs2,y-prediction.voltage)),iterations=iteration)
 end
 (;ions=q,valid=false,residual=sqrt(mean(abs2,y-prediction.voltage)),iterations=p.settings.max_iterations)
end

"""
    reco_ion_signal(response, time, calibration, settings) -> ReconstructedIons

Estimate each pad's baseline from pre-current ADC samples, form voltage averages
and reconstruct signed collection-bucket ions. Select nonlinear reconstruction
from the baseline-corrected measured averages. Propagate calibration-domain,
endpoint-code failures through valid flags and NaN estimates. Nonlinear pads use
the direct integrated inverse; the iterative fit is retained for validation.
No Monte Carlo charge, pad model choice or simulated baseline is read.
"""
function reco_ion_signal(r::IonADCResponse,t::IonReadoutTime,c::IonReconstructionCalibration,s::IonReconstructionSettings)
 reconstruct_ions(r,prepare_ion_reconstruction(r.time,t,c,s))
end

"""Apply an event's prepared inverse to all pads; returned arrays are owned."""
function reconstruct_ions(r::IonADCResponse,p::IonReconstructionPlan)
 size(r.counts,1)==p.adc_time.sample_count||throw(DimensionMismatch("ADC vector length"))
 nc=size(r.counts,2);nb=Int(p.time.bucket_count);m=size(p.response,1)
 base=zeros(Float32,nc);Y=zeros(m,nc)
 lsb=(Float64(p.calibration.adc.maximum_V)-p.calibration.adc.minimum_V)/2^p.calibration.adc.bits
 for j in 1:nc
  zero=mean(view(r.counts,p.groups[1],j))
  base[j]=p.calibration.adc.minimum_V+(zero+.5)*lsb
  for g in 2:length(p.groups);Y[g-1,j]=(zero-mean(view(r.counts,p.groups[g],j)))*lsb;end
 end
 Q=p.inverse*Y;nonlinear=falses(nc);converged=trues(nc);limit=falses(nc);rms=zeros(Float32,nc)
 for j in 1:nc
  limit[j]=any(v->v==0||v==2^p.calibration.adc.bits-1,view(r.counts,:,j))
  nonlinear[j]=p.calibration.linear_response.baseline_V-maximum(view(Y,:,j))<p.settings.nonlinear_below_V
  if limit[j]
   converged[j]=false;rms[j]=NaN32
  elseif nonlinear[j]
   zero=mean(view(r.counts,p.groups[1],j))
   drop=(zero .-Float64.(view(r.counts,:,j))).*lsb
   fit=direct_ion_inverse(p,drop)
   converged[j]=fit.valid;Q[:,j]=fit.ions;rms[j]=fit.residual
  else
   rms[j]=sqrt(mean(abs2,view(Y,:,j)-p.response*view(Q,:,j)))
  end
  !converged[j]&&(Q[:,j].=NaN)
 end
 valid=repeat(reshape(converged,1,:),nb,1)
 ReconstructedIons(r.event_id,r.channels,p.time,Float32.(Q),BitMatrix(valid),IonReconstructionChecks(base,nonlinear,converged,limit,rms))
end
