module IonCSANonlinear
using LinearAlgebra
export NonlinearCSA, nonlinear_waveform, fit_nonlinear_csa, sampled_response
const ELEMENTARY_CHARGE_C = 1.602176634e-19

"""
    NonlinearCSA(knots_V, gain_V_per_e, recovery_V_s; baseline_V, step_s=1e-4)

Own the tables defining a one-state CSA response. The state is voltage excursion
`x = baseline_V - output_V`, with evolution `dx/dt = gain(x)*I/e - recovery(x)`.
Knots increase strictly from zero. Gain is positive in V per elementary input
charge; recovery is nonnegative in V/s and vanishes at zero. Linear interpolation
defines both tables. `step_s` bounds RK4 integration intervals in seconds.

The model describes the noiseless voltage. The existing downstream calibration,
bucket reconstruction and noise covariance remain separate. Inputs are copied
into owned Float64 vectors. Invalid tables or nonfinite parameters raise
`ArgumentError`. A simulated excursion beyond the tables requests fallback.
"""
struct NonlinearCSA
    knots_V::Vector{Float64}
    gain_V_per_e::Vector{Float64}
    recovery_V_s::Vector{Float64}
    baseline_V::Float64
    step_s::Float64
    function NonlinearCSA(knots,gain,recovery;baseline_V,step_s=1e-4)
        k,g,r=Float64.(knots),Float64.(gain),Float64.(recovery)
        length(k)>=2 && length(k)==length(g)==length(r) || throw(ArgumentError("Equal tables with at least two knots required"))
        all(isfinite,k) && first(k)==0 && all(>(0),diff(k)) || throw(ArgumentError("Increasing excursion knots must start at zero"))
        all(x->isfinite(x)&&x>0,g) && all(x->isfinite(x)&&x>=0,r) && first(r)==0 || throw(ArgumentError("Positive gain and nonnegative recovery with zero baseline required"))
        isfinite(baseline_V) && isfinite(step_s) && step_s>0 || throw(ArgumentError("Finite baseline and positive integration step required"))
        new(copy(k),copy(g),copy(r),Float64(baseline_V),Float64(step_s))
    end
end

function bracket(knots,x)
    j=clamp(searchsortedlast(knots,x),1,length(knots)-1)
    j,(x-knots[j])/(knots[j+1]-knots[j])
end
function rate(m,x,current)
    # Negative RK intermediate excursions at equilibrium use the equilibrium law.
    j,a=bracket(m.knots_V,max(0.,x))
    ((1-a)*m.gain_V_per_e[j]+a*m.gain_V_per_e[j+1])*current/ELEMENTARY_CHARGE_C-
        ((1-a)*m.recovery_V_s[j]+a*m.recovery_V_s[j+1])
end

"""
    nonlinear_waveform(model, time_s, current_A; sample_time_s=time_s)

Propagate a piecewise-linear, nonnegative physical current from baseline.
Input times increase strictly and cover all requested sorted sample times.
Each RK4 step ends at a current knot, requested sample or the maximum timestep.
Returns owned `time_s`, `voltage_V`, and `voltage_integral_V_s` arrays together
with `valid`, `reason`, `maximum_excursion_V` and `steps`. The integral is of
absolute output voltage from the first current time, using the same RK stages.
An excursion outside the calibrated range returns `valid=false`; callers must
use SPICE for that whole history. No input array is mutated.
"""
function nonlinear_waveform(m::NonlinearCSA,t::AbstractVector,i::AbstractVector;
                            sample_time_s=t)
    length(t)==length(i)>=2 && all(isfinite,t) && all(>(0),diff(t)) || throw(ArgumentError("Increasing equal-length current arrays required"))
    all(x->isfinite(x)&&x>=0,i) || throw(ArgumentError("Finite nonnegative current required"))
    ts=Float64.(sample_time_s)
    !isempty(ts) && all(isfinite,ts) && all(>=(0),diff(ts)) && first(ts)>=first(t) && last(ts)<=last(t) || throw(ArgumentError("Sorted samples must lie inside current support"))
    volts=fill(NaN,length(ts));integrals=copy(volts);x=0.;area=0.;now=Float64(first(t));j=1;steps=0;maxx=0.
    for (n,target) in enumerate(ts)
        while now<target
            while j<length(t)-1 && t[j+1]<=now;j+=1;end
            stop=min(target,Float64(t[j+1]),now+m.step_s);h=stop-now
            h>0 || error("Integration failed to advance")
            slope=(i[j+1]-i[j])/(t[j+1]-t[j]);u=i[j]+slope*(now-t[j]);um=u+slope*h/2;ue=u+slope*h
            k1=rate(m,x,max(0.,u));x2=x+h*k1/2;k2=rate(m,x2,max(0.,um))
            x3=x+h*k2/2;k3=rate(m,x3,max(0.,um));x4=x+h*k3;k4=rate(m,x4,max(0.,ue))
            xn=x+h*(k1+2k2+2k3+k4)/6
            maxx=max(maxx,x2,x3,x4,xn)
            if !isfinite(xn) || maxx>last(m.knots_V) || xn< -1e-9
                return (;time_s=ts,voltage_V=volts,voltage_integral_V_s=integrals,
                    valid=false,reason="excursion_outside_calibration",maximum_excursion_V=maxx,steps)
            end
            area+=h*(m.baseline_V-(x+2x2+2x3+x4)/6)
            x=max(0.,xn);now=stop;steps+=1
        end
        volts[n]=m.baseline_V-x;integrals[n]=area
    end
    (;time_s=ts,voltage_V=volts,voltage_integral_V_s=integrals,
      valid=true,reason="within_calibration",maximum_excursion_V=maxx,steps)
end

"""
    fit_nonlinear_csa(histories, knots_V; baseline_V, step_s=1e-4,
                      ridge=1e-8, gain_prior=5e-5, tau_prior_s=20.)

Fit a continuous piecewise-linear evolution law from named tuples containing
`time_s`, `current_A`, `voltage_V`. The trapezoidal midpoint voltage increment
equation is linear in knot gains and recovery rates. Time-weighted normal
equations share neighboring knots; a weak diagonal prior regularizes unsupported
directions. Recovery at equilibrium is fixed to zero. Returns `model`, fitted
coefficient vectors and accumulated support time per knot. Nonphysical fitted
coefficients raise `ArgumentError`; scripts must inspect support and choose an
adequately covered domain. Model selection and event splits belong to scripts.
"""
function fit_nonlinear_csa(histories,knots;baseline_V,step_s=1e-4,ridge=1e-8,
                           gain_prior=5e-5,tau_prior_s=20.)
    k=Float64.(knots);n=length(k);N=2n-1
    normal=zeros(N,N);rhs=zeros(N);support=zeros(n);points=0
    for h in histories
        t,i,v=h.time_s,h.current_A,h.voltage_V
        for a in 1:length(t)-1
            dt=t[a+1]-t[a];dt>=1e-6 || continue
            x=baseline_V-(v[a]+v[a+1])/2
            0<=x<=last(k) || continue
            j,w=bracket(k,x);u=(i[a]+i[a+1])/(2ELEMENTARY_CHARGE_C*1e5)
            y=(v[a]-v[a+1])/dt
            indices=Int[j,j+1];values=Float64[(1-w)*u,w*u]
            if j>1;push!(indices,n+j-1);push!(values,-(1-w));end
            push!(indices,n+j);push!(values,-w)
            for (b,ib) in enumerate(indices)
                rhs[ib]+=dt*values[b]*y
                for (c,ic) in enumerate(indices);normal[ib,ic]+=dt*values[b]*values[c];end
            end
            support[j]+=dt*(1-w);support[j+1]+=dt*w;points+=1
        end
    end
    prior=vcat(fill(gain_prior*1e5,n),k[2:end]./tau_prior_s)
    all(>(0),support) || throw(ArgumentError("Every calibration knot requires observed support"))
    lambda=ridge*max(maximum(diag(normal)),1.)
    coeff=(Symmetric(normal)+lambda*I)\(rhs+lambda*prior)
    gain=coeff[1:n]./1e5;recovery=vcat(0.,coeff[n+1:end])
    model=NonlinearCSA(k,gain,recovery;baseline_V,step_s)
    (;model,support_time_s=support,points,gain_V_per_e=gain,recovery_V_s=recovery)
end

"""
    sampled_response(model, current_times_s, current_A, edges_s;
                     averaging_s, alpha, gain_V_per_ion)

Produce voltage samples and the existing discharge-corrected ion estimator at
acquisition edges. Windows centered before current acquisition use baseline;
windows intersecting acquisition integrate baseline plus the propagated voltage.
The current domain must extend through the last averaging window. Returns edge
voltages, average voltages, filtered estimates and the propagation diagnostics.
The downstream three-frame inverse and noise remain separate.
"""
function sampled_response(m,t,i,edges;averaging_s,alpha,gain_V_per_ion)
    averaging_s>0 && 0<alpha<=1 && gain_V_per_ion>0 || throw(ArgumentError("Invalid estimator calibration"))
    left=edges.-averaging_s/2;right=edges.+averaging_s/2
    points=sort!(unique(max.(first(t),vcat(edges,left,right))))
    wave=nonlinear_waveform(m,t,i;sample_time_s=points)
    !wave.valid && return (;valid=false,wave)
    voltage(x)=x<first(t) ? m.baseline_V : wave.voltage_V[searchsortedfirst(points,x)]
    primitive(x)=x<first(t) ? (x-first(t))*m.baseline_V : wave.voltage_integral_V_s[searchsortedfirst(points,x)]
    averages=[(primitive(b)-primitive(a))/averaging_s for (a,b) in zip(left,right)]
    delta=averages.-m.baseline_V
    estimated=-(delta[2:end].-alpha.*delta[1:end-1])./gain_V_per_ion
    (;valid=true,voltage_V=voltage.(edges),average_V=averages,estimated_ions=estimated,wave)
end
end
