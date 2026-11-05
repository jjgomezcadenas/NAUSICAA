"""Compact conditional-count SiPM response.

Included inside SiPMResponse. State totals are R=sum(r_i) and S=sum(r_i^2).
Finite-cell occupancy moments condition on the sampled avalanche count. A
bounded beta draw matches charge mean/variance; recovery follows that same draw.
Dense intervals use a midpoint split with a TOML-configurable number of steps.
"""
struct StatisticalOptions
    substeps::Int
    sparse_trigger_limit::Int
    function StatisticalOptions(substeps::Integer=128,sparse_trigger_limit::Integer=16)
        substeps>0 && sparse_trigger_limit>=0 || throw(ArgumentError("Positive substeps and nonnegative sparse limit required"))
        new(substeps,sparse_trigger_limit)
    end
end

"""Recover aggregate first and second cell moments for a duration measured in tau.

R and S are sums over M cells. This affine update is exact between triggers.
"""
function recover_totals(R,S,M,duration_over_tau)
    d=exp(-duration_over_tau)
    deficit=M-R;second_deficit=M-2R+S
    nextR=M-d*deficit
    nextS=M-2d*deficit+d*d*second_deficit
    (nextR,clamp(nextS,nextR^2/M,nextR))
end

"""Conditional charge moments for k simultaneous triggers assigned uniformly to M cells.

A cell hit at least once contributes its available charge once. Its occupancy
probability is p=1-(1-1/M)^k. The covariance of two occupancy indicators is
(1-2/M)^k-(1-1/M)^(2k). Inputs R=sum(r_i), S=sum(r_i^2) describe the incoming
state. Returns mean charge, charge variance and expected squared charge sum.
This calculation preserves the finite trigger count, including k=1.
"""
function occupancy_moments(R::Real,S::Real,M::Integer,k::Integer)
    M>0 && k>=0 || throw(ArgumentError("Invalid cell/trigger count"))
    k==0 && return (mean=0.,variance=0.,squared_charge=0.,survival=1.)
    M==1 && return (mean=Float64(R),variance=0.,squared_charge=Float64(S),survival=0.)
    k==1 && return (mean=R/M,variance=max(0.,S/M-(R/M)^2),squared_charge=S/M,survival=1-1/M)
    la=k*log1p(-1/Float64(M));a=exp(la);p=-expm1(la)
    lb=M==2 ? -Inf : k*log1p(-2/Float64(M));b=exp(lb)
    cov=exp(2la)*expm1(lb-2la)
    occupied_variance=max(0.,M*a*p+M*(M-1)*cov)
    variance=max(0.,(a-b)*max(0.,S-R*R/M)+(R/M)^2*occupied_variance)
    (mean=p*R,variance,squared_charge=p*S,survival=a)
end

"""Draw a bounded amplitude with specified first two moments.

Uses a beta distribution on [lo,hi]; a zero variance is deterministic and the
maximum variance is a two-point distribution at the bounds. Returns the sample
and whether the requested variance exceeded its bounded feasibility limit.
"""
function bounded_moment_draw(rng,mu,variance,lo,hi)
    hi>=lo || throw(ArgumentError("Invalid amplitude bounds"))
    mu=clamp(Float64(mu),lo,hi);v=max(0.,Float64(variance));span=hi-lo
    maxvar=(mu-lo)*(hi-mu);limited=v>maxvar+1e-10*max(1.,maxvar)
    v=min(v,maxvar)
    (span==0 || v<=1e-14*max(1.,mu^2)) && return (mu,limited)
    p=(mu-lo)/span
    if maxvar-v<=1e-12*max(1.,maxvar)
        return (rand(rng)<p ? hi : lo,limited)
    end
    concentration=maxvar/v-1
    (lo+span*rand(rng,Beta(p*concentration,(1-p)*concentration)),limited)
end

"""Simulate a channel using conditional count moments and a compact recovery state.

Returns the microscopic simulator's count/charge arrays plus numerical
closure diagnostics. Binomial PDE and Poisson dark counts are retained.
Initialization samples the stationary dark recovery of M cells once, then
retains only R and S. Up to sparse_trigger_limit avalanches in a bucket use
individual uniform times and moment-matched single-cell charges. Denser buckets
partition their fixed counts multinomially into substeps and collect at each
substep midpoint. The simultaneous-occupancy mean and variance are exact for
that discretization. The outgoing R subtracts the sampled charge; S subtracts
its exact square for a single trigger and a charge-correlated squared-charge sum for
a batch using Q*S/R, whose expectation is the exact occupancy loss p*S. Optical/dark partitioning uses random-label first/second moments.

The beta sampling law approximates higher moments. State feasibility repairs
and variance limiting are counted explicitly; they are validation diagnostics.
The input vector is preserved and all returned arrays are owned.
"""
function simulate_channel_statistical(rng::AbstractRNG,incident::AbstractVector{<:Integer},
                                      dt::Real,sensor::Sensor;
                                      options::StatisticalOptions=StatisticalOptions())
    isfinite(dt) && dt>0 && all(>=(0),incident) || throw(ArgumentError("Require positive duration and nonnegative counts"))
    n=length(incident);optical=zeros(Int64,n);dark=zeros(Int64,n)
    oq=zeros(n);dq=zeros(n);M=sensor.microcells;tau=sensor.recovery_s;q0=charge_per_pe_C(sensor)
    R=Float64(M);S=Float64(M)
    if tau>0 && sensor.dark_rate_Hz>0
        r=rand(rng,Beta(1.,sensor.dark_rate_Hz*tau/M),M)
        R=sum(r);S=sum(abs2,r)
    end
    repairs=0;limited=0;label_limits=0;maxrepair=0.;batches=0
    for bucket in eachindex(incident)
        ns=rand(rng,Binomial(incident[bucket],sensor.pde));nd=rand(rng,Poisson(sensor.dark_rate_Hz*dt))
        optical[bucket]=ns;dark[bucket]=nd;N=ns+nd
        if tau==0
            oq[bucket]=ns*q0;dq[bucket]=nd*q0;continue
        end
        if N==0
            R,S=recover_totals(R,S,M,dt/tau);continue
        end
        if N<=options.sparse_trigger_limit
            arrivals=sort!([(rand(rng)*dt,i>ns) for i in 1:N];by=first)
            previous=0.
            for (time,is_dark) in arrivals
                R,S=recover_totals(R,S,M,(time-previous)/tau)
                moments=occupancy_moments(R,S,M,1)
                cell_mean=R/M;cell_spread=sqrt(max(0.,(M-1)*(S/M-cell_mean^2)))
                low=clamp(cell_mean-cell_spread,0.,1.)
                high=min(R,min(1.,cell_mean+cell_spread))
                charge,hitlimit=bounded_moment_draw(rng,moments.mean,moments.variance,low,high)
                limited+=hitlimit;R-=charge;S-=charge^2
                fixed=clamp(S,R^2/M,R);correction=abs(fixed-S)
                repairs+=correction>1e-9;maxrepair=max(maxrepair,correction);S=fixed
                if is_dark;dq[bucket]+=charge*q0;else;oq[bucket]+=charge*q0;end
                previous=time;batches+=1
            end
            R,S=recover_totals(R,S,M,(dt-previous)/tau)
        else
            h=dt/options.substeps;remaining_s=ns;remaining_d=nd
            for step in 1:options.substeps
                left=options.substeps-step+1
                ks=left==1 ? remaining_s : rand(rng,Binomial(remaining_s,1/left))
                kd=left==1 ? remaining_d : rand(rng,Binomial(remaining_d,1/left))
                remaining_s-=ks;remaining_d-=kd;k=ks+kd
                R,S=recover_totals(R,S,M,h/(2tau))
                if k>0
                    moments=occupancy_moments(R,S,M,k)
                    cell_mean=R/M;cell_spread=sqrt(max(0.,(M-1)*(S/M-cell_mean^2)))
                    low=clamp(cell_mean-cell_spread,0.,1.)
                    high=min(R,k*min(1.,cell_mean+cell_spread))
                    charge,hitlimit=bounded_moment_draw(rng,moments.mean,moments.variance,low,high)
                    limited+=hitlimit
                    squared=k==1 ? charge^2 : (R>0 ? charge*S/R : 0.)
                    if kd==0
                        signal=charge
                    elseif ks==0
                        signal=0.
                    else
                        fraction=ks/k
                        # At fixed sum Q, k amplitudes in [0,1] have sum of squares
                        # at most floor(Q) + frac(Q)^2.
                        max_squared=floor(charge)+(charge-floor(charge))^2
                        label_squared=clamp(squared,charge^2/k,max(charge^2/k,max_squared))
                        variance=fraction*(1-fraction)*max(0.,k*label_squared-charge^2)/(k-1)
                        signal,hitlimit=bounded_moment_draw(rng,fraction*charge,variance,max(0.,charge-kd),min(charge,Float64(ks)))
                        label_limits+=hitlimit
                    end
                    oq[bucket]+=signal*q0;dq[bucket]+=(charge-signal)*q0
                    R-=charge;S-=squared
                    fixed=clamp(S,R^2/M,R);correction=abs(fixed-S)
                    repairs+=correction>1e-9;maxrepair=max(maxrepair,correction);S=fixed;batches+=1
                end
                R,S=recover_totals(R,S,M,h/(2tau))
            end
        end
    end
    (;incident_photons=Int64.(incident),detected_pe=optical,dark_pulses=dark,
      optical_charge_C=oq,dark_charge_C=dq,charge_C=oq+dq,
      state_repairs=repairs,maximum_state_repair=maxrepair,variance_limits=limited+label_limits,charge_variance_limits=limited,
      label_variance_limits=label_limits,batches)
end

"""
    simulate_channel_illuminated(rng, incident, dt, sensor, illuminated_cells; options)

Statistical response of a uniformly illuminated subset of the physical SiPM.
`illuminated_cells` cells receive all optical photons; the remaining cells
receive dark triggers. The full-device primary dark rate is partitioned by cell
count. Each population uses the same gain and recovery law and independent RNG
substreams, preserving stationary dark initialization and recurrent avalanches.

Returns the standard detected-count and charge histories, separate optical and
dark charge, and `expected_distinct_optical_cells`: the conditional expectation
M*(1-(1-1/M)^N) for N sampled optical triggers in a bucket. This diagnostic
counts each illuminated cell once and does not describe charge-equivalent PE.
`detected_pe` counts optical triggers, including repeated firing of a cell.
"""
function simulate_channel_illuminated(rng::AbstractRNG,incident::AbstractVector{<:Integer},
    dt::Real,sensor::Sensor,illuminated_cells::Integer;
    options::StatisticalOptions=StatisticalOptions())
    M=Int(illuminated_cells);total=sensor.microcells
    0<M<=total || throw(ArgumentError("Illuminated cells must be within the device"))
    lit=Sensor(sensor.pde,M,sensor.gain,sensor.recovery_s,sensor.dark_rate_Hz*M/total)
    a=simulate_channel_statistical(Xoshiro(rand(rng,UInt64)),incident,dt,lit;options)
    distinct=[M==1 ? Float64(n>0) : -M*expm1(n*log1p(-1/Float64(M))) for n in a.detected_pe]
    M==total && return merge(a,(;expected_distinct_optical_cells=distinct))
    unlit=Sensor(0.,total-M,sensor.gain,sensor.recovery_s,sensor.dark_rate_Hz*(total-M)/total)
    b=simulate_channel_statistical(Xoshiro(rand(rng,UInt64)),zeros(Int64,length(incident)),dt,unlit;options)
    (;incident_photons=a.incident_photons,detected_pe=a.detected_pe,
      dark_pulses=a.dark_pulses+b.dark_pulses,optical_charge_C=a.optical_charge_C,
      dark_charge_C=a.dark_charge_C+b.dark_charge_C,charge_C=a.charge_C+b.charge_C,
      expected_distinct_optical_cells=distinct,
      state_repairs=a.state_repairs+b.state_repairs,
      maximum_state_repair=max(a.maximum_state_repair,b.maximum_state_repair),
      variance_limits=a.variance_limits+b.variance_limits,
      charge_variance_limits=a.charge_variance_limits+b.charge_variance_limits,
      label_variance_limits=a.label_variance_limits+b.label_variance_limits,batches=a.batches+b.batches)
end
