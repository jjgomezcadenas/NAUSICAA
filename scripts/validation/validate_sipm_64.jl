using Test
include("../production/run_framed_reconstruction.jl")
const P64=ER.SiPMResponse
"""Validate 64 substeps against individual-microcell simulation at the active
operating point. Both start fully recovered. Use common sampled optical counts
and ASIC noise to compare charge moments and clipping with reduced sampling error.
Dark avalanches and cell occupancy are independently sampled in both methods."""
function validate_64()
 BLAS.set_num_threads(1)
 ec=TOML.parsefile(joinpath(ROOT,"config/readout/bb0nu_electronics_pilot_100.toml"))
 s=prepare_electronics_run(ec).sipm;m=s.illuminated_microcells;total=s.sipm.microcells
 # PDE is sampled once outside each sensor call so both see the same optical counts.
 lit=P64.Sensor(1.,m,s.sipm.gain,s.sipm.recovery_s,s.sipm.dark_rate_Hz*m/total)
 unlit=P64.Sensor(0.,total-m,s.sipm.gain,s.sipm.recovery_s,s.sipm.dark_rate_Hz*(total-m)/total)
 options=P64.StatisticalOptions(64,16;initial_state=:fully_recovered)
 dt=s.asic.bucket_width_s;q0=P64.charge_per_pe_C(s.sipm)
 limit=s.asic.current_limit_A*dt/q0
 patterns=[("dark",zeros(Int64,8)),("weak",fill(Int64(100),8)),("steady",fill(Int64(3000),8)),
  ("clipping_edge",Int64[12000,13000,13500,13700,14000,14500,18000,22000]),
  ("isolated",Int64[0,0,0,22000,0,0,0,0]),
  ("rising_falling",Int64[0,1000,6000,14000,22000,14000,6000,0])]
 n=1000;rows=NamedTuple[];started=time();counter=0
 for (name,photons) in patterns
  counter+=1
  exact=zeros(8,n);approx=zeros(8,n);clip_exact=falses(8,n);clip_approx=falses(8,n)
  optical=Xoshiro(98231+counter);micro=Xoshiro(12004+counter);stat=Xoshiro(9342+counter);asic=Xoshiro(8003+counter)
  for repeat in 1:n
   detected=[rand(optical,P64.Binomial(v,s.sipm.pde)) for v in photons]
   a=P64.simulate_channel(micro,detected,dt,lit;initial_state=:fully_recovered)
   b=P64.simulate_channel(micro,zeros(Int64,8),dt,unlit;initial_state=:fully_recovered)
   c=P64.simulate_channel_statistical(stat,detected,dt,lit;options)
   d=P64.simulate_channel_statistical(stat,zeros(Int64,8),dt,unlit;options)
   exact[:,repeat]=(a.charge_C+b.charge_C)./q0
   approx[:,repeat]=(c.charge_C+d.charge_C)./q0
   noise=s.asic.noise_rms_pe.*randn(asic,8)
   clip_exact[:,repeat]=(exact[:,repeat]+noise).>limit
   clip_approx[:,repeat]=(approx[:,repeat]+noise).>limit
  end
  for j in 1:8
   x=exact[j,:];y=approx[j,:];mx=mean(x);my=mean(y);vx=var(x;corrected=false);vy=var(y;corrected=false)
   mean_se=std(y-x)/sqrt(n)
   variance_se=std((y.-my).^2 .-(x.-mx).^2)/sqrt(n)
   cx=Float64.(clip_exact[j,:]);cy=Float64.(clip_approx[j,:]);clip_se=std(cy-cx)/sqrt(n)
   mtol=max(.05,.005*mx,4mean_se);vtol=max(.01,.05*vx,4variance_se);ctol=max(.005,4clip_se)
   push!(rows,(pattern=name,bucket=j,photons=photons[j],reference_mean_pe=mx,mean_pe=my,mean_difference_pe=my-mx,mean_se_pe=mean_se,
    mean_tolerance_pe=mtol,reference_variance_pe2=vx,variance_pe2=vy,variance_difference_pe2=vy-vx,variance_se_pe2=variance_se,variance_tolerance_pe2=vtol,
    reference_clipping=mean(cx),clipping=mean(cy),clipping_difference=mean(cy-cx),clipping_se=clip_se,clipping_tolerance=ctol,
    mean_pass=abs(my-mx)<=mtol,variance_pass=abs(vy-vx)<=vtol,clipping_pass=abs(mean(cy-cx))<=ctol))
  end
  println("Completed $name ($n repetitions)");flush(stdout)
 end
 dc=TOML.parsefile(joinpath(ROOT,ec["detector_card"]));out=joinpath(ROOT,"results/readout",dc["run_label"],"electronics_sipm_substeps64")
 mkpath(out);target=joinpath(out,"accuracy.h5");isfile(target)&&error("Preserve existing validation")
 tmp=target*".tmp";budget=4*1024^2;GuardedFiles.space_check(out,budget,1024^3)
 try
  h5open(tmp,"w") do f
   f["schema"]="itaca.sipm64_accuracy.v1";f["repetitions"]=n;f["substeps"]=64;f["initial_state"]="fully_recovered"
   f["total_s"]=time()-started;columns(create_group(f,"cases"),rows)
   f["summary/mean_failures"]=count(r->!r.mean_pass,rows);f["summary/variance_failures"]=count(r->!r.variance_pass,rows);f["summary/clipping_failures"]=count(r->!r.clipping_pass,rows)
   bright=filter(r->r.photons>=3000,rows)
   f["summary/maximum_bright_mean_difference_percent"]=100maximum(abs(r.mean_difference_pe)/r.reference_mean_pe for r in bright)
   f["summary/maximum_bright_variance_difference_percent"]=100maximum(abs(r.variance_difference_pe2)/r.reference_variance_pe2 for r in bright)
   f["summary/maximum_clipping_difference_pp"]=100maximum(abs(r.clipping_difference) for r in rows)
   f["sipm_card"]=read(joinpath(ROOT,ec["sipm_card"]),String);f["asic_card"]=read(joinpath(ROOT,ec["asic_card"]),String)
   f["code_files"]=[relpath(@__FILE__,ROOT),"src/ion_signal/sipm_response.jl","src/ion_signal/sipm_statistical.jl"]
   f["code_sha256"]=[bytes2hex(open(sha256,p)) for p in [@__FILE__,joinpath(ROOT,"src/ion_signal/sipm_response.jl"),joinpath(ROOT,"src/ion_signal/sipm_statistical.jl")]]
   flush(f);GuardedFiles.check_size(tmp,budget)
  end
  mv(tmp,target)
 catch;isfile(tmp)&&rm(tmp);rethrow();end
 @testset "64 substeps against individual microcells" begin
  for r in rows
   @test r.mean_pass
   @test r.variance_pass
   @test r.clipping_pass
  end
 end
 println(target)
end
validate_64()
