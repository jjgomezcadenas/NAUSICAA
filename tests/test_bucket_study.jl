using Test,HDF5
include(joinpath(@__DIR__,"../src/ion_signal/bucket_study.jl"))
using .BucketStudy
const B=BucketStudy
@testset "Uniform finite-gap input and bucket estimator" begin
 a=.00025;g=.005;v=3.;q=B.E*(1-B.phi(g,a));pre=.002
 @test B.phi(0.,a)==1
 @test B.phi(g,a)>0
 # Independent Ramo integration of the derivative, including finite-path loss.
 nodes,weights=B.gauss(16)
 W(z)=4a^2/(pi*(z^2+a^2)*sqrt(z^2+2a^2))
 zedges=range(0.,g,length=501)
 charge=0.
 for (lo,hi) in zip(zedges[1:end-1],zedges[2:end])
  mid=(lo+hi)/2;half=(hi-lo)/2
  charge+=B.E*half*sum(w*W(mid+half*n) for (n,w) in zip(nodes,weights))
 end
 @test charge≈q rtol=1e-11
 for frame in (.01,.02)
  @test B.current(-pre,frame,a,g,v)==0
  @test B.current(frame+.001,frame,a,g,v)==0
  @test B.current(frame/2,frame,a,g,v)≈q/frame rtol=1e-12
  t,i=B.pwl(frame,a,g,v,pre,3frame,5e-6)
  @test B.integral(t,i)≈q rtol=1e-10
  tt,ii=B.pwl(frame,a,g,v,pre,3frame,2.5e-6;train=true)
  @test B.integral(tt,ii)≈2q rtol=1e-10
  @test minimum(i)>=0
  # Exact exponential baseline cancels; constant-current response follows CSA ODE.
  tau=.3;alpha=exp(-frame/tau);T=.0003
  ts=collect(range(-pre,3frame,length=20001))
  vs=exp.(-ts./tau)
  @test abs(B.AveragingScan.waveform_gain(ts,vs,T,frame,alpha))<1e-10
  # Band-limited white noise agrees with independent time-domain overlap result.
  f=10. .^ range(-3.,6.,length=1801);asd=ones(length(f))
  sigma=B.AveragingScan.noise_rms(f,asd,T,frame,alpha)
  @test sigma^2≈(1+alpha^2)/(2T) rtol=.001
 end
 @test B.current(.005,.01,a,g,v)≈2B.current(.01,.02,a,g,v) rtol=1e-12
end
@testset "Historical input and completed bucket results" begin
 root=normpath(joinpath(@__DIR__,".."))
 lines=readlines(joinpath(root,"results/ltspice/decks/uniform_r2.pwl"))
 old=reduce(hcat,[[parse(Float64,x) for x in split(line)] for line in lines])
 t,i=B.pwl(.01,.00025,.005,3.,.002,.03,2.5e-6)
 @test length(t)==size(old,2)
 @test maximum(abs.(t.-old[1,:]))<1e-14
 @test maximum(abs.(i.-old[2,:]))/maximum(i)<1e-10
 for ms in (10,20)
  path=joinpath(root,"results/ltspice/bucket_scan/$(ms)ms/analysis.h5")
  if isfile(path)
   h5open(path,"r") do h
    at=attributes(h);n=read(at["verified_N_snr5"])
    counts=read(h["counts"]);snr=read(h["snr"])
    @test snr[findfirst(==(n-1),counts)]<5<snr[findfirst(==(n),counts)]
    @test maximum(read(h["finite_count_refinement_error"]))<1e-3
    @test read(at["charge_error"])<2e-5
    @test abs(sum(read(h["device_variance_fraction"]))-1)<1e-3
   end
  end
 end
end
