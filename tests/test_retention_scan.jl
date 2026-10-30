using Test,HDF5
root=normpath(joinpath(@__DIR__,".."))
include(joinpath(root,"src/ion_signal/retention_scan.jl"))
const R=RetentionScan

@testset "Long-retention numerical references" begin
 for tau in (8.6,20.,30.)
  t=collect(0.:.0005:6.);v=-.003 .*exp.(-t./tau)
  fitted,error,_,_,_=R.D.fit_tail(t,v,.05,5.99,.002)
  @test fitted≈tau rtol=1e-9
  @test error<1e-10
  # Independent covariance reference for centred averages of exponential noise.
  for frame in (.02,.04)
   width=.005;alpha=exp(-frame/tau);sigma=.001
   f=10.0.^range(-7.,6.;length=5201)
   asd=sqrt.(4sigma^2*tau./(1 .+(2pi.*f.*tau).^2))
   actual=R.A.noise_rms(f,asd,width,frame,alpha)
   expected=setprecision(256) do
    q=BigFloat(width)/BigFloat(tau);rho=exp(-BigFloat(frame)/BigFloat(tau))
    variance0=2BigFloat(sigma)^2*(q+expm1(-q))/q^2
    covariance=BigFloat(sigma)^2*rho*(2sinh(q/2)/q)^2
    Float64(sqrt((1+rho^2)*variance0-2rho*covariance))
   end
   @test actual≈expected rtol=1e-4
  end
 end
 mktempdir(joinpath(root,"notes/physics/work")) do dir
  source=joinpath(dir,"source.h5");destination=joinpath(dir,"copy.h5")
  h5open(source,"w") do h
   h["case/20ms/ENC_e_rms"]=[45.];attributes(h["case"])["tau_s"]=30.
  end
  h5open(source,"r") do h
   h5open(destination,"w") do out;HDF5.copy_object(h["case"],out,"case");end
  end
  h5open(destination,"r") do h
   @test read(h["case/20ms/ENC_e_rms"])==[45.]
   @test read(attributes(h["case"])["tau_s"])==30.
  end
 end
end
