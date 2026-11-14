# Compare the original sampler and prepared sampler with the same saved spectrum.
include("../production/run_framed_detector_response.jl")
include("../../src/electronics/ElectronicsResponse.jl")
using .ElectronicsResponse
include("../production/framed_electronics_settings.jl")
function reference_noise(rng,n,s)
 c=s.covariance_V2;m=length(c);eig=real.(ER.rfft(c))
 z=ER.rfft(randn(rng,m));z .*=sqrt.(max.(eig,0.))
 Float32.(ER.irfft(z,m)[1:n])
end
function benchmark_noise()
 card=TOML.parsefile(joinpath(ROOT,"config/readout/bb0nu_electronics_pilot.toml"))
 c=TOML.parsefile(joinpath(ROOT,card["detector_card"]));base=joinpath(ROOT,"results/readout",c["run_label"])
 reference=joinpath(base,"electronics_noise_reference/controls.h5")
 n,pads=h5open(reference,"r") do f
  ids=read(f["events/event_number"]);j=findfirst(==(3),ids)
  (Int64(read(f["events/ion_samples"])[j]),read(f["events/ion_channels"])[j])
 end
 s=electronics_settings(card,IonReadoutTime(0f0,20f0,Int32(25));max_samples=n).ion.csa.noise
 prepared=ER.prepare_csa_noise(n,s);rng=Xoshiro(77)
 reference_noise(rng,n,s);ER.sample_csa_noise!(rng,prepared)
 before=zeros(3);after=zeros(3);bbytes=zeros(Int64,3);abytes=zeros(Int64,3)
 for j in 1:3
  GC.gc();b=@timed for _ in 1:pads;reference_noise(rng,n,s);end
  GC.gc();a=@timed for _ in 1:pads;ER.sample_csa_noise!(rng,prepared);end
  before[j]=b.time;after[j]=a.time;bbytes[j]=b.bytes;abytes[j]=a.bytes
 end
 path=joinpath(base,"electronics/noise_benchmark.h5");tmp=path*".tmp";budget=1024^2
 GuardedFiles.space_check(dirname(path),budget,1024^3)
 try
  h5open(tmp,"w") do f
   for (key,value) in (("reference_s",before),("prepared_s",after),("reference_bytes",bbytes),("prepared_bytes",abytes),
    ("reference_median_s",median(before)),("prepared_median_s",median(after)),("speedup",median(before)/median(after)),
    ("samples",n),("pads",pads),("reference_period",length(s.covariance_V2)),("prepared_period",prepared.period),
    ("code_sha256",bytes2hex(open(sha256,joinpath(ROOT,"src/electronics/ion_signal/noise.jl")))))
    f[key]=value
   end
   flush(f);GuardedFiles.check_size(tmp,budget)
  end
  mv(tmp,path;force=true)
 catch
  isfile(tmp)&&rm(tmp);rethrow()
 end
 println("Noise: ",median(before)," -> ",median(after)," s; speedup ",median(before)/median(after))
 println("Allocated bytes per ",pads," pads: ",median(bbytes)," -> ",median(abytes))
end
benchmark_noise()
