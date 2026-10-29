module BucketNoiseCovariance
using HDF5,LinearAlgebra,SHA,Test
include("averaging_scan.jl")
"""Covariance of processed frame differences, integrated from a one-sided PSD."""
function covariance(f,asd,width,frame,alpha,maxlag;order=8)
 b=[k/sqrt(4k^2-1) for k in 1:order-1];eg=eigen(SymTridiagonal(zeros(order),b))
 nodes=eg.values;weights=2 .* eg.vectors[1,:].^2
 edges=sort!(unique(vcat(f,collect(first(f):(.5/((maxlag+1)*frame+width)):last(f)))))
 lf=log.(f);lp=2 .* log.(asd);result=zeros(maxlag+1)
 for j in 1:length(edges)-1
  half=(edges[j+1]-edges[j])/2;mid=(edges[j+1]+edges[j])/2
  for (node,weight) in zip(nodes,weights)
   x=mid+half*node;lx=log(x);k=clamp(searchsortedlast(f,x),1,length(f)-1)
   psd=exp(lp[k]+(lp[k+1]-lp[k])*(lx-lf[k])/(lf[k+1]-lf[k]))
   base=half*weight*psd*sinc(x*width)^2*((1-alpha)^2+4alpha*sinpi(x*frame)^2)
   for lag in 0:maxlag;result[lag+1]+=base*cospi(2x*frame*lag);end
  end
 end
 result
end
function analyze(root)
 output=joinpath(root,"results/ltspice/bucket_noise_covariance");mkpath(output)
 h5open(joinpath(output,"analysis.h5"),"w") do out
  attributes(out)["source_sha256"]=bytes2hex(sha256(read(@__FILE__)))
  attributes(out)["assumptions"]="Stationary DC operating-point noise; independent local pad circuits; centred 300 us voltage averages and discharge-corrected frame differences; same square-pad circuit benchmark as source studies; 0.1 Hz to 1 MHz. Gaussian noise approximation for Monte Carlo."
  @testset "Processed bucket-noise covariance" begin
   for (name,path,reference,group) in (("20ms","results/ltspice/discharge_calibration/20ms/raw.h5","results/ltspice/discharge_calibration/20ms/paper_readout.h5",""),("40ms","results/ltspice/pitch_4mm_40ms/raw.h5","results/ltspice/pitch_4mm_40ms/analysis.h5",""))
    h5open(joinpath(root,reference),"r") do ref
     a=attributes(ref);frame=read(a["frame_s"]);width=read(a["averaging_s"]);alpha=read(a["alpha"]);expected=read(a["noise_rms_V"])
     h5open(joinpath(root,path),"r") do raw
      f=read(raw["noise_400/frequency"]);asd=read(raw["noise_400/V(onoise)"])
      fc=read(raw["noise_200/frequency"]);ac=read(raw["noise_200/V(onoise)"])
      c=covariance(f,asd,width,frame,alpha,5)
      cc=covariance(fc,ac,width,frame,alpha,5)
      cq=covariance(f,asd,width,frame,alpha,5;order=12)
      spectral=maximum(abs.(cc.-c))/c[1];quadrature=maximum(abs.(cq.-c))/c[1]
      @test sqrt(c[1])≈expected rtol=1e-6
      @test spectral<1e-4
      @test quadrature<1e-6
      @test all(abs.(c).<=c[1]*(1+1e-10))
      matrix=[c[abs(i-j)+1] for i in 1:6,j in 1:6]
      @test minimum(eigvals(Symmetric(matrix)))>0
      g=create_group(out,name);at=attributes(g)
      for (k,v) in ("frame_s"=>frame,"averaging_s"=>width,"alpha"=>alpha,"noise_rms_V"=>sqrt(c[1]),"ENC_e_rms"=>read(a["ENC_e_rms"]),"spectral_error"=>spectral,"quadrature_error"=>quadrature,"source"=>path,"source_sha256"=>bytes2hex(sha256(read(joinpath(root,path)))),"reference"=>reference,"reference_sha256"=>bytes2hex(sha256(read(joinpath(root,reference)))))
       at[k]=v
      end
      g["lag_frames"]=collect(0:5);g["covariance_V2"]=c;g["correlation"]=c./c[1]
      g["covariance_e2"]=c./read(a["gain_V_per_ion"])^2
      println(name," correlations=",c./c[1])
     end
    end
   end
  end
 end
end
end
