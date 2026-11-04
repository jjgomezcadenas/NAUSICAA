using Test,HDF5,Random,Statistics,LinearAlgebra
root=normpath(joinpath(@__DIR__,".."))
include(joinpath(root,"src/ion_signal/ion_csa_response.jl"))
const C=IonCSAResponse
@testset "Ion response contracts and independent numerical references" begin
 E=C.EA
 cells=E.SparseCells(Int64[0,0,1],Int64[0,0,-1],Int64[-2,0,-1],Int64[5,9,2])
 grid=E.GridSpec(3f0,3f0,(0f0,0f0,0f0),(1f0,2f0),Int8(1))
 history=E.IonPadHistory(E.EventPlacement(Int64(81),50f0,100f0),grid,cells,Int64(16))
 before=copy(cells.count);d=C.dense_history(history,.03,Dict("quiet_before_frames"=>2,"recovery_tail_s"=>.06))
 @test d.firstbucket==-4
 @test d.counts[3,1]==5 && d.counts[5,1]==9 && d.counts[4,2]==2
 @test d.counts[4,1]==0
 @test sum(d.counts)==16
 @test cells.count==before
 @test d.edges[1]≈-.12 && all(diff(d.edges).≈.03)
 @test C.S.canonical(Dict("a"=>1,"b"=>2))==C.S.canonical(Dict("b"=>2,"a"=>1))
 seed=C.S.event_seed("test",Int64(81),"pad:0,0")
 @test seed==C.S.event_seed("test",Int64(81),"pad:0,0")
 @test seed!=C.S.event_seed("test",Int64(82),"pad:0,0")
 cv=0.3.^collect(0:63)
 a=C.noise_draw(cv,32,seed);b=C.noise_draw(cv,32,seed)
 @test a==b
 @test any(a.<0)
 L=cholesky(Symmetric([cv[abs(i-j)+1] for i in 1:32,j in 1:32])).L
 draws=L*randn(Xoshiro(seed),32,20000)
 @test abs(mean(abs2,draws)-1)<.02
 @test abs(mean(draws[1:end-1,:].*draws[2:end,:])-.3)<.02
 # Independently integrate the moving average of exponential-covariance noise.
 sigma=.001;tau=.2;width=.005;frame=.03;alpha=exp(-frame/20.)
 f=10.0.^range(-6,5;length=4201);asd=sqrt.(4sigma^2*tau./(1 .+(2pi.*f.*tau).^2))
 actual=C.covariance(f,asd,width,frame,alpha,4)
 expected=setprecision(256) do
  q=BigFloat(width)/BigFloat(tau);dt=BigFloat(frame)/BigFloat(tau);ss=BigFloat(sigma)^2;aa=BigFloat(alpha)
  av(l)=l==0 ? 2ss*(q+expm1(-q))/q^2 : ss*exp(-abs(l)*dt)*(2sinh(q/2)/q)^2
  [Float64((1+aa^2)*av(l)-aa*(av(l-1)+av(l+1))) for l in 0:4]
 end
 @test maximum(abs.(actual.-expected))/expected[1]<1e-4
 mktempdir(joinpath(root,"notes/physics/work")) do dir
  p=joinpath(dir,"scalar.h5")
  h5open(p,"w") do h;h["complete"]=false;write(h["complete"],true);end
  @test h5open(h->Bool(read(h["complete"])),p,"r")
 end
end
@testset "Complete packet sum includes neighbouring frames" begin
 alpha=.98;avg=[.2,4.,4.1];gain=avg[2]-alpha*avg[1]
 k=(;alpha,avg,raw=avg,baseline=.8)
 counts=zeros(Int,10);counts[4]=1
 _,_,signal=C.P.response(counts,k)
 taps=[avg[1],avg[2]-alpha*avg[1],avg[3]-alpha*avg[2]]./gain
 @test sum(signal)/gain≈sum(taps) rtol=1e-12
 @test signal[4]/gain≈1 rtol=1e-12
 @test sum(taps)>1
end
