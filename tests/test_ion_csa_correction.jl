using Test,LinearAlgebra,Random,Statistics,HDF5
root=normpath(joinpath(@__DIR__,".."))
include(joinpath(root,"src/ion_signal/ion_csa_correction.jl"))
const C=IonCSACorrection
@testset "Three-frame inverse, orientation and covariance" begin
 taps=[.02335853,1.,.02028386];n=24
 cv=54.6^2 .*(-.26).^collect(0:n-1)
 op=C.operator(taps,cv,n)
 reference=zeros(n,n)
 for i in 1:n,j in 1:n
  reference[i,j]=i==j ? 1. : j==i+1 ? taps[1] : j==i-1 ? taps[3] : 0.
 end
 @test Matrix(op.M)==reference
 for k in (1,2,n÷2,n-1,n)
  truth=zeros(n);truth[k]=100
  signal=reference*truth;recovered=op.F\signal
  @test recovered≈truth atol=1e-10
  @test sum(recovered)≈100 atol=1e-10
 end
 truth=vcat([0.,0.],[200.,200.,100.,400.,10.,10.,0.,5000.,200.,0.],zeros(n-12))
 @test op.F\(reference*truth)≈truth atol=1e-10
 pulse=zeros(n);pulse[12]=100;signal=op.M*pulse
 @test signal[11]≈100taps[1]
 @test signal[12]≈100
 @test signal[13]≈100taps[3]
 explicit=inv(reference)*op.old*inv(reference)'
 @test op.covariance≈explicit rtol=1e-12
 @test op.sigma.^2≈diag(explicit) rtol=1e-12
 @test isposdef(Symmetric(op.covariance))
 @test op.M*op.covariance*op.M'≈op.old rtol=1e-12
 old=cholesky(Symmetric(op.old)).L*randn(Xoshiro(991),n,30000)
 corrected=op.F\old;standard=corrected./op.sigma
 @test abs(mean(abs2,standard)-1)<.01
 @test abs(mean(standard))<.01
 @test any(corrected.<0)
 @test_throws ArgumentError C.response_matrix([.7,1.,.7],n)
 @test_throws ArgumentError C.operator(taps,cv,30)
 mktempdir(joinpath(root,"work")) do dir
  h5open(joinpath(dir,"replace.h5"),"w") do h
   h["x"]=[1.,2.];C.replace(h,"x",[3.,4.]);@test read(h["x"])==[3.,4.]
  end
 end
end

@testset "Finite-sample uncertainty of correlated noise diagnostics" begin
 rho=.3;npads=7
 covariance=[4. 6rho;6rho 9.]
 v=C.noise_statistic_variances(covariance,npads)
 # Independent Gaussian two-variable moment identities, including covariance.
 @test v.sum≈npads*(2+2rho)
 @test v.squares≈npads*4*(1+rho^2)
 @test v.lag≈npads*(1+rho^2)
 n=5;independent=C.noise_statistic_variances(Matrix{Float64}(I,n,n),npads)
 @test independent.sum==npads*n
 @test independent.squares==2npads*n
 @test independent.lag==npads*(n-1)
 @test_throws ArgumentError C.noise_statistic_variances(ones(1,1),1)
 @test_throws ArgumentError C.noise_statistic_variances(Matrix{Float64}(I,2,2),0)
end
