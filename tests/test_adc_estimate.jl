using Test,Statistics,HDF5
include(joinpath(@__DIR__,"../src/ion_signal/adc_estimate.jl"))
@testset "Quantization noise as an independent two-reading distribution" begin
 n=128;errors=((collect(1:n).-0.5)./n.-0.5)
 for alpha in (0.,.99185,1.)
  distribution=[x-alpha*y for x in errors for y in errors]
  numerical=sqrt(sum(abs2,distribution)/length(distribution))
  analytic=ADCEstimate.quantized_difference_rms(256.,8,alpha)
  @test isapprox(numerical,analytic;rtol=1e-4)
 end
 @test isapprox(ADCEstimate.quantized_difference_rms(2.,12,1.)/ADCEstimate.quantized_difference_rms(2.,12,0.),sqrt(2);rtol=1e-14)
 root=normpath(joinpath(@__DIR__,".."))
 h5open(joinpath(root,"results/ltspice/adc_estimate/analysis.h5"),"r") do h
  @test read(h["populations/ions_per_frame"])==[100.,600.,4000.]
  @test read(attributes(h["span_2V"])["integer_bits_for_ADC_noise_budget"])==12
  @test read(attributes(h["span_1V"])["integer_bits_for_ADC_noise_budget"])==11
 end
end
