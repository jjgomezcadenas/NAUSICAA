using Test,TOML
root=normpath(joinpath(@__DIR__,".."))
include(joinpath(root,"src/ion_signal/retention_compression.jl"))
const PC=RetentionCompression.PadCompression

function context(card)
 c=TOML.parsefile(joinpath(root,card));r=TOML.parsefile(joinpath(root,c["readout_card"]))
 PC.Context(root,c,r,"","test",Dict{String,Any}(),Dict{String,Any}())
end

@testset "Study-specific baseline policy preserves the kernel" begin
 card="config/pad_history_viewer_direct.toml"
 old=context(card);new=context(card)
 # Invoke the unchanged default method directly, then compare the typed study method.
 expected=invoke(PC.kernel,Tuple{Any,Any},old,"p2_bb0nu")
 new.card["baseline_validation_max_error_V"]=1e-7
 actual=PC.kernel(new,"p2_bb0nu")
 for field in (:gain,:noise,:tau,:alpha,:baseline,:enc,:baseline_error)
  @test getproperty(actual,field)==getproperty(expected,field)
 end
 for (directory,name) in (("retention_scan","length_1600um"),("retention_scan_extension","length_4000um"))
  original="results/ltspice/$directory/cards/$(name)_viewer.toml"
  strict=context(original)
  @test_throws AssertionError invoke(PC.kernel,Tuple{Any,Any},strict,"p2_bb0nu")
  validated=context("results/ltspice/$directory/cards/$(name)_validated_worker_viewer.toml")
  for sample in ("p2_bb0nu","p4_bb0nu")
   k=PC.kernel(validated,sample)
   @test k.baseline_error<1e-7
   @test k.baseline_error/k.gain<.001
  end
 end
end
