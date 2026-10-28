@testset "Docstring examples" begin
    for object in (CoordinateFrame, EventInfo, EventHits, center_hits)
        binding = Base.Docs.Binding(EventAnalysis, nameof(object))
        signature = object === center_hits ? Tuple{EventHits,EventInfo} : Union{}
        documentation = Base.Docs.docstr(binding, signature)
        @test documentation isa Base.Docs.DocStr
        source = join(documentation.text)
        examples = collect(eachmatch(r"```julia\n(.*?)```"s, source))
        @test !isempty(examples)
        for example in examples
            sandbox = Module(gensym(:DocExample))
            Core.eval(sandbox, :(using Main.EventAnalysis))
            @test Core.eval(sandbox, Meta.parseall(example.captures[1])) === nothing
        end
    end
end

@testset "Adapted function examples" begin
    for file in (STEP4 ? ("reconstruction.jl", "paths.jl", "production.jl", "diffusion.jl", "collection.jl", "ensemble.jl") : STEP3 ? ("reconstruction.jl", "paths.jl", "production.jl", "diffusion.jl", "collection.jl") : STEP2 ? ("reconstruction.jl", "paths.jl", "production.jl", "diffusion.jl") : ("reconstruction.jl", "paths.jl"))
        source = read(joinpath(ROOT,"src/event_analysis",file),String)
        for example in eachmatch(r"```julia\n(.*?)```"s, source)
            sandbox = Module(gensym(:AdaptedDocExample))
            Core.eval(sandbox, :(using Main.EventAnalysis))
            @test Core.eval(sandbox, Meta.parseall(example.captures[1])) === nothing
        end
    end
end
