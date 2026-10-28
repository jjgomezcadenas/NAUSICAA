using DataFrames, Graphs

@testset "MC graph inspection adapter" begin
    s=MCGraphSettings()
    @test resolved_mc_distance(s)===1.25f0
    @test resolved_mc_distance(MCGraphSettings(voxel_width_mm=1))===2.5f0
    @test resolved_mc_distance(MCGraphSettings(voxel_width_mm=1,distance_mode="fixed_mm"))===1.25f0
    for kw in ((;method="bad"),(;distance_mode="bad"),(;voxel_width_mm=0),
               (;energy_threshold_keV=-1),(;max_distance_mm=Inf),(;k=0))
        @test_throws ArgumentError MCGraphSettings(;kw...)
    end
    function compare_components(v;threshold=1f0)
        original=deepcopy(v)
        c=build_mc_components(v;method="KDT",max_distance_mm=1.25f0,energy_threshold_keV=threshold,k=10)
        tracks=make_tracks(v;method="KDT",max_distance_mm=1.25f0,energy_threshold_kev=threshold,k=10)
        @test v==original
        @test length(c)==length(tracks)
        for (i,(component,t)) in enumerate(zip(c,tracks))
            w=walk_track_from_extremes(t);raw=get_raw_path(t,w.path_indices)
            @test component.component_id==i
            @test component.nodes.x_mm==t.voxels.x
            @test component.nodes.y_mm==t.voxels.y
            @test component.nodes.z_mm==t.voxels.z
            @test component.nodes.energy_MeV==t.voxels.energy
            @test Set(zip(component.edge_src,component.edge_dst))==Set((src(e),dst(e)) for e in edges(t.graph))
            @test component.path_indices==w.path_indices
            @test component.path_s_mm==raw.s
            @test component.length_mm==w.total_length
            @test component.energy_MeV==sum(t.voxels.energy)
            @test component.confidence==w.confidence
            @test isempty(component.path_indices) || component.endpoint_indices==(first(w.path_indices),last(w.path_indices))
        end
        return c
    end
    v=DataFrame(x=Float32[0,1,100,101,300],y=zeros(Float32,5),z=zeros(Float32,5),
                energy=Float32[.001,.002,.001,.001,.0005])
    c=compare_components(v)
    @test length(c)==2 # Exact threshold retained, subthreshold isolated bin omitted.
    @test c[1].length_mm===1f0
    c[1].nodes.x_mm[1]=-50f0
    @test v.x[1]===0f0
    @test isempty(compare_components(v;threshold=100f0))
    single=only(compare_components(v[1:1,:]))
    @test single.endpoint_indices==(1,1) && single.path_s_mm==Float32[0]
    @test isempty(single.edge_src) && single.length_mm===0f0
    for file in ("0nubb_hitdf.h5","electron_0nubb_energy_hitdf.h5")
        path=joinpath(ROOT,"events",file)
        if isfile(path)
            hits,_=read_event(path,Int64(2))
            compare_components(voxelize_event(to_petit_dataframe(hits),s.voxel_width_mm))
        end
    end
    for example in eachmatch(r"```julia\n(.*?)```"s,read(joinpath(ROOT,"src/event_analysis/mc_inspection.jl"),String))
        sandbox=Module(gensym(:InspectionExample));Core.eval(sandbox,:(using Main.EventAnalysis))
        @test Core.eval(sandbox,Meta.parseall(example.captures[1]))===nothing
    end
end
