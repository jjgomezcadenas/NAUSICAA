using Test, DataFrames, HDF5
include(joinpath(@__DIR__, "../src/ion_signal/event_analysis.jl"))
using .EventAnalysis

function sample_frame()
    # Gapped IDs, different event energies, repeated hit IDs in different particles.
    DataFrame(event_id=Int64[0,0,7,7],particle_id=Int32[1,2,1,2],hit_id=Int32[0,0,0,0],
              x_mm=Float32[10,14,-30,-26],y_mm=Float32[4,4,20,20],
              z_mm=Float32[-5,-1,100,104],energy_MeV=Float32[1,3,2,2])
end
sample_vertices() = DataFrame(event_id=Int64[0,7],vertex_x_mm=Float32[10,-30],
                              vertex_y_mm=Float32[4,20],vertex_z_mm=Float32[-5,100])

@testset "Ion expectations and generation vertices" begin
    df = sample_frame()
    v = sample_vertices()
    s = event_statistics(df,v; w_eV=20)
    @test s.energy == [4,4]
    @test s.expected_ions == [200000,200000]
    @test s.ions_per_hit == [100000,100000]
    @test s.cx == [10,-30]
    @test s.cy == [4,20]
    @test s.cz == [-5,100]
    @test s.extent == [4,0,4]
    @test event_statistics(df,v; w_eV=40).expected_ions == s.expected_ions/2
    translated = copy(df)
    translated.x_mm .+= 100; translated.y_mm .-= 30; translated.z_mm .+= 50
    shifted_vertices = copy(v)
    shifted_vertices.vertex_x_mm .+= 100; shifted_vertices.vertex_y_mm .-=30; shifted_vertices.vertex_z_mm .+=50
    t = event_statistics(translated,shifted_vertices; w_eV=20)
    @test t.extent == s.extent
    @test t.cx == s.cx .+ 100
    @test_throws ArgumentError event_statistics(df,v; w_eV=0)
    @test_throws ArgumentError event_statistics(df[1:0,:],v)
    bad = copy(df); bad.energy_MeV[1:2] .= 0
    @test_throws ErrorException event_statistics(bad,v)
    bad = copy(df); bad.energy_MeV[1] = -1
    @test_throws ErrorException event_statistics(bad,v)
    bad = copy(df); bad.x_mm[1] = NaN32
    @test_throws ErrorException event_statistics(bad,v)
    @test_throws ErrorException event_statistics(df[[1,3,2,4],:],v)
    @test_throws ErrorException event_statistics(df,v[1:1,:])
    distant = copy(df); distant.x_mm[2] += 1000
    @test event_statistics(distant,v).cx == s.cx
    @test event_statistics(df,v[[2,1],:]).cx == s.cx
end

@testset "Edges, orientation, and unconditional ensemble means" begin
    edges = [-4.,-2.,0.,2.,4.]
    @test bin_index(-4,edges) == 1
    @test bin_index(-2,edges) == 2
    @test bin_index(0,edges) == 3
    @test bin_index(4,edges) == 4
    @test_throws DomainError bin_index(nextfloat(4.),edges)
    @test_throws DomainError bin_index(NaN,edges)
    grid = (x_edges_mm=edges,y_edges_mm=[-2.,0.,2.],z_edges_mm=edges)
    event = (x_mm=[-4.,4.],y_mm=[1.,-1.],z_mm=[0.,4.],ions=[0.25,0.75])
    hist = project_event(event,grid)
    @test size(hist.xy) == (4,2)
    @test hist.xy[1,2] == .25
    @test hist.xy[4,1] == .75
    @test hist.z == [0,0,.25,.75]
    @test sum(hist.xy) == sum(hist.z) == 1
    df = sample_frame(); df.energy_MeV[3:4] .= 1 # Unequal energies: no per-event normalization.
    v = sample_vertices()
    s = event_statistics(df,v; w_eV=20)
    grid = common_grid([s]; pitch_mm=2,z_bin_mm=2)
    @test all(diff(grid.x_edges_mm) .== 2)
    @test last(grid.x_edges_mm) > maximum(s.extent)
    r = analyze_sample(df,s,grid; w_eV=20,example_event_id=7)
    @test sum(r.xy) ≈ 150000
    @test sum(r.z) ≈ 150000
    @test sum(r.per_hit_ions) ≈ 300000
    @test sum(r.example_hist.xy) ≈ 100000
    # Compare to a simple dense, per-event implementation independent of the reducer.
    reference = [project_event(EventAnalysis.centered_event(df,s,j,20),grid) for j in 1:2]
    @test r.xy ≈ (reference[1].xy+reference[2].xy)/2
    @test r.z ≈ (reference[1].z+reference[2].z)/2
    moved = copy(df); moved.x_mm .+= 1000; moved.z_mm .-= 500
    vm = copy(v); vm.vertex_x_mm .+=1000; vm.vertex_z_mm .-=500
    sm = event_statistics(moved,vm; w_eV=20)
    rm = analyze_sample(moved,sm,grid; w_eV=20,example_event_id=7)
    @test rm.xy == r.xy
    @test rm.z == r.z
    @test_throws KeyError analyze_sample(df,s,grid; example_event_id=5)
end

@testset "HDF5 analysis serialization" begin
    mktempdir(joinpath(@__DIR__,"../work"); prefix="analysis-test-") do dir
        df = sample_frame(); s = event_statistics(df,sample_vertices(); w_eV=20)
        # Non-square matrix detects accidental x/y transposition across languages.
        grid = (x_edges_mm=[-4.,-2.,0.,2.,4.],y_edges_mm=[-2.,0.,2.],z_edges_mm=[-4.,-2.,0.,2.,4.])
        r = analyze_sample(df,s,grid; w_eV=20)
        input = joinpath(dir,"input.h5")
        h5open(input,"w") do f
            f["metadata/source_file"] = "synthetic.next.h5"
            f["metadata/source_sha256"] = "synthetic"
        end
        params = (w_eV=20.,pitch_mm=2.,frame_ms=20.,drift_speed_mm_ms=.1,z_bin_mm=2.,
                  margin_fraction=.1,extra_bins=1,tolerance=1e-11,w_status="Synthetic test",reference_key="Serra:2014zda")
        output = joinpath(dir,"analysis.h5")
        root = normpath(joinpath(@__DIR__,".."))
        write_analysis(output,input,Dict("id"=>"test","label"=>"Test"),s,grid,r,params,
                       joinpath(root,"config/ion_tracks.toml"),joinpath(root,"config/detector.toml"))
        h5open(output,"r") do f
            @test read(f["mean/xy"]) == r.xy
            @test read(f["mean/z"]) == r.z
            @test read(attributes(f["mean/xy"])["axis_order_on_disk"]) == "y,x"
            @test read(f["bins/relative_time_edges_ms"]) == grid.z_edges_mm/.1
            @test read(f["hits/expected_ions"]) == [50000,150000,100000,100000]
            @test read(f["events/event_id"]) == [0,7]
            @test read(f["summary/energy_MeV"]) == [4,0,4,4]
            @test read(f["summary/pooled_mean_ions_per_hit"]) == 100000
        end
    end
end

@testset "Vertex-relative cuts preserve raw totals and account for discarded ions" begin
    df = sample_frame(); v = sample_vertices()
    df.x_mm[3:4] .+= 100
    s = event_statistics(df,v; w_eV=20)
    selected = select_hits(df,s; xy_halfwidth_mm=5)
    @test selected.selected == [true,true,false,false]
    @test s.energy == [4,4]
    @test s.cx == [10,-30]
    grid = common_grid([selected])
    raw_grid = common_grid([s])
    r = analyze_sample(df,s,grid; w_eV=20,xy_halfwidth_mm=5,diagnostic_grid=raw_grid)
    @test r.selected_energy == [4,0]
    @test r.selected_hits == [2,0]
    @test r.selected_ions == [200000,0]
    @test sum(r.xy) ≈ 100000 # Both events remain in the ensemble denominator.
    @test sum(r.z) ≈ 100000
    @test sum(r.diagnostic_x_ions) ≈ 200000
    @test sum(r.diagnostic_y_hits) ≈ 2
    @test sum(r.per_hit_ions[r.selected])+sum(r.per_hit_ions[.!r.selected]) ≈ sum(s.expected_ions)
    @test maximum(r.errors) <= 1e-11
end
