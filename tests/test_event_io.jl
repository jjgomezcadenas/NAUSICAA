using Test, HDF5, DataFrames, SHA
include(joinpath(@__DIR__, "../src/ion_signal/event_io.jl"))
using .EventIO
include(joinpath(@__DIR__, "../src/event_analysis/EventAnalysis.jl"))
using .EventAnalysis: read_event, source_frame

# Fixed-size string compounds exercise the actual NEXUS HDF5 representation.
struct FixtureHit
    event_id::Int64
    x::Float32
    y::Float32
    z::Float32
    time::Float32
    energy::Float32
    label::NTuple{32,UInt8}
    particle_id::Int32
    hit_id::Int32
end
struct FixtureConfig
    param_key::NTuple{64,UInt8}
    param_value::NTuple{64,UInt8}
end
struct FixtureParticle
    event_id::Int64
    particle_id::Int32
    primary::Int8
    initial_x::Float32
    initial_y::Float32
    initial_z::Float32
end
fixed(s, n) = ntuple(i -> i <= ncodeunits(s) ? codeunit(s, i) : UInt8(0), n)
function compound_type(::Type{T}) where T
    d = HDF5.Datatype(HDF5.API.h5t_create(HDF5.API.H5T_COMPOUND, sizeof(T)))
    for (i, name) in enumerate(fieldnames(T))
        F = fieldtype(T, i)
        t = if F <: Tuple
            st = HDF5.Datatype(HDF5.API.h5t_copy(HDF5.API.H5T_C_S1))
            HDF5.API.h5t_set_size(st, sizeof(F))
            st
        else
            datatype(F)
        end
        HDF5.API.h5t_insert(d, string(name), fieldoffset(T, i), t)
        close(t)
    end
    return d
end
function write_compound(parent, name, rows::Vector{T}) where T
    t = compound_type(T)
    d = create_dataset(parent, name, t, (length(rows),))
    isempty(rows) || HDF5.write_dataset(d, t, rows)
    close(d); close(t)
end
hit(event, particle, id; energy=0.1f0, label="NEXT_PRECDR_XENON", x=1234.5f0, time=0.01f0) =
    FixtureHit(event, x, -20.25f0, -1402.4622f0, time, energy, fixed(label, 32), particle, id)
function fixture(path, rows; particles=nothing)
    h5open(path, "w") do f
        mc = create_group(f, "MC")
        write_compound(mc, "hits", rows)
        cfg = [FixtureConfig(fixed("/PhysicsList/RegisterPhysics",64), fixed(value,64))
               for value in ("G4EmStandardPhysics_option4", "G4DecayPhysics")]
        write_compound(mc, "configuration", cfg)
        if particles === nothing
            particles = FixtureParticle[FixtureParticle(id,p,1,100f0,200f0,-300f0)
                                        for id in unique([r.event_id for r in rows]) for p in 1:2]
        end
        write_compound(mc,"particles",particles)
    end
    return path
end

mkpath(joinpath(@__DIR__, "../work"))
mktempdir(joinpath(@__DIR__, "../work"); prefix="event-tests-") do dir
    @testset "Lossless NEXUS extraction and event index" begin
        rows = [hit(0, 1, 0), hit(0, 1, 1; energy=nextfloat(0.1f0)),
                hit(0, 2, 0; energy=1f-10), hit(7, 1, 0; energy=0f0),
                hit(7, 2, 0; energy=2.1f0)]
        input = fixture(joinpath(dir, "sample.next.h5"), rows)
        output = default_output(input)
        original_hash = open(sha256, input)
        @test basename(output) == "sample_hitdf.h5"
        @test basename(default_output("sample.h5")) == "sample_hitdf.h5"
        withenv("NAUSICAA_RAW_DIRECTORY"=>joinpath(dir,"archive")) do
            @test default_output(joinpath(dir,"archive","sample.next.h5")) ==
                joinpath(EventIO.ROOT,"events","sample_hitdf.h5")
            @test default_output(input) == output
        end
        r = convert_hits(input; block_rows=2, chunk_rows=2)
        @test (r.n_hits, r.n_events) == (5, 2)
        @test open(sha256, input) == original_hash
        allhits = read_hitdf(output)
        @test names(allhits) == ["event_id", "particle_id", "hit_id", "x_mm", "y_mm", "z_mm", "energy_MeV", "time_ns"]
        @test eltype.(eachcol(allhits)) == [Int64, Int32, Int32, Float32, Float32, Float32, Float32, Float32]
        @test isequal(allhits.time_ns, [r.time for r in rows])
        @test allhits.x_mm == fill(1234.5f0, 5)
        @test isequal(allhits.energy_MeV, [r.energy for r in rows])
        @test isequal(read_nexus_hits(input; block_rows=1), allhits)
        @test isequal(read_nexus_hits(input; event_id=0, block_rows=2), read_hitdf(output; event_id=0))
        @test read_hitdf(output; event_id=7).particle_id == Int32[1,2]
        @test_throws KeyError read_hitdf(output; event_id=1)
        @test_throws KeyError read_nexus_hits(input; event_id=1)
        @test_throws ArgumentError convert_hits(input)
        @test_throws ArgumentError convert_hits(input; output=input, overwrite=true)
        @test_throws ArgumentError convert_hits(input; block_rows=0)
        @test_throws ArgumentError convert_hits(input; compression=10)
        h5open(output, "r") do f
            @test read(f["events/first_row"]) == [1,4]
            @test read(f["events/n_hits"]) == [3,2]
            @test read(f["events/energy_MeV"]) == [foldl(+, (r.energy for r in rows[1:3]); init=0f0), 2.1f0]
            @test eltype(f["events/energy_MeV"]) == Float32
            @test read_primary_ids(f; event_id=0) == Int32[1,2]
            @test read(f["metadata/label"]) == "NEXT_PRECDR_XENON"
            @test length(unique(read(f["metadata/configuration/param_key"]))) == 1
            @test read(f["metadata/configuration/param_value"]) == ["G4EmStandardPhysics_option4", "G4DecayPhysics"]
            @test read(f["metadata/source_file"]) == "sample.next.h5"
            @test read_hitdf(f; event_id=0).hit_id == [0,1,0]
        end
        @test verify_hitdf(input, output; block_rows=3).n_events == 2
        metadata = read_event_metadata(output)
        @test metadata.vertex_x_mm == Float32[100,100]
        @test metadata.vertex_y_mm == Float32[200,200]
        @test metadata.vertex_z_mm == Float32[-300,-300]
        @test metadata.n_primaries == Int32[2,2]
        # Changing a hit, index entry or event energy must be detected.
        for (column, replacement) in (("hits/x_mm", 1f0), ("events/first_row", Int64(2)),
                                       ("events/energy_MeV", 1f0), ("events/vertex_x_mm", 2f0),
                                       ("hits/time_ns", 20f0), ("events/primary_particle_ids", Int32(9)),
                                       ("events/first_primary_row", Int64(2)))
            convert_hits(input; overwrite=true, block_rows=1)
            h5open(output, "r+") do f
                f[column][1:1] = [replacement]
            end
            @test_throws ErrorException verify_hitdf(input, output)
        end
        convert_hits(input; overwrite=true)
        @test isequal(read_hitdf(output), allhits)
    end
    @testset "Primary vertices must exist and coincide" begin
        for (name,particles) in (
            ("missing",FixtureParticle[]),
            ("secondary_only",[FixtureParticle(0,1,0,0,0,0)]),
            ("different",[FixtureParticle(0,1,1,0,0,0),FixtureParticle(0,2,1,10,0,0)]),
            ("nonfinite_vertex",[FixtureParticle(0,1,1,NaN32,0,0)]))
            input=fixture(joinpath(dir,name*".h5"),[hit(0,1,0)]; particles)
            @test_throws ErrorException convert_hits(input)
            @test !isfile(default_output(input))
        end
        input=fixture(joinpath(dir,"secondary_vertex.h5"),[hit(0,1,0)];
                      particles=[FixtureParticle(0,1,1,1,2,3),FixtureParticle(0,2,0,999,999,999)])
        output=convert_hits(input).output
        @test read_event_metadata(output).vertex_x_mm == Float32[1]
        @test read_event_metadata(output).n_primaries == Int32[1]
    end
    @testset "Reject ambiguous or invalid input without replacing an output" begin
        good = fixture(joinpath(dir, "good.h5"), [hit(0,1,0)])
        out = convert_hits(good).output
        good_hash = open(sha256, out)
        for (name, rows) in (
            ("labels", [hit(0,1,0), hit(0,1,1; label="BUFFER")]),
            ("duplicate", [hit(0,1,0), hit(0,1,0)]),
            ("interleaved", [hit(0,1,0), hit(7,1,0), hit(0,1,1)]),
            ("negative", [hit(0,1,0; energy=-1f0)]),
            ("nonfinite_time", [hit(0,1,0; time=Inf32)]),
            ("nonfinite", [hit(0,1,0; x=NaN32)]))
            bad = fixture(joinpath(dir, name*".h5"), rows)
            @test_throws ErrorException convert_hits(bad; output=out, overwrite=true, block_rows=1)
            @test open(sha256, out) == good_hash
            @test !ispath(default_output(bad))
        end
    end
    @testset "Typed event reader and schema transition" begin
        rows = [hit(42,7,0; time=2f0), hit(42,19,0; time=0f0), hit(103,8,0; time=3f0)]
        particles = [FixtureParticle(42,19,1,1,2,3), FixtureParticle(42,7,1,1,2,3),
                     FixtureParticle(42,33,0,99,99,99), FixtureParticle(103,8,1,4,5,6)]
        input = fixture(joinpath(dir,"actual_ids.h5"), rows; particles)
        output = convert_hits(input).output
        hits, info = read_event(output, Int64(42))
        @test hits.frame == source_frame
        @test hits.time_ns == Float32[2,0]
        @test info.primary_particle_ids == Int32[7,19]
        @test info.vertex_mm == (1f0,2f0,3f0)
        h5open(output,"r") do f
            electron, electron_info = read_event(f, Int64(103))
            @test electron_info.primary_particle_ids == Int32[8]
            @test electron.time_ns == Float32[3]
            hits.x_mm[1] = 0f0
            @test read_event(f,Int64(42))[1].x_mm[1] == 1234.5f0
            @test_throws KeyError read_event(f,Int64(99))
            @test read(f["events/first_primary_row"]) == [1,3]
        end
        h5open(output,"r+") do f
            f["events/first_primary_row"][1:1] = Int64[100]
            @test_throws ErrorException read_event(f,Int64(42))
        end
        # Old schemas remain usable by the existing dataframe reader.
        for schema in ("nausicaa.hitdf.v1", "nausicaa.hitdf.v2")
            convert_hits(input; overwrite=true)
            h5open(output,"r+") do f
                HDF5.delete_attribute(f,"schema")
                attributes(f)["schema"] = schema
                delete_object(f,"hits/time_ns")
                delete_object(f,"events/primary_particle_ids")
                delete_object(f,"events/first_primary_row")
                @test nrow(read_hitdf(f)) == 3
                @test !hasproperty(read_hitdf(f), :time_ns)
                err = try read_event(f,Int64(42)); nothing catch e; e end
                @test err isa ErrorException
                @test occursin("extract_hits.jl --overwrite", sprint(showerror,err))
            end
        end
        bad = fixture(joinpath(dir,"duplicate_primary.h5"), [hit(1,7,0)];
            particles=[FixtureParticle(1,7,1,0,0,0),FixtureParticle(1,7,1,0,0,0)])
        @test_throws ErrorException convert_hits(bad)
    end
    @testset "Empty hit sample" begin
        input = fixture(joinpath(dir, "empty.h5"), FixtureHit[])
        r = convert_hits(input)
        @test (r.n_hits, r.n_events) == (0,0)
        @test isempty(read_hitdf(r.output))
        @test isempty(read_nexus_hits(input))
        @test_throws KeyError read_hitdf(r.output; event_id=0)
    end
end
