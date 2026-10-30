using Test, Random, Statistics, HDF5

@testset "Production moment sampling" begin
    r=Xoshiro(17)
    low=[sample_integer_moments(r,0.2,0.15) for _ in 1:20_000]
    @test all(n->n in (0,1),low)
    @test abs(mean(low)-0.2)<0.01
    @test abs(var(low)-0.16)<0.01
    draws=[sample_integer_moments(r,20.3,0.15) for _ in 1:30_000]
    @test abs(mean(draws)-20.3)<0.06
    @test abs(var(draws)-3.045)<0.12
    floor_draws=[sample_integer_moments(r,1.3,0.15) for _ in 1:10_000]
    @test all(n->n in (1,2),floor_draws)
    @test abs(var(floor_draws)-0.21)<0.015
    @test_throws ArgumentError sample_integer_moments(r,1.0,1.0)
end

@testset "Shared pair and attachment draws" begin
    hits=EventHits(Int64(51),source_frame,Float32[0,2],Float32[0,0],
        Float32[0,2],Float32[0.00000219,0.000219],zeros(Float32,2),
        Int32[1,1],Int32[1,2])
    a=sample_fano_pairs(Xoshiro(42),hits,21.9f0,0.15f0)
    b=sample_fano_pairs(Xoshiro(42),hits,21.9f0,0.15f0)
    @test a==b
    @test a[1] in (0,1)
    @test 0<=a[2]<=20
    fano=[sample_fano_pairs(Xoshiro(i),hits,21.9f0,0.15f0)[2]
          for i in 1:8000]
    @test abs(mean(fano)-10)<0.05
    @test abs(var(fano)-1.5)<0.1
    @test sample_attached_electrons(Xoshiro(2),a,0f0,15f0)==a
    attached=sample_attached_electrons(Xoshiro(2),a,1.5f0,15f0)
    @test all(attached.<=a)
    samples=[sample_attached_electrons(Xoshiro(i),Int64[1000],1.5f0,15f0)[1]
             for i in 1:1000]
    p=exp(-0.1)
    @test abs(mean(samples)-1000p)<1.1
    @test abs(var(samples)-1000p*(1-p))<8
    @test_throws ArgumentError sample_fano_pairs(Xoshiro(1),hits,21.9f0,1f0)
    @test_throws ArgumentError sample_attached_electrons(Xoshiro(1),a,1f0,0f0)
end

@testset "EL and optical moments" begin
    r=Xoshiro(99)
    emitted=[sample_el_photons(r,2,100.2f0,0.1f0) for _ in 1:25_000]
    @test abs(mean(emitted)-200.4)<0.08
    @test abs(var(emitted)-20.04)<0.9
    @test sample_el_photons(r,0,100f0,0.1f0)==0
    incident=[sample_incident_photons(r,100,0.2f0) for _ in 1:20_000]
    @test abs(mean(incident)-20)<0.15
    @test abs(var(incident)-16)<0.55
    @test_throws ArgumentError sample_el_photons(r,1,100f0,1f0)
    @test_throws ArgumentError sample_incident_photons(r,100,1.1f0)
end

@testset "Sparse histories and HDF5 round trip" begin
    placement=EventPlacement(51,100f0,50f0)
    ion_grid=GridSpec(2f0,2f0,(0f0,0f0,0f0),(0f0,0f0),Int8(1))
    photon_grid=GridSpec(10f0,2f0,(0f0,0f0,2.1f0),(0f0,0f0),Int8(-1))
    hits=EventHits(Int64(51),source_frame,Float32[0,2.1],Float32[0,0],
        Float32[0,2.1],Float32[0.001,0.001],zeros(Float32,2),
        Int32[1,1],Int32[1,2])
    ion=collect_ion_pad_history(Xoshiro(4),hits,Int64[1,2],placement,ion_grid;
        sigma_t_mm=0f0,sigma_l_mm=0f0)
    @test ion.initial_pairs==3
    @test ion.cells.pixel_ix==Int64[0,1]
    @test ion.cells.bucket==Int64[0,1]
    @test ion.cells.count==Int64[1,2]
    replay=collect_ion_pad_history(Xoshiro(4),hits,Int64[1,2],
        placement,ion_grid;sigma_t_mm=0f0,sigma_l_mm=0f0)
    @test ion.cells.count==replay.cells.count
    @test ion.cells.pixel_ix==replay.cells.pixel_ix
    photon=collect_photon_sipm_history(Xoshiro(7),hits,Int64[1,2],
        placement,photon_grid;sigma_t_mm=0f0,sigma_l_mm=0f0,
        yield_per_e=100f0,J=0.1f0,acceptance=1f0,
        gas_model=:xe_he15,initial_pairs=3)
    @test photon.initial_pairs==3
    @test photon.surviving_electrons==3
    @test sum(photon.cells.count)==photon.emitted_photons
    @test sort(photon.cells.bucket)==Int64[0,1]
    @test_throws ArgumentError GridSpec(2f0,2f0,(0f0,0f0,0f0),
        (2f0,0f0),Int8(1))
    @test_throws ArgumentError SparseCells(Int64[1,1],Int64[0,0],
        Int64[0,0],Int64[1,2])
    mutable_counts=Int64[3]
    owned=SparseCells(Int64[0],Int64[0],Int64[0],mutable_counts)
    mutable_counts[1]=99
    @test owned.count==Int64[3]
    @test_throws ArgumentError IonPadHistory(placement,ion_grid,ion.cells,2)
    distant=EventHits(Int64(51),source_frame,Float32[0,1000],Float32[0,0],
        Float32[0,1000],Float32[0.001,0.001],zeros(Float32,2),
        Int32[1,1],Int32[1,2])
    extended=collect_ion_pad_history(Xoshiro(4),distant,Int64[1,2],
        placement,ion_grid;sigma_t_mm=0f0,sigma_l_mm=0f0)
    @test sum(extended.cells.count)==3
    @test maximum(extended.cells.pixel_ix)>300
    mktemp() do path,io
        close(io)
        h5open(path,"w") do f
            write_ion_history(create_group(f,"ion"),ion)
            write_photon_history(create_group(f,"photon"),photon)
        end
        h5open(path,"r") do f
            restored_ion=read_ion_history(f["ion"])
            restored_photon=read_photon_history(f["photon"])
            @test restored_ion.initial_pairs==ion.initial_pairs
            @test restored_ion.cells.count==ion.cells.count
            @test restored_ion.grid.phase_mm==ion.grid.phase_mm
            @test restored_photon.gas_model==:xe_he15
            @test restored_photon.cells.count==photon.cells.count
            @test restored_photon.emitted_photons==photon.emitted_photons
        end
        h5open(path,"r+") do f
            write(attributes(f["ion"])["schema"],"invalid")
        end
        h5open(path,"r") do f
            @test_throws ErrorException read_ion_history(f["ion"])
        end
    end
    @test_throws ArgumentError collect_photon_sipm_history(Xoshiro(7),hits,
        Int64[1,2],placement,photon_grid;sigma_t_mm=0f0,sigma_l_mm=0f0,
        yield_per_e=100f0,J=0.1f0,acceptance=1f0,
        gas_model=:xe_he15,initial_pairs=2)
    empty_cells=SparseCells(Int64[],Int64[],Int64[],Int64[])
    zero_history=IonPadHistory(placement,ion_grid,empty_cells,0)
    mktemp() do path,io
        close(io)
        h5open(path,"w") do f
            write_ion_history(create_group(f,"zero"),zero_history)
        end
        h5open(path,"r") do f
            @test isempty(read_ion_history(f["zero"]).cells.count)
        end
    end
end
