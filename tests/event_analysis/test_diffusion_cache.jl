include(joinpath(ROOT,"src/event_analysis/batch/pixel_diffusion_cache.jl"))
@testset "Diffusion cache independent of collection grid" begin
    mktempdir(joinpath(ROOT,"work");prefix="cache-test-") do dir
        source=joinpath(dir,"source.h5");h5open(f->(f["test"]=[1]),source,"w")
        c=Dict("transport"=>Dict("ionisation_energy_eV"=>21.9,"min_cm"=>10.,"max_cm"=>150.,"nbins"=>100,"nsigma"=>3.),
            "run"=>Dict("seed"=>2026,"diffusion_cache_directory"=>dir,"n_events"=>1),
            "pads"=>Dict("pitch_mm"=>2.,"z_bin_width_mm"=>2.))
        det=Dict("gas"=>Dict("temperature_K"=>293.15),"bias"=>Dict("bulk_field_V_m"=>20000.),"geometry"=>Dict("pitch_m"=>.002))
        cache,path=open_pixel_diffusion_cache(c,source,"toy",det,ROOT)
        @test read_cached_diffusion(cache,7)===nothing
        bins=DataFrame(event_id=Int64[7,7],x=Float32[.5,1.5],y=Float32[.5,.5],z=Float32[.5,.5],energy=Float32[.001,.002],ions=Int64[10,20])
        d=(;bins,sampled_ions=Int64(32),retained_ions=Int64(30),outside_ions=Int64(2),
            x_edges_mm=Float32[0,1,2],y_edges_mm=Float32[0,1],z_edges_mm=Float32[0,1])
        save_cached_diffusion(cache,7,d,80f0,1.4f0,UInt64(123))
        result=read_cached_diffusion(cache,7)
        @test result.d.bins==bins && result.d.x_edges_mm==d.x_edges_mm
        @test result.seed==123 && result.d.sampled_ions==32 && result.d.retained_ions==30
        create_group(cache,"events/8")
        @test read_cached_diffusion(cache,8)===nothing
        save_cached_diffusion(cache,8,d,80f0,1.4f0,UInt64(123))
        @test sum(read_cached_diffusion(cache,8).d.bins.ions)==30
        close(cache)
        grid=deepcopy(c);grid["pads"]["pitch_mm"]=3.;grid["pads"]["z_bin_width_mm"]=3.;grid["run"]["n_events"]=8000
        modified=deepcopy(det);modified["geometry"]["pitch_m"]=.003
        cache,same=open_pixel_diffusion_cache(grid,source,"toy",modified,ROOT)
        @test same==path && read_cached_diffusion(cache,7).d.bins==bins
        close(cache)
        changed=deepcopy(c);changed["transport"]["nbins"]=200
        cache,different=open_pixel_diffusion_cache(changed,source,"toy",det,ROOT)
        @test different!=path && read_cached_diffusion(cache,7)===nothing
        close(cache)
    end
end
