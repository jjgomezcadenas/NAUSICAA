using StatsBase

@testset "MC energy display maps" begin
    h=EventHits(1,source_frame,Float32[-1,-0.5,0,0.5,1,1],Float32[0,0,1,1,2,2],
        fill(3f0,6),Float32[0,.1,.2,.3,.4,.5],zeros(Float32,6),ones(Int32,6),Int32.(1:6))
    before=copy(h.energy_MeV)
    for nbins in (2,100)
        m=bin_mc_energy(h;nbins)
        @test m.hit_count==6
        @test m.z_edges_mm==Float32[2.75,3.25]
        @test m.x_edges_mm==Float32.(collect(only(fit(Histogram,h.x_mm;nbins).edges)))
        for (axes,energy,counts) in (((m.x_edges_mm,m.y_edges_mm),m.xy_energy_MeV,m.xy_hit_count),
                                   ((m.x_edges_mm,m.z_edges_mm),m.xz_energy_MeV,m.xz_hit_count),
                                   ((m.y_edges_mm,m.z_edges_mm),m.yz_energy_MeV,m.yz_hit_count))
            @test size(energy)==(length(axes[1])-1,length(axes[2])-1)
            @test sum(counts)==6
            @test isapprox(Float64(sum(energy)),sum(Float64.(h.energy_MeV));rtol=2e-6)
        end
        # Independent exhaustive membership check against original hit coordinates.
        reference=zeros(Float64,size(m.xy_energy_MeV))
        for r in eachindex(h.energy_MeV)
            i=findfirst(j->m.x_edges_mm[j]<=h.x_mm[r] &&
                (h.x_mm[r]<m.x_edges_mm[j+1] || j==length(m.x_edges_mm)-1),1:length(m.x_edges_mm)-1)
            j=findfirst(j->m.y_edges_mm[j]<=h.y_mm[r] &&
                (h.y_mm[r]<m.y_edges_mm[j+1] || j==length(m.y_edges_mm)-1),1:length(m.y_edges_mm)-1)
            reference[i,j]+=Float64(h.energy_MeV[r])
        end
        @test isapprox(Float64.(m.xy_energy_MeV),reference;rtol=2e-6)
        @test sum(m.occupied.hit_count)==6
        @test sum(m.occupied.energy_MeV)≈sum(h.energy_MeV)
        @test h.energy_MeV==before
    end
    @test_throws ArgumentError bin_mc_energy(h;nbins=0)
    @test EventAnalysis._histogram_index(1f0,Float32[0,.5,1])==2
    # Real events: original energy versus every saved projection at two resolutions.
    for sample in ("0nubb","electron_0nubb_energy")
        hits,_=read_event(joinpath(ROOT,"events",sample*"_hitdf.h5"),Int64(0))
        expected=sum(Float64.(hits.energy_MeV))
        for nbins in (100,200)
            m=bin_mc_energy(hits;nbins)
            for array in (m.xy_energy_MeV,m.xz_energy_MeV,m.yz_energy_MeV,m.occupied.energy_MeV)
                @test isapprox(sum(Float64.(array)),expected;rtol=2e-6)
            end
        end
    end
end
