using HDF5, Random
include(joinpath(@__DIR__,"../../src/event_analysis/EventAnalysis.jl"))
using .EventAnalysis

function key_counts(f,id)
    g=f["events/$(id)"]
    keys=zip(read(g["pad_ix"]),read(g["pad_iy"]),read(g["slice_index"]))
    return Dict(Tuple(k)=>v for (k,v) in zip(keys,read(g["signal_count"])))
end

function audit(sample,id)
    cp="results/three_clouds_8000/$(sample)_ion.h5"
    rp="results/grid_readout_8000/$(sample)_ion.h5"
    h5open(cp,"r") do cf
        h5open(rp,"r") do rf
            c=cf["events/$(id)"];r=rf["events/$(id)"]
            source=read(cf["metadata/source_path"])
            hits,_=h5open(f->read_event(f,Int64(id)),source,"r")
            ions=ion_counts(hits,read(cf["metadata/w_eV"])).sampled
            seed=parse(UInt64,read(c["diffusion_seed"]))
            rng=Xoshiro(seed)
            sx=read(c["sigma_t_mm"]);sz=read(c["sigma_l_mm"])
            diffusion_edges=Tuple(read(c["$(a)_edges_mm"]) for a in "xyz")
            # Derive the old grid lattice from its saved reference and phase.
            d=DataFrames.DataFrame(event_id=fill(Int64(id),length(c["ions"])),
                x=read(c["x"]),y=read(c["y"]),z=read(c["z"]),
                energy=Float32.(read(c["ions"])).*read(cf["metadata/w_eV"])./1f6,
                ions=read(c["ions"]))
            axes=pad_axes(d;reference_mm=(read(r["reference_x_mm"]),
                read(r["reference_y_mm"]),read(r["reference_z_mm"])),
                pitch_mm=read(rf["metadata/pitch_mm"]),
                z_bin_width_mm=read(rf["metadata/bucket_mm"]),
                grid_offset_mm=(read(r["grid_offset_x_mm"]),read(r["grid_offset_y_mm"])),
                minimum_xy_span_mm=160f0)
            @assert axes.x_edges_mm==read(r["x_edges_mm"])
            @assert axes.y_edges_mm==read(r["y_edges_mm"])
            grid=(read(r["x_edges_mm"]),read(r["y_edges_mm"]),axes.z_edges_mm)
            direct=Dict{NTuple{3,Int64},Int64}()
            replay_old=Dict{NTuple{3,Int64},Int64}()
            retained=0
            outside_readout=0
            reassigned=0
            for i in eachindex(ions)
                for _ in 1:ions[i]
                    x=hits.x_mm[i]+randn(rng,Float32)*sx
                    y=hits.y_mm[i]+randn(rng,Float32)*sx
                    z=hits.z_mm[i]+randn(rng,Float32)*sz
                    q=(x,y,z)
                    bin=Tuple(EventAnalysis._histogram_index(q[j],diffusion_edges[j]) for j in 1:3)
                    minimum(bin)>0 || continue
                    key=Tuple(Int64(begin
                        idx=EventAnalysis._histogram_index(q[j],grid[j])
                        idx>0 ? idx : floor(Int,(q[j]-first(grid[j]))/(grid[j][2]-grid[j][1]))+1
                    end) for j in 1:3)
                    any(key[j]<1 || key[j]>length(grid[j])-1 for j in 1:3) && (outside_readout+=1)
                    centres=Tuple(diffusion_edges[j][bin[j]]/2f0+
                        diffusion_edges[j][bin[j]+1]/2f0 for j in 1:3)
                    old_key=Tuple(Int64(EventAnalysis._histogram_index(centres[j],grid[j])) for j in 1:3)
                    minimum(old_key)>0 || error("Histogram centre outside original readout grid")
                    replay_old[old_key]=get(replay_old,old_key,0)+1
                    key!=old_key && (reassigned+=1)
                    direct[key]=get(direct,key,0)+1
                    retained+=1
                end
            end
            old=key_counts(rf,id)
            @assert retained==read(c["retained_ions"])==sum(values(old))
            @assert replay_old==old
            if isfile("results/grid_readout_direct_8000/$(sample)_ion.h5")
                h5open("results/grid_readout_direct_8000/$(sample)_ion.h5","r") do nf
                    ng=nf["events/$(id)"]
                    new_edges=(read(ng["x_edges_mm"]),read(ng["y_edges_mm"]),read(ng["z_edges_mm"]))
                    shift=ntuple(j->round(Int64,(first(new_edges[j])-first(grid[j]))/2f0),3)
                    new=Dict((k[1]+shift[1],k[2]+shift[2],k[3]+shift[3])=>v
                             for (k,v) in key_counts(nf,id))
                    @assert direct==new
                end
            end
            l1=sum(abs(get(direct,k,0)-get(old,k,0)) for k in union(keys(direct),keys(old)))
            axis_tv=Float64[]
            for axis in 1:3
                a=Dict{Int64,Int64}();b=Dict{Int64,Int64}()
                for (k,v) in direct;a[k[axis]]=get(a,k[axis],0)+v;end
                for (k,v) in old;b[k[axis]]=get(b,k[axis],0)+v;end
                push!(axis_tv,sum(abs(get(a,k,0)-get(b,k,0)) for k in union(keys(a),keys(b)))/(2retained))
            end
            println(sample," event ",id,": retained=",retained,
                " lower_bound_moved_fraction=",l1/(2retained),
                " actual_cell_reassignment_fraction=",reassigned/retained,
                " direct_cells=",length(direct)," old_cells=",length(old),
                " z_width_mm=",median(diff(diffusion_edges[3])),
                " outside_old_grid=",outside_readout,
                " marginal_TV_xyz=",axis_tv)
        end
    end
end

using DataFrames, Statistics
for id in (0,3790)
    audit("electrons",id)
end
