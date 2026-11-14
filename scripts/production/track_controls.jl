"""Sum hit energy into XY, XZ and YZ cells for presentation."""
function hit_projections(h,pitch,z_width;z_centres=nothing)
 x=sort(unique(h.x_mm));y=sort(unique(h.y_mm));z=z_centres===nothing ? sort(unique(h.z_mm)) : sort(Float32.(z_centres))
 # Channel/time coordinates cover a regular grid; include empty intervening cells.
 grid(v,step)=Float32.(range(first(v),last(v);length=round(Int,(last(v)-first(v))/step)+1))
 x=isempty(x) ? Float32[0] : grid(x,pitch);y=isempty(y) ? Float32[0] : grid(y,pitch)
 # A final incomplete electron bucket can have a different centre.
 z_edges=length(z)==1 ? Float32[z[1]-z_width/2,z[1]+z_width/2] :
  Float32[vcat(z[1]-(z[2]-z[1])/2,(z[1:end-1]+z[2:end])/2,z[end]+(z[end]-z[end-1])/2)...]
 xy=zeros(Float32,length(x),length(y));xz=zeros(Float32,length(x),length(z));yz=zeros(Float32,length(y),length(z))
 zmap=Dict(v=>i for (i,v) in enumerate(z))
 for k in eachindex(h.charge)
  i=round(Int,(h.x_mm[k]-first(x))/pitch)+1;j=round(Int,(h.y_mm[k]-first(y))/pitch)+1
  xy[i,j]+=h.energy_keV[k];xz[i,zmap[h.z_mm[k]]]+=h.energy_keV[k];yz[j,zmap[h.z_mm[k]]]+=h.energy_keV[k]
 end
 Dict("x_edges_mm"=>vcat(x.-pitch/2,last(x)+pitch/2),"y_edges_mm"=>vcat(y.-pitch/2,last(y)+pitch/2),
 "z_edges_mm"=>z_edges,"xy_keV"=>xy,"xz_keV"=>xz,"yz_keV"=>yz,"hit_count"=>length(h.charge))
end

"""Combine histograms with common edges; moments include values beyond the plotted range."""
function combine_histograms(histograms)
 n=sum(h["entries"] for h in histograms)
 mu=sum(h["entries"]*h["mean"] for h in histograms)/n
 second=sum(h["entries"]*(h["std"]^2+h["mean"]^2) for h in histograms)/n
 Dict("edges"=>first(histograms)["edges"],"counts"=>sum(h["counts"] for h in histograms),"entries"=>n,
 "mean"=>mu,"std"=>sqrt(max(0.,second-mu^2)),"underflow"=>sum(h["underflow"] for h in histograms),"overflow"=>sum(h["overflow"] for h in histograms))
end

"""Measure hit-cut effects using input counts only to label signal and noise.
Energy loss is the expected input signal arriving in readout cells, with the
same mean attachment correction and calibration as reconstructed energy.
"""
function hit_cut_statistics(values,valid,input,scale,threshold)
 good=valid .& isfinite.(values)
 positive=good .& (values.>0)
 kept=positive .& (Float64.(values).*scale.>threshold)
 signal=input.>0;noise=.!signal
 input_energy(mask)=sum(Float64,input[mask];init=0.)*scale
 (;positive_hits=count(positive),kept_hits=count(kept),
 noise_hits_before=count(positive.&noise),noise_hits_after=count(kept.&noise),
 signal_hits_before=count(positive.&signal),signal_hits_after=count(kept.&signal),
 input_energy_valid_keV=input_energy(good),input_energy_before_keV=input_energy(positive),
 input_energy_after_keV=input_energy(kept),input_energy_invalid_keV=input_energy(.!good))
end

"""Rebuild the existing kNN graph before and after the hit cut.
Distance is edge_sigma_scale × the carrier's transverse diffusion sigma.
Signal-only diagnostics separate removed noise from lost signal connections.
Endpoints and length use the existing MST-diameter calculation on the
highest-energy component; endpoint pairs are matched without ordering.
"""
function hit_cut_topology(values,valid,input,channels,hits,scale,threshold,sigma_t,k,edge_sigma_scale;compute_endpoints=true)
 selected=findall(valid .& isfinite.(values) .& (values.>0));n=length(selected)
 kept=[Float64(values[q])*scale>threshold for q in selected]
 signal=[input[q]>0 for q in selected]
 distance=Float32(edge_sigma_scale*sigma_t)
 function describe(mask;endpoints=false)
  indices=findall(mask)
  df=EventAnalysis.DataFrame(x=hits.x_mm[indices],y=hits.y_mm[indices],z=hits.z_mm[indices],energy=hits.energy_keV[indices]./1000f0)
  tracks=make_tracks(df;method="kNN",k=Int(k),max_distance_mm=distance,energy_threshold_kev=0f0)
  isempty(tracks)&&return (components=0,largest_fraction=0.,endpoints=fill(NaN32,3,2),length_mm=NaN)
  energies=[sum(Float64,t.voxels.energy) for t in tracks]
  track=tracks[argmax(energies)];ends=fill(NaN32,3,2);length_mm=NaN
  if endpoints
   a,b,_,_,length_mm=find_extremes_mst_diameter(track,EventAnalysis.extract_coords(track))
   d=track.voxels
   ends[:,1]=Float32[d.x[a],d.y[a],d.z[a]];ends[:,2]=Float32[d.x[b],d.y[b],d.z[b]]
  end
  (components=length(tracks),largest_fraction=maximum(energies)/max(sum(energies),eps()),endpoints=ends,length_mm)
 end
 before=describe(trues(n));after=describe(kept)
 sb=describe(signal;endpoints=compute_endpoints);sa=describe(signal.&kept;endpoints=compute_endpoints)
 separation(a,b)=sqrt(sum(abs2,a-b))
 direct=(separation(sb.endpoints[:,1],sa.endpoints[:,1]),separation(sb.endpoints[:,2],sa.endpoints[:,2]))
 swapped=(separation(sb.endpoints[:,1],sa.endpoints[:,2]),separation(sb.endpoints[:,2],sa.endpoints[:,1]))
 moves=sum(direct)<=sum(swapped) ? direct : swapped
 Dict("sigma_t_mm"=>Float64(sigma_t),"max_distance_mm"=>Float64(distance),
 "components_before"=>before.components,"components_after"=>after.components,
 "signal_components_before"=>sb.components,"signal_components_after"=>sa.components,
 "signal_largest_fraction_before"=>sb.largest_fraction,"signal_largest_fraction_after"=>sa.largest_fraction,
 "length_before_mm"=>sb.length_mm,"length_after_mm"=>sa.length_mm,
 "endpoint_mean_shift_mm"=>sum(moves)/2,"endpoint_max_shift_mm"=>maximum(moves))
end
