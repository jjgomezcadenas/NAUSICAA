"""Track selection and track features from the hits of one track of one event.

The functions take the hit positions (mm) and energies (keV) as vectors and
the settings as keywords, so the ITACA and iItaca chains use the same code with
their own cards. Method and thresholds: notes/html/event_analysis/track_selection.html.

Load `src/event_analysis/EventAnalysis.jl` first: the graph of the main track
and its endpoint energies come from there.
"""
module TrackSelection
using ..EventAnalysis

export cell_centred, point_components, track_selection, track_features

"""Coordinates shifted so that hits on a lattice of pitch `width` sit at the centres of the cells of that width.

Hits on such a lattice at an arbitrary position can lie on a cell boundary,
where rounding splits one layer between two cells. When every coordinate is
within 2% of the pitch of the lattice of the first one, all are shifted by the
same amount; otherwise (a readout pitch different from the cell width) they
are returned unchanged.
"""
function cell_centred(x::AbstractVector{Float32},width::Float32)
    isempty(x)&&return x
    steps=(Float64.(x).-Float64(x[1]))./width;all(s->abs(s-round(s))<0.02,steps)||return x
    Float32.(Float64.(x).-(mod(Float64(x[1]),Float64(width))-width/2))
end

"""Label of the connected component of each point: points closer than `radius` are connected. Labels run from 1."""
function point_components(x,y,z,radius::Real)
    n=length(x);parent=collect(1:n);r2=Float64(radius)^2
    root(i)=(while parent[i]!=i;parent[i]=parent[parent[i]];i=parent[i];end;i)
    for i in 1:n-1, j in i+1:n
        (Float64(x[i])-x[j])^2+(Float64(y[i])-y[j])^2+(Float64(z[i])-z[j])^2<r2 || continue
        a=root(i);b=root(j);a==b||(parent[max(a,b)]=min(a,b))
    end
    names=Dict{Int,Int}();[get!(names,root(i),length(names)+1) for i in 1:n]
end

"""
    track_selection(x_mm, y_mm, z_mm, energy_keV; hit_cut_keV, segment_connection_mm,
        segment_min_hits, deposit_min_keV, max_segments, connection_mm)
        -> (segments, deposits, deposit_keV, segments_joined, selected)

Track selection from the components of the hits on the readout lattice.

The hits above `hit_cut_keV` are connected, at their own positions, with the
radius `segment_connection_mm` (touching cells). A component with at least
`segment_min_hits` hits is a track segment. A smaller component with at least
`deposit_min_keV` is a detached deposit: the K X-ray of a photoelectric
absorption, a bremsstrahlung photon. Smaller components are ignored. The event
is selected when it has no detached deposit and either one segment, or at most
`max_segments` segments that belong to one component with the radius
`connection_mm`: a track with a short gap. Segments that stay apart are
separate electrons.

Returns the number of segments and of detached deposits, the largest deposit
energy (keV), whether the segments join, and the selection.
"""
function track_selection(x_mm,y_mm,z_mm,energy_keV;hit_cut_keV::Real,segment_connection_mm::Real,segment_min_hits::Integer,
                         deposit_min_keV::Real,max_segments::Integer,connection_mm::Real)
    keep=findall(>(Float32(hit_cut_keV)),energy_keV);none=(segments=0,deposits=0,deposit_keV=0f0,segments_joined=false,selected=false)
    isempty(keep)&&return none
    x=x_mm[keep];y=y_mm[keep];z=z_mm[keep];e=energy_keV[keep]
    label=point_components(x,y,z,segment_connection_mm);n=maximum(label);hits=zeros(Int,n);energy=zeros(Float32,n)
    for (i,l) in enumerate(label);hits[l]+=1;energy[l]+=e[i];end
    segment=hits.>=Int(segment_min_hits);deposit=.!segment .& (energy.>=Float32(deposit_min_keV))
    wide=point_components(x,y,z,connection_mm);joined=length(unique(wide[i] for i in eachindex(label) if segment[label[i]]))<=1
    segments=count(segment);deposits=count(deposit)
    (segments=segments,deposits=deposits,deposit_keV=deposits==0 ? 0f0 : maximum(energy[deposit]),segments_joined=joined&&segments>0,
     selected=deposits==0&&segments>=1&&(segments==1||(segments<=Int(max_segments)&&joined)))
end

"""
    track_features(x_mm, y_mm, z_mm, energy_keV; hit_cut_keV, voxel_width_mm, connection_mm, blob_radius_mm)

Graph features of one track of one event, from its hits.

Hits above `hit_cut_keV` are voxelised, in cells centred on the readout lattice
where its pitch equals `voxel_width_mm`, and connected with the radius
`connection_mm`. The main track is the component with most voxels. Its length
is the length of the tree path between its two endpoints; its extensions are
the ranges of its voxel centres plus one voxel width. E1 >= E2 are the energies
in spheres of radius `blob_radius_mm` around the two endpoints.

Returns `hits`, `components`, `main_keV`, `second_keV`, `main_voxels`,
`length_mm`, `extent_x_mm`, `extent_y_mm`, `extent_z_mm`, `E1_keV`, `E2_keV`.
"""
function track_features(x_mm::AbstractVector{Float32},y_mm::AbstractVector{Float32},z_mm::AbstractVector{Float32},energy_keV::AbstractVector{Float32};
                        hit_cut_keV::Real,voxel_width_mm::Real,connection_mm::Real,blob_radius_mm::Real)
    cut=Float32(hit_cut_keV);width=Float32(voxel_width_mm);radius=Float32(connection_mm);blob=Float64(blob_radius_mm)
    keep=energy_keV .> cut
    nothing_found=(hits=count(keep),components=0,main_keV=0f0,second_keV=0f0,main_voxels=0,length_mm=NaN32,
        extent_x_mm=NaN32,extent_y_mm=NaN32,extent_z_mm=NaN32,E1_keV=NaN32,E2_keV=NaN32)
    any(keep)||return nothing_found
    df=EventAnalysis.DataFrame(x=cell_centred(x_mm[keep],width),y=cell_centred(y_mm[keep],width),z=cell_centred(z_mm[keep],width),energy=energy_keV[keep]./1000f0)
    tracks=make_tracks(voxelize_event(df,width);method="KDT",max_distance_mm=radius,energy_threshold_kev=0f0)
    isempty(tracks)&&return nothing_found
    energies=sort(Float32[1000sum(t.voxels.energy) for t in tracks];rev=true)
    track=tracks[argmax([size(t.voxels,1) for t in tracks])];v=track.voxels
    u,w,path,_,_=EventAnalysis.find_extremes_mst_diameter(track,EventAnalysis.extract_coords(track))
    len=sum(sqrt((v.x[path[i+1]]-v.x[path[i]])^2+(v.y[path[i+1]]-v.y[path[i]])^2+(v.z[path[i+1]]-v.z[path[i]])^2) for i in 1:length(path)-1;init=0.0)
    b=EventAnalysis.energy_in_spheres_around_extremes(track,(extremes=(v[u,:],v[w,:]),),blob)
    e1,e2=1000 .*(max(b.blob1_energy,b.blob2_energy),min(b.blob1_energy,b.blob2_energy))
    extent(c)=Float32(maximum(c)-minimum(c)+width)
    (hits=count(keep),components=length(tracks),main_keV=1000f0*Float32(sum(v.energy)),second_keV=length(energies)<2 ? 0f0 : energies[2],
        main_voxels=size(v,1),length_mm=Float32(len),extent_x_mm=extent(v.x),extent_y_mm=extent(v.y),extent_z_mm=extent(v.z),
        E1_keV=Float32(e1),E2_keV=Float32(e2))
end

end # module
