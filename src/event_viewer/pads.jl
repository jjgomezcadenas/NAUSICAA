"""Collect the full transported ion cloud on the detector's regular hole/pad grid.

Use the retained carrier positions' ion-weighted x/y centre and minimum z as reference, with one
reproducible uniform transverse phase in [0,pitch) per event realisation. The
phase uses a separate RNG stream from diffusion and drift. Read pitch and ion
speed from the saved transport detector card, explicitly converting SI to mm
and mm/ms. Replay the same diffusion seed to determine the reference and to
assign each carrier directly to its physical pad/time bucket.
Return the sparse pad-by-time history, axes, projections, phase, reference and
SI charge per ion. Every retained ion is collected, including low-count bins.
The geometric approximation models ideal routing and collection. Reconstruction
A/B/C decisions are separate downstream observations.
"""
function collect_pads(t)
    c=CONFIG["pads"];det=TOML.parse(t.detector_toml)
    pitch=Float32(det["geometry"]["pitch_m"])*1000f0
    speed=Float32(det["gas"]["bulk_speed_m_s"]) # 1 m/s = 1 mm/ms.
    seed=parse(UInt64,t.seed)
    diffusion_seed=seed⊻0x9e3779b97f4a7c15
    counts=ion_counts(t.hits,t.ionisation_energy_eV).sampled
    reference=diffused_carrier_reference(Xoshiro(diffusion_seed),t.hits,counts;
        sigma_t_mm=t.sigma_t_mm,sigma_l_mm=t.sigma_l_mm,
        nsigma=Float32(CONFIG["transport"]["nsigma"])).reference_mm
    phase=sample_grid_offset_mm(Xoshiro(seed⊻0xd1b54a32d192ed03),pitch)
    history=diffuse_to_grid(Xoshiro(diffusion_seed),t.hits,counts;
        sigma_t_mm=t.sigma_t_mm,sigma_l_mm=t.sigma_l_mm,
        nsigma=Float32(CONFIG["transport"]["nsigma"]),reference_mm=reference,
        grid_offset_mm=phase,pitch_mm=pitch,
        bucket_width_mm=Float32(c["z_bin_width_mm"]),
        drift_speed_mm_ms=speed,drift_direction=Int(c["drift_direction"]),
        minimum_xy_span_mm=Float32(c["minimum_xy_span_mm"]))
    first_edge_indices=(round(Int64,(first(history.x_edges_mm)-reference[1]-phase[1])/pitch),
        round(Int64,(first(history.y_edges_mm)-reference[2]-phase[2])/pitch),
        round(Int64,(first(history.z_edges_mm)-reference[3])/Float32(c["z_bin_width_mm"])))
    axes=(;reference_mm=reference,grid_offset_mm=phase,first_edge_indices)
    projections=pad_projections(history)
    @assert history.total_ions==t.d.retained_ions
    return (;axes,history,projections,pitch_mm=pitch,drift_speed_mm_ms=speed,
        z_bin_width_mm=Float32(c["z_bin_width_mm"]),charge_per_ion_C=1.602176634e-19)
end

"""Save pad charge totals and complete sparse pad/time histories with integer ledgers.

Display XY sums every time slice; X-time/Y-time sum the opposite pad coordinate.
The 3D cells show individual pad/time counts. Axes are source x/y in mm and
relative arrival time in ms. Save source-z edges and the direction convention
alongside the time axis. All reductions use Int64; physical signal charge is
ion count times the positive elementary charge. No collection loss is introduced.
"""
function save_pads(f,p)
    b="pads/";h=p.history;a=p.axes;q=p.projections
    for key in (:pitch_mm,:z_bin_width_mm,:drift_speed_mm_ms,:charge_per_ion_C)
        f[b*string(key)]=getproperty(p,key)
    end
    f[b*"reference_mm"]=collect(a.reference_mm);f[b*"grid_offset_mm"]=collect(a.grid_offset_mm)
    f[b*"first_edge_indices"]=collect(a.first_edge_indices)
    f[b*"total_ions"]=h.total_ions;f[b*"total_charge_C"]=h.total_ions*p.charge_per_ion_C
    f[b*"total_charge_fC"]=h.total_ions*p.charge_per_ion_C*1e15
    f[b*"drift_direction"]=h.drift_direction;f[b*"time_origin_z_mm"]=h.time_origin_z_mm
    f[b*"z_edges_mm"]=h.z_edges_mm;f[b*"time_edges_ms"]=h.time_edges_ms
    save_columns(create_group(f,b*"history"),h.counts)
    f[b*"history/charge_fC"]=h.counts.ions.*(p.charge_per_ion_C*1e15)
    f[b*"xy_ions"]=q.xy_ions;f[b*"slice_ions"]=q.slice_ions
    f[b*"xy_charge_C"]=q.xy_ions.*p.charge_per_ion_C
    f[b*"display/axes/x_edges_mm"]=h.x_edges_mm
    f[b*"display/axes/y_edges_mm"]=h.y_edges_mm
    f[b*"display/axes/z_edges_ms"]=h.time_edges_ms
    nx,ny,nt=length(h.x_edges_mm)-1,length(h.y_edges_mm)-1,length(h.time_edges_ms)-1
    xt=zeros(Int64,nx,nt);yt=zeros(Int64,ny,nt)
    occupied=EventAnalysis.DataFrame(ix=Int64[],iy=Int64[],iz=Int64[],ions=Int64[],hit_count=Int64[])
    for r in eachrow(h.counts)
        xt[r.pad_ix,r.slice_index]+=r.ions;yt[r.pad_iy,r.slice_index]+=r.ions
        push!(occupied,(r.pad_ix,r.pad_iy,r.slice_index,r.ions,r.ions))
    end
    for (plane,counts) in (("xy",q.xy_ions),("xz",xt),("yz",yt))
        f[b*"display/projections/"*plane*"_ions"]=counts
        f[b*"display/projections/"*plane*"_charge_fC"]=counts.*(p.charge_per_ion_C*1e15)
        f[b*"display/projections/"*plane*"_hit_count"]=counts
    end
    save_columns(create_group(f,b*"display/occupied"),occupied)
    f[b*"display/occupied/charge_fC"]=occupied.ions.*(p.charge_per_ion_C*1e15)
end
