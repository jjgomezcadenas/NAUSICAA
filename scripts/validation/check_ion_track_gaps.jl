# Check broken itracks against the graph containing every edge within 6σₜ.
include("../production/run_framed_reconstruction.jl")
function check_gaps(scanpath)
 BLAS.set_num_threads(1)
 ids=h5open(joinpath(scanpath,"knn_scan_summary.h5"),"r") do f;read(f["ions/k_10/broken_event_numbers"]);end
 rc=TOML.parsefile(joinpath(ROOT,"config/readout/bb0nu_tracks_ion4mm_20ms_500_cut1keV.toml"))
 ec=TOML.parsefile(joinpath(ROOT,rc["electronics_card"]));dc=TOML.parsefile(joinpath(ROOT,ec["detector_card"]))
 prepared=prepare_electronics_run(ec);workspace=ER.CSAWorkspace();c=rc["ions"]
 settings=IonReconstructionSettings(c["baseline_ms"],c["average_ms"],c["nonlinear_below_V"],c["max_iterations"],c["convergence_ions"])
 positions,frames=h5open(joinpath(ROOT,rc["selection_controls"]),"r") do f
  e=f["events"];allids=read(e["event_number"]);rows=[findfirst(==(id),allids) for id in ids]
  ps=[EventPlacement(id,read(e["ion_drift_cm"])[j],read(e["electron_drift_cm"])[j]) for (id,j) in zip(ids,rows)]
  fs=[EventFrame(mc_vertex_frame,ntuple(k->read(e["frame_$("xyz"[k])_min"])[j],3),ntuple(k->read(e["frame_$("xyz"[k])_max"])[j],3)) for j in rows]
  ps,fs
 end
 rows=NamedTuple[]
 h5open(joinpath(ROOT,dc["hits"]),"r") do input
  index=Dict(id=>j for (j,id) in enumerate(read(input["events/event_id"])))
  for (ordinal,id) in enumerate(ids)
   p=positions[ordinal];frame=frames[ordinal];h,info=read_mc_event(input,id;event_row=index[id]);h=hits_in_frame(center_hits(h,info),frame)
   ds=event_settings(dc,p);rng(tag)=Xoshiro(seedfor(dc["seed"],id,tag));erng(tag)=Xoshiro(seedfor(ec["seed"],id,tag))
   pairs=sample_fano_pairs(rng("pairs"),h,Float32(dc["W_eV"]),Float32(dc["fano_factor"]))
   ions=sim_ion_detector_response(rng("ions"),h,pairs,p,frame,ds;max_vector_bytes=dc["max_vector_bytes"])
   s=electronics_settings(prepared,ions.ions.time)
   adc,_=ER.sim_ion_electronics_prepared(erng("csa_noise"),erng("adc_noise"),ions,s.ion,workspace)
   cal=IonReconstructionCalibration(s.ion.adc,s.ion.csa.current_pulse,s.ion.csa.linear_response,s.ion.csa.nonlinear_response,s.ion.csa.noise)
   plan=SR.prepare_ion_reconstruction(adc.time,ions.ions.time,cal,settings);ri=SR.reconstruct_ions(adc,plan)
   hits=reconstruct_hits(ri,p,ds;keV_per_count=Float64(ds.W_eV)/1000,minimum_energy_keV=1.)
   selected=findall(ri.valid .& isfinite.(ri.ions) .& (Float64.(ri.ions).*Float64(ds.W_eV)./1000 .>1.))
   keep=ions.ions.counts[selected].>0
   df=EventAnalysis.DataFrame(x=hits.x_mm[keep],y=hits.y_mm[keep],z=hits.z_mm[keep],energy=hits.energy_keV[keep]./1000f0)
   distance=6f0*ds.ion_transport.sigma_t_mm
   knn=make_tracks(df;method="kNN",k=30,max_distance_mm=distance,energy_threshold_kev=0f0)
   radius=make_tracks(df;method="KDT",max_distance_mm=distance,energy_threshold_kev=0f0)
   smallest=Inf
   for i in 1:length(radius),j in i+1:length(radius)
    a=radius[i].voxels;b=radius[j].voxels
    for m in 1:size(a,1),n in 1:size(b,1)
     d=sqrt(Float64(a.x[m]-b.x[n])^2+Float64(a.y[m]-b.y[n])^2+Float64(a.z[m]-b.z[n])^2)
     smallest=min(smallest,d)
    end
   end
   push!(rows,(event_number=id,max_distance_mm=Float64(distance),minimum_component_gap_mm=smallest,knn_components=length(knn),radius_components=length(radius)))
   println(last(rows));flush(stdout)
  end
 end
 target=joinpath(scanpath,"ion_gap_check.h5");isfile(target)&&error("Preserve existing check")
 budget=1024^2;GuardedFiles.space_check(scanpath,budget,1024^3);tmp=target*".tmp"
 try
  h5open(tmp,"w") do f
   f["schema"]="itaca.ion_gap_check.v1";columns(create_group(f,"events"),rows)
   f["card"]=read(joinpath(ROOT,"config/readout/bb0nu_tracks_ion4mm_20ms_500_cut1keV.toml"),String)
   flush(f);GuardedFiles.check_size(tmp,budget)
  end
  mv(tmp,target)
 catch;isfile(tmp)&&rm(tmp);rethrow();end
end
check_gaps(abspath(ARGS[1]))
