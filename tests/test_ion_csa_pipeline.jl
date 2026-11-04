using Test, HDF5, TOML, LinearAlgebra
root=normpath(joinpath(@__DIR__,".."))
include(joinpath(root,"src/ion_signal/ion_csa_pipeline.jl"))
const C=IonCSAProduction

function fixture(dir)
    base=joinpath(dir,"bb0nu","realisation")
    source=joinpath(base,"detector","ions","detector_id","ions.h5")
    parent=joinpath(base,"response","ions","detector_id","corrected_id","ion_response.h5")
    mkpath(dirname(source));mkpath(dirname(parent))
    placement=(reference_mm=Float32[0,0,0],phase_mm=Float32[.1,.2],ion_drift_cm=Float32(20),electron_drift_cm=Float32(130),drift_direction=Int8(1),bucket_width_mm=Float32(3))
    responsevalues=TOML.parsefile(joinpath(root,"config/specs/response/ion_csa/csa_v1.toml"))
    calibrationdir=joinpath(dir,"calibration");mkpath(calibrationdir)
    deck=joinpath(dir,"source.cir");write(deck,"fixture source deck\n")
    model=joinpath(dir,"corrected.inc");write(model,"fixture transistor model\n")
    calibration_card=C.write_toml(joinpath(dir,"calibration.toml"),Dict("source_deck"=>deck,"model"=>model))
    inputs=joinpath(calibrationdir,"inputs.h5")
    h5open(inputs,"w") do h;h["fixture"]=true;end
    design=TOML.parsefile(joinpath(root,responsevalues["design_card"]))
    h5open(joinpath(calibrationdir,"analysis.h5"),"w") do h
        attributes(h)["config_sha256"]=C.sha(calibration_card)
        g=create_group(h,responsevalues["variant"]);a=attributes(g)
        a["source_sha256"]=C.sha(deck);a["model_sha256"]=C.sha(model)
        a["inputs_sha256"]=C.sha(inputs)
        for (key,value) in (("length_m",design["feedback"]["length_m"]),
                            ("width_m",design["feedback"]["width_m"]),
                            ("gate_V",design["feedback"]["gate_V"]),
                            ("m4_gate_V",design["bias"]["m4_gate_V"]))
            a[key]=value
        end
    end
    responsevalues["calibration_directory"]=calibrationdir
    responsevalues["calibration_card"]=calibration_card
    responsepath=C.write_toml(joinpath(dir,"csa.toml"),responsevalues)
    responsecard=read(responsepath,String)
    reconstruction=read(joinpath(root,"config/specs/response/ion_csa/reconstruction_v1.toml"),String)
    h5open(source,"w") do h
        attributes(h)["schema"]="nausicaa.ion_pad_file.pilot.v1"
        m=create_group(h,"metadata");m["complete"]=true;m["sample"]="bb0nu";m["event_count"]=1
        m["branch_card_toml"]="pad_pitch_mm = 3\nbucket_width_ms = 30\n"
        g=create_group(h,"events/5");g["pixel_ix"]=Int64[-50,1000];g["pixel_iy"]=Int64[0,700]
        g["bucket"]=Int64[-1,0];g["count"]=Int64[2,3];g["initial_pairs"]=Int64(5);g["pair_checksum"]="shared_primary_5"
        for (key,value) in pairs(placement);g[string(key)]=value;end
    end
    h5open(parent,"w") do h
        attributes(h)["schema"]="nausicaa.ion_csa_response.pilot.v2"
        m=create_group(h,"metadata")
        values=(complete=true,signal_unit="ion_charges",time_unit="s",matrix_axes_julia="frame,pad",matrix_axes_hdf5="pad,frame",
            reconstruction_method="three_frame_inverse",source_path=source,source_sha256=C.sha(source),sample="bb0nu",source_event_count=1,
            selected_event_ids=Int64[5],event_count=1,frame_s=.03,pitch_mm=3,input_ion_total=5,response_id="corrected_id",
            response_card_toml=responsecard,reconstruction_card_toml=reconstruction,parent_code_hash="fixture_backend_code",
            physical_model_id=C.S.content_hash(join((C.S.canonical(responsevalues),C.sha(joinpath(calibrationdir,"analysis.h5")),"fixture_backend_code"),"|")))
        for (k,v) in pairs(values);m[string(k)]=v;end
        g=create_group(h,"events/5");g["event_id"]=Int64(5);g["pixel_ix"]=Int64[-50,1000];g["pixel_iy"]=Int64[0,700]
        g["first_bucket"]=Int64(-2);g["edge_time_s"]=collect(-2:2).*.03
        n=Int64[0 0;2 0;0 3;0 0];z=Float32[-1 -2;1 0;0 1;-2 -1]
        g["counts"]=n;g["noiseless_estimated_ions"]=Float32.(n);g["noise_ions"]=z;g["recovered_ions"]=Float32.(n).+z
        g["noise_sigma_ions"]=ones(4);g["model_used"]=Int8[0,1];g["noise_seed"]=UInt64[8,9]
        g["initial_pairs"]=Int64(5);g["pair_checksum"]="shared_primary_5"
        for (k,v) in pairs(placement);g[string(k)]=v;end
        cov=create_group(h,"covariance_by_frames/4");cov["covariance_e2"]=Matrix{Float64}(I,4,4)
    end
    index=joinpath(dir,"index.h5")
    h5open(index,"w") do h;h["bb0nu/ions/3/path"]=source;end
    study=joinpath(dir,"corrected_study");mkpath(study)
    h5open(joinpath(study,"summary.h5"),"w") do h
        g=create_group(h,"bb0nu_p3");attributes(g)["response_path"]=parent;attributes(g)["response_sha256"]=C.sha(parent)
    end
    card=joinpath(dir,"run.toml")
    values=Dict("schema"=>"nausicaa.ion_csa_pipeline_run.v1","response_card"=>responsepath,
        "reconstruction_card"=>"config/specs/response/ion_csa/reconstruction_v1.toml","ion_index"=>index,"samples"=>["bb0nu"],
        "events_per_sample"=>1,"selection"=>"evenly_spaced_sorted_event_ids","channel_selection"=>C.CHANNEL_POLICY,
        "absent_channel_policy"=>C.ABSENT_POLICY,"negative_estimate_clipping"=>false,"measured_zero_suppression"=>false,
        "results_directory"=>joinpath(dir,"published_runs"),"branches"=>[Dict("pitch_mm"=>3,"corrected_study"=>study)])
    C.write_toml(card,values)
    (;source,parent,index,study,card,values,deck,model,calibration_card,inputs)
end

@testset "Compact reconstructed matrices and production nodes" begin
 mktempdir(joinpath(root,"notes/physics/work")) do dir
    f=fixture(dir)
    e=C.read_event(f.parent,5)
    @test size(e.reconstructed_ions)==(4,2)
    @test e.pixel_ix==[-50,1000] # distant deposit adds one pixel, no spatial rectangle
    @test e.first_bucket==-2 && e.edge_time_s≈[-.06,-.03,0,.03,.06]
    @test e.reference_mm==Float32[0,0,0] && e.phase_mm==Float32[.1,.2] && e.bucket_width_mm==3
    @test e.input_ions[:,1]==[0,2,0,0]
    @test e.input_ions[:,2]==[0,0,3,0]
    @test e.reconstructed_ions[1,1]==-1 && e.reconstructed_ions[4,2]==-1
    @test e.reconstructed_ions==e.noiseless_ions.+e.noise_ions
    e.input_ions[2,1]=99
    @test C.read_event(f.parent,5).input_ions[2,1]==2 # owned data
    @test_throws ArgumentError C.read_event(f.parent,99)
    @test C.read_covariance(f.parent,4)==Matrix{Float64}(I,4,4)
    @test_throws ArgumentError C.read_covariance(f.parent,3)
    ledger=C.validate_response(f.parent)
    @test ledger.input_ion_total==5 && ledger.event_count==1 && ledger.pads==2
    @test ledger.zero_input_buckets==6 && ledger.negative_estimates==4
    manifest=C.run(root,f.card);m=TOML.parsefile(manifest);output=only(m["responses"])["path"]
    @test m["complete"] && m["channel_selection"]==C.CHANNEL_POLICY
    @test isfile(output) && output!=f.parent
    @test C.validate_response(output).input_ion_total==5
    @test C.read_event(output,5).reconstructed_ions==C.read_event(f.parent,5).reconstructed_ions
    h5open(output,"r") do h
        @test C.schema(h)==C.PRODUCTION_SCHEMA
        @test C.str(h["metadata"],"absent_channel_policy")==C.ABSENT_POLICY
        @test !read(h["metadata/measured_zero_suppression"]) && !read(h["metadata/negative_estimate_clipping"])
    end
    hash=C.sha(output)
    @test C.run(root,f.card)==manifest && C.sha(output)==hash
    @test isfile(joinpath(dirname(dirname(dirname(dirname(f.source)))),"runs",m["graph_id"],"ion_manifest.toml"))
    # Comments and key order do not change the parsed run identity.
    original=C.settings(root,f.card).identity
    open(f.card,"a") do io;write(io,"\n# presentation-only comment\n");end
    @test C.settings(root,f.card).identity==original
    # A changed source deck or transistor model invalidates the cached calibration.
    for resource in (f.deck,f.model,f.calibration_card)
        body=read(resource,String)
        write(resource,body*"# changed\n")
        @test_throws ArgumentError C.settings(root,f.card)
        write(resource,body)
    end
    backup=joinpath(dir,"inputs_backup.h5");cp(f.inputs,backup)
    h5open(f.inputs,"r+") do h;h["changed"]=true;end
    @test_throws ArgumentError C.settings(root,f.card)
    cp(backup,f.inputs;force=true)
    # A mismatched acquisition policy fails before publishing.
    bad=copy(f.values);bad["channel_selection"]="rectangle";C.write_toml(f.card,bad)
    @test_throws ArgumentError C.settings(root,f.card)
    C.write_toml(f.card,f.values)
    # Saved-study provenance must match the response hash recorded at reduction.
    h5open(joinpath(f.study,"summary.h5"),"r+") do h
        g=h["bb0nu_p3"];HDF5.delete_attribute(g,"response_sha256");attributes(g)["response_sha256"]="changed"
    end
    @test_throws ArgumentError C.study_responses(C.settings(root,f.card),f.study)
    # Calibration contents change the execution identity with unchanged cards.
    before=C.settings(root,f.card).identity
    calibration=joinpath(dir,"calibration/analysis.h5")
    h5open(calibration,"r+") do h;h["changed"]=true;end
    @test C.settings(root,f.card).identity!=before
    @test_throws ArgumentError C.check_selection(C.settings(root,f.card),f.parent,ledger,3)
    # Changes to individual sparse cell placement are detected with unchanged totals.
    broken=joinpath(dir,"broken.h5");cp(f.parent,broken)
    h5open(broken,"r+") do h;write(h["events/5/counts"],Int64[0 0;0 0;2 3;0 0]);end
    @test_throws ArgumentError C.validate_response(broken)
    # Changed source contents fail the recorded source checksum.
    h5open(f.source,"r+") do h;h["changed"]=true;end
    @test_throws ArgumentError C.validate_response(f.parent)
 end
end
