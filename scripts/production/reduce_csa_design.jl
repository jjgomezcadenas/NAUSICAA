# Gather saved CSA design and pilot results for the technical note.
# Check matched inputs and copy the saved summaries for presentation.

using HDF5,SHA,TOML,Statistics,Test
const ROOT=normpath(joinpath(@__DIR__,"../.."))
sha(p)=open(io->bytes2hex(sha256(io)),p)
attrs(g)=Dict(k=>read(attributes(g)[k]) for k in keys(attributes(g)))
function run(card=joinpath(ROOT,"config/specs/runs/csa_design_note.toml"))
 c=TOML.parsefile(card);@assert c["schema"]=="nausicaa.csa_design_note.v1"
 sources=String[card,joinpath(ROOT,c["design_card"]),joinpath(ROOT,c["response_card"]),joinpath(ROOT,c["reconstruction_card"]),joinpath(ROOT,c["parameters"]),joinpath(ROOT,c["retention"]),joinpath(ROOT,c["calibration"]),@__FILE__]
 studies=Dict{Int,String}()
 for p in c["pilots"]
  index=joinpath(ROOT,p["index"]);push!(sources,index)
  @assert isfile(index) "Missing completed pilot index: $index"
  studies[p["pitch_mm"]]=joinpath(ROOT,TOML.parsefile(index)["results"])
  append!(sources,[joinpath(studies[p["pitch_mm"]],f) for f in ("summary.h5","manifest.toml","independent_validation.json")])
 end
 outdir=joinpath(ROOT,c["results_directory"]);mkpath(outdir);dest=joinpath(outdir,"summary.h5")
 references=Dict{String,Tuple{Vector{Int64},Vector{Int64}}}();checks=0
 h5open(dest*".tmp","w") do out
  attributes(out)["schema"]="nausicaa.csa_design_summary.v1"
  for (name,key) in (("design","design_card"),("response","response_card"),("reconstruction","reconstruction_card"))
   out[name*"_card_toml"]=read(joinpath(ROOT,c[key]),String)
  end
  # Preserve the simulation package values; derive supply power in Julia.
  params=read(joinpath(ROOT,c["parameters"]),String)
  out["parameters_json"]=params
  idd=parse(Float64,match(r"\"I\(Vdd\)\"\s*:\s*([-+0-9.eE]+)",params).captures[1])
  supply=parse(Float64,match(r"\"V\(vdda\)\"\s*:\s*([-+0-9.eE]+)",params).captures[1])
  out["supply_current_A"]=abs(idd);out["supply_power_W"]=abs(idd)*supply
  retention=create_group(out,"retention")
  h5open(joinpath(ROOT,c["retention"]),"r") do h
   for name in keys(h);HDF5.copy_object(h[name],retention,name);end
  end
  pilots=create_group(out,"pilots")
  @testset "CSA note matched pilot inputs and summaries" begin
   for pitch in sort(collect(keys(studies)))
    study=studies[pitch];pgroup=create_group(pilots,string(pitch));attributes(pgroup)["study_directory"]=relpath(study,ROOT)
    h5open(joinpath(study,"summary.h5"),"r") do h
     for name in read(h["groups"])
      g=h[name];a=attrs(g);sample=String(a["sample"]);ids=Int64.(read(g["event_id"]));ions=Int64.(read(g["initial_ions"]))
      @test length(ids)==100;checks+=1
      @test a["pitch_mm"]==pitch;checks+=1
      @test a["frame_s"]≈pitch*.01;checks+=1
      if haskey(references,sample)
       @test references[sample][1]==ids;checks+=1
       @test references[sample][2]==ions;checks+=1
      else
       references[sample]=(ids,ions)
      end
      HDF5.copy_object(g,pgroup,sample);destg=pgroup[sample]
      aa=attributes(destg)
      for (key,input) in (("noiseless","noiseless_charge_bias_fraction"),("noisy","noisy_charge_bias_fraction"))
       v=read(g[input]);aa[key*"_bias_std_fraction"]=std(v);aa[key*"_bias_median_fraction"]=median(v)
      end
      @test sum(read(g["positive_cells"]))==a["positive_cells"];checks+=1
      @test sum(ions)==a["input_ion_total"];checks+=1
     end
     HDF5.copy_object(h["validation"],pgroup,"validation")
     HDF5.copy_object(h["calibration"],pgroup,"calibration")
    end
   end
  end
  attributes(out)["validation_checks"]=checks
  out["source_path"]=relpath.(sources,ROOT);out["source_sha256"]=sha.(sources)
  attributes(out)["complete"]=true
 end
 mv(dest*".tmp",dest;force=true)
 open(joinpath(outdir,"manifest.toml"),"w") do io
  TOML.print(io,Dict("schema"=>"nausicaa.csa_design_summary.v1","summary"=>relpath(dest,ROOT),"summary_sha256"=>sha(dest),"card"=>relpath(card,ROOT),"source_paths"=>relpath.(sources,ROOT),"source_sha256"=>sha.(sources),"validation_checks"=>checks))
 end
 println(dest)
end
run(isempty(ARGS) ? joinpath(ROOT,"config/specs/runs/csa_design_note.toml") : abspath(ARGS[1]))
