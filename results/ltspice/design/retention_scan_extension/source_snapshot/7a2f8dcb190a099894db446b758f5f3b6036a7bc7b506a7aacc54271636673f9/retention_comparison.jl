module RetentionComparison
using HDF5,TOML,SHA,Test
sha(p)=open(io->bytes2hex(sha256(io)),p)

"""Combine completed retention studies while preserving each source's calibration and provenance."""
function combine(root,card)
 c=TOML.parsefile(joinpath(root,card));outdir=joinpath(root,c["output"]);mkpath(outdir)
 write(joinpath(outdir,".gitignore"),"*.h5\n*.h5.tmp\n")
 @testset "Combined retention scan provenance and completeness" begin
 for filename in ("analysis.h5","presentation.h5")
  names=String[]
  h5open(joinpath(outdir,filename),"w") do out
   at=attributes(out);at["config_sha256"]=sha(joinpath(root,card));at["combiner_sha256"]=sha(@__FILE__)
   for source in c["scans"]
    file=joinpath(root,source,filename)
    h5open(file,"r") do h
     groupnames=filter(k->h[k] isa HDF5.Group,collect(keys(h)))
     for name in groupnames
      @test !(name in names);push!(names,name)
      HDF5.copy_object(h[name],out,name)
      attributes(out[name])["combined_source"]=relpath(file,root)
      attributes(out[name])["combined_source_sha256"]=sha(file)
      for frame in ("20ms","40ms")
       @test read(out[name*"/"*frame*"/ENC_e_rms"])[1]>0
      end
     end
    end
   end
   order=sortperm(names;by=n->read(attributes(out[n])["length_m"]))
   out["variant_names"]=names[order]
  end
 end
 end
end
end
