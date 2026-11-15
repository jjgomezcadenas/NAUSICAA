using Test,HDF5
include("../../scripts/production/merge_hit_batches.jl")
@testset "Chaining the hit files of batches" begin
 mktempdir() do dir
  settings=HitOutputSettings(1f0,1f0,4,4,8*1024^2)
  meta(label,source;code="c1")=HitFileMetadata(label,"bb0nu",Xe_85_He_15,Dict("run"=>"card"),Dict("source"=>source,"selection"=>"sel-"*source,"code"=>code))
  frame=EventFrame(mc_vertex_frame,(-10f0,-10f0,-10f0),(10f0,10f0,10f0));e=EventEnergy(100f0,200f0,2f0,3f0,0,0)
  hits(n,v)=SavedHits(fill(Float32(v),n),zeros(Float32,n),ones(Float32,n),fill(Float32(v+2),n))
  where(z)=HitEventLocation((1f0,2f0,Float32(z)),frame)
  function batch(name,label,source,events;code="c1",located=true)
   path=joinpath(dir,name);w=open_hit_output(path,meta(label,source;code),settings;located)
   for (id,n) in events;write_hit_event(w,HitEvent(HitEventInfo(id,70f0,60f0,frame,1f0,2f0),e,hits(n,id),hits(2n,id));location=located ? where(id) : nothing);end
   finish_hit_output(w)
  end
  a=batch("a.h5","run/batches/b00000","s0",[(Int64(0),3),(Int64(1),0)]);b=batch("b.h5","run/batches/b00001","s1",[(Int64(100),2)])
  target=joinpath(dir,"out/hits.h5")
  @test merge_hit_batches([a,b],target)==target
  @test !isfile(target*".tmp")
  h5open(target,"r") do f
   @test read(f["schema"],String)=="itaca.reconstructed_hits.v2"
   @test read(f["events/event_number"])==[0,1,100]
   @test read(f["events/vertex_z_mm"])==Float32[0,1,100]
   @test read(f["itrack/offsets"])==[0,3,3,5] && read(f["etrack/offsets"])==[0,6,6,10]
   @test read(f["itrack/x_mm"])==Float32[0,0,0,100,100] && read(f["etrack/energy_keV"])==Float32[2,2,2,2,2,2,102,102,102,102]
   @test read(f["batches/run_label"])==["run/batches/b00000","run/batches/b00001"] && read(f["batches/events"])==[2,1]
   @test read(f["batches/source_sha256"])==["s0","s1"]
   @test read(f["configuration/checksums/code"])=="c1"
  end
  r=open_hit_input(target)                                   # the chained file is a valid hit file
  @test read_hit_event(r,Int64(100)).itrack.x_mm==Float32[100,100] && read_hit_location(r,Int64(1)).vertex_mm==(1f0,2f0,1f0)
  close_hit_input(r)
  @test_throws ArgumentError merge_hit_batches([a,b],target)                                   # no overwrite
  @test_throws ErrorException merge_hit_batches([a,a],joinpath(dir,"dup.h5"))                   # duplicate event numbers
  c=batch("c.h5","run/batches/b00002","s2",[(Int64(200),1)];code="c2")
  @test_throws ErrorException merge_hit_batches([a,c],joinpath(dir,"code.h5"))                  # different code checksum
  d=batch("d.h5","run/batches/b00003","s3",[(Int64(300),1)];located=false)
  @test_throws ErrorException merge_hit_batches([a,d],joinpath(dir,"schema.h5"))                # different schema
  @test !isfile(joinpath(dir,"dup.h5")) && !isfile(joinpath(dir,"dup.h5.tmp"))
  @test_throws ArgumentError merge_hit_batches(String[],joinpath(dir,"none.h5"))
 end
end
