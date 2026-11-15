using Test,HDF5
@testset "Compact hit output" begin
 mktempdir() do dir
  settings(level=4,budget=8*1024^2)=HitOutputSettings(1f0,1f0,level,4,budget)
  meta=HitFileMetadata("test_ion4mm_el4mm","bb0nu",Xe_85_He_15,Dict("test"=>"seed=42"),Dict("source"=>"test"))
  frame=EventFrame(mc_vertex_frame,(-10f0,-10f0,-10f0),(10f0,10f0,10f0))
  h=SavedHits(Float32[1,2,3,4,5],zeros(Float32,5),ones(Float32,5),Float32[2,3,4,5,6])
  empty=SavedHits(Float32[],Float32[],Float32[],Float32[])
  e=EventEnergy(-3f0,12f0,-.066f0,4f0,0,0)
  event(id,ih=h,eh=empty)=HitEvent(HitEventInfo(id,70f0,80f0,frame,1f0,8f0),e,ih,eh)
  for level in (0,4)
   path=joinpath(dir,"hits$level.h5");w=open_hit_output(path,meta,settings(level))
   write_hit_event(w,event(7));write_hit_event(w,event(91,empty,h))
   @test_throws ArgumentError write_hit_event(w,event(7))
   @test finish_hit_output(w)==path
   @test !isfile(path*".tmp")
   @test_throws ArgumentError open_hit_output(path,meta,settings())
   r=open_hit_input(path)
   for id in Int64[91,7]
    actual=read_hit_event(r,id);expected=event(id,id==7 ? h : empty,id==7 ? empty : h)
    @test actual.info.event_number==id
    @test actual.info.frame.lower_mm==frame.lower_mm
    @test actual.info.ion_sigma_t_mm==1f0
    for key in fieldnames(EventEnergy);@test isequal(getfield(actual.energy,key),getfield(e,key));end
    for arm in (:itrack,:etrack),key in fieldnames(SavedHits)
     @test getfield(getfield(actual,arm),key)==getfield(getfield(expected,arm),key)
    end
   end
   @test_throws KeyError read_hit_event(r,Int64(12));close_hit_input(r)
   h5open(path,"r") do f
    @test Set(keys(f))==Set(["schema","events","itrack","etrack","configuration"])
    @test read(f["itrack/offsets"])==[0,5,5]
    @test read(f["etrack/offsets"])==[0,0,5]
   end
  end
  # Event location (schema v2): vertex in detector mm and the fixed frame
  cube=EventFrame(mc_vertex_frame,(-90f0,-100f0,-110f0),(110f0,100f0,90f0));where=HitEventLocation((12f0,-34f0,560f0),cube)
  path=joinpath(dir,"located.h5");w=open_hit_output(path,meta,settings();located=true)
  @test_throws ArgumentError write_hit_event(w,event(7))                       # a located file needs the location
  write_hit_event(w,event(7);location=where);write_hit_event(w,event(8,empty,h);location=HitEventLocation((0f0,0f0,125f0),cube))
  @test finish_hit_output(w)==path
  r=open_hit_input(path)
  @test read_hit_location(r,Int64(7)).vertex_mm==(12f0,-34f0,560f0)
  @test read_hit_location(r,Int64(8)).vertex_mm==(0f0,0f0,125f0)
  @test read_hit_location(r,Int64(7)).fixed_frame.lower_mm==cube.lower_mm && read_hit_location(r,Int64(7)).fixed_frame.upper_mm==cube.upper_mm
  @test read_hit_event(r,Int64(7)).itrack.energy_keV==h.energy_keV            # the version 1 content is unchanged
  @test_throws KeyError read_hit_location(r,Int64(12));close_hit_input(r)
  h5open(path,"r") do f
   @test read(f["schema"],String)=="itaca.reconstructed_hits.v2"
   @test read(f["events/vertex_z_mm"])==Float32[560,125]
   @test read(f["events/fixed_frame_x_min_mm"])==Float32[-90,-90]
  end
  plain=joinpath(dir,"plain.h5");w=open_hit_output(plain,meta,settings())
  @test_throws ArgumentError write_hit_event(w,event(7);location=where)        # a version 1 file takes no location
  write_hit_event(w,event(7));finish_hit_output(w)
  r=open_hit_input(plain);@test_throws ArgumentError read_hit_location(r,Int64(7));close_hit_input(r)
  h5open(plain,"r") do f;@test read(f["schema"],String)=="itaca.reconstructed_hits.v1"&&!haskey(f["events"],"vertex_z_mm");end
  @test_throws ArgumentError HitOutputSettings(.5f0,1f0,4,4,1024)
  path=joinpath(dir,"failed.h5");w=open_hit_output(path,meta,settings())
  boundary=SavedHits(Float32[0],Float32[0],Float32[0],Float32[1])
  @test_throws ArgumentError write_hit_event(w,event(1,boundary))
  abort_hit_output(w);@test !isfile(path)&&!isfile(path*".tmp")
  @test_throws ErrorException open_hit_output(path,meta,settings(4,1))
  @test !isfile(path*".tmp")
  w=open_hit_output(path,meta,settings());write_hit_event(w,event(2));w.settings=settings(4,filesize(w.temporary_path))
  big=SavedHits(rand(Float32,10000),rand(Float32,10000),rand(Float32,10000),2f0.+rand(Float32,10000))
  @test_throws ErrorException write_hit_event(w,event(3,big));abort_hit_output(w)
  @test !isfile(path*".tmp")
 end
end
