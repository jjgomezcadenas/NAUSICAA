using TOML
include(joinpath(@__DIR__,"../ion_signal/pad_compression.jl"))
include(joinpath(@__DIR__,"standard_events.jl"))
const ROOT=normpath(joinpath(@__DIR__,"../.."))
const WORK=abspath(ARGS[1])
const C=PadCompression.Context(ROOT,abspath(ARGS[2]))
println("READY\t",C.directory);flush(stdout)
for id in eachline(stdin)
 request=joinpath(WORK,id*".request.toml");response=joinpath(WORK,id*".response.toml")
 result=try
  r=TOML.parsefile(request);op=r["op"];path=""
  if op=="plan"
   path=PadCompression.plan(C)
  elseif op=="standard"
   path=PadHistoryNavigation.standard(C,r["sample"])
  elseif op=="aggregate"
   path=PadCompression.aggregate(C)
  elseif op in ("events","pads")
   path=PadCompression.navigation(C,op,r["sample"],get(r,"event_id",0),joinpath(WORK,id*".h5"))
  elseif op in ("history","prepare","reduce")
   fields=(C,r["sample"],r["event_id"],r["pixel_x"],r["pixel_y"])
   if op=="history";path=PadCompression.history(fields...)
   elseif op=="prepare";path=PadCompression.prepare(fields...;refine=get(r,"refine",false))
   else;path=PadCompression.reduce(fields...;refine=get(r,"refine",false));end
  else;error("Unknown operation $op");end
  Dict("status"=>"ok","path"=>path)
 catch e
  showerror(stderr,e,catch_backtrace());println(stderr)
  Dict("status"=>"error","message"=>sprint(showerror,e))
 end
 open(response,"w") do io;TOML.print(io,result);end
 println(id);flush(stdout)
end
