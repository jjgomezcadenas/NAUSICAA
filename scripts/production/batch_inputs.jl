# Input files of a run card and the selection of one input batch.

"""Path of a run-card input: `\$NAUSICAA_DATA` is this machine's data folder
(environment variable), `~` the home folder; other paths are relative to root."""
function input_path(root::AbstractString,s::AbstractString)::String
    if occursin("\$NAUSICAA_DATA",s)
        haskey(ENV,"NAUSICAA_DATA") || error("Set NAUSICAA_DATA to the folder holding the input files of $s")
        s=replace(s,"\$NAUSICAA_DATA"=>ENV["NAUSICAA_DATA"])
    end
    joinpath(root,expanduser(s))
end

"""Tag of an input batch: b00000, b00001, ..."""
batch_tag(batch::Integer)::String=batch>=0 ? "b"*lpad(batch,5,'0') : throw(ArgumentError("Batch numbers start at 0"))

"""Restrict a resolved run to one input batch.

`{batch}` in `hits` and `particles` becomes the batch tag, and the results go
under `<run_label>/batches/<tag>/`. Event numbers are unique across the
batches of a sample, so the per-event seeds are unchanged.
"""
function select_batch!(run::AbstractDict,batch::Integer)
    tag=batch_tag(batch)
    for k in ("hits","particles")
        occursin("{batch}",run[k]) || error("Run card key $k has no {batch} to select a batch")
        run[k]=replace(run[k],"{batch}"=>tag)
    end
    run["run_label"]=joinpath(run["run_label"],"batches",tag)
    run
end

"""Read `--batch=N` from the arguments; return the batch (or nothing) and the other arguments."""
function batch_argument(args::Vector{String})
    k=findall(startswith("--batch="),args)
    length(k)<=1 || error("Give --batch once")
    isempty(k) ? (nothing,args) : (parse(Int,split(args[k[1]],"=")[2]),args[setdiff(eachindex(args),k)])
end
