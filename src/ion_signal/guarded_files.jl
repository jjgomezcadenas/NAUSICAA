"""Bounded, verified file copies and temporary publication transactions."""
module GuardedFiles
using SHA
export checked_copy, staged_copy, check_size
const BLOCK_BYTES=8*1024^2
const RESERVE_BYTES=1024^3
sha(path)=open(io->bytes2hex(sha256(io)),path)

"""Require a regular temporary file no larger than its explicit byte budget."""
function check_size(path,maximum_bytes::Integer)
    isfile(path) && !islink(path) || throw(ArgumentError("Expected a regular temporary file: $path"))
    actual=filesize(path)
    actual<=maximum_bytes || error("Temporary file exceeds size budget: $actual > $maximum_bytes bytes: $path")
    actual
end

function space_check(directory,remaining,reserve)
    available=diskstat(directory).available
    Int128(available)>=Int128(remaining)+reserve ||
        error("Insufficient free space: need $remaining bytes plus $reserve reserve; available $available")
end

"""Copy exactly the source's initial length with bounded writes and SHA-256 validation.

The caller owns the empty destination. Eight-MiB writes avoid the runtime's
32-bit sendfile transfer counter. Reject shrinking/growing sources, short writes,
content mismatches and insufficient space. Return the verified source digest.
"""
function copy_bytes(source,destination;expected_sha256=nothing,block_bytes=BLOCK_BYTES,
                    reserve_bytes=RESERVE_BYTES)
    0<block_bytes<=BLOCK_BYTES || throw(ArgumentError("Copy block must be in 1:$BLOCK_BYTES bytes"))
    reserve_bytes>=0 || throw(ArgumentError("Negative free-space reserve"))
    before=stat(source);n=before.size
    space_check(dirname(destination),n,reserve_bytes)
    expected=expected_sha256===nothing ? sha(source) : expected_sha256
    buffer=Vector{UInt8}(undef,min(block_bytes,max(1,n)))
    open(source,"r") do input
      open(destination,"w") do output
        copied=0
        while copied<n
            space_check(dirname(destination),n-copied,reserve_bytes)
            wanted=min(length(buffer),n-copied)
            got=readbytes!(input,buffer,wanted)
            got==wanted || error("Source ended before its recorded size")
            write(output,view(buffer,1:got))==got || error("Short file write")
            flush(output)
            copied+=got
            check_size(destination,n)==copied || error("Copy made unexpected progress")
        end
        eof(input) || error("Source grew during copy")
      end
    end
    after=stat(source)
    (before.size,before.mtime,before.inode)==(after.size,after.mtime,after.inode) ||
        error("Source changed during copy")
    check_size(destination,n)==n || error("Copy length differs from source")
    sha(destination)==expected || error("Copy SHA-256 differs from source")
    expected
end

"""
    staged_copy(update!, source, target; expected_sha256=nothing,
                growth_percent=10, reserve_bytes=1GiB, force=false)

Create a unique temporary file beside `target`, copy and verify the immutable
source, then call `update!(temporary, check!)` for metadata edits and validation.
The callback calls `check!()` after each flushed write. The temporary's maximum
size is the original size plus `growth_percent` percent (10% by default); check
again before publishing by a same-directory rename. Check free disk space before
copying and during each block. The default leaves a 1 GiB reserve.

On an exception, including cancellation, delete only this invocation's temporary
file. Existing targets remain intact until successful publication. Abrupt process
termination can leave a uniquely named temporary; it is never reused implicitly.
`force=false` rejects existing targets. Symbolic-link targets and self-copies are
rejected. This function protects copy/publication steps; newly generated datasets
need a byte budget derived from their output shape rather than their input size.
"""
function staged_copy(update!::Function,source::AbstractString,target::AbstractString;
        expected_sha256=nothing,growth_percent::Integer=10,reserve_bytes::Integer=RESERVE_BYTES,
        block_bytes::Integer=BLOCK_BYTES,force::Bool=false)
    0<=growth_percent<=10 || throw(ArgumentError("Growth allowance must be between 0 and 10 percent"))
    isfile(source) || throw(ArgumentError("Copy source is absent: $source"))
    abspath(source)!=abspath(target) || throw(ArgumentError("Source and target coincide"))
    islink(target) && throw(ArgumentError("Target is a symbolic link"))
    ispath(target) && (!force || !isfile(target)) && throw(ArgumentError("Target already exists: $target"))
    ispath(target) && stat(source).device==stat(target).device && stat(source).inode==stat(target).inode &&
        throw(ArgumentError("Source and target refer to the same file"))
    n=filesize(source);maximum_bytes=n+fld(Int128(n)*growth_percent,100)
    mkpath(dirname(abspath(target)))
    space_check(dirname(abspath(target)),maximum_bytes,reserve_bytes)
    temporary,handle=mktemp(dirname(abspath(target));cleanup=false);close(handle)
    try
        copy_bytes(source,temporary;expected_sha256,block_bytes,reserve_bytes)
        check!()=check_size(temporary,maximum_bytes)
        update!(temporary,check!)
        check!()
        !force && ispath(target) && error("Target appeared during publication: $target")
        Base.rename(temporary,target)
    finally
        ispath(temporary) && rm(temporary;force=true)
    end
    target
end

"""Atomically copy a file with exact byte length and SHA-256; preserve old target on failure."""
function checked_copy(source,target;kwargs...)
    staged_copy((temporary,check!)->nothing,source,target;growth_percent=0,kwargs...)
end
end
