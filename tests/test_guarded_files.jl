using Test,SHA
include("../src/ion_signal/guarded_files.jl")
const G=GuardedFiles
@testset "Bounded copies and atomic publication" begin
 mktempdir() do dir
    src=joinpath(dir,"source");dst=joinpath(dir,"target")
    bytes=UInt8.(mod.(0:199999,251));write(src,bytes)
    digest=bytes2hex(sha256(bytes))
    @test G.checked_copy(src,dst;expected_sha256=digest,block_bytes=1024)==dst
    @test read(dst)==bytes
    @test_throws ArgumentError G.checked_copy(src,dst)
    before=sort(readdir(dir))
    @test_throws ErrorException G.checked_copy(src,dst;force=true,expected_sha256="bad")
    @test read(dst)==bytes && sort(readdir(dir))==before
    @test_throws ErrorException G.checked_copy(src,dst;force=true,reserve_bytes=typemax(Int64))
    @test read(dst)==bytes && sort(readdir(dir))==before
    @test_throws ArgumentError G.checked_copy(src,src;force=true)
    @test_throws ArgumentError G.checked_copy(src,dst;force=true,block_bytes=2^32)
    @test sort(readdir(dir))==before
    # Simulated failed metadata writer: cross the limit then check immediately.
    temp=Ref("")
    @test_throws ErrorException G.staged_copy(src,dst;force=true) do p,check!
        temp[]=p
        open(p,"a") do f;write(f,zeros(UInt8,20001));end
        check!()
    end
    @test !ispath(temp[]) && read(dst)==bytes && sort(readdir(dir))==before
    @test_throws ErrorException G.staged_copy(src,dst;force=true) do p,check!
        error("injected validation failure")
    end
    @test sort(readdir(dir))==before && read(dst)==bytes
    @test_throws InterruptException G.staged_copy(src,dst;force=true) do p,check!
        throw(InterruptException())
    end
    @test sort(readdir(dir))==before
    G.staged_copy(src,dst;force=true) do p,check!
        open(p,"a") do f;write(f,zeros(UInt8,20000));end
        @test check!()==220000
    end
    @test filesize(dst)==220000 && read(src)==bytes
    empty=joinpath(dir,"empty");write(empty,"")
    G.checked_copy(empty,dst;force=true)
    @test filesize(dst)==0
    # Exact 10% check for a sparse oversize file without allocating its contents.
    open(dst,"w") do f;truncate(f,220001);end
    @test_throws ErrorException G.check_size(dst,220000)
 end
end
