"""Deterministic geometries shared by the pinned Petit reference and adapted tests."""
function reconstruction_fixtures(T=Float32)
    table(x,y,z; energy=ones(length(x))) = DataFrame(event_id=fill(Int64(42),length(x)),
        x=T.(x), y=T.(y), z=T.(z), energy=T.(energy))
    xs = [0,0.7,1.8,3.1,4.6,6.3]
    line = table(xs,zeros(6),zeros(6); energy=[4,1,1,1,1,5])
    curve = table([0,0.6,1.5,2.5,3.2,3.7,3.8], [0,0.1,0.4,1,1.9,3,4.3], zeros(7))
    branch = table([0,1,2,3,4,2,2], [0,0,0,0,0,1,2], zeros(7))
    dense = table([i+0.037j+0.009k for i in 0:3 for j in 0:2 for k in 0:1],
                  [j+0.051i for i in 0:3 for j in 0:2 for k in 0:1],
                  [k+0.027i for i in 0:3 for j in 0:2 for k in 0:1];
                  energy=[1+0.13i+0.07j+0.11k for i in 0:3 for j in 0:2 for k in 0:1])
    lattice = table([i for i in 0:4 for j in 0:3], [j for i in 0:4 for j in 0:3], zeros(20))
    return Dict(
        "line" => (data=line, cutoff=T(1.8), k=2, tied=false),
        "curve" => (data=curve, cutoff=T(1.5), k=2, tied=false),
        "branch" => (data=branch, cutoff=T(1.01), k=4, tied=true),
        "dense" => (data=dense, cutoff=T(2.4), k=10, tied=false),
        "lattice" => (data=lattice, cutoff=T(1.5), k=4, tied=true),
        "separated" => (data=table([0,0.5,1.3,10,10.8],zeros(5),zeros(5)), cutoff=T(1.4), k=2,tied=false))
end

"""Encode graph edges as parent-table row pairs for reference comparison."""
function parent_edges(tracks)
    result = Tuple{Int,Int}[]
    for t in tracks, e in edges(t.graph)
        rows = only(t.components)
        push!(result,minmax(rows[src(e)],rows[dst(e)]))
    end
    sort!(result)
    return isempty(result) ? zeros(Int,2,0) : hcat(([a,b] for (a,b) in result)...)
end
