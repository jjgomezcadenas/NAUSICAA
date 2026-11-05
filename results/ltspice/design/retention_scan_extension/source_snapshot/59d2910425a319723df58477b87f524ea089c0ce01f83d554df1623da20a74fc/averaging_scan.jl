module AveragingScan
using LinearAlgebra
export waveform_gain, noise_rms, read_columns

function read_columns(path, wanted)
    text = endswith(path, ".gz") ? read(`gzip -dc $path`, String) : read(path, String)
    lines = split(strip(text), '\n'); names = split(first(lines), ',')
    indices = [findfirst(==(name), names) for name in wanted]
    @assert all(!isnothing, indices)
    columns = [Float64[] for _ in wanted]
    for line in lines[2:end]
        cells = split(line, ',')
        for (col, j) in zip(columns, indices)
            push!(col, parse(Float64, cells[j]))
        end
    end
    return columns
end

"""Exact integration of the saved piecewise-linear transient over centered windows."""
function waveform_gain(t, v, width, frame, alpha)
    @assert width > 0 && all(diff(t) .> 0)
    integral = vcat(0., cumsum(diff(t) .* (v[1:end-1] .+ v[2:end]) ./ 2))
    function primitive(x)
        @assert first(t) <= x <= last(t)
        j = min(searchsortedlast(t, x), length(t)-1); dx = x-t[j]
        integral[j] + v[j]*dx + (v[j+1]-v[j])/(t[j+1]-t[j])*dx^2/2
    end
    average(center) = (primitive(center+width/2)-primitive(center-width/2))/width
    -(average(frame)-alpha*average(0.))
end

"""One-sided stationary output PSD, logarithmic interpolation, Gauss integration.
Panel edges include spectral samples and half-period steps of the fastest filter
oscillation. The transfer includes both centered averaging and frame subtraction.
"""
function noise_rms(f, asd, width, frame, alpha; order=8)
    @assert all(diff(f) .> 0) && all(asd .> 0)
    b = [k/sqrt(4k^2-1) for k in 1:order-1]
    eig = eigen(SymTridiagonal(zeros(order), b))
    nodes = eig.values; weights = 2 .* eig.vectors[1,:].^2
    edges = sort!(unique(vcat(f, collect(first(f):(.5/(frame+width)):last(f)))))
    logf = log.(f); logpsd = 2 .* log.(asd); variance = 0.
    for j in 1:length(edges)-1
        half=(edges[j+1]-edges[j])/2; mid=(edges[j+1]+edges[j])/2
        for (node, weight) in zip(nodes, weights)
            x=mid+half*node; lx=log(x)
            k=clamp(searchsortedlast(f,x),1,length(f)-1)
            psd=exp(logpsd[k]+(logpsd[k+1]-logpsd[k])*(lx-logf[k])/(logf[k+1]-logf[k]))
            transfer2=sinc(x*width)^2*((1-alpha)^2+4alpha*sinpi(x*frame)^2)
            variance += half*weight*psd*transfer2
        end
    end
    sqrt(variance)
end
end
