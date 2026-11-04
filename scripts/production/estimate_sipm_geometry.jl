#!/usr/bin/env julia
using HDF5,TOML,SHA
const ROOT=normpath(joinpath(@__DIR__,"../.."))

"""Depth-averaged fraction of isotropic light incident on a centred square SiPM.

Emission is uniform along the hole axis. All lengths are mm. The square's
corners must lie within the hole radius, making the emitting-hole aperture
unrestrictive for rays to this aligned sensor at nonnegative gap. The solid
angle at distance d is 4atan(a²/(d*sqrt(d²+2a²))), with half-side a. Midpoint
integration resolves its depth average; two resolutions check convergence.
"""
function square_acceptance(side,depth,gap;steps=100000)
    a=side/2
    sum(atan(a*a,(gap+(i-.5)*depth/steps)*sqrt((gap+(i-.5)*depth/steps)^2+2a*a))/pi
        for i in 1:steps)/steps
end

"""Exact aligned circular-SiPM acceptance for radius equal to the hole radius."""
function circle_acceptance(radius,depth,gap)
    (depth+hypot(gap,radius)-hypot(gap+depth,radius))/(2depth)
end

"""Compare the current circular sensor with smaller squares at zero and 10 mm gap.

Store geometrical probabilities, incident-light ratios and microcell occupancy
ratios. These describe the sensor aligned with one emitting hole. A complete
array also collects contributions from neighbouring holes.
"""
function estimate_geometry()
    ep=joinpath(ROOT,"config/specs/detector/electrons/elcc_250ns_ep2.toml")
    sp=joinpath(ROOT,"config/specs/response/sipm/sipm_statistical_v1.toml")
    e=TOML.parsefile(ep);s=TOML.parsefile(sp)
    radius=e["hole_diameter_mm"]/2;depth=e["hole_depth_mm"]
    isapprox(e["sipm_diameter_mm"],2radius) || error("Circular reference must cover the exit")
    out=joinpath(ROOT,"results/sipm_geometry_estimate");mkpath(out)
    write(joinpath(out,".gitignore"),"*.h5\n")
    h5open(joinpath(out,"estimate.h5"),"w") do f
        meta=create_group(f,"metadata");meta["complete"]=false
        meta["electron_card"]=read(ep,String);meta["sensor_card"]=read(sp,String)
        meta["source"]=read(@__FILE__,String);meta["source_sha256"]=bytes2hex(sha256(read(@__FILE__)))
        meta["hole_radius_mm"]=radius;meta["hole_depth_mm"]=depth
        baseline=circle_acceptance(radius,depth,0.)
        pitch_mm=s["microcell_pitch_m"]*1000;circle_area=pi*radius^2
        for gap in (0.,10.)
            reference=circle_acceptance(radius,depth,gap)
            g=create_group(f,"gap_$(Int(gap))");g["gap_mm"]=gap
            g["circle_acceptance_percent"]=100reference
            for side in (3.,2.,1.)
                side/sqrt(2)<=radius || error("Square corners extend beyond aperture")
                p=square_acceptance(side,depth,gap)
                coarse=square_acceptance(side,depth,gap;steps=50000)
                abs(p-coarse)<1e-9 || error("Depth quadrature convergence")
                row=create_group(g,"square_$(Int(side))")
                for (key,value) in pairs((;side_mm=side,acceptance_percent=100p,
                    relative_to_circle_percent=100p/reference,relative_to_original_percent=100p/baseline,
                    original_reduction_factor=baseline/p,area_ratio=side^2/circle_area,
                    microcells=round(Int,side^2/pitch_mm^2),
                    average_occupancy_ratio=(p/reference)/(side^2/circle_area),quadrature_error=abs(p-coarse)))
                    row[string(key)]=value
                end
                println("gap=$gap mm, square=$side mm: A=$(100p)%, relative=$(100p/reference)%")
            end
        end
        write(meta["complete"],true)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__;estimate_geometry();end
