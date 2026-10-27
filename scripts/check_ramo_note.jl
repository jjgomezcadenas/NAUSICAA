# Independent numerical checks for the explanation of manuscript Eq. 6.
using HDF5,SHA,TOML,LinearAlgebra
root=normpath(joinpath(@__DIR__,".."));cfgpath=joinpath(root,"config/detector.toml")
cfg=TOML.parsefile(cfgpath);geo=cfg["geometry"];bias=cfg["bias"];gas=cfg["gas"]
e=1.602176634e-19;a=geo["pad_width_m"]/2;g=geo["gap_m"]
mu=gas["bulk_speed_m_s"]/bias["bulk_field_V_m"]
field=(bias["ic_voltage_V"]-bias["pad_voltage_V"])/g;v=mu*field
phi(z,a)=2/pi*atan(a^2,z*sqrt(z^2+2a^2))
weight(z,a)=4a^2/(pi*(z^2+a^2)*sqrt(z^2+2a^2))
function quadrature(n)
    eig=eigen(SymTridiagonal(zeros(n),[k/sqrt(4k^2-1) for k in 1:n-1]))
    eig.values,2eig.vectors[1,:].^2
end
function integrate(f,l,r,n)
    x,w=quadrature(n);sum(w.*f.((r+l)/2 .+(r-l)/2 .*x))*(r-l)/2
end
qfraction=1-phi(g,a)
q128=integrate(z->weight(z,a),0.,g,128)
q256=integrate(z->weight(z,a),0.,g,256)
@assert abs(q128-qfraction)<1e-11 && abs(q128-q256)<1e-11
# Check the 2D Poisson-kernel integral independently of the solid-angle formula.
x,w=quadrature(64)
phi_integral=sum(w[j]*w[k]*a^3/(a^2+(a*x[j])^2+(a*x[k])^2)^1.5 for j in eachindex(x),k in eachindex(x))/(2pi)
@assert abs(phi_integral-phi(a,a))<1e-12
zvalues=a.*[.1,.5,1.,2.,10.,20.]
err=maximum(abs.([-(phi(z+a*1e-5,a)-phi(z-a*1e-5,a))/(2a*1e-5)/weight(z,a)-1 for z in zvalues]))
@assert err<1e-7
@assert maximum(abs.([2weight(2z,2a)/weight(z,a)-1 for z in zvalues]))<1e-14
@assert phi(0.,a)==1 && phi(1e6a,a)<1e-11
peak=e*v*weight(0.,a); closed=2sqrt(2)/pi*e*v/a
@assert abs(peak/closed-1)<1e-14
z90=a*sqrt(1/sin(pi*.1/2)-1)
@assert abs(1-phi(z90,a)-.9)<1e-14
# The schematic's illustrative parallel-plate example, kept distinct from the pad.
rogelio_length=1.25;rogelio_speed=10.;rogelio_count=1e4
frame=.01
out=joinpath(root,"results/ramo_check");mkpath(out)
h5open(joinpath(out,"values.h5"),"w") do h
    for (key,value) in ("pad_side_m"=>2a,"pad_half_side_m"=>a,"gap_m"=>g,"charge_C"=>e,
        "mobility_m2_V_s"=>mu,"gap_field_V_m"=>field,"gap_speed_m_s"=>v,"peak_current_A"=>peak,
        "weight_at_pad_per_m"=>weight(0.,a),"gap_transit_s"=>g/v,"gap_weighting_potential"=>phi(g,a),
        "gap_charge_fraction"=>qfraction,"final_90percent_full_charge_s"=>z90/v,
        "frame_s"=>frame,"uniform_plateau_per_ion_A"=>e*qfraction/frame,
        "rogelio_distance_m"=>rogelio_length,"rogelio_speed_m_s"=>rogelio_speed,"rogelio_ion_count"=>rogelio_count,
        "rogelio_charge_for_10ke_C"=>e*rogelio_count,"rogelio_parallel_plate_current_A"=>e*rogelio_count*rogelio_speed/rogelio_length,
        "ions_in_16fC"=>16e-15/e,"charge_integral_error"=>abs(q128-qfraction),
        "quadrature_convergence_error"=>abs(q128-q256),"poisson_integral_error"=>abs(phi_integral-phi(a,a)),"derivative_relative_error"=>err)
        h[key]=value
    end
    files=[cfgpath,@__FILE__,joinpath(root,"notes/physics/nausica.tex"),joinpath(root,"LTSpice/Spice2/CSAfolded180nm.asc")]
    h["source_paths"]=[relpath(p,root) for p in files]
    h["source_sha256"]=[bytes2hex(sha256(read(p))) for p in files]
    attributes(h)["assumptions"]="Square 2a by 2a in infinite grounded plane; homogeneous upper half-space; pad weighting potential 1, rest of plane 0, zero at infinity; no IFG or surface dielectrics in weighting solution; on-axis constant-speed positive singly charged ion; finite gap trajectory without renormalization."
end
println("a=$(a*1e3) mm; v=$v m/s; peak=$(peak*1e15) fA; Q/e=$qfraction; T90=$(z90/v*1e6) us")
println("10 ke=$(rogelio_count*e*1e15) fC; planar current=$(e*rogelio_count*rogelio_speed/rogelio_length*1e15) fA")
