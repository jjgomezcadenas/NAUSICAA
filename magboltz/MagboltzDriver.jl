module MagboltzDriver

using HDF5, SHA, TOML
include("../src/ion_signal/guarded_files.jl")
using .GuardedFiles

export GasPoint, read_point, density_factor, effective_pressure_bar,
       number_density_m3, reduced_field_Td, input_card, parse_output,
       diffusion_units, run_point, save_table

const K_B_J_K = 1.380649e-23
const BAR_TO_PA = 1e5
const BAR_TO_TORR = 750.061683
const THERMAL_COLLISIONS_PER_NMAX = 40_000_000
const XE_ENRICHED_MOLAR_MASS_G_MOL = 135.778115758
const HE_MOLAR_MASS_G_MOL = 4.0026
const NUMBER = raw"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DdEe][+-]?\d+)?"

"""One Magboltz point. Pressure is the physical chamber pressure. `density_mode`
selects plain ideal-gas pressure (`:plain`), the handoff's partial-number
formula (`:partial_number`), or the number density implied by the Geant4
material density and mass fractions (`:geant4_material`). The selected density
is expressed as an equivalent ideal-gas pressure in the Magboltz input card.
With thermal motion, the downloaded source multiplies `nmax` by 40 million.
"""
Base.@kwdef struct GasPoint
    gas::Symbol = :Xe_85_He_15
    temperature_K::Float64 = 290.0
    pressure_bar::Float64 = 17.5
    field_V_cm::Float64 = 300.0
    density_mode::Symbol = :plain
    nmax::Int = 2
    energy_limit_eV::Float64 = 0.0
    penning::Bool = false
    thermal_motion::Bool = true
end

function validate(p::GasPoint)
    p.gas in (:Xe_85_He_15, :Xe_100) || throw(ArgumentError("Unknown gas"))
    p.density_mode in (:plain, :partial_number, :geant4_material) || throw(ArgumentError("Unknown density mode"))
    all(isfinite, (p.temperature_K,p.pressure_bar,p.field_V_cm,p.energy_limit_eV)) || throw(ArgumentError("Nonfinite input"))
    p.temperature_K>273.15 && p.pressure_bar>0 && p.field_V_cm>0 || throw(ArgumentError("Invalid temperature, pressure or field"))
    p.nmax>0 && p.energy_limit_eV>=0 || throw(ArgumentError("Invalid collisions or energy limit"))
    nothing
end

"""Read one point from a TOML card with explicit units."""
function read_point(path::AbstractString)::GasPoint
    c=TOML.parsefile(path)
    p=GasPoint(gas=Symbol(c["gas"]),temperature_K=Float64(c["temperature_K"]),
        pressure_bar=Float64(c["pressure_bar"]),field_V_cm=Float64(c["field_V_cm"]),
        density_mode=Symbol(c["density_mode"]),nmax=Int(c["nmax"]),
        energy_limit_eV=Float64(get(c,"energy_limit_eV",0.0)),
        penning=Bool(get(c,"penning",false)),thermal_motion=Bool(get(c,"thermal_motion",true)))
    validate(p);p
end

"""Ratio of gas number density to ideal density at the physical pressure.
Both corrections use the Xe compressibility in `ItacaChamber.cc` and ideal He.
The Geant4 material mode follows its mixture mass density and nominal mass
fractions. The enriched Xe molar mass uses the listed isotope mass numbers;
the resulting density factor is insensitive to the sub-percent isotope-mass
differences.
"""
function density_factor(p::GasPoint)::Float64
    validate(p)
    p.density_mode==:plain && return 1.0
    tr=p.temperature_K/289.733
    pr=p.pressure_bar/58.420
    z_xe=1+(0.083-0.422/tr^1.6)*pr/tr
    z_xe>0 || throw(ArgumentError("Nonpositive Xe compressibility"))
    p.gas==:Xe_100 && return 1/z_xe
    p.density_mode==:partial_number && return 0.85/z_xe+0.15
    xe_mass=0.85*XE_ENRICHED_MOLAR_MASS_G_MOL
    he_mass=0.15*HE_MOLAR_MASS_G_MOL
    (xe_mass/z_xe+he_mass)/(xe_mass+he_mass)
end

effective_pressure_bar(p::GasPoint)::Float64=p.pressure_bar*density_factor(p)
number_density_m3(p::GasPoint)::Float64=effective_pressure_bar(p)*BAR_TO_PA/(K_B_J_K*p.temperature_K)
reduced_field_Td(p::GasPoint)::Float64=(100*p.field_V_cm/number_density_m3(p))/1e-21

"""Make the five free-format input cards read by this Fortran source.
Gas IDs Xe=7 and He-4=3. The final zero-gas card terminates the run.
"""
function input_card(p::GasPoint)::String
    validate(p)
    mixture=p.gas==:Xe_85_He_15
    gas_ids=mixture ? "7 3 0 0 0 0" : "7 0 0 0 0 0"
    fractions=mixture ? "85 15 0 0 0 0" : "100 0 0 0 0 0"
    join(("$(mixture ? 2 : 1) $(p.nmax) $(Int(p.penning)) $(Int(p.thermal_motion)) $(p.energy_limit_eV)",
          gas_ids,
          "$fractions $(p.temperature_K-273.15) $(effective_pressure_bar(p)*BAR_TO_TORR)",
          "$(p.field_V_cm) 0 0",
          "0 0 0 0 0"),"\n")*"\n"
end

parse_fortran_number(s::AbstractString)=parse(Float64,replace(s,'D'=>'E','d'=>'E'))

function captures(text::AbstractString,pattern::Regex,label::AbstractString)
    m=match(pattern,text)
    m===nothing && throw(ArgumentError("Missing Magboltz $label"))
    map(parse_fortran_number,m.captures)
end

function diffusion_output(text::AbstractString,label::String,other::String)
    parts=split(text,label;limit=2)
    length(parts)==2 || throw(ArgumentError("Missing $label"))
    block=split(parts[2],other;limit=2)[1]
    d,err=captures(block,Regex(raw"=\s*("*NUMBER*raw")\s*\+\-\s*("*NUMBER*raw")\s*%"),label)
    spatial,spatial_err=captures(block,Regex(raw"=\s*("*NUMBER*raw")\s+MICRONS/CENTIMETER\*\*0\.5\s*\+\-\s*("*NUMBER*raw")\s*%"),label*" spatial")
    (diffusion_cm2_s=d,diffusion_error_percent=err,spatial_um_sqrtcm=spatial,spatial_error_percent=spatial_err)
end

"""Parse the zero-magnetic-field drift velocity and diffusion with their
reported statistical errors. Velocity in microns/ns is numerically mm/us.
"""
function parse_output(text::AbstractString)
    velocity,velocity_err=captures(text,Regex(raw"Z DRIFT VELOCITY\s*=\s*("*NUMBER*raw")\s+MICRONS/NANOSECOND\s*\+\-\s*("*NUMBER*raw")\s*%"),"Z drift velocity")
    transverse=diffusion_output(text,"TRANSVERSE DIFFUSION","LONGITUDINAL DIFFUSION")
    longitudinal=diffusion_output(text,"LONGITUDINAL DIFFUSION","IONISATION RATE")
    m=match(Regex(raw"NUMBER OF COLLISIONS IN FINAL ENERGY BIN\s*=\s*("*NUMBER*raw")"),text)
    last_bin=m===nothing ? NaN : parse_fortran_number(m.captures[1])
    (;velocity_mm_us=velocity,velocity_error_percent=velocity_err,
      transverse_diffusion_cm2_s=transverse.diffusion_cm2_s,
      transverse_diffusion_error_percent=transverse.diffusion_error_percent,
      transverse_um_sqrtcm=transverse.spatial_um_sqrtcm,
      transverse_spatial_error_percent=transverse.spatial_error_percent,
      longitudinal_diffusion_cm2_s=longitudinal.diffusion_cm2_s,
      longitudinal_diffusion_error_percent=longitudinal.diffusion_error_percent,
      longitudinal_um_sqrtcm=longitudinal.spatial_um_sqrtcm,
      longitudinal_spatial_error_percent=longitudinal.spatial_error_percent,
      final_energy_bin_collisions=last_bin)
end

"""Convert Magboltz's microns/sqrt(cm) to mm/sqrt(m) and to the card
coefficient in mm*sqrt(bar)/sqrt(cm). The card coefficient uses the physical
pressure, including when Magboltz receives an equivalent ideal pressure.
"""
function diffusion_units(spatial_um_sqrtcm::Real,physical_pressure_bar::Real)
    spatial_um_sqrtcm>=0 && physical_pressure_bar>0 || throw(ArgumentError("Invalid diffusion or pressure"))
    (mm_sqrtm=Float64(spatial_um_sqrtcm)*0.01,
     card_mm_sqrtbar_sqrtcm=Float64(spatial_um_sqrtcm)*0.001*sqrt(physical_pressure_bar))
end

"""Run one point with its input and raw output retained in `directory`.
The caller selects an already compiled Fortran binary and an empty directory.
"""
function run_point(p::GasPoint,binary::AbstractString,directory::AbstractString)
    validate(p)
    isfile(binary) || throw(ArgumentError("Magboltz binary missing"))
    ispath(directory) && throw(ArgumentError("Run directory exists"))
    mkpath(directory)
    input=joinpath(directory,"input.card");output=joinpath(directory,"magboltz.out");error=joinpath(directory,"magboltz.err")
    write(input,input_card(p))
    command=Cmd(Cmd([abspath(binary)]);dir=abspath(directory))
    elapsed=@elapsed run(pipeline(command;stdin=input,stdout=output,stderr=error))
    result=parse_output(read(output,String))
    (;point=p,result,elapsed_s=elapsed,input,output,error)
end

"""Save one or more completed points and source identities to HDF5."""
function save_table(path::AbstractString,rows::AbstractVector;source_path::AbstractString,binary_path::AbstractString)
    isempty(rows) && throw(ArgumentError("No Magboltz points"))
    ispath(path) && throw(ArgumentError("Result exists"))
    mkpath(dirname(path))
    tmp=path*".tmp"
    ispath(tmp) && throw(ArgumentError("Temporary result exists"))
    raw_outputs=[read(r.output,String) for r in rows]
    raw_errors=[read(r.error,String) for r in rows]
    cards=[read(r.input,String) for r in rows]
    budget=256*1024+3sum(ncodeunits,raw_outputs)+3sum(ncodeunits,raw_errors)+3sum(ncodeunits,cards)+2048length(rows)
    GuardedFiles.space_check(dirname(path),budget,1024^3)
    try
      h5open(tmp,"w") do f
        f["schema"]="itaca.magboltz.transport.v1"
        f["source_url"]="https://magboltz.web.cern.ch/magboltz/magboltz-11.21.f"
        f["source_filename_version"]="11.21"
        f["source_header_version"]="11.20"
        f["source_sha256"]=bytes2hex(open(sha256,source_path))
        f["binary_sha256"]=bytes2hex(open(sha256,binary_path))
        f["driver_sha256"]=bytes2hex(open(sha256,@__FILE__))
        f["gas"]=string.([row.point.gas for row in rows])
        f["density_mode"]=string.([row.point.density_mode for row in rows])
        for (name,getter) in (("temperature_K",r->r.point.temperature_K),
            ("pressure_bar",r->r.point.pressure_bar),("effective_pressure_bar",r->effective_pressure_bar(r.point)),
            ("field_V_cm",r->r.point.field_V_cm),("reduced_field_Td",r->reduced_field_Td(r.point)),
            ("nmax",r->r.point.nmax),("elapsed_s",r->r.elapsed_s))
            f[name]=[getter(row) for row in rows]
        end
        f["collision_target"]=[r.point.nmax*THERMAL_COLLISIONS_PER_NMAX for r in rows]
        create_group(f,"results")
        for name in keys(rows[1].result)
            f["results/$(name)"]=[getfield(row.result,name) for row in rows]
        end
        f["results/velocity_error_mm_us"]=[r.result.velocity_mm_us*r.result.velocity_error_percent/100 for r in rows]
        f["results/drift_time_150cm_ms"]=[1.5/r.result.velocity_mm_us for r in rows]
        for (label,field) in (("transverse","transverse_um_sqrtcm"),("longitudinal","longitudinal_um_sqrtcm"))
            f["results/$(label)_mm_sqrtm"]=[diffusion_units(getfield(r.result,Symbol(field)),r.point.pressure_bar).mm_sqrtm for r in rows]
            f["results/$(label)_card_mm_sqrtbar_sqrtcm"]=[diffusion_units(getfield(r.result,Symbol(field)),r.point.pressure_bar).card_mm_sqrtbar_sqrtcm for r in rows]
            f["results/$(label)_card_error_mm_sqrtbar_sqrtcm"]=[diffusion_units(getfield(r.result,Symbol(field)),r.point.pressure_bar).card_mm_sqrtbar_sqrtcm*getfield(r.result,Symbol(label*"_spatial_error_percent"))/100 for r in rows]
        end
        f["input_cards"]=cards
        f["raw_outputs"]=raw_outputs
        f["raw_errors"]=raw_errors
      end
      GuardedFiles.check_size(tmp,budget)
      Base.rename(tmp,path)
    catch
      isfile(tmp) && rm(tmp)
      rethrow()
    end
    path
end

end
