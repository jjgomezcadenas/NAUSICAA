"""
Radiogenic background of the ITACA detectors: the inner copper shield (ICS) in
which the 214Bi and 208Tl decays are generated.

A shielding card (`config/background/*.toml`, schema
`nausicaa.background_shielding.v1`) names the detector family and variant that
give the fiducial gas, and holds the copper (material, density, thickness,
specific activities, depth sampling length) and the xenon gap and passive
layer on each face. `ShieldGeometry` is the copper around that gas, in cm and
in the detector coordinates of the Geant4 geometry (z = 0 to the fiducial
length along the axis): the cavity of inner radius `inner_cm` between the
copper disks at `zmin_cm` and `zmax_cm`, with the copper `thickness_cm` on
every face.

The Geant4 region `ICS` (ItacaChamber.cc) draws the decay vertex in three
steps: the piece (barrel, cathode disk at z < zmin, anode disk at z > zmax)
with probability proportional to its volume, the depth from the inner face
with density exp(-depth/lambda) truncated at the thickness (uniform when
lambda = 0), and a uniform position on the surface at that depth.
`sample_vertex` reproduces those steps. The uniform distribution in the
copper is restored by `weight`, A(depth) / (V_piece p(depth)), where A is the
area of the surface at that depth and p the depth density: the mean of the
weights over the sampled vertices is one, and the weighted sum of the events
in the energy window over the number of decays generated estimates the
fraction of uniform decays that reach the window.
"""
module Background

using TOML
using SHA
using Random

export ShieldCard, load_shielding, ShieldGeometry, shield_geometry, cavity_length_cm, outer_cm,
       volume_barrel_cm3, volume_disk_cm3, volume_cm3, copper_mass_kg, activity_Bq, decays_per_year,
       depth_density, sample_vertex, piece, depth_cm, weight, SECONDS_PER_YEAR

const SCHEMA = "nausicaa.background_shielding.v1"
const FAMILIES = ("iitaca", "eitaca")
const FACES = ("barrel", "anode", "cathode")
const ICS_KEYS = ("material", "density_g_cm3", "thickness_cm", "bi214_Bq_kg", "tl208_Bq_kg", "sampling_length_cm")
"Julian year, as in the activity budget of ITACA revisited."
const SECONDS_PER_YEAR = 365.25 * 86400

"""A shielding card: path, text, SHA-256 of the text, and the parsed values."""
struct ShieldCard
    path::String
    text::String
    sha256::String
    card::Dict{String,Any}
end

"""
    load_shielding(path)::ShieldCard

Read a shielding card and require its schema, the family, the variant path,
the saved window, the `[ics]` values and the gap, layer material and layer
thickness of the three faces.
"""
function load_shielding(path::AbstractString)::ShieldCard
    text = read(path, String); d = TOML.parse(text)
    get(d, "schema", "") == SCHEMA || throw(ArgumentError("Expected schema $SCHEMA in $path"))
    for k in ("name", "family", "variant", "saved_window", "ics", FACES...)
        haskey(d, k) || throw(ArgumentError("Missing $k in $path"))
    end
    d["family"] in FAMILIES || throw(ArgumentError("Unknown family $(d["family"]) in $path"))
    d["saved_window"] in ("gamma", "general") || throw(ArgumentError("saved_window must be gamma or general in $path"))
    for k in ICS_KEYS
        haskey(d["ics"], k) || throw(ArgumentError("Missing ics.$k in $path"))
    end
    d["ics"]["thickness_cm"] > 0 || throw(ArgumentError("ics.thickness_cm must be positive in $path"))
    d["ics"]["sampling_length_cm"] >= 0 || throw(ArgumentError("ics.sampling_length_cm must not be negative in $path"))
    for f in FACES, k in ("gap_cm", "layer_material", "layer_thickness_cm")
        haskey(d[f], k) || throw(ArgumentError("Missing $f.$k in $path"))
        k == "layer_material" || d[f][k] >= 0 || throw(ArgumentError("$f.$k must not be negative in $path"))
    end
    ShieldCard(String(path), text, bytes2hex(sha256(text)), d)
end

"""The copper around the fiducial gas, in cm (see the module documentation)."""
struct ShieldGeometry
    radius_cm::Float64            # fiducial gas radius
    length_cm::Float64            # fiducial gas length, z = 0 to length
    inner_cm::Float64             # inner radius of the copper barrel
    zmin_cm::Float64              # inner face of the cathode disk (z < 0 side)
    zmax_cm::Float64              # inner face of the anode disk
    thickness_cm::Float64         # copper on every face
    sampling_length_cm::Float64   # depth sampling length lambda; 0 = uniform
    density_g_cm3::Float64        # copper density
end

"""
    shield_geometry(card, radius_cm, length_cm)::ShieldGeometry

The geometry of a card around a fiducial gas of the given radius and length,
as ItacaChamber::ConstructShielding builds it.
"""
function shield_geometry(card::ShieldCard, radius_cm::Real, length_cm::Real)::ShieldGeometry
    c = card.card; b = c["barrel"]; a = c["anode"]; k = c["cathode"]; ics = c["ics"]
    ShieldGeometry(radius_cm, length_cm, radius_cm + b["gap_cm"] + b["layer_thickness_cm"],
                   -(k["gap_cm"] + k["layer_thickness_cm"]), length_cm + a["gap_cm"] + a["layer_thickness_cm"],
                   ics["thickness_cm"], ics["sampling_length_cm"], ics["density_g_cm3"])
end

cavity_length_cm(g::ShieldGeometry) = g.zmax_cm - g.zmin_cm
outer_cm(g::ShieldGeometry) = g.inner_cm + g.thickness_cm
volume_barrel_cm3(g::ShieldGeometry) = π * (outer_cm(g)^2 - g.inner_cm^2) * cavity_length_cm(g)
volume_disk_cm3(g::ShieldGeometry) = π * outer_cm(g)^2 * g.thickness_cm
volume_cm3(g::ShieldGeometry) = volume_barrel_cm3(g) + 2 * volume_disk_cm3(g)
copper_mass_kg(g::ShieldGeometry) = volume_cm3(g) * g.density_g_cm3 / 1000
"Activity of the copper in Bq from its specific activity in Bq/kg."
activity_Bq(g::ShieldGeometry, Bq_kg::Real) = copper_mass_kg(g) * Bq_kg
decays_per_year(g::ShieldGeometry, Bq_kg::Real) = activity_Bq(g, Bq_kg) * SECONDS_PER_YEAR

"""Density of the sampled depth at `depth_cm`: exponential truncated at the thickness, or uniform."""
function depth_density(g::ShieldGeometry, depth_cm::Real)
    0 <= depth_cm <= g.thickness_cm || return 0.0
    λ = g.sampling_length_cm; t = g.thickness_cm
    λ > 0 ? exp(-depth_cm / λ) / (λ * (1 - exp(-t / λ))) : 1 / t
end

"""
    sample_vertex(g, rng) -> (x_cm, y_cm, z_cm)

One decay vertex with the steps of region ICS of ItacaChamber::GenerateVertex,
in the same order of the random numbers: depth, piece, azimuth, then the
radius (disks) or the axial position (barrel).
"""
function sample_vertex(g::ShieldGeometry, rng::AbstractRNG)
    t = g.thickness_cm; λ = g.sampling_length_cm; inner = g.inner_cm; outer = outer_cm(g); L = cavity_length_cm(g)
    v_barrel = volume_barrel_cm3(g); v_disk = volume_disk_cm3(g)
    depth = λ > 0 ? -λ * log(1 - rand(rng) * (1 - exp(-t / λ))) : t * rand(rng)
    u = rand(rng) * (v_barrel + 2 * v_disk)
    φ = 2π * rand(rng)
    if u < v_barrel
        r = inner + depth; z = g.zmin_cm + L * rand(rng)
    else
        r = outer * sqrt(rand(rng)); z = u < v_barrel + v_disk ? g.zmin_cm - depth : g.zmax_cm + depth
    end
    (r * cos(φ), r * sin(φ), z)
end

"""
    piece(g, x_cm, y_cm, z_cm) -> :barrel | :cathode | :anode

The piece of the copper holding a point; an error for a point outside the copper.
"""
function piece(g::ShieldGeometry, x::Real, y::Real, z::Real; tolerance_cm = 1e-4)
    r = hypot(x, y); t = g.thickness_cm; out = outer_cm(g)
    r <= out + tolerance_cm || throw(ArgumentError("Point outside the copper: r = $r cm"))
    if z < g.zmin_cm
        z >= g.zmin_cm - t - tolerance_cm || throw(ArgumentError("Point beyond the cathode disk: z = $z cm"))
        return :cathode
    elseif z > g.zmax_cm
        z <= g.zmax_cm + t + tolerance_cm || throw(ArgumentError("Point beyond the anode disk: z = $z cm"))
        return :anode
    end
    r >= g.inner_cm - tolerance_cm || throw(ArgumentError("Point inside the cavity: r = $r cm, z = $z cm"))
    :barrel
end

"""Depth of a point from the inner face of its piece, clamped to [0, thickness]."""
function depth_cm(g::ShieldGeometry, x::Real, y::Real, z::Real)
    p = piece(g, x, y, z)
    d = p == :barrel ? hypot(x, y) - g.inner_cm : p == :cathode ? g.zmin_cm - z : z - g.zmax_cm
    clamp(d, 0.0, g.thickness_cm)
end

"""
    weight(g, x_cm, y_cm, z_cm)

Weight of a vertex drawn by `sample_vertex` (or region ICS) that restores the
uniform distribution in the copper: A(depth) / (V_piece p(depth)), with A the
area of the surface at that depth (2 pi (inner + depth) L_cavity for the
barrel, pi outer^2 for a disk) and V_piece the volume of the piece. Its mean
over the sampled vertices is one.
"""
function weight(g::ShieldGeometry, x::Real, y::Real, z::Real)
    p = piece(g, x, y, z); d = depth_cm(g, x, y, z)
    area = p == :barrel ? 2π * (g.inner_cm + d) * cavity_length_cm(g) : π * outer_cm(g)^2
    volume = p == :barrel ? volume_barrel_cm3(g) : volume_disk_cm3(g)
    area / (volume * depth_density(g, d))
end

end # module
