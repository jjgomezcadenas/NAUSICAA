"""
    EventPlacement(event_id::Int64, ion_drift_cm::Float32,
                   electron_drift_cm::Float32)

One event location shared by the ion, pure-Xe electron and Xe--He electron
branches. `event_id` is the source ID within its sample; both distances are
nonnegative Float32 cm. The run boundary verifies that their sum is the
configured chamber length. This immutable value preserves its inputs.
"""
struct EventPlacement
    event_id::Int64
    ion_drift_cm::Float32
    electron_drift_cm::Float32
    function EventPlacement(event_id::Int64,ion_drift_cm::Float32,electron_drift_cm::Float32)
        all(x->isfinite(x)&&x>=0f0,(ion_drift_cm,electron_drift_cm)) ||
            throw(ArgumentError("Drift lengths must be finite and nonnegative cm"))
        new(event_id,ion_drift_cm,electron_drift_cm)
    end
end

"""
    GridSpec(pitch_mm::Float32, bucket_width_mm::Float32,
             reference_mm::NTuple{3,Float32}, phase_mm::NTuple{2,Float32},
             drift_direction::Int8)

Physical transverse pixel lattice and longitudinal bucket, with every length
in mm. `reference_mm` is the common source-frame edge reference; the two
independent transverse phases lie in `[0,pitch_mm)`. The phase is fixed for an
event and detector branch. `drift_direction=+1` gives later ion buckets at
larger source z; `-1` gives later electron buckets at smaller source z. Bucket
indices may be signed because diffusion may cross the reference edge.
"""
struct GridSpec
    pitch_mm::Float32
    bucket_width_mm::Float32
    reference_mm::NTuple{3,Float32}
    phase_mm::NTuple{2,Float32}
    drift_direction::Int8
    function GridSpec(pitch_mm::Float32,bucket_width_mm::Float32,
            reference_mm::NTuple{3,Float32},phase_mm::NTuple{2,Float32},
            drift_direction::Int8)
        isfinite(pitch_mm)&&pitch_mm>0f0&&isfinite(bucket_width_mm)&&bucket_width_mm>0f0 ||
            throw(ArgumentError("Grid pitch and bucket width must be finite positive mm"))
        all(isfinite,reference_mm)&&all(isfinite,phase_mm)&&
            all(p->0f0<=p<pitch_mm,phase_mm) ||
            throw(ArgumentError("Grid reference must be finite and phase in [0,pitch)"))
        drift_direction in (Int8(-1),Int8(1)) ||
            throw(ArgumentError("Grid drift direction must be ±1"))
        new(pitch_mm,bucket_width_mm,reference_mm,phase_mm,drift_direction)
    end
end

"""
    SparseCells(pixel_ix, pixel_iy, bucket, count)

Owned column vectors for occupied physical pixel/time cells. Every vector is
`Vector{Int64}`; cell indices are signed and `count` is a strictly positive
number of ions or incident photons. Keys `(pixel_ix,pixel_iy,bucket)` are
unique and sorted lexicographically. Zero-signal events use four empty vectors.
Construction copies all input arrays; their later mutation cannot change the
history. The enclosing history identifies the count unit.
"""
struct SparseCells
    pixel_ix::Vector{Int64}
    pixel_iy::Vector{Int64}
    bucket::Vector{Int64}
    count::Vector{Int64}
    function SparseCells(pixel_ix::Vector{Int64},pixel_iy::Vector{Int64},
            bucket::Vector{Int64},count::Vector{Int64})
        n=length(pixel_ix)
        all(v->length(v)==n,(pixel_iy,bucket,count)) ||
            throw(DimensionMismatch("Sparse cell columns must have equal lengths"))
        all(>(0),count) || throw(ArgumentError("Sparse cell counts must be positive"))
        keys=collect(zip(pixel_ix,pixel_iy,bucket))
        issorted(keys)&&allunique(keys) ||
            throw(ArgumentError("Sparse cell keys must be unique and sorted"))
        new(copy(pixel_ix),copy(pixel_iy),copy(bucket),copy(count))
    end
end

"""
    IonPadHistory(placement, grid, cells, initial_pairs)

One event's ion counts by cathode pad and time bucket. `grid` supplies mm
geometry and phase; the ion frame duration is recorded in the detector card.
`initial_pairs::Int64` is the shared primary pair count. Every ion reaches
one cell in the unbounded-grid pilot, so `sum(cells.count)==initial_pairs` is
enforced. The HDF5 writer saves this structure in a versioned event group.
"""
struct IonPadHistory
    placement::EventPlacement
    grid::GridSpec
    cells::SparseCells
    initial_pairs::Int64
    function IonPadHistory(placement::EventPlacement,grid::GridSpec,
            cells::SparseCells,initial_pairs::Int64)
        initial_pairs>=0&&foldl(Base.Checked.checked_add,cells.count;init=Int64(0))==initial_pairs ||
            throw(ArgumentError("Ion cell sum must equal initial pair count"))
        new(placement,grid,cells,initial_pairs)
    end
end

"""
    PhotonSiPMHistory(placement, gas_model, grid, cells, initial_pairs,
                      surviving_electrons, emitted_photons)

Incident VUV photons by anode SiPM and time bucket for `gas_model` equal to
`:pure_xe` or `:xe_he15`. The initial pair count is shared with the ion branch.
Attachment leaves `surviving_electrons`; EL produces `emitted_photons`; the
sum of sparse cell counts is the incident-photon total and cannot exceed the
emitted total. The electron bucket duration is recorded in the detector card.
SiPM conversion to photoelectrons belongs to a later response stage.
"""
struct PhotonSiPMHistory
    placement::EventPlacement
    gas_model::Symbol
    grid::GridSpec
    cells::SparseCells
    initial_pairs::Int64
    surviving_electrons::Int64
    emitted_photons::Int64
    function PhotonSiPMHistory(placement::EventPlacement,gas_model::Symbol,
            grid::GridSpec,cells::SparseCells,initial_pairs::Int64,
            surviving_electrons::Int64,emitted_photons::Int64)
        gas_model in (:pure_xe,:xe_he15) || throw(ArgumentError("Unknown gas model"))
        0<=surviving_electrons<=initial_pairs && emitted_photons>=0 &&
            foldl(Base.Checked.checked_add,cells.count;init=Int64(0))<=emitted_photons ||
            throw(ArgumentError("Photon/electron ledger is inconsistent"))
        new(placement,gas_model,grid,cells,initial_pairs,
            surviving_electrons,emitted_photons)
    end
end

"""
    sample_integer_moments(rng, mean::Float64, variance_ratio::Float64)::Int64

Sample a nonnegative integer with target mean `mean` and target variance
`variance_ratio*mean`, for `0 ≤ variance_ratio < 1`. The adjacent-integer
law has the minimum possible variance `f*(1-f)`, where `f` is the fractional
part of `mean`. Mix that law with `Binomial(ceil(mean/(1-ratio)), mean/n)`
to reach the target variance when feasible. Below the integer variance floor,
use the adjacent law. `mean=0` returns zero. Inputs are preserved; the caller
owns the RNG. Invalid or unrepresentable parameters raise an error.
"""
function sample_integer_moments(rng::AbstractRNG,mean::Float64,variance_ratio::Float64)::Int64
    isfinite(mean)&&0<=mean<typemax(Int64) &&
        isfinite(variance_ratio)&&0<=variance_ratio<1 ||
        throw(ArgumentError("Require representable mean ≥ 0 and variance ratio in [0,1)"))
    mean==0 && return Int64(0)
    lo=floor(Int64,mean);f=mean-Float64(lo)
    vfloor=f*(1-f);target=variance_ratio*mean
    adjacent()=lo+Int64(rand(rng)<f)
    target<=vfloor+eps(mean) && return adjacent()
    nfloat=ceil(mean/(1-variance_ratio))
    nfloat<typemax(Int64) || throw(OverflowError("Binomial trial count overflows Int64"))
    n=Int64(nfloat)
    vbin=mean*(1-mean/Float64(n))
    mix=clamp((target-vfloor)/(vbin-vfloor),0.0,1.0)
    return rand(rng)<mix ? Int64(rand(rng,Binomial(n,mean/Float64(n)))) : adjacent()
end

"""
    sample_fano_pairs(rng, hits::EventHits, W_eV::Float32, F::Float32)::Vector{Int64}

Sample one initial electron–ion pair count for each ordered source hit. The
per-hit mean is `1e6*energy_MeV/W_eV`; for means below one, the integer draw
is Bernoulli(mean). For larger means, use `sample_integer_moments` with target
variance `F*mean`. The same resulting vector supplies both carrier branches.
`F` is a sub-Poisson variance ratio in [0,1). Preserve hits and consume only
the supplied RNG. Invalid inputs and sum overflow raise errors.
"""
function sample_fano_pairs(rng::AbstractRNG,hits::EventHits,
        W_eV::Float32,F::Float32)::Vector{Int64}
    hits.frame==source_frame || throw(ArgumentError("Pair sampling requires source-frame hits"))
    _validate_hit_vectors(hits.x_mm,hits.y_mm,hits.z_mm,hits.energy_MeV,
                          hits.time_ns,hits.particle_id,hits.hit_id)
    _check_w(W_eV)
    isfinite(F)&&0f0<=F<1f0 || throw(ArgumentError("F must be in [0,1)"))
    counts=Vector{Int64}(undef,length(hits.energy_MeV))
    for i in eachindex(counts)
        mean=Float64(hits.energy_MeV[i])*1e6/Float64(W_eV)
        counts[i]=sample_integer_moments(rng,mean,Float64(F))
    end
    foldl(Base.Checked.checked_add,counts;init=Int64(0))
    return counts
end

"""
    sample_attached_electrons(rng, pairs, drift_ms, lifetime_ms)::Vector{Int64}

For each hit draw surviving electrons from Binomial(pairs[i], exp(-drift/lifetime)).
The drift and lifetime are finite nonnegative and positive milliseconds,
respectively. Preserve the common pair vector; zero drift returns a copy.
"""
function sample_attached_electrons(rng::AbstractRNG,pairs::Vector{Int64},
        drift_ms::Float32,lifetime_ms::Float32)::Vector{Int64}
    all(>=(0),pairs)&&isfinite(drift_ms)&&drift_ms>=0f0&&
        isfinite(lifetime_ms)&&lifetime_ms>0f0 ||
        throw(ArgumentError("Invalid attachment inputs"))
    p=exp(-Float64(drift_ms)/Float64(lifetime_ms))
    return Int64[rand(rng,Binomial(n,p)) for n in pairs]
end

"""
    sample_el_photons(rng, electrons, yield_per_e, J)::Int64

Draw emitted VUV photons from the conditional moments mean `electrons*yield`
and variance `electrons*J*yield`. `J` is the Oliveira EL variance ratio,
supplied by the run model for the selected E/P. The integer variance floor is
handled by `sample_integer_moments`. Require finite nonnegative yield and
`0 ≤ J < 1`; return zero for zero electrons.
"""
function sample_el_photons(rng::AbstractRNG,electrons::Int64,
        yield_per_e::Float32,J::Float32)::Int64
    electrons>=0&&isfinite(yield_per_e)&&yield_per_e>=0f0&&
        isfinite(J)&&0f0<=J<1f0 || throw(ArgumentError("Invalid EL inputs"))
    mean=Float64(electrons)*Float64(yield_per_e)
    return sample_integer_moments(rng,mean,Float64(J))
end

"""
    sample_incident_photons(rng, emitted::Int64, acceptance::Float32)::Int64

Draw the number of emitted VUV photons that reach the exit-plane SiPM directly.
With `emitted` fixed, the result is binomial with mean `emitted*acceptance`
and variance `emitted*acceptance*(1-acceptance)`. The acceptance comes from
`axial_light_acceptance` for absorbing hole walls and an exit-covering SiPM.
The output is an incident-photon count, before PDE and electronics. Preserve
the inputs; reject negative counts or probabilities outside [0,1].
"""
function sample_incident_photons(rng::AbstractRNG,emitted::Int64,
        acceptance::Float32)::Int64
    emitted>=0&&isfinite(acceptance)&&0f0<=acceptance<=1f0 ||
        throw(ArgumentError("Invalid optical inputs"))
    return Int64(rand(rng,Binomial(emitted,Float64(acceptance))))
end

"""Convert an accumulated integer-count dictionary into owned sorted sparse cells."""
function _sparse_cells(sums::Dict{NTuple{3,Int64},Int64})::SparseCells
    keys_sorted=sort!(collect(keys(sums)))
    return SparseCells([k[1] for k in keys_sorted],[k[2] for k in keys_sorted],
                       [k[3] for k in keys_sorted],[sums[k] for k in keys_sorted])
end

"""Sample each carrier position into its physical grid cell, without clipping."""
function _sample_carriers_to_cells(rng::AbstractRNG,hits::EventHits,
        carriers::Vector{Int64},grid::GridSpec,
        sigma_t_mm::Float32,sigma_l_mm::Float32)::Dict{NTuple{3,Int64},Int64}
    hits.frame==source_frame || throw(ArgumentError("Collection requires source-frame hits"))
    _validate_hit_vectors(hits.x_mm,hits.y_mm,hits.z_mm,hits.energy_MeV,
                          hits.time_ns,hits.particle_id,hits.hit_id)
    length(carriers)==length(hits.x_mm) || throw(DimensionMismatch("One count per hit"))
    all(>=(0),carriers) && all(x->isfinite(x)&&x>=0f0,(sigma_t_mm,sigma_l_mm)) ||
        throw(ArgumentError("Invalid carrier counts or widths"))
    bases=(grid.reference_mm[1]+grid.phase_mm[1],
           grid.reference_mm[2]+grid.phase_mm[2],grid.reference_mm[3])
    sums=Dict{NTuple{3,Int64},Int64}()
    for i in eachindex(carriers),_ in 1:carriers[i]
        x=hits.x_mm[i]+randn(rng,Float32)*sigma_t_mm
        y=hits.y_mm[i]+randn(rng,Float32)*sigma_t_mm
        z=hits.z_mm[i]+randn(rng,Float32)*sigma_l_mm
        all(isfinite,(x,y,z)) || throw(ArgumentError("Diffused coordinate is nonfinite"))
        key=(_direct_lattice_index(x,bases[1],grid.pitch_mm),
             _direct_lattice_index(y,bases[2],grid.pitch_mm),
             _direct_lattice_index(grid.drift_direction==1 ? z : 2bases[3]-z,
                                   bases[3],grid.bucket_width_mm))
        sums[key]=Base.Checked.checked_add(get(sums,key,Int64(0)),Int64(1))
    end
    return sums
end

"""
    collect_ion_pad_history(rng, hits, pairs, placement, grid; sigma_t_mm, sigma_l_mm)

Diffuse every initial ion from a source-frame hit directly into one pad and
one longitudinal bucket. Return `IonPadHistory` with an exact event-level
charge ledger. The source arrays, pair counts and grid are preserved.
"""
function collect_ion_pad_history(rng::AbstractRNG,hits::EventHits,
        pairs::Vector{Int64},placement::EventPlacement,grid::GridSpec;
        sigma_t_mm::Float32,sigma_l_mm::Float32)::IonPadHistory
    hits.event_id==placement.event_id || throw(ArgumentError("Event ID mismatch"))
    grid.drift_direction==1 || throw(ArgumentError("Ion bucket direction must be +1"))
    sums=_sample_carriers_to_cells(rng,hits,pairs,grid,sigma_t_mm,sigma_l_mm)
    total=foldl(Base.Checked.checked_add,pairs;init=Int64(0))
    return IonPadHistory(placement,grid,_sparse_cells(sums),total)
end

"""
    collect_photon_sipm_history(rng, hits, surviving, placement, grid;
        sigma_t_mm, sigma_l_mm, yield_per_e, J, acceptance, gas_model,
        initial_pairs)::PhotonSiPMHistory

Diffuse surviving electrons into EL holes and time buckets. For each occupied
cell, sample emitted VUV photons using `J`, then incident photons using the
analytic geometric acceptance. Assign each cell's photons to the parent
electrons' arrival bucket, matching the approved first timing approximation.
Return sparse incident-photon counts with an
exact emitted/incident ledger. A single event RNG drives diffusion and EL;
the caller supplies a separate stream for each gas branch.
"""
function collect_photon_sipm_history(rng::AbstractRNG,hits::EventHits,
        surviving::Vector{Int64},placement::EventPlacement,grid::GridSpec;
        sigma_t_mm::Float32,sigma_l_mm::Float32,
        yield_per_e::Float32,J::Float32,acceptance::Float32,
        gas_model::Symbol,initial_pairs::Int64)::PhotonSiPMHistory
    hits.event_id==placement.event_id || throw(ArgumentError("Event ID mismatch"))
    grid.drift_direction==-1 || throw(ArgumentError("Electron bucket direction must be -1"))
    gas_model in (:pure_xe,:xe_he15) || throw(ArgumentError("Unknown gas model"))
    0<=foldl(Base.Checked.checked_add,surviving;init=Int64(0))<=initial_pairs ||
        throw(ArgumentError("Survivors exceed initial pairs"))
    electron_cells=_sample_carriers_to_cells(rng,hits,surviving,grid,sigma_t_mm,sigma_l_mm)
    photons=Dict{NTuple{3,Int64},Int64}();emitted_total=Int64(0)
    for key in sort!(collect(keys(electron_cells)))
        emitted=sample_el_photons(rng,electron_cells[key],yield_per_e,J)
        incident=sample_incident_photons(rng,emitted,acceptance)
        emitted_total=Base.Checked.checked_add(emitted_total,emitted)
        incident>0 && (photons[key]=incident)
    end
    return PhotonSiPMHistory(placement,gas_model,grid,_sparse_cells(photons),
        initial_pairs,foldl(Base.Checked.checked_add,surviving;init=Int64(0)),
        emitted_total)
end

const ION_HISTORY_SCHEMA="nausicaa.ion_pad_history.v1"
const PHOTON_HISTORY_SCHEMA="nausicaa.photon_sipm_history.v1"

"""Write common event history fields and compressed sparse vectors to HDF5."""
function _write_history_fields(g::HDF5.Group,placement::EventPlacement,
        grid::GridSpec,cells::SparseCells)
    g["event_id"]=placement.event_id
    g["ion_drift_cm"]=placement.ion_drift_cm
    g["electron_drift_cm"]=placement.electron_drift_cm
    g["pitch_mm"]=grid.pitch_mm
    g["bucket_width_mm"]=grid.bucket_width_mm
    g["reference_mm"]=collect(grid.reference_mm)
    g["phase_mm"]=collect(grid.phase_mm)
    g["drift_direction"]=grid.drift_direction
    for (name,values) in (("pixel_ix",cells.pixel_ix),("pixel_iy",cells.pixel_iy),
                          ("bucket",cells.bucket),("count",cells.count))
        if isempty(values)
            g[name]=values
        else
            d=create_dataset(g,name,Int64,(length(values),);
                chunk=(min(length(values),4096),),shuffle=true,deflate=3)
            write(d,values)
        end
    end
    return nothing
end

"""Read and validate common event history fields from HDF5."""
function _read_history_fields(g::HDF5.Group)
    placement=EventPlacement(read(g["event_id"]),read(g["ion_drift_cm"]),
                             read(g["electron_drift_cm"]))
    reference=Tuple(read(g["reference_mm"]))::NTuple{3,Float32}
    phase=Tuple(read(g["phase_mm"]))::NTuple{2,Float32}
    grid=GridSpec(read(g["pitch_mm"]),read(g["bucket_width_mm"]),
                  reference,phase,read(g["drift_direction"]))
    cells=SparseCells(read(g["pixel_ix"]),read(g["pixel_iy"]),
                      read(g["bucket"]),read(g["count"]))
    return placement,grid,cells
end

"""
    write_ion_history(group::HDF5.Group, history::IonPadHistory)::Nothing

Write one versioned event group with the placement, physical grid, exact
initial-pair ledger and four compressed sparse cell columns. Empty cell
vectors are stored as empty HDF5 datasets. The parent file stores the source,
configuration, seed, code identity and completion marker. The caller owns and
closes the group; this function preserves the input history.
"""
function write_ion_history(g::HDF5.Group,h::IonPadHistory)::Nothing
    attributes(g)["schema"]=ION_HISTORY_SCHEMA
    _write_history_fields(g,h.placement,h.grid,h.cells)
    g["initial_pairs"]=h.initial_pairs
    return nothing
end

"""
    read_ion_history(group::HDF5.Group)::IonPadHistory

Read an owned ion history from a versioned event group. Reconstructors check
finite geometry, sorted unique sparse keys and exact charge conservation.
An incompatible schema or missing HDF5 field raises an error. The group stays
open under the caller's ownership.
"""
function read_ion_history(g::HDF5.Group)::IonPadHistory
    read(attributes(g)["schema"])==ION_HISTORY_SCHEMA || error("Ion history schema mismatch")
    placement,grid,cells=_read_history_fields(g)
    return IonPadHistory(placement,grid,cells,read(g["initial_pairs"]))
end

"""
    write_photon_history(group::HDF5.Group, history::PhotonSiPMHistory)::Nothing

Write one versioned gas-specific event group. Besides the common placement,
grid and sparse columns, record initial pairs, surviving electrons, emitted
photons and gas-model identity. Sparse counts are VUV photons incident on
SiPMs; the parent file records provenance and completion. Preserve the input
history and leave the HDF5 handle open for the caller.
"""
function write_photon_history(g::HDF5.Group,h::PhotonSiPMHistory)::Nothing
    attributes(g)["schema"]=PHOTON_HISTORY_SCHEMA
    _write_history_fields(g,h.placement,h.grid,h.cells)
    g["gas_model"]=String(h.gas_model)
    g["initial_pairs"]=h.initial_pairs
    g["surviving_electrons"]=h.surviving_electrons
    g["emitted_photons"]=h.emitted_photons
    return nothing
end

"""
    read_photon_history(group::HDF5.Group)::PhotonSiPMHistory

Read an owned gas-specific photon history from HDF5. Validate the schema,
grid, sparse keys and inequalities `surviving_electrons ≤ initial_pairs` and
`sum(cells.count) ≤ emitted_photons`. The caller owns and closes the group.
"""
function read_photon_history(g::HDF5.Group)::PhotonSiPMHistory
    read(attributes(g)["schema"])==PHOTON_HISTORY_SCHEMA || error("Photon history schema mismatch")
    placement,grid,cells=_read_history_fields(g)
    return PhotonSiPMHistory(placement,Symbol(read(g["gas_model"])),grid,cells,
        read(g["initial_pairs"]),read(g["surviving_electrons"]),
        read(g["emitted_photons"]))
end
