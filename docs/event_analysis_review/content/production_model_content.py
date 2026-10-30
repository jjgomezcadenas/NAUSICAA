"""Production detector-stage contracts for the event-analysis review guide."""

from html import escape


def production_page(section, table, code):
    flow = section('Production boundary', '''
<p>One source event produces one shared integer ionisation realisation. Its ions
reach the cathode pad grid; its electrons drift the complementary distance to
the anode. The electron branch is evaluated twice, for pure Xe and Xe–15% He.
The outputs are three independent, sparse HDF5 histories: ions per pad and
time bucket, and incident VUV photons per SiPM and time bucket for each gas.
The SiPM and electronics stage consumes those histories later.</p>
<p>The existing voxel, graph, backbone and truth-path implementation remains the
reconstruction core. This chapter specifies the new production functions and
data structures. Physical formulas and their provenance are in
<a href="../production_formulas.pdf">production_formulas.pdf</a>; the
detector/run dependency graph is in
<a href="../detector_response_pipeline_plan.md">the pipeline plan</a>.</p>
''', 'production-boundary')

    types = section('Concrete data types', table(
        ['Type', 'Fields and units', 'Invariant and role'], [
            ('<code>EventPlacement</code>',
             '<code>event_id::Int64</code>; <code>ion_drift_cm::Float32</code>; '
             '<code>electron_drift_cm::Float32</code>',
             'The two lengths sum to the run card’s chamber length. One placement is shared by all three branches.'),
            ('<code>GridSpec</code>',
             '<code>pitch_mm::Float32</code>; <code>bucket_width_mm::Float32</code>; '
             '<code>reference_mm::NTuple{3,Float32}</code>; '
             '<code>phase_mm::NTuple{2,Float32}</code>; '
             '<code>drift_direction::Int8</code>',
             'A physical lattice. The pilot uses the source-energy-weighted x/y centre, minimum source z for ions and maximum source z for electrons as references. The transverse phase is constant across an event and lies in [0,pitch). Bucket indices increase with arrival time.'),
            ('<code>SparseCells</code>',
             '<code>pixel_ix::Vector{Int64}</code>; <code>pixel_iy::Vector{Int64}</code>; '
             '<code>bucket::Vector{Int64}</code>; <code>count::Vector{Int64}</code>',
             'Equal-length owned arrays, sorted by cell key; unique keys, strictly positive counts. One row is one occupied pixel/time cell.'),
            ('<code>IonPadHistory</code>',
             '<code>placement::EventPlacement</code>; <code>grid::GridSpec</code>; '
             '<code>cells::SparseCells</code>; <code>initial_pairs::Int64</code>',
             'For this unbounded-grid pilot, sum(cells.count) equals initial_pairs exactly.'),
            ('<code>PhotonSiPMHistory</code>',
             '<code>placement::EventPlacement</code>; <code>gas_model::Symbol</code>; '
             '<code>grid::GridSpec</code>; <code>cells::SparseCells</code>; '
             '<code>initial_pairs::Int64</code>; <code>surviving_electrons::Int64</code>; '
             '<code>emitted_photons::Int64</code>',
             'sum(cells.count) is the incident-photon total; it cannot exceed emitted_photons. Gas identity is explicit.'),
        ]), 'production-types')

    funcs = section('Function interfaces and sampling contracts', table(
        ['Function', 'Input → output', 'Contract and required test'], [
            ('<code>sample_fano_pairs</code>',
             '<code>(rng, hits::EventHits, W_eV::Float32, F::Float32) → Vector{Int64}</code>',
             'One count per hit; a single draw seeds both branches. Bernoulli(μ) for μ&lt;1. Above that, the integer sampler targets mean μ and variance Fμ where feasible; the minimum integer variance is {μ}(1−{μ}). Test low occupancy, moments, replay and charge equality.'),
            ('<code>sample_attached_electrons</code>',
             '<code>(rng, pairs::Vector{Int64}, drift_ms::Float32, lifetime_ms::Float32) → Vector{Int64}</code>',
             'Independent binomial survival with p=exp(−drift_ms/lifetime_ms). Ion counts retain their original values. Test zero drift, mean, variance and complementary lengths.'),
            ('<code>sample_el_photons</code>',
             '<code>(rng, electrons::Int64, yield_per_e::Float32, J::Float32) → Int64</code>',
             'Conditional mean nₑY and variance nₑJY. Approximate Oliveira Figure 6 readings supplied by the user are J=0.10 at E/P=1.5, 0.05 at 2 and 0.01 at 3 kV/cm/bar. Intermediate fields require an explicit J. Test moments, zero input and overflow.'),
            ('<code>sample_incident_photons</code>',
             '<code>(rng, emitted::Int64, acceptance::Float32) → Int64</code>',
             'Binomial geometric collection. The optical acceptance is the existing analytic function. Test bounds, mean and variance.'),
            ('<code>collect_ion_pad_history</code>',
             '<code>(rng, hits, pairs, placement, grid; σT_mm, σL_mm) → IonPadHistory</code>',
             'Diffuse every ion directly into a physical pad and bucket; no display histogram or clipping window. Test exact sum, edge assignment and seed replay.'),
            ('<code>collect_photon_sipm_history</code>',
             '<code>(rng, hits, surviving, placement, grid; σT_mm, σL_mm, Y, J, A, gas_model, initial_pairs) → PhotonSiPMHistory</code>',
             'Diffuse electrons to holes, then sample emitted and incident photons in each occupied electron cell. Test gas labels, time order, optical ledger and cross-branch pair identity.'),
            ('<code>write/read_*_history</code>',
             '<code>(HDF5.Group, history) → nothing/history</code>',
             'Versioned, compressed sparse vectors with exact round trip and schema validation. Test malformed and incomplete groups.'),
        ]), 'production-functions')

    rng = section('Integer fluctuations and random streams', '''
<p>For an integer random variable of mean μ, the lowest attainable variance is
f(1−f), with f the fractional part of μ. The pilot sampler uses the adjacent
integer law at that lower bound and mixes it with a binomial law to attain the
requested sub-Poisson variance when possible. The same mechanism serves Fano
pair counts and the conditional EL count. It records when the bound overrides
the requested variance; tests cover the low-occupancy regime.</p>
<p>Each event has separate deterministic streams for placement, pairs, ion
diffusion, electron attachment, each gas’s electron diffusion and EL light,
and each grid phase. Streams are derived from a run seed, source-content
identity, event ID and stage label. Replaying a branch reads the same ordered
source hits and reproduces the common pair vector. A per-event checksum of that
vector is stored in all three outputs. Worker scheduling and changes to either
grid cannot change the pair draw.</p>
''', 'production-rng')

    files = section('Files and script boundary', '''
<p>The terminal command <code>python3.13 scripts/production/run_pilot.py</code>
presents the parameter form, saves five TOML cards under
<code>config/specs/generated/</code> and launches the Julia pilot. Use
<code>--defaults</code> for the complete default run or <code>--dry-run</code>
to inspect the resolved cards. The script-side resolver reads TOML, validates
cross-card constraints, computes canonical identities and chooses output
paths. It passes typed values to <code>EventAnalysis</code>; physics functions
do not parse TOML. The pilot produces <code>realisation.h5</code>, <code>ions.h5</code>,
<code>photons_pure_xe.h5</code> and <code>photons_xe_he15.h5</code> for a small
event subset. Each HDF5 file records schema, gas, units, exact card values,
source hash, seeds, event IDs, pair checksum and charge/light ledgers. Files
are written to temporary paths and marked complete after the event ledgers are
saved.</p>
<p>Each event group stores the four sparse vectors, the physical lattice and
its phase, placement, and ledger fields. Ion and photon files have different
schema identifiers. The run command generates run-specific formulas TeX/PDF
and a numerical HDF5 validation summary. The checked default pilot has a
<a href="../../latex_generated/production_pilot/26f5b7a9092738ed/formulas.pdf">generated formulas and results PDF</a>.</p>
''', 'production-files')

    pilot = section('Pilot acceptance', table(['Check', 'Acceptance'], [
        ('Common realisation', 'Identical pair checksums in ion, pure-Xe and Xe–He files for every event; counts independent of branch order and worker count.'),
        ('Ion ledger', 'Every sampled ion appears in exactly one pad/bucket cell; sum of sparse counts equals initial pairs.'),
        ('Electron/photon ledger', 'Surviving electrons ≤ initial pairs; incident photons ≤ emitted photons; all sparse counts positive and unique.'),
        ('Statistical tests', 'Independent numerical checks of Bernoulli, Fano and J moments, attachment and geometric binomial collection.'),
        ('HDF5 replay', 'Exact type/array/metadata round trip and rejection of incompatible schema or incomplete files.'),
        ('Runtime', 'Benchmark a short run for both samples and three streams; record wall time and peak memory before selecting a full-run event count.'),
    ]), 'production-pilot')

    return dict(id='production-model', title='8. Production detector model',
                summary='Shared pair sampling, branch-specific transport, sparse histories and pilot validation.',
                status='Pilot implemented · J anchor readings recorded',
                body=flow + types + funcs + rng + files + pilot)
