"""Historical group-two proposal, superseded by the Petit-first assessment.

Retained as the source of the earlier review discussion. The active builder
uses petit_review_content.py; new interfaces follow the prerequisite review.
"""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[3]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

import re


def voxel_pages(plan, table, code):
    declarations = {}
    for name in ('VoxelGrid', 'VoxelCloud'):
        match = re.search(rf'^struct {name}\n.*?^end$', plan, re.M | re.S)
        if not match:
            raise ValueError(f'Missing proposed {name} declaration')
        declarations[name] = match.group()
    match = re.search(r'^voxelize_hits\(.*$', plan, re.M)
    if not match:
        raise ValueError('Missing proposed voxelize_hits signature')
    declarations['voxelize_hits'] = match.group()

    example_hits = table(['Input hit j', 'x, y, z · mm', 'Energy · MeV', 'Cell key'], [
        ('1', '(0.25, 0.25, 0.25)', '0.001', '(0, 0, 0)'),
        ('2', '(-0.25, 0.25, 0.25)', '0.002', '(-1, 0, 0)'),
        ('3', '(0.75, 0.25, 0.25)', '0.003', '(0, 0, 0)'),
        ('4', '(1.0, 0.25, 0.25)', '0.004', '(1, 0, 0)'),
    ])
    example_cloud = table(['Cloud index v', 'Cell key', 'Centre · mm', 'Energy · MeV'], [
        ('1', '(-1, 0, 0)', '(-0.5, 0.5, 0.5)', '0.002'),
        ('2', '(0, 0, 0)', '(0.5, 0.5, 0.5)', '0.004'),
        ('3', '(1, 0, 0)', '(1.5, 0.5, 0.5)', '0.004'),
    ])
    example_caption = ('<p class="caption">Illustrative 1 mm cells with origin '
                       '(0, 0, 0). Energies are displayed as decimal values; '
                       'stored sums use Float32.</p>')

    pages = [
        dict(id='voxel-grid', name='VoxelGrid', kind='Data type', number='04',
             summary='Cubic cell size and grid origin.',
             intro='Define the spatial grid used to aggregate hits before building the track graph. The size sets the analysis resolution; the origin sets the grid alignment relative to the trace.',
             declaration=declarations['VoxelGrid'],
             fields=table(['Field', 'Julia type / units', 'Proposed meaning'], [
                 ('<code>size_mm</code>', '<code>Float32</code> · mm', 'Finite, positive side length of each cubic cell.'),
                 ('<code>origin_mm</code>', '<code>NTuple{3,Float32}</code> · mm', 'Finite (x, y, z) corner of cell (0, 0, 0), expressed in the input hit frame.'),
             ]),
             example=code('grid = VoxelGrid(1f0, (0f0, 0f0, 0f0))') + table(
                 ['x · mm', 'x cell index', 'Cell interval · mm'], [
                     ('−0.25', '−1', '[−1, 0)'), ('0', '0', '[0, 1)'),
                     ('0.75', '0', '[0, 1)'), ('1', '1', '[1, 2)'),
                 ]) + '<p class="caption">The same assignment rule applies independently to y and z. The 1 mm size is an example for reviewing the interface.</p>',
             contract='<ul><li>Immutable scalar and tuple fields.</li><li>Cell index on x: <code>floor(Int, (x_mm - origin_mm[1]) / size_mm)</code>.</li><li>Cells include their lower boundary and exclude their upper boundary.</li><li>Invalid size or origin raises <code>ArgumentError</code>.</li></ul>',
             tests=[('Size and origin validation', 'Accept finite positive sizes and finite origins; reject zero/negative size, NaN and Inf.'),
                    ('Known cell assignments', 'Check positive and negative positions, exact boundaries and adjacent Float32 values in voxelize_hits.'),
                    ('Grid translation', 'Translate hits and grid together; preserve cell assignments within the agreed precision convention.')],
             decisions=[('Resolution setting', 'Choose voxel size through the backbone stability study. Sensor pitch remains a separate detector setting.'),
                        ('Origin convention', 'Use a cell corner as the origin. An origin of (0, 0, 0) places the true vertex on a grid corner for vertex-centred hits; origin shifts support the planned grid-phase study.')],
             docstring='''    VoxelGrid(size_mm::Float32, origin_mm::NTuple{3,Float32})

Specify an axis-aligned cubic grid for hit aggregation.

# Arguments
- `size_mm`: finite positive side length in mm.
- `origin_mm`: finite (x, y, z) corner of cell (0, 0, 0), in mm and in
  the same coordinate frame as the input EventHits.

# Cell convention
The integer key on each axis is floor((coordinate - origin) / size).
Each interval includes its lower boundary and excludes its upper boundary.
Keys can be negative. The geometric centre is origin + (key + 0.5) * size.

# Validation and ownership
The constructor stores immutable Float32 values. Throw ArgumentError for
nonfinite values or a nonpositive size.

# Example
```julia
grid = VoxelGrid(1f0, (0f0, 0f0, 0f0))
@assert grid.size_mm == 1f0
@assert grid.origin_mm == (0f0, 0f0, 0f0)
```
'''),
        dict(id='voxel-cloud', name='VoxelCloud', kind='Data type', number='05',
             summary='Occupied grid cells and their summed energies.',
             intro='Represent one event as occupied grid cells. Each voxel has a geometric centre and the sum of its deposited energies.',
             declaration=declarations['VoxelCloud'],
             fields=table(['Field', 'Julia type / units', 'Proposed meaning'], [
                 ('<code>event_id</code>', '<code>Int64</code>', 'Source event ID, preserved from EventHits.'),
                 ('<code>frame</code>', '<code>CoordinateFrame</code>', 'Coordinate system shared by the centres and grid origin; preserved from EventHits.'),
                 ('<code>grid</code>', '<code>VoxelGrid</code>', 'Grid size and origin used to form these voxels.'),
                 ('<code>x_mm</code>, <code>y_mm</code>, <code>z_mm</code>', '<code>Vector{Float32}</code> · mm', 'V geometric cell centres, in ascending lexicographic (ix, iy, iz) cell-key order.'),
                 ('<code>energy_MeV</code>', '<code>Vector{Float32}</code> · MeV', 'V finite, nonnegative energy sums, one per occupied cell.'),
             ]),
             example=example_hits + example_cloud + example_caption,
             contract='<ul><li>Centre and energy vectors have the same positive length V.</li><li>Each occupied cell is represented once, in canonical cell-key order.</li><li>Each occupied cell contains at least one hit. Zero-energy occupied cells remain represented.</li><li>Construction copies the arrays; the stored grid specification is immutable.</li><li>Shape mismatches raise <code>DimensionMismatch</code>; invalid values, grid centres or ordering raise <code>ArgumentError</code>.</li></ul>',
             tests=[('Four-hit worked example', 'Three cells with the stated centres and energy sums.'),
                    ('Centre geometry and ordering', 'Validate finite grid centres, unique cells and lexicographic cell-key order, with Float32 centre rounding accounted for.'),
                    ('Field and shape contracts', 'Check lengths, finite values, energy sign and concrete field types.'),
                    ('Array ownership', 'Changing constructor inputs preserves the constructed cloud.')],
             decisions=[('Coordinate representation', 'Use separate x_mm, y_mm and z_mm vectors, following EventHits. Carry frame and grid so the centre coordinates have an explicit interpretation.'),
                        ('Original hits', 'Keep EventHits alongside VoxelCloud for the later hit-to-backbone distance calculation.')],
             docstring='''    VoxelCloud(event_id::Int64, frame::CoordinateFrame, grid::VoxelGrid,
               x_mm::Vector{Float32}, y_mm::Vector{Float32},
               z_mm::Vector{Float32}, energy_MeV::Vector{Float32})

Store the occupied grid cells of one event and their deposited energy sums.

# Fields
`event_id` and `frame` identify the event and its coordinate system.
`grid` records the grid size and origin in that frame. Position vectors give
the geometric centres of V cells in mm. `energy_MeV` gives their deposited
energy sums. All four physical vectors have length V and Float32 elements.

Cells follow ascending lexicographic (ix, iy, iz) order. Each cell appears once.
Each occupied cell contains at least one hit, including cells with zero energy.

# Validation and ownership
V is positive. Centres are finite grid-centre positions at Float32
precision; energies are finite and nonnegative.
Construction copies each array separately. Throw DimensionMismatch for unequal
centre/energy lengths; throw ArgumentError for invalid values or cell ordering.
Agreement with the originating hits is tested in voxelize_hits.

# Example
For 1 mm cells at origin (0, 0, 0), hits at x=[0.25, -0.25, 0.75, 1.0],
all with y=z=0.25 mm, occupy cells with x centres [-0.5, 0.5, 1.5] mm.
'''),
        dict(id='voxelize', name='voxelize_hits', kind='Pure function', number='06',
             summary='Assign hits to cells and sum their deposited energies.',
             intro='Bin a hit cloud on a 3D grid, summing deposited energies in each occupied cell. The pipeline supplies vertex-centred hits and an explicit grid specification.',
             declaration=declarations['voxelize_hits'],
             fields=table(['Argument / return', 'Julia type', 'Proposed contract'], [
                 ('<code>hits</code>', '<code>EventHits</code>', 'One valid event with its original hit order.'),
                 ('<code>grid</code>', '<code>VoxelGrid</code>', 'Grid definition in hits.frame.'),
                 ('Return', '<code>VoxelCloud</code>', 'Same event ID and frame; grid definition, cell centres and energy sums.'),
             ]),
             example=example_hits + code('''hits = EventHits(Int64(42), vertex_frame,
    Float32[0.25, -0.25, 0.75, 1.0],
    Float32[0.25, 0.25, 0.25, 0.25],
    Float32[0.25, 0.25, 0.25, 0.25],
    Float32[0.001, 0.002, 0.003, 0.004],
    Float32[0, 0, 0, 0], Int32[1, 1, 1, 1], Int32[0, 1, 2, 3])
cloud = voxelize_hits(hits, VoxelGrid(1f0, (0f0, 0f0, 0f0)))''') + example_cloud + example_caption,
             contract='<ol><li>Validate the event and grid.</li><li>Compute each hit’s integer (ix, iy, iz) cell key.</li><li>Sort occupied keys lexicographically to assign voxel indices.</li><li>Sum each cell’s energy in Float32, visiting hits in (particle_id, hit_id) order.</li><li>Calculate the geometric cell centres.</li><li>Return independently owned arrays, preserving the input values.</li></ol>',
             tests=[('Cell assignment and boundaries', 'Known 3D cells, negative coordinates, exact edges and adjacent Float32 values; repeat with nonzero origins and non-unit sizes.'),
                    ('Complete binning', 'Every hit contributes once to its expected cell; each output cell is occupied, including zero-energy cells.'),
                    ('Energy conservation', 'Compare each voxel and event total against independently grouped Float64 sums; report absolute errors and their fraction of 0.0125 MeV.'),
                    ('Permutation reproducibility', 'Reordering hits preserves cell order and energy sums.'),
                    ('Frame and translations', 'Preserve event ID/frame. A common hit/grid translation preserves cell keys and energies at the agreed precision.'),
                    ('Purity and validation', 'Check input snapshots and output ownership; reject invalid mutable hit contents and unrepresentable cell keys or centres.')],
             decisions=[('Energy accumulation', 'Propose Float32 sums with canonical hit-identity ordering. Validate the error on fixtures and representative NEXUS events against the 12.5 keV resolution scale before accepting this reduction.'),
                        ('Numerical cell boundaries', 'Review the lower-inclusive rule using stored Float32 inputs. Boundary tests will determine the arithmetic needed to honour it for shifted grids and general voxel sizes.'),
                        ('Graph threshold', 'Voxelisation includes every hit and occupied cell. GraphSpec will supply the subsequent voxel-energy threshold.')],
             docstring='''    voxelize_hits(hits::EventHits, grid::VoxelGrid) -> VoxelCloud

Bin one event's hits into cubic cells and sum their deposited energies.

# Inputs and output
`hits` supplies finite Float32 positions in mm, energies in MeV and unique
(particle_id, hit_id) identities. `grid` defines the grid in hits.frame.
Return a VoxelCloud carrying hits.event_id, hits.frame and grid.

# Assignment and ordering
On each axis, the cell key is floor((coordinate - origin) / size). Lower cell
boundaries are inclusive and upper boundaries exclusive. Sort occupied keys
lexicographically. Store geometric cell centres as Float32 positions in mm.

# Energy reduction
Every hit contributes its energy once, including zero-energy deposits.
Within a cell, visit hits in (particle_id, hit_id) order and accumulate in
Float32. This fixes summation order under input permutations. Validate errors
against independent Float64 reference sums before accepting this reduction.

# Ownership and errors
Inputs retain their values; output arrays are independently owned.
Validate current hit-array contents on entry. Throw DimensionMismatch for
unequal hit-vector lengths. Throw ArgumentError for invalid values or
identities, or cell keys/centres outside their representable range.

# Worked example
With 1 mm cells at origin (0, 0, 0), x=[0.25, -0.25, 0.75, 1.0] mm and
y=z=0.25 mm, return x centres [-0.5, 0.5, 1.5] mm.
For energies [0.001, 0.002, 0.003, 0.004] MeV, the cell sums are approximately
[0.002, 0.004, 0.004] MeV, with Float32 rounding.
'''),
    ]
    for page in pages:
        page['proposed'] = True
    return pages
