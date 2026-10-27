# Pixel area budget

Reproduce from the repository root:

```sh
julia --project=. tests/test_pixel_area.jl
julia --project=. scripts/build_pixel_area.jl
python scripts/present_pixel_area.py
```

Configuration: config/pixel_area.toml. Device dimensions are parsed directly
from LTSpice/variants/m4_x4_m3_x4/noise.cir, including multiplicity. Pitch and
pad width come from config/detector.toml. Julia saves numerical results,
rectangle coordinates and source hashes to analysis.h5. Python presents them.

The square 2 mm pitch gives 4 mm2 per pixel. The stricter
1 mm2 budget includes the pad and 50 um clearance, transistor envelopes,
a 50 x 50 um feedback/passive reserve, and explicit unimplemented allocations
for buffer/selection, control/routing and shared services. ADCs are shared
outside the pixel. No circuitry is placed under the pad keepout.

Channel area: 706.647 um2. CSA footprint with 10 um device
margins: 0.0752 mm2. Total nominal budget: 0.8752 mm2;
with 20 um device margins: 0.968 mm2. Readout allocations
are planning reserves, not measured circuit footprints or completed layouts.
Public GF180 rules supply a scale comparison, not DRC rules for the supplied
180 nm transistor card. Routing, extraction, die edges, ESD and packaging
will be evaluated with the actual process and array design.

Topmetal-S complete-chip areas are 4.72181 mm2 (Yang 2024) and
3.9 mm2 (Liang 2024). A complete Yang-sized common block
amortized over 64 pixels would be 0.0737782 mm2/pixel.
This is a reference-area comparison; it does not define the NAUSICA common
circuit. Two 0.4 mm edge bands geometrically accommodate the
12.8 mm2 common-area allocation in a 16 mm
square example array while preserving the 2 mm pad pitch.

Report and figures: latex_generated/pixel_area/. Source generated without
compiling or modifying the manuscript. Tests check multiplicity, physical
areas, containment, overlap rejection and dimensional scaling.
