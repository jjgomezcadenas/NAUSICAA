# Event display: zoom to the main track — handoff

Date: 2026-11-15. Branch: `new-config`. Requested by the user; written by the
agent of the 2026-11-14/15 session.

## Goal

Add a button to the event viewer that zooms every view into the main track.
The zoomed view must contain the main track in full. It may leave gamma
deposits outside. A second press returns to the full frame.

User's words: "the frame is huge and often details of the track are tiny. We
should add a button that zooms into the track (even if it leaves a gamma
outside), containing the main track in full."

## Rules that apply (read `AGENTS.md` first)

- No production run, cluster job or long computation without the user's OK.
  Regenerating a report from a saved hit file (seconds) is part of this task.
- No destructive step and no push without the user's OK. Commit only when the
  user says "commit".
- Julia calculates and saves numbers; Python makes the page. The viewer's
  JavaScript draws what the saved data contain; it does not reconstruct.
- Follow `notes/html/event_analysis/plots_viewer.html` and preserve its
  reviewed wording when extending it.
- Do not change the fast trigger, the hit file or the run cards.

## The viewer today

| Item | Path |
|---|---|
| Presentation data (Julia): histograms, projections, graph paths, blobs | `scripts/production/prepare_hit_presentation.jl`, `hit_projections` in `scripts/production/track_controls.jl` |
| Page and viewer (Python + inline JavaScript) | `scripts/production/present_tracks.py` |
| Viewer check without a browser | `scripts/validation/check_hit_viewer.js` (`node scripts/validation/check_hit_viewer.js <presentation folder>`) |
| Documentation | `notes/html/event_analysis/plots_viewer.html` |
| Reports to regenerate | `results/readout/bb0nu_itaca_pilot_10x100/vh1_hits/presentation/` and `single_electrons_itaca_pilot_10x100/...` (linked from `notes/html/index.html`, row 7) |

How a report is made (seconds each; `NAUSICAA_DATA` is not needed):

```sh
JULIA_DEPOT_PATH="$PWD/work/julia-depot:$HOME/.julia" julia --startup-file=no --project=. \
  scripts/production/prepare_hit_presentation.jl results/readout/<run_label>/vh1_hits/hits.h5
work/event-viewer-venv/bin/python scripts/production/present_tracks.py \
  results/readout/<run_label>/vh1_hits/presentation/controls.h5
```

`prepare_hit_presentation.jl` refuses to overwrite an existing
`presentation/controls.h5` ("Preserve existing presentation data"). To
regenerate, move the old `presentation` folder aside (for example into
`work/`) or pass another output folder as third argument; do not delete it
without the user's OK.

What the viewer shows, per displayed event (one in 50, `example_every` in
`[track_analysis]`):

- Two tabs, itrack and etrack. Each has three energy heatmaps (XY, XZ, YZ) and
  a 3D view with the graph paths and the two endpoint spheres.
- Since 2026-11-15, for hit files of the fixed-frame selection (schema
  `itaca.reconstructed_hits.v2`), every view covers the whole fixed frame: a
  cube of 160 mm, the same for both tracks. `prepare_hit_presentation.jl`
  builds the X, Y and Z cell centres over the cube (`covering_lattice`,
  `extend_centres`) and passes them to `hit_projections`. This was done so
  that all views and events share one scale; before, Z covered the cube and X,
  Y only the occupied cells, which distorted the views.
- In the JavaScript, the axis ranges are computed once per event from the cell
  edges of both tracks plus 12.5% on each side (`ranges` in `receiveEvent`),
  and all heatmap axes are `fixedrange: true`, `constrain: 'domain'`, with
  `scaleanchor` for equal scale. Zoom and pan are disabled on purpose; the 3D
  view only rotates, and a `plotly_relayout` handler restores the camera
  distance.
- The event data are separate files, `viewer/event_<id>.js`, loaded by a
  script tag; each calls `receiveEvent(data)`.

Positions are vertex-relative mm. The itrack may hold a few noise hits
anywhere in the cube (about 4 per event at the 3 keV hit cut; the display cut
is `ion_hit_cut_keV` of `[track_analysis]`, 5 keV in the pilot cards).

## What is needed

1. **The box of the main track, per track and event, saved by the Julia
   step.** The viewer has no such box today. Two candidate definitions; the
   user has to choose, so ask:
   - reconstructed: the bounds of the voxels of the largest graph component of
     that track (`tracks[1]` in `prepare_hit_presentation.jl`; the components
     are already computed there). Available for any hit file, and it is what
     the reconstruction calls the main track. If the track is broken it shows
     only the largest piece.
   - Monte Carlo: the main-track box of the fast trigger (primary electrons
     and their electron descendants with no gamma in the chain). The hit file
     stores `frame_*` (main track, contained gamma deposits and a 20 mm
     margin, limited to the cube), not the main-track box itself, so this
     would need a new column from the production.
   The session agent's suggestion is the reconstructed box, with one common
   zoom box per event: the union of the itrack and etrack boxes, so both tabs
   keep the same scale.
2. **A margin** around that box (suggest 10 mm or 10% of the largest side,
   whichever is larger) and **equal scale in X, Y and Z**: the zoom region
   should be a cube whose side is the largest side of the box plus margins,
   centred on the box and limited to the frame.
3. **A button** in the viewer's control row (next to "Show connecting paths"),
   "Zoom to track" / "Full frame". It switches the `ranges` used by the three
   heatmaps and the 3D scene of both tabs. The state is kept when the event
   changes. Default: ask the user (full frame today).
4. **The heatmaps are not recomputed**: only the axis ranges change, so cells
   outside the zoom region are simply not shown. Colour scale: keep the
   current per-track maximum.
5. **Checks**: extend `scripts/validation/check_hit_viewer.js` to press the
   button and verify that the ranges of all eight plots change to the zoom
   region, that the region contains the saved main-track box, and that a second
   press restores the full-frame ranges. Add a Julia test for the box
   calculation if it is a new function.
6. **Documentation**: one or two items in `plots_viewer.html`.
7. **Regenerate** the two pilot reports and check events with a small track in
   a corner of the cube (events near the detector boundary have the track
   off-centre by up to 57 mm).

## Things to watch

- `present_tracks.py` writes the selected events' data with `.tolist()`; the
  matrices are indexed `[second axis][first axis]` as Plotly expects. Do not
  transpose them.
- `hit_projections` looks hit Z values up in the list of Z centres by exact
  Float32 equality. Do not recompute Z centres with a different formula.
- Reports of the earlier productions (schema v1, no fixed frame) must still be
  produced as before; the button should work there too if the box is the
  reconstructed one.
- The viewer has no build step and must keep working from `file://`.

## State of the tree

`new-config` has uncommitted changes by the session agent (fixed frame moved
inside the detector, EL exclusion cut, full-frame display, pilot cards at
160 mm, index link). Do not revert them; if they are still uncommitted when
you start, ask the user before touching the same files
(`prepare_hit_presentation.jl`, `present_tracks.py`, `track_controls.jl`).

## Report back

To the user, with: the definition of the main-track box used, the zoom rule,
a list of the events checked, the test results, and the regenerated report
links.
