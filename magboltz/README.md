# Xe–He electron transport with Magboltz

The driver in this directory prepares input for CERN's Fortran Magboltz,
parses electron drift and diffusion, and writes a numerical HDF5 table.
`present_magboltz.py` makes plots and an HTML report from that table.

The source downloaded from CERN on 2026-11-15 is
`https://magboltz.web.cern.ch/magboltz/magboltz-11.21.f`.
Its SHA-256 is
`efa3244aa7c033872864997ee005b50d3732a59080e3957d13be9bf6d9569cbc`.
The CERN page labels it 11.21; the file header labels it 11.20. Both
identifiers belong in result provenance. The source has a copyright notice
and no explicit redistribution terms, so `source/` remains local.

Build on macOS:

```sh
mkdir -p magboltz/build
gfortran -O2 -std=legacy -ffixed-line-length-none \
  -fallow-argument-mismatch magboltz/source/magboltz-11.21.f \
  -o magboltz/build/magboltz-11.21
```

The input format is taken from this source, especially `SUBROUTINE READIN`.
Gas identifiers are Xe = 7 and He-4 = 3. Temperature is supplied in Celsius,
pressure in Torr, field in V/cm, and the gas fractions in percent. The pilot
card uses `IPEN=0`, `ITHRM=1`, `EFINAL=0`, and zero magnetic field.
For `ITHRM=1`, `SETUPT` multiplies the `NMAX` input by 40 million; the pilot's
`NMAX=2` requests 80 million collisions.

The physical pressure and the ideal-gas-equivalent pressure are both recorded.
The Xe–15% He correction uses the number density implied by the Geant4
material implementation, as selected by the user. The plain-pressure pilots
are separate comparison points.

The approved workflow and reference data are in
`docs/handoffs/magboltz_xe_he.md`. Results belong in
`results/transport/magboltz_xe_he_290K/`.

The first plain-pressure pilot completed in 18.07 s with 80 million
collisions. Its reported drift-velocity error is 7.12%; 4,654 collisions
landed in the final energy bin, above the source's recommended limit of 500.
`pilot_refined.toml` prepares a 400-million-collision check with a 1.5 eV
energy limit. It completed in 76.75 s: velocity 0.4152 mm/µs ± 3.20%,
transverse diffusion 1.9261 mm/√m ± 4.84%, and longitudinal diffusion
2.7859 mm/√m ± 5.69%. The final energy bin contained zero collisions.

The corresponding plain-pressure points at 15 and 16 bar used the same gas,
temperature, field, energy limit and 400-million-collision setting:

| Pressure (bar) | Drift (mm/µs) | Transverse (mm/√m) | Longitudinal (mm/√m) | Time (s) |
|---:|---:|---:|---:|---:|
| 15 | 0.6236 ± 2.75% | 2.2296 ± 4.49% | 3.6766 ± 5.55% | 81.58 |
| 16 | 0.5322 ± 3.68% | 2.1146 ± 4.06% | 3.4569 ± 4.16% | 80.75 |
| 17.5 | 0.4152 ± 3.20% | 1.9261 ± 4.84% | 2.7859 ± 5.69% | 76.75 |

All three runs have zero collisions in the final energy bin.

The Geant4-material density correction gives:

| Physical pressure (bar) | Equivalent ideal pressure (bar) | Drift (mm/µs) | Transverse (mm/√m) | Longitudinal (mm/√m) | Time (s) |
|---:|---:|---:|---:|---:|---:|
| 15 | 16.4184 | 0.5167 ± 4.41% | 1.9905 ± 4.66% | 3.6071 ± 8.64% | 81.54 |
| 16 | 17.6241 | 0.4219 ± 3.38% | 1.9137 ± 3.89% | 2.9137 ± 8.11% | 80.86 |
| 17.5 | 19.4617 | 0.3442 ± 4.55% | 1.8195 ± 4.71% | 2.6577 ± 11.42% | 80.17 |

These runs also have zero collisions in the final energy bin. The mixture
number density follows `ItacaChamber.cc`; the Xe isotope molar mass uses the
same abundances and their integer mass numbers for the pressure conversion.

At 350 V/cm, with the same correction and collision settings:

| Physical pressure (bar) | Drift (mm/µs) | Transverse (mm/√m) | Longitudinal (mm/√m) | Time (s) |
|---:|---:|---:|---:|---:|
| 15 | 0.7057 ± 3.46% | 2.2867 ± 5.28% | 3.8611 ± 5.65% | 82.16 |
| 16 | 0.5684 ± 2.67% | 2.0115 ± 2.46% | 3.4227 ± 4.61% | 80.79 |
| 17.5 | 0.4807 ± 4.99% | 1.7951 ± 5.12% | 3.0758 ± 10.12% | 76.51 |

All three 350 V/cm runs have zero collisions in the final energy bin.
