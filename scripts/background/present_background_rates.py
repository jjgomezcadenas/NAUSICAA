#!/usr/bin/env python3
"""Page of the radiogenic background of the copper shield (prepare_background_rates.jl).

  python scripts/background/present_background_rates.py [results/background]

Reads rates.h5 and writes index.html in the same folder, with the tables of
the ITACA revisited budget: copper activity, window fractions, efficiencies
after the window, and the background per year and per tonne-year of iItaca and
eITACA. Nothing is recalculated.
"""
import sys
from pathlib import Path
import h5py
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'iitaca'))
from present_iitaca_detector import STYLE, table, text

ROOT = Path(__file__).resolve().parents[2]
DET = ('iitaca', 'eitaca'); ISO = (('bi214', '²¹⁴Bi'), ('tl208', '²⁰⁸Tl'))


def sig(x, n=3):
    return f'{x:.{n}g}'


def pm(x, e, n=2):
    return f'{x:.{n}g} ± {e:.{n - 1 if n > 1 else 1}g}'


def pct(x, e=None):
    return f'{100 * x:.2f} %' if e is None else f'{100 * x:.2f} ± {100 * e:.2f} %'


def main(folder):
    out = Path(folder)
    with h5py.File(out / 'rates.h5') as f:
        d = {}
        for k in DET:
            g = f['detectors'][k]; d[k] = {n: (text(g[n][()]) if g[n].dtype.kind in 'SO' else g[n][()]) for n in g if isinstance(g[n], h5py.Dataset)}
            for iso, _ in ISO:
                h = g[iso]; d[k][iso] = {}
                for grp in ('window', 'selection'):
                    d[k][iso][grp] = {n: (text(h[grp][n][()]) if h[grp][n].dtype.kind in 'SO' else h[grp][n][()]) for n in h[grp]}
                for n in h:
                    if isinstance(h[n], h5py.Dataset): d[k][iso][n] = h[n][()]
        ref = {n: f['reference/itaca_revisited'][n][()] for n in f['reference/itaca_revisited']}
        roi = text(f['roi'][()]); count_file = text(f['tl208_count_toml'][()])
    L = {k: d[k]['label'] for k in DET}

    budget = table(['Detector', 'Shield', 'Copper', 'Mass (t)', 'Xenon (t)', '²¹⁴Bi (mBq)', '²⁰⁸Tl (mBq)', '²¹⁴Bi decays/yr', '²⁰⁸Tl decays/yr'],
                   [[L[k], d[k]['shield'], f"{d[k]['copper_thickness_cm']:.0f} cm, inner radius {d[k]['copper_inner_radius_cm']:.1f} cm, cavity {d[k]['copper_cavity_length_cm']:.1f} cm",
                     f"{d[k]['copper_mass_kg'] / 1000:.1f}", f"{d[k]['xenon_mass_kg'] / 1000:.2f}",
                     sig(1000 * d[k]['bi214']['activity_Bq']), sig(1000 * d[k]['tl208']['activity_Bq']),
                     sig(d[k]['bi214']['decays_per_year']), sig(d[k]['tl208']['decays_per_year'])] for k in DET])
    activities = table(['Material', '²¹⁴Bi (Bq/kg)', '²⁰⁸Tl (Bq/kg)', 'Source'],
                       [['Copper (ICS)', sig(d['iitaca']['bi214']['specific_activity_Bq_kg'], 2), sig(d['iitaca']['tl208']['specific_activity_Bq_kg'], 2), 'ITACA revisited, Table "Specific activities" (NEXT radiopurity measurements)']])
    wrows = []
    for k in DET:
        for iso, name in ISO:
            w = d[k][iso]['window']; cl = [text(c) for c in w['classes']]; fr = w['fraction']; er = w['fraction_error']; ct = w['counts']
            i1, im = cl.index('gamma'), cl.index('multiple')
            wrows.append([L[k], name, f"{w['generated_decays']:.0e}", int(w['saved_events']), f"{fr[i1]:.3g} ± {er[i1]:.2g} ({ct[i1]})",
                          f"{fr[im]:.2g} ({ct[im]})", f"{ct[cl.index('electron')]}", pm(d[k][iso]['window']['events_per_year'], d[k][iso]['window']['events_per_year_error'])])
    window = table(['Detector', 'Isotope', 'Decays', 'Saved in the window', 'One entering gamma, f_w (events)', 'Several entering particles, rejected (events)', 'One entering electron (events)', 'Events/yr in the window, one gamma'], wrows)
    erows = []
    for k in DET:
        for iso, name in ISO:
            s = d[k][iso]['selection']
            erows.append([L[k], name, int(s['saved_window_events']), int(s['track_selected_roi']), pct(s['eps_track'], s['eps_track_error']),
                          f"{pct(s['eps_net'], s['eps_net_error'])} ({int(s['net_events'])} of {int(s['net_sample'])})", f"{s['eps_after_window']:.3g} ± {s['eps_after_window_error']:.2g}", s['note']])
    eff = table(['Detector', 'Isotope', 'Saved window events', 'Track-selected in the ROI', 'ε_track', 'ε_net at 80 % signal', 'ε after the window', 'Source'], erows)
    srows = [[name] + [pm(d[k][iso]['events_per_year'], d[k][iso]['events_per_year_error']) for k in DET] + [f"{d['eitaca'][iso]['events_per_year'] / d['iitaca'][iso]['events_per_year']:.2f}"] for iso, name in ISO]
    srows.append(['<b>Total</b>'] + [f"<b>{pm(d[k]['total_events_per_year'], d[k]['total_events_per_year_error'])}</b>" for k in DET] + [f"<b>{d['eitaca']['total_events_per_year'] / d['iitaca']['total_events_per_year']:.2f}</b>"])
    per_year = table(['Isotope'] + [f'{L[k]} (events/yr)' for k in DET] + ['Ratio eITACA / iItaca'], srows)
    trows = [[name] + [pm(d[k][iso]['events_per_tonne_year'], d[k][iso]['events_per_tonne_year_error']) for k in DET] for iso, name in ISO]
    trows.append(['<b>Total</b>'] + [f"<b>{pm(d[k]['total_events_per_tonne_year'], d[k]['total_events_per_tonne_year_error'])}</b>" for k in DET])
    per_tonne = table(['Isotope'] + [f'{L[k]} (events / t yr)' for k in DET], trows)
    ir = ref['total_after_roi_per_tonne_year']; irs = ref['bi214_single_track_per_year']; irt = ref['tl208_single_track_per_year']
    up = '../../'
    page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Radiogenic background</title><style>{STYLE}</style></head><body><header><a href="{up}notes/html/index.html">Documentation index</a> · <a href="{up}notes/html/background/index.html">Method</a> · <a href="{up}notes/html/iitaca/index.html">iItaca</a> · <a href="{up}notes/html/eitaca/index.html">eITACA</a><h1>²¹⁴Bi and ²⁰⁸Tl from the copper shield in the ROI · iItaca and eITACA</h1><p>Decays generated in the inner copper shield, followed through the materials and the gas; events in the 2420–2490 keV window, the track selection, the network at 80 % signal and the ROI ({roi}). Format of the ITACA revisited budget.</p></header><main>
<section class="panel"><h2>Radioactivity budget</h2>{activities}{budget}<p>Both detectors: fiducial gas of radius 2.0 m and length 1.5 m at 16 bar; the copper is 20 cm thick on every face, around the xenon cavity and the passive layers of the shielding cards (<code>config/background/</code>). The pressure vessel (Ti Gr5, about 10 t, 9 Bq of ²¹⁴Bi) is outside the copper and not simulated: behind 20 cm it adds at most 5 % of the ICS ²¹⁴Bi and 1 % of its ²⁰⁸Tl (E₂ attenuation, oblique paths included). The other components of ITACA revisited (fibre barrel, field-cage rings, SiPM boards, EL structure, cathode mesh; for iItaca also the light guides and MARS) are not simulated; in that paper they add 5 % (²¹⁴Bi) and 23 % (²⁰⁸Tl, two thirds from 270 kg of HDPE light guides at 8 µBq/kg) to the ICS.</p></section>
<section class="panel"><h2>Window</h2><p>Fraction of the decays with 2420–2490 keV deposited in the fiducial gas, weighted to the uniform distribution in the copper (depth sampled as e<sup>−z/2.9 cm</sup>), by the particles that entered the gas. Events with several entering particles leave two separate deposits and are rejected.</p>{window}</section>
<section class="panel"><h2>After the window</h2><p>Per saved window event of the gamma productions of each detector: the track selection in the ROI, then the network at 80 % signal. The topology of a window event does not depend on where the gamma came from, so these efficiencies apply to the copper decays.</p>{eff}</section>
<section class="panel"><h2>Background</h2>{per_year}{per_tonne}<p>Errors: window fraction (weights), track selection and network (binomial), in quadrature. The specific activities carry no error; the rates scale with them. ITACA revisited, for comparison (all components, 15 cm of copper at normal-incidence attenuation, 2400–2500 keV, resolution {ref['resolution_fwhm_percent']:.1f} % FWHM, 1 t): after the single-track cut ²¹⁴Bi {irs[0]:.2f} (ion) / {irs[1]:.2f} (electron) and ²⁰⁸Tl {irt[0]:.2f} / {irt[1]:.1f} events per year; after topology and the ROI {ir[0]:.2f} events per tonne-year with the ion track, {ir[1]:.1f} for a conventional pure-xenon HPXeEL.</p></section>
<section class="panel"><h2>Files</h2><p><code>results/background/rates.h5</code> (prepare_background_rates.jl) from the pilot window tables <code>results/background/&lt;shield&gt;/pilot/&lt;isotope&gt;/window.h5</code>, the iItaca gamma-window results, the eITACA results and the ²⁰⁸Tl window count <code>{count_file}</code>.</p></section>
</main></body></html>'''
    (out / 'index.html').write_text(page); print(out / 'index.html')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else ROOT / 'results/background')
