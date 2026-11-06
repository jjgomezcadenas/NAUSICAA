#!/usr/bin/env python3
"""Present fresh-output CSA timing measurements and the saved scaling estimate."""
from pathlib import Path
import sys,tomllib
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts/documentation'))
from build_report import build_report

def main():
    path=Path(tomllib.loads((ROOT/'results/ion_csa_nonlinear/scaling/latest.toml').read_text())['summary'])
    data=tomllib.loads(path.read_text());rows=[r for r in data['measurements'] if not r['warmup']]
    old=tomllib.loads((ROOT/'results/elm_matched_production/run_timings.toml').read_text())['csa_wall_s']
    rate=data['seconds_per_event'];projected=data['estimated_8000_event_work_s']
    out=ROOT/'latex_generated/ion_readout/ion_csa_nonlinear/scaling';out.mkdir(parents=True,exist_ok=True)
    table='\n'.join(f"{r['events']} & {r['work_without_presentation_s']:.3f} & {r['calibration_setup_s']:.3f} & {r['candidate_histories']} & {r['spice_fallback_histories']} "+r'\\' for r in rows)
    narrative=(f"Fresh outputs, after a 10-event warm-up: {rows[0]['work_without_presentation_s']:.3f} s for {rows[0]['events']} events and "
        f"{rows[1]['work_without_presentation_s']:.3f} s for {rows[1]['events']} events. "
        f"The incremental cost is {rate:.4f} s/event. The projected event-processing cost for 8,000 events is "
        f"{projected:.1f} s ({projected/60:.2f} min), excluding one-time calibration, process loading and presentation. "
        f"The previous linear/SPICE 8,000-event job took {old:.1f} s ({old/60:.2f} min), or {old/8000:.4f} s/event including its fixed overhead.")
    methods=("The benchmark uses bb0nu at 4 mm / 40 ms. Each run uses a recorded fresh noise seed to force a fresh response identity; "
        "absence of filtered response files is asserted before running. Source ion histories and all physical/electronics parameters are unchanged. "
        "Both measured runs use the same warmed Julia process. Covariance calibration is computed and timed separately before event processing; "
        "presentation is timed and subtracted. Processing includes linear screening, nonlinear evaluation, noise generation, reconstruction, "
        "summaries, output writing, publication and consistency checks. The two-point slope removes remaining fixed processing overhead. "
        "No SPICE controls or fallbacks were needed. The 8,000-event value is an extrapolation, subject to event-complexity and disk-throughput variation.")
    tex=r'''\documentclass[11pt]{article}
\usepackage[a4paper,margin=23mm]{geometry}
\usepackage{booktabs,hyperref}
\title{ITACA CSA: per-event production scaling}
\author{Generated from saved timing measurements}\date{2026-11-06}
\begin{document}\maketitle
'''+narrative+r'''
\begin{center}\begin{tabular}{rrrrr}\toprule
Events & Processing (s) & Calibration (s) & Nonlinear pads & SPICE\\\midrule
'''+table+r'''
\bottomrule\end{tabular}\end{center}
'''+methods+f"\n\nThe fitted residual fixed cost is {data['fixed_warm_cost_s']:.3f} s. "
    tex+=r"The earlier 28.46 s pilot reused a filtered response and is excluded from this scaling analysis."+"\n\n"
    tex+=r"Reproduce with \texttt{scripts/production/benchmark\_ion\_csa\_scaling.jl}; the selected timing artifact is recorded in \path{results/ion_csa_nonlinear/scaling/latest.toml}."+r'\end{document}'
    (out/'report.tex').write_text(tex)
    htmlrows=''.join(f"<tr><td>{r['events']}</td><td>{r['work_without_presentation_s']:.3f}</td><td>{r['calibration_setup_s']:.3f}</td></tr>" for r in rows)
    (out/'index.html').write_text(f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>CSA per-event scaling</title>
<style>body{{font:17px system-ui;max-width:950px;margin:40px auto;line-height:1.5}}table{{border-collapse:collapse}}th,td{{padding:10px;border:1px solid #ccd}}</style>
<h1>CSA per-event production scaling</h1><p>{narrative}</p>
<table><tr><th>Events</th><th>Event processing (s)</th><th>Separate calibration (s)</th></tr>{htmlrows}</table>
<p>{methods}</p><p>The earlier 28.46 s run reused a filtered response; it is excluded from this analysis.</p>
<p><a href="report.pdf">PDF report</a></p></html>''')
    build_report(out/'report.tex');print(out/'index.html')

if __name__=='__main__':main()
