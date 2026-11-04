#!/usr/bin/env python3
"""Present saved statistical-versus-microscopic response validation."""
from pathlib import Path
import html,shutil,subprocess,tomllib,os
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]

def plot_label(name):
    labels={'bright':'20k','dark':'Dark','few':'16 triggers','moderate':'1k',
            'transient':'Transient','very_bright':'100k','pde_moderate':'PDE / 1k','pde_bright':'PDE / 100k'}
    return labels.get(name,name.replace('ep2p5','2.5').replace('ep3','3.0').replace('_pure_xe',' / Xe').replace('_xe_he15',' / Xe-He').replace('_independent','\nnew seeds'))

def main():
    """Render numerical comparisons computed in Julia to HTML, PNG and LaTeX/PDF."""
    out=ROOT/'results/sipm_statistical_validation';latex=ROOT/'latex_generated/sipm_statistical_validation';latex.mkdir(parents=True,exist_ok=True)
    rows=[];figures=[];independent='';diagnostics=''
    pilots=[]
    for path in (ROOT/'results/sipm_statistical_pilot').glob('*/manifest.toml'):
        manifest=tomllib.loads(path.read_text())
        if manifest.get('complete') and str(manifest.get('run_card','')).endswith('sipm_bb0nu_statistical_pilot.toml'):
            pilots.append(path)
    pilot=max(pilots,key=lambda p:p.stat().st_mtime) if pilots else None
    pilot_status='The two-event pilot completed for both gases, including the ASIC stage. ' if pilot else 'The two-event pilot card is ready. '
    pilot_link=f' · <a href="{html.escape(os.path.relpath(pilot,out))}">Two-event pilot manifest</a>' if pilot else ''

    with h5py.File(out/'validation.h5') as f:
        if not f['metadata/complete'][()]:raise ValueError('Incomplete benchmark')
        summary=f['default_summary']
        headline=(f"Default: {summary['substeps'][()]} intervals per 250 ns bucket. "
            f"Maximum saved-event mean-charge difference: {summary['max_event_mean_bias_percent'][()]:.4f}%. "
            f"Measured sensor speedup: {summary['min_event_speedup'][()]:.1f} to {summary['max_event_speedup'][()]:.1f} times. "
            +pilot_status+"Full production remains stopped.")
        diagnostics=(f"At 128 intervals, {summary['state_repairs'][()]:,} second-state-moment repairs occurred in "
            f"{summary['statistical_batches'][()]:,} aggregate draws; the largest correction to S was "
            f"{summary['maximum_state_repair'][()]:.4g}. Variance-limit adjustments: {summary['variance_limits'][()]}. "
            "State repairs adjust the dimensionless second recovery sum S. ")
        for category in ('synthetic','events'):
            names=list(f[category]);fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
            for step in (32,64,128):
                keys=[f'{category}_{name}_steps_{step}' for name in names]
                valid=[(i,k) for i,k in enumerate(keys) if k in f['comparison']]
                if not valid:continue
                x=[i for i,k in valid];groups=[f['comparison/'+k] for i,k in valid]
                for ax,field,error,title in [
                    (axes[0,0],'total_mean_bias_percent','total_mean_bias_se_percent','Mean total-charge bias [%]'),
                    (axes[0,1],'total_variance_ratio','total_variance_ratio_se','Total-charge variance / microscopic'),
                    (axes[1,0],'speedup',None,'Runtime speedup'),
                    (axes[1,1],'state_repairs',None,'State feasibility repairs')]:
                    y=[float(g[field][()]) for g in groups]
                    errors=[float(g[error][()]) for g in groups] if error else None
                    ax.errorbar(x,y,yerr=errors,marker='o',capsize=3,label=f'{step} intervals')
                    ax.set_title(title);ax.set_xticks(range(len(names)),[plot_label(n) for n in names],rotation=25,ha='right');ax.grid(alpha=.2);ax.legend(fontsize=8)
            axes[0,0].axhline(0,color='gray',lw=.8);axes[0,1].axhline(1,color='gray',lw=.8)
            if category=='events':
                for i,name in enumerate(names):
                    g=f['comparison/'+category+'_'+name+'_steps_128']
                    if 'controlled_variance_ratio' in g:
                        axes[0,1].errorbar([i],[g['controlled_variance_ratio'][()]],yerr=[g['controlled_variance_ratio_se'][()]],
                            fmt='s',color='forestgreen',capsize=4,label='128 + count control')
                axes[0,1].legend(fontsize=8)
            axes[1,0].set_yscale('log');axes[1,0].axhline(1,color='gray',lw=.8)
            fig.suptitle('Synthetic channel ensembles' if category=='synthetic' else 'Fixed-photon event replicas · both gases and E/P settings')
            filename=category+'_comparison.png';fig.savefig(out/filename,dpi=150);plt.close(fig);shutil.copy2(out/filename,latex/filename);figures.append(filename)
        for key in f['comparison']:
            g=f['comparison/'+key]
            rows.append((key,float(g['total_mean_bias_percent'][()]),float(g['total_variance_ratio'][()]),float(g['speedup'][()]),
                bool(g['mean_pass'][()]),bool(g['variance_pass'][()]),
                None if key.startswith('events_') else bool(g['covariance_pass'][()]),
                int(g['state_repairs'][()]),int(g['variance_limits'][()])))
            if 'reference_detected_variance_ratio' in g:
                independent=(f"Independent-seed event check: charge variance ratio = {g['total_variance_ratio'][()]:.3f} "
                    f"+/- {g['total_variance_ratio_se'][()]:.3f} (one standard error). "
                    f"Detected-count variance / analytical binomial variance: microscopic "
                    f"{g['reference_detected_variance_ratio'][()]:.3f} +/- {g['reference_detected_variance_ratio_se'][()]:.3f}; "
                    f"statistical {g['statistical_detected_variance_ratio'][()]:.3f} +/- {g['statistical_detected_variance_ratio_se'][()]:.3f}. "
                    f"Using the known binomial second moment as a control variate gives charge-variance ratio "
                    f"{g['controlled_variance_ratio'][()]:.3f} +/- {g['controlled_variance_ratio_se'][()]:.3f}. "
                    "The variance estimator subtracts 0.5 times the difference between the sampled squared detected-count deviation "
                    "from its known mean and its analytical variance. That subtraction has expectation zero; the fixed coefficient "
                    "is the same for both models. Standard errors come from the per-replica estimator contributions.")
    explanation=('The statistical model preserves binomial PDE and Poisson dark-count sampling. Fixed counts are partitioned into short intervals; '
        'finite-cell occupancy gives conditional charge moments. A bounded beta draw supplies charge fluctuations, and the same charge updates recovery. '
        'The second recovery moment uses a closure preserving its conditional expected loss. Sparse buckets use individual arrival times and aggregate-state draws. '
        'The physical setup has 7,854 cells, 73.4 ns recovery and 250 ns buckets. Synthetic optical cases use PDE = 1 and zero dark rate to isolate response fluctuations; '
        'the two additional PDE cases use PDE = 0.30 and 10,000 replicas. '
        'Real-event replicas use PDE = 0.30 and the production dark rate, with the same fixed photon input in each replica. '
        'The first bb0nu event of each of the four photon files (both gases at E/P = 2.5 and 3) is repeated 256 times. '
        'The E/P = 2.5 pure-Xe case is also checked with 512 independent new seeds. '
        'Other synthetic cases use 2,000 replicas, except sparse pulses (20,000). '
        'The card selects 128 intervals; 32 and 64 measure discretization sensitivity. '
        'Error bars show one Monte Carlo standard error. Pilot review budgets are 0.05% mean discrepancy and 5% variance discrepancy, '
        'or four Monte Carlo standard errors when larger. Covariance uses four standard errors plus 5% of the neighbouring standard-deviation product. '
        'Synthetic table checks evaluate each bucket and its adjacent covariance. Real-event checks evaluate total-event charge; covariance is marked N/A. '
        'Plots show total-channel or total-event charge, conditional on fixed photon input. This spread is the sensor contribution, not the full event energy resolution. '
        'Sensor-computation speedups use median time per real-event replica; synthetic speedups use the complete ensemble runtime. '
        'Feasibility repairs and bounded-variance limiting are counted explicitly.')
    headers=['Case','Mean bias %','Variance ratio','Speedup','Means','Variances','Covariances','State repairs','Variance limits']
    def val(v):return 'N/A' if v is None else ('pass' if v else 'FAIL') if isinstance(v,bool) else f'{v:.4g}' if isinstance(v,float) else str(v)
    table=''.join('<tr>'+''.join('<td>'+html.escape(val(v))+'</td>' for v in row)+'</tr>' for row in rows)
    page='<!doctype html><meta charset="utf-8"><title>SiPM statistical model validation</title><style>body{font:16px system-ui;max-width:1400px;margin:35px auto;padding:0 20px;color:#17334a}p{line-height:1.6}td,th{padding:7px;border-bottom:1px solid #ddd;text-align:right}img{width:100%}details{padding:12px;background:#f1f5f8}summary{cursor:pointer}</style><h1>SiPM statistical model · validation pilot</h1><p>'+headline+'</p><p><a href="../../latex_generated/sipm_statistical_validation/validation.pdf">PDF report</a> · <a href="../../latex_generated/sipm_response/sipm_electronics_model.pdf">Model derivation</a>'+pilot_link+'</p><details><summary>Methods and comparison criteria</summary><p>'+explanation+'</p></details><h2>Variance follow-up</h2><p>'+independent+'</p>'+''.join(f'<img src="{x}">' for x in figures)+'<h2>Numerical diagnostics</h2><p>'+diagnostics+'</p><table><tr>'+''.join(f'<th>{h}</th>' for h in headers)+'</tr>'+table+'</table>'
    (out/'validation.html').write_text(page)
    tex=r'''\documentclass[10pt]{article}
\usepackage[a4paper,margin=17mm]{geometry}
\usepackage{graphicx,booktabs,hyperref}
\begin{document}\title{SiPM statistical response: validation pilot}\author{NAUSICAA}\date{}\maketitle
'''+('\n\n'.join([headline,explanation,independent,diagnostics])).replace('%',r'\%')+'\n'
    for filename in figures:tex+='\\clearpage\\noindent\\includegraphics[width=\\textwidth]{'+filename+'}\n'
    tex+='\\clearpage\\scriptsize\\noindent\\resizebox{\\textwidth}{!}{\\begin{tabular}{lrrrrrrrr}\\toprule\nCase & Bias (\\%) & Var. ratio & Speedup & Mean & Variance & Covariance & Repairs & Limits \\\\ \\midrule\n'
    for row in rows:tex+=' & '.join(val(v).replace('_',r'\_') for v in row)+r' \\'+'\n'
    tex+='\\bottomrule\\end{tabular}}\\end{document}\n';(latex/'validation.tex').write_text(tex)
    (latex/'.gitignore').write_text('*.aux\n*.log\n*.out\n*.fls\n*.fdb_latexmk\n')
    subprocess.run(['latexmk','-pdf','-interaction=nonstopmode','-halt-on-error','validation.tex'],cwd=latex,check=True,stdout=subprocess.DEVNULL)
    print(out/'validation.html')
if __name__=='__main__':main()
