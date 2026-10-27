#!/usr/bin/env python3
"""Present the final circuit and calibrated 20 ms Julia/HDF5 results.

Generates manuscript macros, tables, figures and source provenance without
rerunning LTspice or repeating the physical reduction.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import tomllib
import numpy as np
import h5py

ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'work/matplotlib-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

OUT=ROOT/'notes/physics/generated'
PREVIEW=ROOT/'notes/physics/work/spice-paper'


def selected_circuit():
    config = tomllib.loads((ROOT / "config/bucket_study.toml").read_text())
    return ROOT / config["circuit"]


def circuit_figure(deck):
    """Draw the seven devices with drain/source connectivity from the saved deck."""
    nodes={};params={}
    for line in deck.read_text().splitlines():
        a=line.split()
        if not a:continue
        if a[0].lower().startswith('m'):
            nodes[a[0]]=a[1:5]
        if a[0].lower() in ('.param','.parameter'):
            k,v=a[1].split('=');params[k]=v
    expected={'M1':['s','in','med','Vdda'],'M2':['Vdda','Vnmos','s','0'],
        'M3':['med','Vnload','0','0'],'M4':['med','Vncasc','out','0'],
        'M5':['vd6','Vp5','out','Vdda'],'M6':['Vdda','Vp6','vd6','Vdda'],
        'Mf':['out','Vdis','in','0']}
    assert nodes==expected
    def scalar(value):
        value=value.strip('{}');value=params.get(value,value)
        suffix={'u':1e-6,'m':1e-3}
        return float(value[:-1])*suffix[value[-1]] if value[-1] in suffix else float(value)
    biases={a[0]:scalar(a[3]) for line in deck.read_text().splitlines()
        if (a:=line.split()) and a[0] in ('Vdd','Vnmos','Vnload','Vp6','Vp5','Vp1','Vf1')}
    assert set(biases)=={'Vdd','Vnmos','Vnload','Vp6','Vp5','Vp1','Vf1'}
    assert all(np.isfinite(v) for v in biases.values())
    devices=[]
    device_rows=[]
    for line in deck.read_text().splitlines():
        a=line.split()
        if not a or a[0] not in expected:continue
        settings=dict(v.split('=') for v in a[6:])
        length=scalar(settings['l'])*1e6;width=scalar(settings['w'])*1e6
        mult=scalar(settings['m']);kind='PMOS' if a[5].endswith('P') else 'NMOS'
        devices.append(dict(name=a[0],type=kind,width_um=width,length_um=length,multiplicity=mult))
        name=r'$M_f$' if a[0]=='Mf' else a[0]
        device_rows.append(' & '.join([name,kind,f'{width:g}',f'{length:g}',f'{mult:g}'])+' '+chr(92)*2)
    table=r'''\begin{tabular}{llrrr}
\toprule
Device & Type & $W$ ($\mu$m) & $L$ ($\mu$m) & $m$ \\
\midrule
'''+ '\n'.join(device_rows)+'\n'+r'\bottomrule'+'\n'+r'\end{tabular}'+'\n'
    (OUT/'spice_devices_table.tex').write_text(f'% Generated from {deck.relative_to(ROOT)}.\n'+table)
    fig,ax=plt.subplots(figsize=(8.4,6.8),layout='constrained')
    ink='#263747';accent='#245e8b'
    def wire(points,**kwargs):ax.plot(*zip(*points),color=ink,lw=1.1,**kwargs)
    def label(x,y,s,**kwargs):ax.text(x,y,s,fontsize=10,color=ink,**kwargs)
    def dot(x,y):ax.plot(x,y,'o',color=ink,ms=3)
    def ground(x,y):
        for j,w in enumerate((.35,.23,.10)):wire([(x-w,y-.12*j),(x+w,y-.12*j)])
    def mos(x,y,name,p=False,horizontal=False,reverse=False):
        # Local drain=(0,+.65), source=(0,-.65), gate=(-.8,0).
        def tr(a,b):return (x+(b if reverse else -b),y+a) if horizontal else (x+a,y+b)
        def line(points):wire([tr(a,b) for a,b in points])
        line([(-.05,.65),(-.05,.4),(.25,.4)])
        line([(.25,-.4),(-.05,-.4),(-.05,-.65)])
        line([(.25,-.4),(.25,.4)])
        line([(-.15,-.42),(-.15,.42)])
        line([(-.8,0),(-.27 if p else -.15,0)])
        if p:ax.add_patch(Circle(tr(-.21,0),.06,fill=False,ec=ink,lw=1))
        if horizontal:
            label(x+.1,y+.5,name,ha='center')
            ax.text(x+(.65 if reverse else -.65),y+.12,'D',fontsize=7,ha='center')
            ax.text(x+(-.65 if reverse else .65),y+.12,'S',fontsize=7,ha='center')
        else:
            label(x+.45,y,name,va='center')
            ax.text(x+.06,y+.62,'D',fontsize=7)
            ax.text(x+.06,y-.66,'S',fontsize=7)
    # Main branches: M2--M1--M3 on the left, M6--M5 on the right.
    for x,y,n,p in [(1.5,6,'M2',False),(1.5,4,'M1',True),(1.5,2,'M3',False),(6,6,'M6',True),(6,4,'M5',True)]:
        mos(x,y,n,p)
    wire([(1.45,6.65),(1.45,7.15),(5.95,7.15),(5.95,6.65)])
    label(3.7,7.3,rf'$V_{{DD}}={biases["Vdd"]:g}$ V',ha='center')
    for x in (1.45,5.95):wire([(x,5.35),(x,4.65)])
    wire([(1.45,3.35),(1.45,2.65)])
    wire([(1.45,1.35),(1.45,.7)]);ground(1.45,.7)
    mos(3.7,3,'M4',horizontal=True)
    # Rotated MOS drain/source at y=2.95.
    wire([(1.45,2.95),(3.05,2.95)])
    wire([(4.35,2.95),(5.95,2.95),(5.95,3.35)])
    wire([(5.95,2.95),(8,2.95)])
    for x in (1.45,5.95):dot(x,2.95)
    label(8.05,2.95,r'$V_{\rm out}$',va='center')
    # Input to M1 gate; Cd1 and signal current enter this node.
    wire([(-1,4),(.7,4)])
    dot(-1,4)
    ax.annotate(r'$i_N(t)$',xy=(-1,4),xytext=(-2.4,4),ha='center',va='center',fontsize=11,color=accent,
        arrowprops={'arrowstyle':'->','color':accent})
    wire([(-1,4),(-1,2.1)]);wire([(-1.35,2.1),(-.65,2.1)])
    wire([(-1.35,1.94),(-.65,1.94)]);wire([(-1,1.94),(-1,.7)]);ground(-1,.7)
    label(-1.5,2.02,r'$C_{\rm pad}$',ha='right',va='center')
    # Bias gates.
    for x,y,key in [(1.5,6,'Vnmos'),(1.5,2,'Vnload'),(6,6,'Vp6'),(6,4,'Vp5')]:
        label(x-.85,y,f'{biases[key]:.4f}'.rstrip('0').rstrip('.')+' V',ha='right',va='center')
    wire([(3.7,2.2),(3.7,1.8)])
    label(3.7,1.65,f'{biases["Vp1"]:.4f}'.rstrip('0').rstrip('.')+' V',ha='center',va='top')
    # Parallel feedback capacitor and Mf. Mf D at output, S at input.
    wire([(-1,4),(-1,8.2),(2.9,8.2)])
    wire([(3.1,8.2),(8,8.2),(8,2.95)])
    wire([(2.9,7.92),(2.9,8.48)]);wire([(3.1,7.92),(3.1,8.48)])
    label(3.35,7.85,r'$C_f$',ha='left')
    mos(3,9.6,'$M_f$',horizontal=True,reverse=True)
    wire([(-1,8.2),(-1,9.55),(2.35,9.55)])
    wire([(3.65,9.55),(8,9.55),(8,8.2)])
    wire([(3,8.8),(3.45,8.8)]);label(3.55,8.8,f'{biases["Vf1"]*1e3:g} mV',ha='left',va='center')
    ax.text(8.7,.6,'NMOS bulk: ground\nPMOS bulk: $V_{DD}$',ha='right',va='bottom',fontsize=9,color=ink)
    ax.set(xlim=(-2.8,9),ylim=(.2,10.4));ax.set_aspect('equal');ax.axis('off')
    fig.savefig(OUT/'spice_circuit.pdf',metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(PREVIEW/'spice_circuit.png',dpi=170)
    plt.close(fig)
    config=tomllib.loads((ROOT/'config/bucket_study.toml').read_text())
    analysis=(ROOT/config['noise_source']).with_name('analysis.h5')
    with h5py.File(analysis) as h:
        operating=h[config['noise_group']].attrs
        output_V=float(operating['dc_output_V'])
        power_W=float(operating['power_W'])
    (OUT/'spice_operating_values.tex').write_text(
        f'% Generated from {analysis.relative_to(ROOT)}, {config["noise_group"]}.\n'
        +r'\newcommand{\SpiceBaselineV}{'+f'{output_V:.3f}'+"}\n"
        +r'\newcommand{\SpicePowerUw}{'+f'{power_W*1e6:.0f}'+"}\n")
    sources=[deck,deck.parent/'corrected.inc',ROOT/'config/bucket_study.toml',analysis,Path(__file__).resolve()]
    provenance={
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        'biases_V':biases,'devices':devices,
        'operating_point':{'dc_output_V':output_V,'quiescent_power_W':power_W},
        'generated_files':['spice_circuit.pdf','spice_devices_table.tex','spice_operating_values.tex'],
    }
    (OUT/'spice_circuit_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')



def discharge_value():
    """Update the correction time from the fitted circuit tail, preserving other macros."""
    config_path=ROOT/'config/discharge_calibration.toml'
    config=tomllib.loads(config_path.read_text())
    source=ROOT/config['output']/'analysis.h5'
    with h5py.File(source) as h:
        tau=float(h.attrs['tau_s'])
    value=f'{tau:.3f}'.rstrip('0').rstrip('.')
    target=OUT/'spice_values.tex'
    lines=target.read_text().splitlines()
    prefix='\\newcommand{\\SpiceCorrectionSeconds}{'
    assert sum(line.startswith(prefix) for line in lines)==1
    target.write_text('\n'.join(prefix+value+'}' if line.startswith(prefix) else line for line in lines)+'\n')
    sources=[source,config_path,Path(__file__).resolve()]
    provenance={'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                'macro':'SpiceCorrectionSeconds','tau_s':tau,'displayed_value':value,
                'method':'Exponential fit to baseline-subtracted, current-free CSA tail'}
    (OUT/'spice_discharge_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    return value

def paper_data_path():
    config=tomllib.loads((ROOT/'config/discharge_calibration.toml').read_text())
    return ROOT/config['output']/'paper_readout.h5'


def display(value, places):
    return f'{value:.{places}f}'.rstrip('0').rstrip('.')


def write_macros(path, values, source):
    path.write_text(f'% Generated from {source}; do not edit.\n'+
                   ''.join('\\newcommand{\\'+key+'}{'+value+'}\n' for key,value in values.items()))


def plot_averaging():
    """Present the saved Julia scan, with the selected point checked against calibration."""
    source=paper_data_path()
    with h5py.File(source) as h:
        a=dict(h.attrs);scan=h['averaging']
        widths=scan['width_s'][:]*1e6;gains=scan['gain_V_per_ion'][:]
        noises=scan['noise_rms_V'][:];enc=scan['ENC_e_rms'][:]
        selected=int(scan.attrs['selected_index'])-1
        errors={key:float(max(scan[key][:])) for key in ('quadrature_error','spectral_error','transient_error')}
    assert np.isclose(gains[selected],a['gain_V_per_ion'],rtol=1e-10)
    assert np.isclose(noises[selected],a['noise_rms_V'],rtol=1e-10)
    OUT.mkdir(parents=True,exist_ok=True);PREVIEW.mkdir(parents=True,exist_ok=True)
    with plt.rc_context({'font.size':11,'axes.spines.top':False,'axes.spines.right':False}):
        fig,ax=plt.subplots(figsize=(6.7,3.8),layout='constrained')
        ax.plot(widths,enc,'o-',color='#245e8b',lw=1.5)
        ax.plot(widths[selected],enc[selected],'o',color='#b46618',ms=8,zorder=3)
        ax.annotate(f'{widths[selected]:g} µs: {enc[selected]:.1f} electrons RMS',
                    xy=(widths[selected],enc[selected]),xytext=(.40,.85),textcoords='axes fraction',
                    arrowprops={'arrowstyle':'->','color':'#b46618'},color='#925317')
        ax.set(xlabel='Averaging time (µs)',ylabel='ENC (electrons RMS)',xticks=widths)
        ax.grid(alpha=.18)
        fig.savefig(OUT/'spice_averaging.pdf',metadata={'CreationDate':None,'ModDate':None})
        fig.savefig(PREVIEW/'spice_averaging.png',dpi=170);plt.close(fig)
    write_macros(OUT/'spice_averaging_values.tex',
                 {'SpiceScanFirstEnc':display(enc[0],1),'SpiceScanLastEnc':display(enc[-1],1)},source.relative_to(ROOT))
    reportdir=ROOT/'latex_generated/calibrated_spice_averaging';reportdir.mkdir(parents=True,exist_ok=True)
    import shutil
    shutil.copyfile(OUT/'spice_averaging.pdf',reportdir/'spice_averaging.pdf')
    rows='\n'.join(f'{w:g} & {display(g*1e6,3)} & {display(n*1e3,3)} & {display(e,2)} '+r'\\' for w,g,n,e in zip(widths,gains,noises,enc))
    report=r'''\documentclass[11pt]{article}
\usepackage[margin=24mm]{geometry}
\usepackage{graphicx,booktabs,amsmath,siunitx}
\begin{document}
\section*{Averaging scan for the calibrated NAUSICAA readout}
The selected 180 nm circuit has M3 dimensions $200\times2\,\mu\mathrm{m}$,
M4 dimensions $2.88\times1.44\,\mu\mathrm{m}$, external input capacitance
1.25 pF, feedback capacitance 3 fF, supply 1.8 V and temperature 293.15 K.
The transistor models include intrinsic capacitances and thermal and flicker noise.
Uniform ion collection lasts @FRAME@ ms; the calibrated discharge time is
@TAU@ s. Voltage averages are centred at the frame boundaries.

The input is the finite-gap current of a central trajectory to a 0.5 mm
square pad in an infinite grounded plane, with a homogeneous gas half-space.
The weighting potential is one on the pad, zero on the surrounding plane,
and zero at infinity. The ion moves at 3 m/s over the 5 mm gap. The source
preserves $e[1-\phi_w(g)]$ per ion.

Julia recalculates the signal gain for every averaging width and integrates
the stationary output noise from 0.1 Hz to 1 MHz with the same averaging and
frame-difference filter. The circuit and raw transients are unchanged.
\begin{center}
\includegraphics[width=.9\linewidth]{spice_averaging.pdf}
\begin{tabular}{rrrr}\toprule
Average ($\mu$s) & Gain ($\mu$V/ion) & Noise (mV RMS) & ENC ($e$ RMS)\\\midrule
@ROWS@
\bottomrule\end{tabular}\end{center}
The selected 300 $\mu$s point reproduces the saved discharge-calibration gain
and noise to a relative tolerance of $10^{-10}$. Checks halve the transient
step, compare 200 and 400 spectral points per decade and increase quadrature
order from 8 to 16. Maximum relative changes are @TERR@, @SERR@ and @QERR@,
respectively.

Reproduction: run \texttt{scripts/build\_averaging\_scan.jl}, then
\texttt{scripts/plot\_averaging\_scan.py}. The Julia result is
\texttt{results/ltspice/discharge\_calibration/20ms/paper\_readout.h5}.
It records hashes of the raw data, calibration, configuration, analysis and filter.
\end{document}
'''
    replacements={'FRAME':f"{a['frame_s']*1e3:g}",'TAU':display(a['tau_s'],6),'ROWS':rows,
                  'TERR':r'\num{'+f"{errors['transient_error']:.3g}"+'}',
                  'SERR':r'\num{'+f"{errors['spectral_error']:.3g}"+'}',
                  'QERR':r'\num{'+f"{errors['quadrature_error']:.3g}"+'}'}
    for key,value in replacements.items():report=report.replace('@'+key+'@',value)
    (reportdir/'report.tex').write_text(report)



def noise_budget_table():
    """Present saved device variance fractions after averaging and frame subtraction."""
    source=paper_data_path()
    with h5py.File(source) as h:
        fractions=dict(zip(h['device_names'].asstr()[:],map(float,h['device_variance_fraction'][:])))
    assert set(fractions)=={'V(m1)','V(m2)','V(m3)','V(m4)','V(m5)','V(m6)','V(mf)'}
    assert np.isclose(sum(fractions.values()),1,rtol=1e-3)
    groups=[('$M_f$',['V(mf)']),('$M_4$',['V(m4)']),
            ('$M_2$',['V(m2)']),('$M_3$',['V(m3)']),
            ('$M_1$, $M_5$, $M_6$',['V(m1)','V(m5)','V(m6)'])]
    bs=chr(92);row_end=bs*2
    rows=[label+' & '+display(100*sum(fractions[key] for key in keys),1)+' '+row_end for label,keys in groups]
    lines=[f'% Generated from {source.relative_to(ROOT)}:device_variance_fraction.',
           bs+'begin{tabular}{lr}',bs+'toprule',
           'Transistor & Contribution ('+bs+'%) '+row_end,bs+'midrule',
           *rows,bs+'bottomrule',bs+'end{tabular}']
    (OUT/'spice_noise_budget.tex').write_text(chr(10).join(lines)+chr(10))
    return fractions


def main():
    """Generate manuscript assets from the calibrated circuit's Julia/HDF5 results."""
    OUT.mkdir(parents=True,exist_ok=True);PREVIEW.mkdir(parents=True,exist_ok=True)
    source=paper_data_path()
    with h5py.File(source) as h:
        a=dict(h.attrs)
        t=h['time_s'][:];i=h['current_A_per_ion'][:];v=h['voltage_V_per_ion'][:]
        sample_times=h['sample_time_s'][:];samples=h['sample_voltage_V_per_ion'][:]
        counts=h['counts'][:];snrs=h['snr'][:]
    cfgpath=ROOT/'config/discharge_calibration.toml';cfg=tomllib.loads(cfgpath.read_text())
    assert a['config_sha256'].decode()==hashlib.sha256(cfgpath.read_bytes()).hexdigest()
    calibration=ROOT/cfg['output']/'analysis.h5'
    assert a['calibration_sha256'].decode()==hashlib.sha256(calibration.read_bytes()).hexdigest()
    dt=float(a['frame_s']);T=float(a['averaging_s']);n=int(a['verified_N_snr5'])
    checks=dict(zip(counts,snrs));assert checks[n-1]<5<checks[n]
    deck=selected_circuit();lines=[line.split() for line in deck.read_text().splitlines() if line.strip()]
    def scalar(value):
        suffix={'f':1e-15,'p':1e-12,'n':1e-9,'u':1e-6,'m':1e-3}
        return float(value[:-1])*suffix[value[-1]] if value[-1] in suffix else float(value)
    def component(name):return next(x for x in lines if x[0]==name)
    pad=scalar(component('Cd1')[3]);cf=scalar(component('Cf')[3]);supply=scalar(component('Vdd')[3])
    temperature=scalar(component('.temp')[1])+273.15
    with plt.rc_context({'font.size':11,'axes.spines.top':False,'axes.spines.right':False}):
        fig,axs=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
        axs[0].plot(t*1e3,i*1e15,color='#245e8b')
        axs[0].set(xlabel='Time (ms)',ylabel='Current per ion (fA)')
        axs[1].plot(t*1e3,v*1e6,color='#245e8b',label='CSA output')
        axs[1].plot(sample_times*1e3,samples*1e6,'o',color='#b46618',label='Averaged readings')
        axs[1].set(xlabel='Time (ms)',ylabel='Voltage per ion (µV)')
        axs[1].legend(frameon=False,loc='upper right',fontsize=9)
        for ax in axs:
            ax.set_xlim(-2,dt*1e3*1.5)
            for edge in sample_times:ax.axvspan((edge-T/2)*1e3,(edge+T/2)*1e3,color='#b46618',alpha=.18)
            ax.grid(alpha=.18)
        fig.savefig(OUT/'spice_ion_response.pdf',metadata={'CreationDate':None,'ModDate':None})
        fig.savefig(PREVIEW/'spice_ion_response.png',dpi=170);plt.close(fig)
    values={
        'SpicePadPf':f'{pad*1e12:g}','SpiceCfFf':f'{cf*1e15:g}',
        'SpiceFrameMs':f'{dt*1e3:g}','SpiceAverageUs':f'{T*1e6:g}',
        'SpiceCorrectionSeconds':display(a['tau_s'],3),'SpiceTemperatureK':f'{temperature:g}',
        'SpiceLowHz':f"{a['frequency_low_Hz']:g}",'SpiceHighMHz':f"{a['frequency_high_Hz']/1e6:g}",
        'SpiceGainUv':display(a['gain_V_per_ion']*1e6,3),'SpiceNoiseMv':display(a['noise_rms_V']*1e3,3),
        'SpiceEncIons':display(a['ENC_e_rms'],2),'SpiceMinimumIons':str(n),'SpiceBelowIons':str(n-1),
        'SpiceBelowSnr':display(checks[n-1],4),'SpiceAboveSnr':display(checks[n],4),
        'SpicePreviousPercent':display(a['previous_frame_residual_fraction']*100,3),
        'SpiceFollowingPercent':display(a['following_frame_fraction']*100,3),'SpiceSupplyV':f'{supply:g}'}
    write_macros(OUT/'spice_values.tex',values,source.relative_to(ROOT))
    circuit_figure(deck);discharge_value();plot_averaging()
    device_variances=noise_budget_table()
    sources=[source,calibration,cfgpath,deck,ROOT/cfg['model'],Path(__file__).resolve()]
    provenance={'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                'frame_s':dt,'averaging_s':T,'correction_tau_s':float(a['tau_s']),
                'gain_V_per_ion':float(a['gain_V_per_ion']),'noise_rms_V':float(a['noise_rms_V']),
                'ENC_e_rms':float(a['ENC_e_rms']),'verified_minimum_ions':n,'macros':values,
                'flicker_variance_fraction':float(a['flicker_variance_fraction']),
                'device_variance_fraction':device_variances}
    (OUT/'spice_paper_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(json.dumps(values,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--circuit-only',action='store_true',help='Generate selected circuit schematic and device table')
    parser.add_argument('--discharge-only',action='store_true',help='Update only the calibrated discharge-time macro')
    args=parser.parse_args()
    if args.discharge_only:discharge_value()
    elif args.circuit_only:
        OUT.mkdir(parents=True,exist_ok=True);PREVIEW.mkdir(parents=True,exist_ok=True)
        circuit_figure(selected_circuit())
    else:main()
