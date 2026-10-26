#!/usr/bin/env python3
"""Generate manuscript assets from the saved 1.25 pF, 300 us SPICE study.

No simulator invocation. Recheck signal/noise calibration against the stored
transients and one-sided output spectra. Scientific definitions live in
ion_signal.shaping. Original schematic models and saved results are unchanged.
"""
from pathlib import Path
import hashlib
import json
import os
import re
import numpy as np
from ion_signal.shaping import box_average, averaged_frame_transfer, filtered_rms

ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'work/matplotlib-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

DATA=ROOT/'results/ltspice/capacitance/pad_1.25pF'
OUT=ROOT/'notes/physics/generated'
PREVIEW=ROOT/'notes/physics/work/spice-paper'


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
    assert biases=={'Vdd':1.8,'Vnmos':1.6,'Vnload':.6,'Vp6':1.1,'Vp5':.8,'Vp1':.62,'Vf1':.01}
    device_rows=[]
    for line in deck.read_text().splitlines():
        a=line.split()
        if not a or a[0] not in expected:continue
        settings=dict(v.split('=') for v in a[6:])
        length=scalar(settings['l'])*1e6;width=scalar(settings['w'])*1e6
        mult=scalar(settings['m']);kind='PMOS' if a[5].endswith('P') else 'NMOS'
        name=r'$M_f$' if a[0]=='Mf' else a[0]
        device_rows.append(' & '.join([name,kind,f'{width:g}',f'{length:g}',f'{mult:g}'])+' '+chr(92)*2)
    table=r'''\begin{tabular}{llrrr}
\toprule
Device & Type & $W$ ($\mu$m) & $L$ ($\mu$m) & $m$ \\
\midrule
'''+ '\n'.join(device_rows)+'\n'+r'\bottomrule'+'\n'+r'\end{tabular}'+'\n'
    (OUT/'spice_devices_table.tex').write_text('% Generated from the simulated transistor deck.\n'+table)
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
    label(3.7,7.3,r'$V_{DD}=1.8$ V',ha='center')
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
    for x,y,v in [(1.5,6,'1.6 V'),(1.5,2,'0.6 V'),(6,6,'1.1 V'),(6,4,'0.8 V')]:
        label(x-.85,y,v,ha='right',va='center')
    wire([(3.7,2.2),(3.7,1.8)]);label(3.7,1.65,'0.62 V',ha='center',va='top')
    # Parallel feedback capacitor and Mf. Mf D at output, S at input.
    wire([(-1,4),(-1,8.2),(2.9,8.2)])
    wire([(3.1,8.2),(8,8.2),(8,2.95)])
    wire([(2.9,7.92),(2.9,8.48)]);wire([(3.1,7.92),(3.1,8.48)])
    label(3.35,7.85,r'$C_f$',ha='left')
    mos(3,9.6,'$M_f$',horizontal=True,reverse=True)
    wire([(-1,8.2),(-1,9.55),(2.35,9.55)])
    wire([(3.65,9.55),(8,9.55),(8,8.2)])
    wire([(3,8.8),(3.45,8.8)]);label(3.55,8.8,'10 mV',ha='left',va='center')
    ax.text(8.7,.6,'NMOS bulk: ground\nPMOS bulk: $V_{DD}$',ha='right',va='bottom',fontsize=9,color=ink)
    ax.set(xlim=(-2.8,9),ylim=(.2,10.4));ax.set_aspect('equal');ax.axis('off')
    fig.savefig(OUT/'spice_circuit.pdf',metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(PREVIEW/'spice_circuit.png',dpi=170)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True,exist_ok=True);PREVIEW.mkdir(parents=True,exist_ok=True)
    report=json.loads((DATA/'study.json').read_text());c=report['conditions'];r=report['ions']
    T=report['averaging_s'];dt=c['frame_s'];tau=c['estimator_nominal_tau_s'];alpha=np.exp(-dt/tau)
    wave=np.loadtxt(DATA/'ions_N1_r2.csv.gz',delimiter=',',skiprows=1)
    t,i,v=wave[:,0],wave[:,1],wave[:,-1]
    samples=box_average(t,v,np.array([0.,dt])+T/2,T)
    gain=float(-(samples[1]-alpha*samples[0]))
    a=np.loadtxt(DATA/'noise_400.csv',delimiter=',',skiprows=1)
    names=report['noise_400_columns'];f=a[:,0];asd=a[:,names.index('V(onoise)')]
    H=lambda x:averaged_frame_transfer(x,T,dt,alpha)
    sigma=filtered_rms(f,asd,H,dt+T)
    assert np.isclose(gain,r['signal_V_per_ion'],rtol=1e-10)
    assert np.isclose(sigma,r['output_noise_rms_V'],rtol=1e-10)
    n=r['verified_minimum_integer_N'];checks={x['ions']:x for x in r['finite_count_checks']}
    assert checks[n-1]['snr']<=5<checks[n]['snr']
    with plt.rc_context({'font.size':11,'axes.spines.top':False,'axes.spines.right':False}):
        fig,axs=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
        axs[0].plot(t*1e3,i*1e15,color='#245e8b')
        axs[0].set(xlabel='Time (ms)',ylabel='Current per ion (fA)',xlim=(-2,15))
        axs[1].plot(t*1e3,v*1e6,color='#245e8b',label='CSA output')
        axs[1].plot(np.array([0,dt])*1e3,samples*1e6,'o',color='#b46618',label='Averaged readings')
        axs[1].set(xlabel='Time (ms)',ylabel='Voltage per ion (µV)',xlim=(-2,15))
        axs[1].legend(frameon=False,loc='upper right',fontsize=9)
        for ax in axs:
            for edge in (0,dt):ax.axvspan((edge-T/2)*1e3,(edge+T/2)*1e3,color='#b46618',alpha=.18)
            ax.grid(alpha=.18)
        fig.savefig(OUT/'spice_ion_response.pdf',metadata={'CreationDate':None,'ModDate':None})
        fig.savefig(PREVIEW/'spice_ion_response.png',dpi=170);plt.close(fig)
        # Full and filtered spectra; filter lobes resolved explicitly at 5 Hz.
        dense=np.unique(np.r_[f,np.arange(f[0],f[-1],5.)])
        ad=np.exp(np.interp(np.log(dense),np.log(f),np.log(asd)))
        fig,ax=plt.subplots(figsize=(7.2,4.0),layout='constrained')
        ax.loglog(f,asd*1e6,label='CSA output',color='#245e8b')
        ax.loglog(dense,ad*np.abs(H(dense))*1e6,label='After the frame filter',color='#b46618',lw=.8)
        ax.set(xlabel='Frequency (Hz)',ylabel=r'Noise amplitude ($\mu$V/$\sqrt{\mathrm{Hz}}$)',xlim=(f[0],f[-1]),ylim=(1e-4,2e4))
        ax.legend(frameon=False);ax.grid(alpha=.18,which='major')
        fig.savefig(OUT/'spice_noise.pdf',metadata={'CreationDate':None,'ModDate':None})
        fig.savefig(PREVIEW/'spice_noise.png',dpi=170);plt.close(fig)
    values={
        'SpicePadPf':f"{report['external_pad_F']*1e12:g}", 'SpiceCfFf':f"{c['feedback_F']*1e15:g}",
        'SpiceFrameMs':f'{dt*1e3:g}','SpiceAverageUs':f'{T*1e6:g}','SpiceCorrectionSeconds':f'{tau:g}',
        'SpiceTemperatureK':f"{c['temperature_K']:g}",'SpiceLowHz':f'{f[0]:g}','SpiceHighMHz':f'{f[-1]*1e-6:g}',
        'SpiceGainUv':f'{gain*1e6:.2f}','SpiceNoiseMv':f'{sigma*1e3:.3f}','SpiceEncIons':f'{sigma/gain:.2f}',
        'SpiceMinimumIons':str(n),'SpiceBelowIons':str(n-1),'SpiceBelowSnr':f"{checks[n-1]['snr']:.3f}",
        'SpiceAboveSnr':f"{checks[n]['snr']:.3f}",
        'SpicePreviousPercent':f"{100*r['neighbour_packet_response_relative']['-0.01']:.2f}",
        'SpiceFollowingPercent':f"{100*r['neighbour_packet_response_relative']['0.01']:.2f}"}
    deck=DATA/f'decks/ions_N{n}_r2.cir'
    supply=re.search(r'^Vdd Vdda 0 (\S+)$',deck.read_text(),re.M)[1]
    values['SpiceSupplyV']=supply
    # Keep the chosen rounding, without padding displayed numbers with zeros.
    values = {key: value.rstrip("0").rstrip(".")
               if "." in value and "e" not in value.lower() else value
               for key, value in values.items()}
    (OUT/'spice_values.tex').write_text('% Generated by scripts/build_spice_paper.py from saved transistor results.\n'+
        ''.join('\\newcommand{\\'+k+'}{'+v+'}\n' for k,v in values.items()))
    circuit_figure(deck)
    sources=[DATA/'study.json',DATA/'ions_N1_r2.csv.gz',DATA/'noise_400.csv',deck,
        ROOT/'src/ion_signal/shaping.py',Path(__file__).resolve()]
    output={'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        'conditions':c,'external_pad_F':report['external_pad_F'],'averaging_s':T,
        'signal_V_per_ion':gain,'noise_rms_V':sigma,'enc_ions_rms':sigma/gain,
        'verified_minimum_ions':n,'macros':values}
    (OUT/'spice_paper_provenance.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:output[k] for k in ('signal_V_per_ion','noise_rms_V','enc_ions_rms','verified_minimum_ions')},indent=2))

if __name__=='__main__':main()
