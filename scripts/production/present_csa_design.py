#!/usr/bin/env python3
"""Present the Julia CSA-note summary as LaTeX tables, macros and figures."""
from pathlib import Path
import argparse, hashlib, json, os, re
import h5py, numpy as np
ROOT=Path(__file__).resolve().parents[2]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'notes/physics/work/csa_design_matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LABELS={'bb0nu':'Double beta','single_electrons':'Single electron'}
COLORS={2:'#5f7e43',3:'#176589',4:'#bd6333'}
def attr(g,k):
 v=g.attrs[k]
 return v.decode() if isinstance(v,bytes) else v

def table(path,header,rows,columns):
 path.write_text('\\begin{tabular}{'+columns+'}\n\\toprule\n'+header+' \\\\\n\\midrule\n'+'\n'.join(' & '.join(row)+r' \\' for row in rows)+'\n\\bottomrule\n\\end{tabular}\n')

def save(fig,out,name):
 fig.tight_layout()
 for ext in ('pdf','png'):fig.savefig(out/(name+'.'+ext),dpi=170,bbox_inches='tight')
 plt.close(fig)

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,default=ROOT/'results/csa_design');a=p.parse_args()
 source=a.results/'summary.h5';out=ROOT/'latex_generated/ion_readout/csa_design';out.mkdir(parents=True,exist_ok=True)
 with h5py.File(source) as h:
  assert h.attrs['complete']
  parameters=json.loads(h['parameters_json'][()].decode());response=__import__('tomllib').loads(h['response_card_toml'][()].decode());design=__import__('tomllib').loads(h['design_card_toml'][()].decode())
  pilots=h['pilots'];pitches=sorted(map(int,pilots.keys()));cal={f'{pitch*10}ms':pilots[str(pitch)]['calibration'][f'{pitch*10}ms'] for pitch in pitches};selected=parameters['circuit']
  macros={'CsaNominalTau':f"{design['nominal_discharge_time_s']:g}",'CsaFittedTau':f"{attr(cal['20ms'],'tau_s'):.2f}",'CsaAverageMs':f"{1e3*response['averaging_s']:g}",'CsaCfFf':f"{1e15*parameters['feedback_capacitance_F']:g}",'CsaCinPf':f"{1e12*parameters['external_input_capacitance_F']:g}",'CsaPadMm':f"{1e3*response['pad_side_m']:g}",'CsaGapMm':f"{1e3*response['gap_m']:g}",'CsaGapSpeed':f"{response['gap_speed_m_s']:g}",'CsaBaseline':f"{selected['baseline_V']:.3f}",'CsaSupply':f"{parameters['voltage_sources']['Vdd']['voltage_V']:g}",'CsaCurrentUa':f"{1e6*h['supply_current_A'][()]:.1f}",'CsaPowerMw':f"{1e3*h['supply_power_W'][()]:.3f}",'CsaTemperature':f"{parameters['temperature_K']:.2f}",'CsaNoiseLow':f"{response['noise_low_Hz']:g}",'CsaNoiseHighMHz':f"{response['noise_high_Hz']/1e6:g}",'CsaRecoveryS':f"{response['recovery_tail_s']:g}",'CsaQuietFrames':f"{response['quiet_before_frames']}",'CsaLinearTolerance':f"{response['linear_error_tolerance_ENC']:g}",'CsaNoteChecks':str(h.attrs['validation_checks'])}
  for pitch in pitches:
   g=pilots[str(pitch)]['bb0nu'];macros['CsaEnc'+{2:'Two',3:'Three',4:'Four'}[pitch]]=f"{attr(g,'ENC_e_rms'):.2f}"
  (out/'values.tex').write_text('% Generated from saved Julia reductions.\n'+'\n'.join('\\newcommand{\\'+k+'}{'+v+'}' for k,v in macros.items())+'\n')
  rows=[]
  for name,t in parameters['transistors'].items():rows.append([('$M_f$' if name=='Mf' else '$M_'+name[1:]+'$'),'PMOS' if t['model'].endswith('P') else 'NMOS',f"{t['width_m']*1e6:g}",f"{t['length_m']*1e6:g}",f"{t['multiplicity']:g}"])
  table(out/'transistors.tex',r'Device & Type & $W$ ($\mu$m) & $L$ ($\mu$m) & Multiplicity',rows,'llrrr')
  rows=[]
  for name,s in parameters['voltage_sources'].items():rows.append([r'\texttt{'+s['positive_node'].replace('_',r'\_')+'}',f"{s['voltage_V']:.4f}".rstrip('0').rstrip('.')])
  table(out/'biases.tex','Node & Voltage (V)',rows,'lr')
  rows=[]
  for pitch in pitches:
   c=cal[f'{pitch*10}ms'];rows.append([str(pitch),str(pitch*10),f"{attr(c,'tau_s'):.2f}",f"{1e6*attr(c,'gain_V_per_ion'):.2f}",f"{attr(c,'precorrection_ENC_e_rms'):.2f}",f"{attr(c,'ENC_e_rms'):.2f}",f"{c['covariance_e2'][1]/c['covariance_e2'][0]:.3f}"])
  table(out/'calibration.tex',r'Pitch (mm) & Frame (ms) & $\tau_f$ (s) & $g$ ($\mu$V/ion) & ENC$_y$ & ENC$_{\widehat N}$ & $r_1$',rows,'rrrrrrr')
  rows=[]
  for pitch in pitches:
   for sample in LABELS:
    g=pilots[str(pitch)][sample];rows.append([str(pitch),LABELS[sample],str(len(g['event_id'])),f"{attr(g,'mean_ions_per_positive_cell'):.1f}",f"{attr(g,'inclusive_mean_snr'):.2f}",f"{attr(g,'event_weighted_mean_snr'):.2f}",str(attr(g,'pads')),str(attr(g,'direct_spice_pads'))])
  table(out/'signals.tex',r'Pitch & Sample & Events & Ions/cell & Cell S/N & Event S/N & Pads & SPICE',rows,'rlrrrrrr')
  rows=[]
  for pitch in pitches:
   for sample in LABELS:
    g=pilots[str(pitch)][sample];rows.append([str(pitch),LABELS[sample],f"{100*g['fraction_events_with_compression'][2]:.1f}",f"{100*attr(g,'mean_compression_deficit_fraction'):.3f}",f"{100*attr(g,'mean_noiseless_charge_bias_fraction'):+.3f}",f"{100*attr(g,'mean_noisy_charge_bias_fraction'):+.3f}",f"{100*attr(g,'noisy_bias_std_fraction'):.2f}"])
  table(out/'charge.tex',r'Pitch & Sample & Local loss events (\%) & Deficit (\%) & Bias (\%) & Noisy bias (\%) & Spread (\%)',rows,'rlrrrrr')
  rows=[]
  for pitch in pitches:
   g=pilots[str(pitch)]['validation'];rows.append([str(pitch),f"{attr(g,'maximum_remaining_control_error_ENC'):.3g}",f"{np.max(g['refinement_error_ENC'][:]):.3g}"])
  table(out/'validation.tex',r'Pitch (mm) & Maximum fast-response error (ENC) & Solver refinement (ENC)',rows,'rrr')
  # All plotted quantities have already been calculated and saved by Julia.
  fig,axes=plt.subplots(1,2,figsize=(10.5,3.6))
  retention=h['retention'];names=sorted(retention,key=lambda k:attr(retention[k],'length_m'))
  for frame in (20,40):
   tau=[attr(retention[n][f'{frame}ms'],'tau_s') for n in names];enc=[attr(retention[n][f'{frame}ms'],'ENC_e_rms') for n in names]
   axes[0].plot(tau,enc,'o-',label=f'{frame} ms frames')
  axes[0].axvline(attr(cal['20ms'],'tau_s'),ls=':',color='0.5');axes[0].set(xlabel='Fitted discharge time (s)',ylabel='Scalar frame ENC (e⁻ RMS)');axes[0].legend(frameon=False)
  r=retention[design['selected_variant']];axes[1].loglog(r['frequency_Hz'][:],r['output_asd_V_sqrtHz'][:],color=COLORS[3]);axes[1].set(xlabel='Frequency (Hz)',ylabel='Output noise ASD (V / √Hz)')
  for ax in axes:ax.grid(alpha=.2)
  save(fig,out,'retention_noise')
  fig,axes=plt.subplots(1,2,figsize=(10.5,3.6),sharex=True,sharey=True)
  for ax,sample in zip(axes,LABELS):
   for pitch in pitches:
    g=pilots[str(pitch)][sample]['snr'];counts=g['count'][:];normalizer=sum(counts)+attr(g,'underflow')+attr(g,'overflow');ax.stairs(counts/normalizer,g['edges'][:],color=COLORS[pitch],label=f'{pitch} mm / {10*pitch} ms')
   ax.set(xscale='log',xlabel='Expected S/N per positive pixel–frame',title=LABELS[sample]);ax.grid(alpha=.2);ax.legend(frameon=False,fontsize=8)
  axes[0].set_ylabel('Fraction of cells / bin');save(fig,out,'pilot_snr')
  fig,axes=plt.subplots(1,2,figsize=(10.5,3.6))
  for pitch in pitches:
   c=cal[f'{pitch*10}ms'];cv=c['covariance_e2'][:];axes[0].plot(np.arange(len(cv)),cv/cv[0],'o-',color=COLORS[pitch],label=f'{pitch*10} ms')
  axes[0].set(xlabel='Frame lag',ylabel='Reconstructed-noise correlation');axes[0].legend(frameon=False)
  for pitch in pitches:
   g=pilots[str(pitch)]['bb0nu/noisy_bias'];counts=g['count'][:];axes[1].stairs(counts/sum(counts),100*g['edges'][:],color=COLORS[pitch],label=f'{pitch} mm')
  axes[1].set(xlabel='Reconstructed-charge bias (%)',ylabel='Fraction of double-beta events / bin');axes[1].legend(frameon=False)
  for ax in axes:ax.grid(alpha=.2)
  save(fig,out,'covariance_charge')
  fig,axes=plt.subplots(3,2,figsize=(10.5,8),sharex='col')
  for col,(pitch,key) in enumerate(((3,'1'),(4,'2'))):
   e=pilots[str(pitch)]['bb0nu/examples'][key];edges=e['edge_time_s'][:];x=(edges[:-1]+edges[1:])/2
   axes[0,col].stairs(e['counts'][:],edges,color=COLORS[pitch]);axes[0,col].set_title(f"{pitch} mm · event {attr(e,'event_id')} · pixel ({attr(e,'pixel_ix')}, {attr(e,'pixel_iy')})",fontsize=10)
   axes[1,col].plot(edges,1e3*e['output_at_edges_V'][:],color='#744b83',label='CSA output');axes[1,col].plot(edges,1e3*e['average_at_edges_V'][:],color='#bd6333',ls='--',label='5 ms mean')
   axes[2,col].scatter(x,e['recovered_ions'][:],s=8,color='0.65',label='With noise');axes[2,col].plot(x,e['noiseless_estimated_ions'][:],color='#bd6333',label='CSA');axes[2,col].plot(x,e['linear_estimated_ions'][:],color='#176589',ls=':',label='Linear')
   axes[2,col].set_xlabel('Event-relative time (s)')
   for ax in axes[:,col]:ax.set_xlim(*e['event_window_s'][:]);ax.grid(alpha=.2)
   axes[1,col].legend(frameon=False,fontsize=8);axes[2,col].legend(frameon=False,fontsize=8)
  axes[0,0].set_ylabel('Input ions / bucket');axes[1,0].set_ylabel('Output voltage (mV)');axes[2,0].set_ylabel('Reconstructed ions / bucket');save(fig,out,'histories')
  sources=[{'path':p.decode() if isinstance(p,bytes) else p,'sha256':s.decode() if isinstance(s,bytes) else s} for p,s in zip(h['source_path'][:],h['source_sha256'][:])]
  sources += [{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in (ROOT/'notes/physics/csa_design.tex',ROOT/'notes/physics/figures/Nausica_Image.jpg',ROOT/'notes/physics/references.bib')]
  sources += [{'path':str(source.relative_to(ROOT)),'sha256':hashlib.sha256(source.read_bytes()).hexdigest()},{'path':str(Path(__file__).relative_to(ROOT)),'sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}]
  (out/'sources.json').write_text(json.dumps({'sources':sources},indent=2)+'\n')
 import shutil
 shutil.copy2(ROOT/'notes/physics/csa_design.tex',out/'csa_design.tex')
 shutil.copy2(ROOT/'notes/physics/figures/Nausica_Image.jpg',out/'Nausica_Image.jpg')
 bib=(ROOT/'notes/physics/references.bib').read_text()
 entries=[]
 for key in ('riegler2016','salazar2021'):
  found=re.search(r'(?ms)^@[^\n{]+\{'+key+r',.*?(?=^@|\Z)',bib)
  assert found,key
  entries.append(found.group().strip())
 (out/'references.bib').write_text('\n\n'.join(entries)+'\n')
 print(out)
if __name__=='__main__':main()
