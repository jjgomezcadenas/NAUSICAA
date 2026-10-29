#!/usr/bin/env python3
"""Reproduce milestone-one numbers, tables and figures from detector.toml.

Run from any directory with PYTHONPATH=<repo>/src, or use `make physics`.
No network access or external event sample is required.
"""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[2]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

from pathlib import Path
import argparse
import hashlib
import json
import platform
import os

ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'work/matplotlib-cache'))
os.environ.setdefault('XDG_CACHE_HOME', str(ROOT/'work/cache'))
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from ion_signal.physics import (Detector, E_CHARGE, K_B, square_phi, square_weight,
                               square_z_at_phi, disk_weight, thermal_spread)


def tex_table(path, columns, headings, rows):
    lines = [r'\begin{tabular}{'+columns+'}', r'\toprule',
             ' & '.join(headings)+r' \\', r'\midrule']
    lines.extend(' & '.join(map(str,row))+r' \\' for row in rows)
    lines += [r'\bottomrule', r'\end{tabular}', '']
    path.write_text('\n'.join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT/'config/detector.toml')
    args = parser.parse_args()
    d = Detector.from_toml(args.config)
    out, gen = ROOT/'results', ROOT/'notes/physics/generated'
    out.mkdir(exist_ok=True)
    gen.mkdir(parents=True,exist_ok=True)
    a = d.pad_width_m/2
    v = d.mobility_m2_V_s*d.gap_field_V_m
    ipk = E_CHARGE*v*square_weight(0,a)
    half_z = brentq(lambda u: square_weight(u*a,a)/square_weight(0,a)-.5,0,10)*a
    full_area, err = quad(lambda u:a*square_weight(a*u,a),0,np.inf,
                          epsabs=1e-12,epsrel=1e-12)
    zcheck=np.geomspace(.02*a,20*a,100)
    h=a*1e-5
    derivative_error = float(np.max(np.abs((square_phi(zcheck-h,a)-square_phi(zcheck+h,a))
                                           /(2*h)/square_weight(zcheck,a)-1)))
    vals = dict(open_fraction=d.open_fraction,
                entrance_threshold_V_m=d.bulk_field_V_m/d.open_fraction,
                flux_capacity_ratio=d.open_fraction*d.hole_field_V_m/d.bulk_field_V_m,
                mobility_m2_V_s=d.mobility_m2_V_s,
                thermal_D_m2_s=d.thermal_D_m2_s,
                gap_speed_m_s=v, peak_A=ipk,
                characteristic_time_s=a/v, half_peak_duration_s=half_z/v,
                final_90pct_duration_s=float(square_z_at_phi(.1,a)/v),
                final_99pct_duration_s=float(square_z_at_phi(.01,a)/v),
                central_90pct_duration_s=float((square_z_at_phi(.05,a)-square_z_at_phi(.95,a))/v),
                rise_10_to_90_charge_s=float((square_z_at_phi(.1,a)-square_z_at_phi(.9,a))/v),
                phi_at_gap=float(square_phi(d.gap_m,a)),
                finite_gap_charge_e=float(1-square_phi(d.gap_m,a)),
                disk_peak_A=float(E_CHARGE*v*disk_weight(0,a)),
                integral_full_model_e=full_area, quadrature_error_estimate=err,
                derivative_max_relative_error=derivative_error)
    regions=[]
    for name,L,E in d.regions():
        speed=d.mobility_m2_V_s*E
        sigma=float(thermal_spread(L,E,d.temperature_K))
        regions.append(dict(name=name,length_m=L,field_V_m=E,speed_m_s=speed,
                            transit_s=L/speed,free_diffusion_sigma_m=sigma,
                            local_time_sigma_s=sigma/speed,
                            reduced_field_kV_cm_bar=E/1e5/d.pressure_bar))
    drift_rows=[]
    for L in (.1,.75,d.bulk_reference_length_m):
        sigma=float(thermal_spread(L,d.bulk_field_V_m,d.temperature_K))
        drift_rows.append(dict(length_m=L,transit_s=L/d.bulk_speed_m_s,
                               sigma_m=sigma,time_sigma_s=sigma/d.bulk_speed_m_s))
    distance_rows=[]
    for zmm in (5,1,.5,.25,.1,.05,0):
        z=zmm*1e-3
        distance_rows.append(dict(z_mm=zmm,phi=float(square_phi(z,a)),
                                  current_fA=float(E_CHARGE*v*square_weight(z,a)*1e15)))
    payload={'model':'central square pad in a grounded infinite plane; constant mobility; thermal diffusion',
             'config_sha256':hashlib.sha256(args.config.read_bytes()).hexdigest(),
             'environment':{'python':platform.python_version(),'numpy':np.__version__,
                            'scipy':scipy.__version__,'matplotlib':matplotlib.__version__},
             'constants':{'e_C':E_CHARGE,'k_B_J_K':K_B},'benchmarks':vals,
             'nominal_regions':regions,'bulk_diffusion':drift_rows,'distance_table':distance_rows}
    (out/'physics_benchmarks.json').write_text(json.dumps(payload,indent=2)+'\n')
    np.savetxt(out/'single_ion_distance_samples.txt',
               np.array([[r['z_mm'],r['phi'],r['current_fA']] for r in distance_rows]),
               header='z_mm  weighting_potential  current_fA_before_collection',fmt='%.12g')
    macros={'OpenFraction':f'{d.open_fraction:.5f}',
            'EntranceThreshold':f'{vals["entrance_threshold_V_m"]/1e5:.3f}',
            'FluxRatio':f'{vals["flux_capacity_ratio"]:.3f}',
            'Mobility':f'{d.mobility_m2_V_s/1e-6:.3f}',
            'DiffusionCoefficient':f'{d.thermal_D_m2_s/1e-7:.3f}',
            'PeakCurrent':f'{ipk*1e15:.3f}',
            'DiskPeakCurrent':f'{vals["disk_peak_A"]*1e15:.3f}',
            'FinalNinetyTime':f'{vals["final_90pct_duration_s"]*1e6:.2f}',
            'GapChargeFraction':f'{vals["finite_gap_charge_e"]:.7f}',
            'GapPotential':f'{vals["phi_at_gap"]:.7f}',
            'GapMissingPercent':f'{vals["phi_at_gap"]*100:.4f}',
            'SquareHalfWidth':f'{a*1e3:.3f}',
            'GapSpeed':f'{v:.2f}',
            'ConfigHash':payload['config_sha256'][:12]}
    # Keep the chosen rounding, without padding displayed numbers with zeros.
    macros = {key: value.rstrip("0").rstrip(".")
               if "." in value and "e" not in value.lower() else value
               for key, value in macros.items()}
    (gen/'values.tex').write_text('% Generated by scripts/benchmarks/build_physics.py; do not edit.\n'+
                                 '\n'.join(chr(92)+'newcommand{'+chr(92)+k+'}{'+value+'}' for k,value in macros.items())+'\n')
    tex_table(gen/'parameters.tex','llp{5.6cm}', ['Parameter','Value','Status'],[
        ('Ion charge',r'$+e$','Molecular '+r'$\mathrm{Xe}_2^+$'+' assumption'),
        ('Pressure',f'{d.pressure_bar:g} bar','Detector design'),
        ('Temperature',f'{d.temperature_K:g} K','Room-temperature benchmark'),
        ('Bulk reference length',f'{d.bulk_reference_length_m:g} m','Reference drift; not track length'),
        ('Bulk drift speed',f'{d.bulk_speed_m_s:g} m/s','Operating assumption, not a mobility measurement'),
        ('Hole diameter / pitch',f'{d.hole_diameter_m*1e3:g} / {d.pitch_m*1e3:g} mm','Circular holes on a square pitch'),
        ('IFG thickness / pad gap',f'{d.ifg_thickness_m*1e3:g} / {d.gap_m*1e3:g} mm','Revised geometry'),
        ('Pad width',f'{d.pad_width_m*1e3:g} mm','Square used for analytical benchmark; outline to confirm')])
    tex_table(gen/'regions.tex','lrrrr', ['Region',r'$E$ (kV/cm)',r'$E/p$ (kV/cm/bar)',
                                         r'$v$ (m/s)',r'$t$ (ms)'],
              [[r['name'],f'{r["field_V_m"]/1e5:g}',f'{r["reduced_field_kV_cm_bar"]:.4f}',
                f'{r["speed_m_s"]:g}',f'{r["transit_s"]*1e3:.3f}'] for r in regions])
    tex_table(gen/'diffusion.tex','rrrr',[r'$L$ (m)',r'$t_d$ (s)',r'$\sigma$ (mm)',r'$\sigma_t$ (ms)'],
              [[f'{r["length_m"]:g}',f'{r["transit_s"]:g}',f'{r["sigma_m"]*1e3:.3f}',
                f'{r["time_sigma_s"]*1e3:.3f}'] for r in drift_rows])
    tex_table(gen/'pulse_metrics.tex','lr',['Quantity','Value'],[
        ('Peak immediately before collection',f'{ipk*1e15:.3f} fA'),
        (r'Characteristic time $a/v$',f'{a/v*1e6:.2f} '+r'$\mu$s'),
        ('Time above half peak',f'{half_z/v*1e6:.2f} '+r'$\mu$s'),
        (r'Final interval containing 90\% of charge',f'{vals["final_90pct_duration_s"]*1e6:.2f} '+r'$\mu$s'),
        (r'Final interval containing 99\% of charge',f'{vals["final_99pct_duration_s"]*1e6:.2f} '+r'$\mu$s'),
        (r'Central 90\% interval (5--95\%)',f'{vals["central_90pct_duration_s"]*1e6:.2f} '+r'$\mu$s'),
        (r'Cumulative charge 10--90\% rise',f'{vals["rise_10_to_90_charge_s"]*1e6:.2f} '+r'$\mu$s'),
        (f'Charge induced from {d.gap_m*1e3:g} mm to the metal',f'{vals["finite_gap_charge_e"]:.7f} '+r'$e$')])
    tex_table(gen/'distance.tex','rrr',[r'$z$ (mm)',r'$\phi_w(z)$',r'$i(z)$ (fA)'],
              [[f'{r["z_mm"]:g}',f'{r["phi"]:.6f}',f'{r["current_fA"]:.6f}'] for r in distance_rows])

    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                         'axes.spines.right':False,'pdf.fonttype':42,'axes.labelcolor':'#243247',
                         'text.color':'#243247','axes.titlesize':11})
    def save(fig,name):
        fig.savefig(gen/(name+'.pdf'),bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None})
        fig.savefig(gen/(name+'.png'),dpi=160,bbox_inches='tight')
        plt.close(fig)

    fig,ax=plt.subplots(figsize=(9,4.7))
    ax.set(xlim=(-2.5,5.2),ylim=(-2,20))
    ax.axis('off')
    # Schematic: metal faces are separate from the unspecified IFG body.
    ax.add_patch(Rectangle((-2,5),4,10,facecolor='#edf0f4',edgecolor='none'))
    ax.add_patch(Rectangle((-.25,5),.5,10,facecolor='white',edgecolor='#9aa6b4',linestyle='--'))
    for z in (5,15):
        for x1,x2 in ((-2,-.25),(.25,2)):
            ax.plot([x1,x2],[z,z],color='#495d75',lw=5)
    ax.plot([-.25,.25],[0,0],color='#16816c',lw=6)
    ax.annotate('',(0,.45),(0,18.5),arrowprops={'arrowstyle':'->','color':'#bc5e28','lw':1.8})
    ax.text(.30,18,'Positive-ion drift',color='#bc5e28',fontsize=10)
    ax.text(-1.85,10,'IFG body\nmaterial TBD',fontsize=10,va='center')
    ax.text(2.15,15,f'Gate: {d.gate_voltage_V/1000:+g} kV\nz = {(d.gap_m+d.ifg_thickness_m)*1000:g} mm',va='center',fontsize=10)
    ax.text(2.15,5,f'IC: {d.ic_voltage_V/1000:+g} kV\nz = {d.gap_m*1000:g} mm',va='center',fontsize=10)
    ax.text(.5,-.1,f'Selected pad: {d.pad_voltage_V:g} V, z = 0',va='center',fontsize=10)
    ax.text(-1.85,18.5,f'Bulk: {d.bulk_field_V_m/100:g} V/cm',fontsize=11)
    ax.text(.55,10,f'{d.ifg_thickness_m*1000:g} mm\nnominal {d.hole_field_V_m/1e5:g} kV/cm',va='center',fontsize=10)
    ax.text(.55,2.5,f'{d.gap_m*1000:g} mm\nnominal {d.gap_field_V_m/1e5:g} kV/cm',va='center',fontsize=10)
    ax.text(-2,-1.55,f'Hole diameter {d.hole_diameter_m*1000:g} mm; pitch {d.pitch_m*1000:g} mm. Schematic dimensions.',fontsize=9)
    save(fig,'geometry')

    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(9,3.5),layout='constrained')
    L=np.linspace(0,d.bulk_reference_length_m,400)
    sigma=thermal_spread(L,d.bulk_field_V_m,d.temperature_K)
    ax1.plot(L,sigma*1e3,color='#1769a1',lw=2)
    ax1.set(xlabel='Bulk drift distance (m)',ylabel='One-coordinate RMS spread (mm)',title='Thermal diffusion benchmark')
    ax2.plot(L,sigma/d.bulk_speed_m_s*1e3,color='#16816c',lw=2)
    ax2.set(xlabel='Bulk drift distance (m)',ylabel='Arrival-time RMS spread (ms)',title='Longitudinal spread mapped to time')
    for ax in (ax1,ax2): ax.grid(alpha=.2)
    save(fig,'transport')

    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(9,3.7),layout='constrained')
    z=np.linspace(0,1.5e-3,1500)
    ax1.plot(z*1e3,square_phi(z,a),color='#1769a1',lw=2)
    ax1.set(xlabel='Height above pad (mm)',ylabel='Weighting potential',ylim=(0,1.05),title='Central square-pad potential')
    t=np.linspace(-.55e-3,0,2000)
    ax2.plot(t*1e3,E_CHARGE*v*square_weight(-v*t,a)*1e15,color='#16816c',lw=2,label=f'Square, {2*a*1000:g} mm side')
    ax2.plot(t*1e3,E_CHARGE*v*disk_weight(-v*t,a)*1e15,color='#bc5e28',lw=1.5,ls='--',label=f'Disk, {2*a*1000:g} mm diameter')
    ax2.plot([0,0,.06],[ipk*1e15,0,0],color='#16816c',lw=2)
    ax2.axvspan(-vals['final_90pct_duration_s']*1e3,0,color='#16816c',alpha=.1)
    ax2.set(xlabel='Time relative to collection (ms)',ylabel='Induced current (fA)',xlim=(-.55,.06),
            ylim=(0,1.1*max(ipk,vals['disk_peak_A'])*1e15),title='Constant-speed central approach')
    ax2.legend(fontsize=8,loc='upper left',frameon=False)
    for ax in (ax1,ax2): ax.grid(alpha=.2)
    save(fig,'single_ion')
    print(json.dumps({'generated':str(gen),'peak_fA':ipk*1e15,
                      'final_90pct_us':vals['final_90pct_duration_s']*1e6,
                      'bulk_sigma_mm':drift_rows[-1]['sigma_m']*1e3,
                      'integral_e':full_area,'derivative_relative_error':derivative_error},indent=2))


if __name__=='__main__':
    main()
