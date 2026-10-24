#!/usr/bin/env python3
"""Reproduce supplied CSA and generate/run matched ion transients and noise.

Native LTspice 17.2.4 batch takes a text deck. Original .asc/.inc are preserved.
Derived decks retain M5's supplied dimensions; unused parameters are not fixed.
"""
from pathlib import Path
import argparse, hashlib, json, subprocess, tomllib
from build_realistic_csa import spice_comparator
import numpy as np
from ion_signal.spice import schematic_netlist, read_raw, sampled_charge
from ion_signal.physics import Detector, E_CHARGE, square_phi
from ion_signal.realistic_csa import uniform_ion_waveform, pwl_current, CSA, response

ROOT=Path(__file__).resolve().parents[1]

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def run(simulator,deck):
    result=subprocess.run([str(simulator),'-b',str(deck.resolve())],cwd=deck.parent,capture_output=True,timeout=120)
    log=deck.with_suffix('.log')
    if log.exists():
        b=log.read_bytes();text=b.decode('utf-16le' if b[1:2]==b'\x00' else 'utf-8',errors='replace')
        deck.with_suffix('.run.txt').write_text('\n'.join(line.rstrip() for line in text.splitlines()).rstrip()+'\n')
    if result.returncode: raise RuntimeError(f'{deck}: simulator exit {result.returncode}; inspect .log')
    return read_raw(deck.with_suffix('.raw'))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--simulator',type=Path,required=True)
    ap.add_argument('--stage',choices=['original','study','all'],default='all')
    ap.add_argument('--out',type=Path,default=ROOT/'results/ltspice')
    args=ap.parse_args(); args.simulator=args.simulator.resolve()
    out=args.out.resolve(); decks=out/'decks'; decks.mkdir(parents=True,exist_ok=True)
    inc=ROOT/'LTSpice/180nm_professional.inc'; asc=ROOT/'LTSpice/CSAfolded180nm.asc'
    # Normalize only the derived library's line endings/trailing whitespace.
    (decks/inc.name).write_text('\n'.join(line.rstrip() for line in inc.read_text().splitlines()).strip()+'\n')
    elements,directives=schematic_netlist(asc)
    # Independent, explicit connectivity expectations from the supplied wiring.
    expected={'M1':['s','in','med','Vdda'],'M2':['Vdda','Vnmos','s','0'],
      'M3':['med','Vnload','0','0'],'M4':['med','Vncasc','out','0'],
      'M5':['vd6','Vp5','out','Vdda'],'M6':['Vdda','Vp6','vd6','Vdda'],
      'Mf':['out','Vdis','in','0'],'Cf':['out','in'],'Cd1':['in','0']}
    for line in elements:
        parts=line.split()
        if parts[0] in expected and parts[1:1+len(expected[parts[0]])]!=expected[parts[0]]:
            raise AssertionError(f'Connectivity mismatch: {line}')
    base='* Derived from supplied CSAfolded180nm.asc; original device sizes and biases\n'+'\n'.join(elements)+'\n'
    original=decks/'original.cir';original.write_text(base+'\n'.join(directives)+'\n.end\n')
    report={'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in (asc,inc)},'simulator':str(args.simulator),
      'assumptions':['Supplied device models and M5 dimensions unchanged.',
        '5 pF external pad capacitance PLUS intrinsic transistor capacitances.',
        'Central square pad in grounded plane; uniform drift in finite 5 mm gap; no renormalization.',
        'Uniform mean arrivals, not a sampled random track.',
        'Paper estimator keeps Cf=3 fF and nominal tau=0.3 s; no fitted gain correction.',
        'Stationary model-default transistor noise; no added resistor noise; model provenance/flicker calibration unconfirmed.']}
    report['code_and_config_sha256']={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__).resolve(),ROOT/'src/ion_signal/spice.py',ROOT/'src/ion_signal/realistic_csa.py',ROOT/'config/detector.toml',ROOT/'config/realistic_csa.toml')}
    report['simulator_sha256']=sha(args.simulator)
    report['simulator_version']='LTspice 17.2.4 for MacOS (native batch)'
    if args.stage in ('original','all'):
        raw_data,raw_header=run(args.simulator,original)
        report['original_default_solver']={'points':len(raw_data['time']),'minimum_output_V':float(raw_data['V(out)'].min())}
        precise=decks/'original_precise.cir'
        precise.write_text(original.read_text().replace('.tran 10','.tran 0 10 0 100u').replace('.end', '.options reltol=1e-7 abstol=1e-22 vntol=1e-11 numdgt=15 plotwinsize=0\n.end'))
        data,header=run(args.simulator,precise)
        t=data['time'];v=data['V(out)'];i=data['I(Is1)']
        report['original']={'baseline_output_V':float(v[0]),'minimum_output_V':float(v.min()),
          'final_output_V':float(v[-1]),'integrated_charge_C':float(np.trapezoid(i,t)),
          'source_charge_ions':float(np.trapezoid(i,t)/E_CHARGE),'raw_header':header}
        np.savetxt(out/'original_transient.csv.gz',np.c_[t,i,data['V(in)'],v][::10],delimiter=',',header='time_s,current_A,input_V,output_V',comments='')
        print(json.dumps(report['original'],indent=2))
    if args.stage in ('study','all'):
        cfg=tomllib.loads((ROOT/'config/realistic_csa.toml').read_text()); csa=CSA(**cfg['electronics'])
        det=Detector.from_toml(ROOT/'config/detector.toml');dt=cfg['readout']['sample_interval_s']
        records=[]
        common=[d for d in directives if not d.lower().startswith(('.tran','.meas'))]
        # Explicit temperature; tight current/voltage tolerances for single-ion signals.
        options='.temp 20\n.options reltol=1e-5 abstol=1e-18 vntol=1e-9 numdgt=15 plotwinsize=0\n'
        def body(source):
            lines=[]
            for line in elements:
                if line.startswith('Is1 '): lines.append(source)
                elif line.startswith(('B1 ','Vion1 ')): continue
                elif line.startswith('Cd1 '): lines.append('Cd1 in 0 5p')
                else: lines.append(line)
            return '* NAUSICA matched benchmark; supplied transistor sizes/biases unchanged\n'+'\n'.join(lines+common)+'\n'+options
        for refinement in (1,2,3):
            step=5e-6/min(refinement,2);offset=.002
            if refinement==3:
                options=options.replace('reltol=1e-5 abstol=1e-18 vntol=1e-9','reltol=1e-6 abstol=1e-19 vntol=1e-10')
            wave=uniform_ion_waveform(det,0.,dt,start_s=-offset,end_s=.03,max_step_s=step)
            t,i=pwl_current(wave)
            np.savetxt(decks/f'uniform_r{refinement}.pwl',np.c_[t,i],fmt='%.17e')
            q=np.trapezoid(i,t); expected_q=E_CHARGE*(1-square_phi(det.gap_m,det.pad_width_m/2))
            if abs(q/expected_q-1)>1e-9: raise AssertionError('Finite-gap PWL charge error')
            counts=(0,1,68,100,300,1000,1574,1575,1577,1578,1579,1580,1600,3000,10000) if refinement<3 else (0,1,1577,1578,1579)
            for n in counts:
                pwl=decks/f'input_N{n}_r{refinement}.pwl'; np.savetxt(pwl,np.c_[t,n*i],fmt='%.17e')
                deck=decks/f'uniform_N{n}_r{refinement}.cir'
                deck.write_text(body(f'Is1 0 in PWL FILE="{pwl.name}"')+
                  f'.tran 0 {t[-1]:.17g} 0 {step:.17g}\n.save V(in) V(out) I(Is1)\n.end\n')
                data,header=run(args.simulator,deck);time=data['time'];v=data['V(out)'];vi=data['V(in)']
                # Matched zero-input run removes bias drift, including solver settling.
                if n==0: baseline=(time.copy(),v.copy());delta=v-v
                else: delta=v-np.interp(time,*baseline)
                qhat=sampled_charge(time,delta,[offset,offset+dt],csa.feedback_capacitance_F,csa.nominal_discharge_time_s)[0]
                cumulative=np.r_[0.,np.cumsum((n*i[1:]+n*i[:-1])/2*np.diff(t))]
                qframe=float(np.diff(np.interp([offset,offset+dt],t,cumulative))[0])
                r=dict(ions=n,refinement=refinement,pwl_step_s=step,maximum_timestep_s=step,
                    baseline_output_V=float(v[0]),input_charge_C=qframe,estimated_charge_C=float(qhat),
                    gain_e_per_ion=float(qhat/E_CHARGE/n) if n else None,
                    recovered_fraction=float(qhat/qframe) if n else None,
                    peak_output_change_V=float(-delta.min()),minimum_output_V=float(v.min()))
                records.append(r);print(json.dumps(r),flush=True)
                np.savetxt(out/f'uniform_N{n}_r{refinement}.csv.gz',np.c_[time-offset,data['I(Is1)'],vi,v,delta],delimiter=',',header='physical_time_s,current_A,input_V,output_V,baseline_subtracted_output_V',comments='')
        report['transients']=records
        # Same physical waveform in the analytical behavioral comparator.
        deck=decks/'behavioral.cir'
        spice_comparator(deck,'uniform_r2.pwl',csa,t[-1],2.5e-6)
        deck.write_text(deck.read_text().replace('.end','.options plotwinsize=0\n.end'))
        bd,bh=run(args.simulator,deck)
        bq=float(sampled_charge(bd['time'],bd['V(out)'],[offset,offset+dt],csa.feedback_capacitance_F,csa.nominal_discharge_time_s)[0])
        pv=response(wave,csa)[:,1]
        pq=float(sampled_charge(wave.edges_s,pv,[0,dt],csa.feedback_capacitance_F,csa.nominal_discharge_time_s)[0])
        report['behavioral_comparison']={'spice_charge_C':bq,'python_charge_C':pq,'relative_difference':abs(bq/pq-1)}
        if abs(bq/pq-1)>1e-4: raise ArithmeticError('Behavioral SPICE/Python mismatch')
        # Full closed-loop transimpedance for future amplifier fits.
        deck=decks/'ac.cir'
        deck.write_text(body('Is1 0 in 0 AC 1')+'.ac dec 400 .1 1Meg\n.save V(in) V(out) I(Is1)\n.end\n')
        ad,ah=run(args.simulator,deck)
        np.savetxt(out/'ac.csv',np.column_stack([ad['frequency'].real,ad['V(in)'].real,ad['V(in)'].imag,ad['V(out)'].real,ad['V(out)'].imag]),delimiter=',',header='frequency_Hz,input_transimpedance_real_ohm,input_transimpedance_imag_ohm,output_transimpedance_real_ohm,output_transimpedance_imag_ohm',comments='')
        deck=decks/'operating_point.cir'
        deck.write_text(body('Is1 0 in 0 AC 1')+'.op\n.end\n')
        od,oh=run(args.simulator,deck)
        report['operating_point']={name:float(value[0]) for name,value in od.items()}
        report['supply_power_W']=-1.8*report['operating_point']['I(Vdd)']

        # The two-sample noise filter oscillates on 1/dt=100 Hz. A dense log grid
        # alone cannot resolve its high-frequency oscillations; analysis resamples
        # PSD with explicit phase integration below (see analyze_ltspice_study).
        for density in (200,400):
            deck=decks/f'noise_{density}.cir'
            deck.write_text(body('Is1 0 in 0 AC 1')+f'.noise V(out) Is1 dec {density} .1 1Meg\n.end\n')
            data,header=run(args.simulator,deck)
            print(header,flush=True)
            names=list(data);f=data[names[0]].real
            np.savetxt(out/f'noise_{density}.csv',np.column_stack([f,*[data[x].real for x in names[1:]]]),delimiter=',',header=','.join(names),comments='')
            report[f'noise_{density}_columns']=names
    (out/'study.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
