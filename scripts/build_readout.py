#!/usr/bin/env python3
"""Generate analytical current/CSA benchmarks, offline HTML and SPICE inputs."""
from pathlib import Path
import hashlib
import json
import re
import tomllib
import platform
import numpy as np
import scipy
from ion_signal.physics import Detector
from ion_signal.readout import make_waveform, csa_response, waveform_metrics

ROOT=Path(__file__).resolve().parents[1]


def plot_indices(wave, voltages):
    """Keep uniform time samples, charge quantiles, and actual voltage peaks."""
    n=len(wave.charge_per_bin_C)
    c=np.cumsum(wave.charge_per_bin_C)
    quantiles=np.searchsorted(c,np.linspace(0,1,501)*c[-1])+1
    points=np.r_[np.linspace(0,n,min(n+1,1001)).astype(int),
                 np.clip(quantiles,0,n),
                 [np.argmin(v) for v in voltages],
                 np.argmax(wave.current_A),0,n]
    return np.unique(points)


def export_spice(wave, path, relative_tolerance):
    """Piecewise-linear current with explicit charge error, no rescaling.

    Bin-centre interpolation plus flat half-bin ends preserves the charge on
    the full grid. Coarser exports are accepted only after charge comparison.
    SPICE time starts at zero; physical time offset is recorded in the JSON.
    Positive current is defined as flowing into the CSA input node.
    """
    t=(wave.edges_s[:-1]+wave.edges_s[1:])/2
    t=np.r_[wave.edges_s[0],t,wave.edges_s[-1]]-wave.edges_s[0]
    i=np.r_[wave.current_A[0],wave.current_A,wave.current_A[-1]]
    stride=max(1,len(t)//5000)
    while True:
        ix=np.unique(np.r_[np.arange(0,len(t),stride),np.argmax(i),len(t)-1])
        q=float(np.trapezoid(i[ix],t[ix]))
        error=(q-wave.total_charge_C)/wave.total_charge_C
        if abs(error)<=relative_tolerance or stride==1:
            break
        stride=max(1,stride//2)
    if abs(error)>relative_tolerance:
        raise RuntimeError('SPICE export failed charge check')
    np.savetxt(path,np.c_[t[ix],i[ix]],fmt='%.15e')
    return dict(file=str(path.relative_to(ROOT)),time_offset_s=float(wave.edges_s[0]),
                exported_charge_C=q,relative_charge_error=error,points=len(ix),
                convention='two columns: time_s current_A; positive current into input; '
                'PWL interpolates sampled current, no charge renormalization')


def main():
    detector_path=ROOT/'config/detector.toml'
    config_path=ROOT/'config/readout.toml'
    config=tomllib.loads(config_path.read_text())
    d=Detector.from_toml(detector_path)
    el=config['electronics']; num=config['numerics']
    cf=el['feedback_capacitance_F']
    out=ROOT/'results/readout'; spice=out/'spice'
    spice.mkdir(parents=True,exist_ok=True)
    web=ROOT/'notes/physics/readout'
    payload=dict(model='Linear CSA driven by finite-gap square-pad currents',
                 detector_sha256=hashlib.sha256(detector_path.read_bytes()).hexdigest(),
                 readout_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest(),
                 environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
                 electronics=el,provenance=config['provenance'],cases=[])
    summaries=[]
    for case in config['cases']:
        kwargs={k:v for k,v in case.items() if k not in ('id','label')}
        wave,meta=make_waveform(d,**kwargs,max_step_s=num['max_step_s'],
                                tail_sigmas=num['gaussian_tail_sigmas'])
        times=el['discharge_times_s']+[np.inf]
        voltages=[csa_response(wave,cf,tau) for tau in times]
        ix=plot_indices(wave,voltages)
        # Current values are bin averages shown at the corresponding left edge.
        currents=np.r_[wave.current_A,0.]
        responses=[]
        for tau,v in zip(times,voltages):
            metrics=[waveform_metrics(wave,v,cf,el['sample_interval_s'],phase)
                     for phase in el['sample_phases_s']]
            responses.append(dict(discharge_time_s=None if np.isinf(tau) else tau,
                                  feedback_resistance_ohm=None if np.isinf(tau) else tau/cf,
                                  voltage_V=v[ix].tolist(),metrics=metrics))
        record=dict(id=case['id'],label=case['label'],metadata=meta,
                    time_s=wave.edges_s[ix].tolist(),current_A=currents[ix].tolist(),
                    peak_current_A=float(np.max(wave.current_A)),responses=responses,
                    spice=export_spice(wave,spice/(case['id']+'.pwl'),num['spice_charge_relative_tolerance']))
        payload['cases'].append(record)
        summaries.append({k:v for k,v in record.items() if k not in ('time_s','current_A','responses')}
                         | {'responses':[{k:v for k,v in r.items() if k!='voltage_V'} for r in responses]})
    (out/'benchmarks.json').write_text(json.dumps({**payload,'cases':summaries},indent=2,allow_nan=False)+'\n')
    source=web/'rogelio_reference.html'
    svg=re.search(r'<svg\b[\s\S]*?</svg>',source.read_text()).group(0)
    template=(ROOT/'scripts/templates/readout.html').read_text()
    html=template.replace('__MODEL_DATA__',json.dumps(payload,separators=(',',':'),allow_nan=False)).replace('__SCHEMATIC__',svg)
    (web/'analytical_model.html').write_text(html)
    print('Created',web/'analytical_model.html')
    print('SPICE input current files:',spice)
    for case in summaries:
        response=case['responses'][3]['metrics'][0]
        print(case['id'], 'tau=0.1 s peak=',response['peak_output_magnitude_V'],
              'V; retention=',response['peak_fraction_of_ideal'])


if __name__=='__main__':
    main()
