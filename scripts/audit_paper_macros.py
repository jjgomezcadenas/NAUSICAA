#!/usr/bin/env python3
"""Check every generated numerical macro loaded by nausica.tex against its source."""
from pathlib import Path
import hashlib,json,re,tomllib
import h5py
ROOT=Path(__file__).resolve().parents[1];NOTE=ROOT/'notes/physics'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def active(s):return '\n'.join(re.split(r'(?<!\\)%',line)[0] for line in s.splitlines())
def main():
    text=active((NOTE/'nausica.tex').read_text())
    files=[NOTE/(x+'.tex') for x in re.findall(r'\\input\{(generated/[^}]+)\}',text) if x.endswith('_values') or x=='generated/values']
    macros={};locations={}
    for file in files:
        for name,value in re.findall(r'\\newcommand\{\\(\w+)\}\{([^{}]+)\}',active(file.read_text())):
            assert name not in macros,f'Duplicate macro {name}'
            macros[name]=value;locations[name]=str(file.relative_to(ROOT))
    checks={}
    def check(name,expected,source):
        value=macros[name]
        if isinstance(expected,str):assert value==expected,(name,value,expected)
        else:
            places=len(value.split('.')[1]) if '.' in value else 0
            tolerance=.50001*10**(-places) if places else 1e-9
            # Integer power is deliberately rounded to whole microwatts.
            if name=='SpicePowerUw':tolerance=.50001
            assert abs(float(value)-float(expected))<=tolerance,(name,value,expected)
        checks[name]={'displayed':value,'expected':expected,'source':source,'macro_file':locations[name]}
    physfile=ROOT/'results/physics_benchmarks.json';phys=json.loads(physfile.read_text());b=phys['benchmarks']
    detfile=ROOT/'config/detector.toml';det=tomllib.loads(detfile.read_text())
    assert phys['config_sha256']==sha(detfile)
    mapping={'OpenFraction':('open_fraction',1),'EntranceThreshold':('entrance_threshold_V_m',1e-5),'FluxRatio':('flux_capacity_ratio',1),'Mobility':('mobility_m2_V_s',1e6),'DiffusionCoefficient':('thermal_D_m2_s',1e7),'PeakCurrent':('peak_A',1e15),'DiskPeakCurrent':('disk_peak_A',1e15),'FinalNinetyTime':('final_90pct_duration_s',1e6),'GapChargeFraction':('finite_gap_charge_e',1),'GapPotential':('phi_at_gap',1),'GapMissingPercent':('phi_at_gap',100),'GapSpeed':('gap_speed_m_s',1)}
    for name,(key,scale) in mapping.items():check(name,b[key]*scale,f'results/physics_benchmarks.json:benchmarks/{key}')
    check('SquareHalfWidth',det['geometry']['pad_width_m']*500,'config/detector.toml:geometry/pad_width_m divided by two, mm')
    check('ConfigHash',sha(detfile)[:12],'config/detector.toml:sha256 prefix')
    readfile=ROOT/'results/readout/benchmarks.json';readout=json.loads(readfile.read_text())
    assert readout['detector_sha256']==sha(detfile)
    readcfg=tomllib.loads((ROOT/'config/readout.toml').read_text())
    check('SignalPadWidth',det['geometry']['pad_width_m']*1000,'config/detector.toml:geometry/pad_width_m, mm')
    check('SignalGap',det['geometry']['gap_m']*1000,'config/detector.toml:geometry/gap_m, mm')
    region=next(x for x in phys['nominal_regions'] if abs(x['length_m']-det['geometry']['gap_m'])<1e-12)
    check('SignalGapField',region['field_V_m']/1e5,'results/physics_benchmarks.json:gap field, kV/cm')
    check('ReadoutCapacitance',readcfg['electronics']['feedback_capacitance_F']*1e15,'config/readout.toml:feedback_capacitance_F, fF')
    case=next(x for x in readout['cases'] if x['id']=='single_ion')
    check('SingleIonIdealVoltage',case['responses'][0]['metrics'][0]['ideal_output_magnitude_V']*1e6,'results/readout/benchmarks.json:single_ion ideal_output_magnitude_V, microvolt')
    cfg=tomllib.loads((ROOT/'config/discharge_calibration.toml').read_text())
    datafile=ROOT/cfg['output']/'paper_readout.h5'
    with h5py.File(datafile) as h:
        a=dict(h.attrs);snrs=dict(zip(h['counts'][:],h['snr'][:]))
        fields={'SpiceFrameMs':('frame_s',1e3),'SpiceAverageUs':('averaging_s',1e6),'SpiceCorrectionSeconds':('tau_s',1),'SpiceLowHz':('frequency_low_Hz',1),'SpiceHighMHz':('frequency_high_Hz',1e-6),'SpiceGainUv':('gain_V_per_ion',1e6),'SpiceNoiseMv':('noise_rms_V',1e3),'SpiceEncIons':('ENC_e_rms',1),'SpiceMinimumIons':('verified_N_snr5',1),'SpicePreviousPercent':('previous_frame_residual_fraction',100),'SpiceFollowingPercent':('following_frame_fraction',100)}
        for name,(key,scale) in fields.items():check(name,float(a[key])*scale,f'{datafile.relative_to(ROOT)}:attributes/{key}')
        n=int(a['verified_N_snr5'])
        check('SpiceBelowIons',n-1,f'{datafile.relative_to(ROOT)}:threshold lower bracket')
        for name,count in [('SpiceBelowSnr',n-1),('SpiceAboveSnr',n)]:check(name,float(snrs[count]),f'{datafile.relative_to(ROOT)}:snr for N={count}')
        assert snrs[n-1]<5<snrs[n]
        enc=h['averaging/ENC_e_rms'][:]
        check('SpiceScanFirstEnc',float(enc[0]),f'{datafile.relative_to(ROOT)}:averaging first point')
        check('SpiceScanLastEnc',float(enc[-1]),f'{datafile.relative_to(ROOT)}:averaging last point')
        assert float(a['flicker_variance_fraction'])>0
    deck=ROOT/cfg['circuit'];lines={x[0]:x for line in deck.read_text().splitlines() if (x:=line.split())}
    def scalar(value):
        scales={'f':1e-15,'p':1e-12,'n':1e-9,'u':1e-6,'m':1e-3}
        return float(value[:-1])*scales[value[-1]] if value[-1] in scales else float(value)
    for name,component,scale in [('SpiceCfFf','Cf',1e15),('SpicePadPf','Cd1',1e12),('SpiceSupplyV','Vdd',1)]:check(name,scalar(lines[component][3])*scale,f'{deck.relative_to(ROOT)}:{component}')
    check('SpiceTemperatureK',scalar(lines['.temp'][1])+273.15,f'{deck.relative_to(ROOT)}:.temp converted to K')
    assert float(macros['ReadoutCapacitance'])==float(macros['SpiceCfFf'])
    circuit=json.loads((NOTE/'generated/spice_circuit_provenance.json').read_text())
    check('SpiceBaselineV',circuit['operating_point']['dc_output_V'],'generated/spice_circuit_provenance.json:operating_point/dc_output_V')
    check('SpicePowerUw',circuit['operating_point']['quiescent_power_W']*1e6,'generated/spice_circuit_provenance.json:operating_point/quiescent_power_W, microwatt')
    assert set(checks)==set(macros),f'Unchecked macros: {set(macros)-set(checks)}'
    body=text[text.index('\\begin{document}'):]+active((NOTE/'appendix_spice_readout.tex').read_text())
    for name in checks:checks[name]['used_in_active_prose']=bool(re.search(r'\\'+name+r'\b',body))
    assert 'zero contribution in the transistor' not in body
    assert not re.search(r'\\SI\{10\}\{ms\}',body)
    payload={'status':'passed','checked_macro_count':len(checks),'macros':checks,
             'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files+[NOTE/'nausica.tex',NOTE/'appendix_spice_readout.tex',datafile,physfile,readfile,Path(__file__).resolve()]}}
    (datafile.parent/'macro_audit.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(f'Checked all {len(checks)} loaded numerical macros; all agree with source calculations.')
if __name__=='__main__':main()
