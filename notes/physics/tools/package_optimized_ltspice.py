#!/usr/bin/env python3
"""Freeze the selected optimized CSA and readout settings for Vicente's review."""
from pathlib import Path
import hashlib,json,re,shutil,subprocess,sys
import h5py
ROOT=Path(__file__).resolve().parents[3]
DEST=ROOT/'LTSpice/optimized_2026-10-29'
RESULT=ROOT/'results/ltspice/feedback_optimization'
SELECTED='length_400um_gate0_rounded'
DECKS=RESULT/'decks'/SELECTED
records=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def copy(source,relative):
 target=DEST/relative;target.parent.mkdir(parents=True,exist_ok=True)
 if target.exists() and target.read_bytes()!=source.read_bytes():raise ValueError(f'Preserving {target}')
 shutil.copyfile(source,target)
 records.append({'file':str(relative),'source':str(source.relative_to(ROOT)),'sha256':sha(source),'bytes':source.stat().st_size})
def drawings():
 sys.path.insert(0,str(ROOT/'src'))
 from ion_signal.spice import schematic_netlist
 original=ROOT/'LTSpice/Spice2/v3/CSAfolded180nm-2_Astronio.asc'
 template='\n'.join(line for line in original.read_text(encoding='latin1').splitlines() if not line.startswith(('TEXT ','DATAFLAG ')))+'\n'
 copy(DECKS/'corrected.inc',Path('schematic/180nm_professional.inc'))
 cases=[('NAUSICA_40ms_400ions','40ms_single_N400_r2'),('NAUSICA_20ms_120ions','20ms_single_N120_r2'),('NAUSICA_40ms_baseline','40ms_single_N0_r2'),('NAUSICA_20ms_baseline','20ms_single_N0_r2'),('NAUSICA_noise','noise_400')]
 for name,reference in cases:
  deck=DECKS/(reference+'.cir');lines=deck.read_text().splitlines()
  elements={x.split()[0]:x for x in lines if x and not x.startswith(('*','.'))}
  def symbol(match):
   block=match.group();instance=re.search(r'^SYMATTR InstName (\S+)',block,re.M)[1]
   if instance in ('Vion1','B1'):return ''
   if instance in ('Is1','Vnload','Vp1','Vf1','Cd1'):
    block=re.sub(r'^SYMATTR Value .*$', 'SYMATTR Value '+elements[instance].split(maxsplit=3)[3],block,flags=re.M)
    block=re.sub(r'^SYMATTR Value2 .*\n?','',block,flags=re.M)
   return block
  text=re.sub(r'^SYMBOL .*?(?=^SYMBOL |\Z)',symbol,template,flags=re.M|re.S)
  text=text.replace('SHEET 1 1220 900\n','SHEET 1 1220 900\nWIRE -688 48 -688 -32\n')
  directives=[x.replace('.lib corrected.inc','.lib 180nm_professional.inc') for x in lines if x.startswith('.') and x!='.end']
  text+=''.join(f'TEXT -888 {-456+24*i} Left 2 !{line}\n' for i,line in enumerate(directives))
  target=DEST/'schematic'/(name+'.asc');target.write_text(text,encoding='ascii')
  actual,commands=schematic_netlist(target)
  norm=lambda xs:{x.split()[0]:' '.join(x.lower().split()) for x in xs}
  assert norm(actual)==norm(elements.values()) and commands==directives,name
  records.append({'file':str(target.relative_to(DEST)),'source':str(original.relative_to(ROOT)),'source_sha256':sha(original),'electrical_reference':str(deck.relative_to(ROOT)),'electrical_reference_sha256':sha(deck),'sha256':sha(target),'bytes':target.stat().st_size,'validation':'Parsed elements and directives equal the archived optimized netlist; model filename localized.'})
def main():
 DEST.mkdir(parents=True,exist_ok=True)
 for deck in sorted(DECKS.glob('*.cir')):
  copy(deck,Path('netlists')/deck.name)
  log=deck.with_suffix('.run.txt')
  if log.exists():copy(log,Path('netlists')/log.name)
  for dep in re.findall(r'PWL FILE="([^"]+)"',deck.read_text()):
   relative=Path('inputs')/Path(dep).name
   if not any(r['file']==str(relative) for r in records):copy(deck.parent/dep,relative)
 copy(DECKS/'corrected.inc',Path('netlists/corrected.inc'))
 drawings()
 for name in ('feedback_optimization.toml','feedback_combined.toml','detector.toml'):copy(ROOT/'config'/name,Path('reference')/name)
 for name in ('report.tex','noise_comparison.pdf','noise_comparison.png','provenance.json'):copy(ROOT/'latex_generated/feedback_optimization'/name,Path('report')/name)
 source=RESULT/'analysis.h5'
 with h5py.File(source) as h:
  g=h[SELECTED];a=dict(g.attrs)
  parameters={'selected_case':SELECTED,'Mf_width_m':float(a['width_m']),'Mf_length_m':float(a['length_m']),'Mf_gate_V':float(a['gate_V']),'M4_gate_V':float(a['m4_gate_V']),'supply_V':1.8,'temperature_K':293.15,'external_input_capacitance_F':1.25e-12,'feedback_capacitance_F':3e-15,'selected_averaging_s':.005,'noise_low_Hz':.1,'noise_high_Hz':1e6,'model_sha256':sha(DECKS/'corrected.inc'),'analysis_source':str(source.relative_to(ROOT)),'analysis_source_sha256':sha(source),'cases':{}}
  for key in ('20ms','40ms'):
   f=g[key];j=int(f.attrs['best_index'])-1
   parameters['cases'][key]={'frame_s':float(f.attrs['frame_s']),'mean_ions':int(f.attrs['mean_ions']),'discharge_time_s':float(f.attrs['tau_s']),'alpha':float(f.attrs['alpha']),**{field:float(f[field][j]) for field in ('gain_V_per_ion','noise_rms_V','ENC_e_rms','mean_signal_V','mean_snr','previous_frame_fraction','following_frame_fraction')}}
  target=DEST/'parameters.h5'
  with h5py.File(target,'w') as saved:
   saved.attrs['source']=str(source.relative_to(ROOT));saved.attrs['source_sha256']=sha(source)
   h.copy(g,saved,name=SELECTED);h.copy(h['reference'],saved,name='reference')
  records.append({'file':'parameters.h5','source':str(source.relative_to(ROOT)),'source_sha256':sha(source),'sha256':sha(target),'bytes':target.stat().st_size,'transformation':'Copy selected circuit and reference groups unchanged.'})
 (DEST/'parameters.json').write_text(json.dumps(parameters,indent=2)+'\n')
 for deck in list(DEST.rglob('*.cir'))+list(DEST.rglob('*.asc')):
  for dep in re.findall(r'(?:\.lib\s+|PWL FILE=)"?([^\s"\n]+)',deck.read_text(),re.I):assert (deck.parent/dep).is_file(),(deck,dep)
 (DEST/'README.md').write_text(README)
 (DEST/'manifest.json').write_text(json.dumps({'selected':SELECTED,'base_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'generator':str(Path(__file__).resolve().relative_to(ROOT)),'generator_sha256':sha(Path(__file__)),'files':records},indent=2)+'\n')
 print(DEST,';',len(records),'source-tracked files; 5 schematic/netlist comparisons and all local dependencies checked.')
README='''# Optimized NAUSICA simulation — 29 October 2026

Snapshot for Vicente Herrero and Rogelio Palomo. The selected circuit is
`length_400um_gate0_rounded` from the feedback/averaging optimization. The paper's
previous circuit remains in the separate `LTSpice/paper_2026-10-28/` snapshot.

## Open and run

- `schematic/NAUSICA_40ms_400ions.asc`: complete editable circuit, 400 ions / 40 ms.
- `schematic/NAUSICA_20ms_120ions.asc`: 120 ions / 20 ms.
- `schematic/NAUSICA_40ms_baseline.asc` and `NAUSICA_20ms_baseline.asc`: matching zero-input references.
- `schematic/NAUSICA_noise.asc`: total output-noise spectrum, 0.1 Hz–1 MHz.
- `schematic/180nm_professional.inc`: exact transistor model card used in these runs.

Keep the entire directory together. Each current source reads a PWL file from
`../inputs/`: first column time in seconds, second column current in amperes.
Open an .asc in LTspice and press Run. Plot V(out), V(in) and I(Is1).
The exact batch netlists, readable run logs and their local model card are in
`netlists/`. Both temporal resolutions and all tested input populations are
included. The graphical circuits are electrically checked against those netlists.

## Selected parameters

Mf W=0.36 um, L=400 um, gate=0 V; M4 gate=0.6164 V. Other transistor dimensions
and biases are those of the selected paper circuit and are explicit in every
netlist. Supply 1.8 V, temperature 20 C, external input capacitance 1.25 pF plus
intrinsic model capacitances, feedback capacitor 3 fF. Model-card thermal and
flicker noise are included. Numerical settings disable waveform compression.

The readout averages the output over centred 5 ms windows at boundaries 20 or
40 ms apart, then takes D_n = -(Vbar_n - alpha*Vbar_(n-1)), with
alpha=exp(-frame/discharge_time). Averaging and differencing are performed in
the Julia analysis, after subtracting the matched zero-input output baseline.
They are applied to the signal and to the noise spectrum consistently.
`parameters.json` and `parameters.h5` contain the fitted discharge times, gains,
noise, ENC, neighbouring-frame responses and source hashes from the saved study.
`report/report.tex` and its scientific figures describe the full scan and checks.

The current uses uniform mean collection times within the specified bucket,
a central trajectory through the final 5 mm at 3 m/s, and the 0.5 mm square-pad
weighting benchmark in an infinite grounded plane with uniform gas permittivity.
The selected pad has unit weighting potential, the surrounding plane and
infinity zero. The finite-path integral e[1-phi_w(g)] is retained. Each isolated
packet begins collection at simulation time 2 ms. Both 2.5 and 5 us maximum
transient steps are included; the noise scan uses 200 and 400 points per decade.

The longer feedback retention will be tested with full sparse pad histories,
including zero-charge buckets and high-charge blob deposits. Those tests and ADC
range/saturation handling are the next study. The netlists here preserve the
already completed isolated-packet calculation, including inputs through 2000 ions.

## Reproduction

From the repository root:

    PYTHONPATH=src python scripts/run_feedback_pipeline.py --julia /path/to/julia --simulator /path/to/LTspice
    python notes/physics/tools/package_optimized_ltspice.py

Raw acquisitions and all scan variants are in `results/ltspice/feedback_optimization/`.
`manifest.json` identifies source files and SHA-256 checksums. The circuit study
used macOS LTspice 17.2.4.
'''
if __name__=='__main__':main()
