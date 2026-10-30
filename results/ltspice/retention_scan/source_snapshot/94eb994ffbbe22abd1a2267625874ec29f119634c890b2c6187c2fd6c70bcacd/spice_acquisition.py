"""Acquire native LTspice waveforms and preserve readable simulator logs."""
import subprocess
from .spice import read_raw

def run(simulator,deck,*,timeout_s=120):
    result=subprocess.run([str(simulator),'-b',str(deck.resolve())],cwd=deck.parent,capture_output=True,timeout=timeout_s)
    log=deck.with_suffix('.log')
    if log.exists():
        b=log.read_bytes();text=b.decode('utf-16le' if b[1:2]==b'\x00' else 'utf-8',errors='replace')
        deck.with_suffix('.run.txt').write_text('\n'.join(line.rstrip() for line in text.splitlines()).rstrip()+'\n')
    if result.returncode: raise RuntimeError(f'{deck}: simulator exit {result.returncode}; inspect .log')
    return read_raw(deck.with_suffix('.raw'))
