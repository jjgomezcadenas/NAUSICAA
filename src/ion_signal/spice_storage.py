"""Preserve circuit text and store acquired waveform groups in HDF5."""
import hashlib
import numpy as np

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def preserve(p, text):
    if p.exists():
        assert p.read_text() == text, f'Refusing overwrite {p}'
    else: p.write_text(text)


def save(h, name, data, **attrs):
    g=h.create_group(name);g.attrs.update(attrs)
    for k,v in data.items():g.create_dataset(k,data=np.asarray(v).real,compression='gzip',compression_opts=9,shuffle=True)
    h.file.flush()
