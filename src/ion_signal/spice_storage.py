"""Preserve circuit text and store acquired waveform groups in HDF5."""
import hashlib
from pathlib import Path
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


def source_manifest_matches(saved, current, root):
    """Compare input hashes through recorded-path relocation maps, retaining saved records.

    Source paths may be repository-relative or absolute. Two manifests describe
    the same inputs when their resolved paths and hashes agree. A changed hash
    or an additional/missing source is a mismatch.
    """
    from .project_paths import resolve_recorded_path
    def canonical(manifest):
        return {(str(resolve_recorded_path(path, root=root).resolve()), digest)
                for path, digest in manifest.items()}
    return canonical(saved) == canonical(current)
