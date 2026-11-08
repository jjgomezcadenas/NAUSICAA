"""Resolve a retired corrected response to its verified published numerical copy."""
from pathlib import Path
import hashlib
import tomllib
import h5py


def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def resolve_response(path, expected_sha256=None):
    """Return the retained HDF5 path, validating a retirement receipt when needed.

    A receipt preserves the original corrected-file hash separately from the
    published-file hash. The published metadata must identify that exact parent
    and declare an unchanged numerical payload. Existing files are checked against
    ``expected_sha256`` when supplied. Missing or inconsistent receipts raise
    ValueError; the original summary and its provenance remain unchanged.
    """
    path = Path(path).resolve()
    if path.is_file():
        if expected_sha256 is not None and sha256(path) != expected_sha256:
            raise ValueError(f'Response checksum mismatch: {path}')
        return path
    receipt = tomllib.loads(path.with_suffix('.canonical.toml').read_text())
    target = (path.parent / receipt['published_relative_path']).resolve()
    if (receipt['schema'] != 'itaca.canonical_response.v1'
            or receipt['original_filename'] != path.name
            or (expected_sha256 is not None and receipt['original_sha256'] != expected_sha256)
            or sha256(target) != receipt['published_sha256']):
        raise ValueError(f'Canonical response receipt mismatch: {path}')
    with h5py.File(target) as h:
        text = lambda k: h['metadata/' + k][()].decode()
        if (h.attrs['schema'].decode() != 'nausicaa.ion_csa_response.v1'
                or not h['metadata/complete'][()]
                or text('parent_response_sha256') != receipt['original_sha256']
                or Path(text('parent_response_path')).name != path.name
                or text('numerical_payload') != 'unchanged_from_corrected_parent'):
            raise ValueError(f'Published response parent mismatch: {path}')
    return target
