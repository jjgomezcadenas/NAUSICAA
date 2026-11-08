"""Publication receipts preserve original provenance through response deduplication."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
import h5py
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ion_signal.canonical_response import resolve_response, sha256

class CanonicalResponseTest(unittest.TestCase):
    def test_retired_parent_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            old=Path(tmp).resolve()/'old.h5';new=Path(tmp).resolve()/'published.h5'
            old.write_bytes(b'original corrected file');oldhash=sha256(old)
            self.assertEqual(resolve_response(old,oldhash),old)
            with self.assertRaises(ValueError):resolve_response(old,'wrong')
            with h5py.File(new,'w') as h:
                h.attrs['schema']=np.bytes_('nausicaa.ion_csa_response.v1')
                h['metadata/complete']=True
                h['metadata/parent_response_path']=np.bytes_(str(old))
                h['metadata/parent_response_sha256']=np.bytes_(oldhash)
                h['metadata/numerical_payload']=np.bytes_('unchanged_from_corrected_parent')
            receipt=dict(schema='itaca.canonical_response.v1',original_filename=old.name,
                original_sha256=oldhash,published_sha256=sha256(new),published_relative_path=new.name)
            def save():old.with_suffix('.canonical.toml').write_text('\n'.join(f'{k} = {json.dumps(v)}' for k,v in receipt.items()))
            save();old.unlink()
            self.assertEqual(resolve_response(old,oldhash),new)
            with self.assertRaises(ValueError):resolve_response(old,'wrong')
            receipt['published_sha256']='wrong';save()
            with self.assertRaises(ValueError):resolve_response(old,oldhash)
            receipt['published_sha256']=sha256(new);save()
            with h5py.File(new,'r+') as h:h['metadata/complete'][...]=False
            receipt['published_sha256']=sha256(new);save()
            with self.assertRaises(ValueError):resolve_response(old,oldhash)

if __name__=='__main__':unittest.main()
