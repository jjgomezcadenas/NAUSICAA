"""Check the HDF5 boundary that prevents redundant SPICE acquisitions."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
import json
import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('csa_launcher',ROOT/'scripts/production/run_ion_csa_pilot.py')
launcher=importlib.util.module_from_spec(spec)
spec.loader.exec_module(launcher)

class NonlinearDispatch(unittest.TestCase):
    def test_absent_control_measurement_is_json_null(self):
        spec=importlib.util.spec_from_file_location('csa_presenter',ROOT/'scripts/production/present_ion_csa_pilot.py')
        presenter=importlib.util.module_from_spec(spec);spec.loader.exec_module(presenter)
        payload={'control_error':presenter.plain(np.float64('nan')),
                 'control_errors':presenter.plain(np.array([],dtype=float))}
        self.assertEqual(json.loads(json.dumps(payload,allow_nan=False)),
                         {'control_error':None,'control_errors':[]})

    def test_julia_bytes_python_strings_and_spice(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'work') as directory:
            job={'folder':directory}
            self.assertFalse(launcher.nonlinear_case(job))
            for value,expected in ((np.bytes_('calibrated_nonlinear'),True),
                                   ('calibrated_nonlinear',True),('spice',False)):
                with h5py.File(Path(directory)/'response.h5','w') as h:
                    h.attrs['response_model']=value
                self.assertEqual(launcher.nonlinear_case(job),expected)
            with h5py.File(Path(directory)/'response.h5','w'):
                pass
            self.assertFalse(launcher.nonlinear_case(job))

if __name__=='__main__':unittest.main()
