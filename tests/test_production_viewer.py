"""Projection conservation and signed response checks for the production viewer."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from cloud_viewer.production import ProductionSource
class Scalar:
    def __init__(self,x):self.x=x
    def __getitem__(self,_):return self.x
class ProjectionTest(unittest.TestCase):
    def test_signed_coordinates_and_charge(self):
        g={k:Scalar(v) for k,v in dict(pitch_mm=3,reference_mm=np.array([10,20,30]),phase_mm=np.array([1,2]),event_id=7,ion_drift_cm=50,electron_drift_cm=100,pair_checksum=b'abc').items()}
        panel=ProductionSource.panel(g,[-1,0],[0,1],[10,11],np.array([-5.,12.]),30,'ion','ions',1,'Test')
        self.assertEqual(panel['edges']['x'],[8.,11.,14.])
        self.assertEqual(panel['edges']['z'],[0.,30.,60.])
        self.assertTrue(panel['signed'])
        for image in panel['projections'].values():self.assertEqual(np.sum(image),7)
        self.assertEqual(panel['cloud3d']['ions'],[-5.,12.])
if __name__=='__main__':unittest.main()
