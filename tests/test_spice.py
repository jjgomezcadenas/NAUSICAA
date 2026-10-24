import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from ion_signal.spice import schematic_netlist, transform_pin, read_raw, sampled_charge, charge_noise

class SpiceTests(unittest.TestCase):
    def test_feedback_mirrored_pin_order(self):
        self.assertEqual(transform_pin((48,0),(144,736),'M270'),(144,688))
        self.assertEqual(transform_pin((0,80),(144,736),'M270'),(64,736))
        self.assertEqual(transform_pin((48,96),(144,736),'M270'),(48,688))

    def test_supplied_connectivity(self):
        root=Path(__file__).resolve().parents[1]
        el,d=schematic_netlist(root/'LTSpice/CSAfolded180nm.asc')
        self.assertIn('Mf out Vdis in 0 180nmN l={LMf} w={WMf} m=1',el)
        self.assertIn('M1 s in med Vdda 180nmP l={Lcs} w={Wcs} m={Mp}',el)
        self.assertIn('Cf out in 3f',el)
        self.assertEqual(len(el),19)

    def test_raw_double_and_compression_sign(self):
        with TemporaryDirectory() as tmp:
            p=Path(tmp)/'x.raw'
            h='Title: Test\nFlags: real double\nNo. Variables: 2\nNo. Points: 3\nVariables:\n0 time time\n1 V(out) voltage\nBinary:\n'
            p.write_bytes(h.encode('utf-16le')+np.array([[0,1],[-.1,2],[.2,3]],dtype='<f8').tobytes())
            d,_=read_raw(p)
            np.testing.assert_array_equal(d['time'],[0,.1,.2])
            np.testing.assert_array_equal(d['V(out)'],[1,2,3])

    def test_constant_current_recovery_with_initial_state(self):
        cf=3e-15;tau=.3;dt=.01;i=1e-15
        t=np.linspace(0,dt,1001);v=.02*np.exp(-t/tau)-i*tau/cf*(1-np.exp(-t/tau))
        q=sampled_charge(t,v,[0,dt],cf,tau)[0]
        self.assertAlmostEqual(q/(i*dt),1,places=12)

    def test_noise_against_analytic_lorentzian_covariance(self):
        # S(f)=4*tau/(1+(2*pi*f*tau)^2): variance=1 and covariance=exp(-dt/tau).
        tau=.3;dt=.01;cf=3e-15;f=np.geomspace(1e-7,1e6,20000)
        asd=np.sqrt(4*tau/(1+(2*np.pi*f*tau)**2));a=np.exp(-dt/tau)
        k=cf*(dt/tau)/(1-a);exact=k*np.sqrt(1-a*a)
        got=charge_noise(f,asd,dt,cf,tau)
        self.assertLess(abs(got/exact-1),2e-5)
        self.assertLess(abs(charge_noise(f,2*asd,dt,cf,tau)/got-2),1e-12)

if __name__=='__main__':unittest.main()
