"""Independent charge, response and convergence checks for the analytical CSA."""
from pathlib import Path
import unittest
import numpy as np
from scipy.integrate import quad
from ion_signal.physics import Detector, E_CHARGE, square_phi, square_weight
from ion_signal.readout import (CurrentWaveform, single_ion_bins, make_waveform,
                               csa_response, arrival_counts, waveform_metrics)

ROOT = Path(__file__).resolve().parents[1]


class ReadoutTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = Detector.from_toml(ROOT/'config/detector.toml')

    def test_exact_finite_gap_charge(self):
        expected = E_CHARGE*(1-square_phi(self.d.gap_m,self.d.pad_width_m/2))
        for step in (50e-6,5e-6,1e-6):
            wave = single_ion_bins(self.d,step)
            self.assertAlmostEqual(wave.total_charge_C/expected,1.,places=13)
            self.assertLess(wave.total_charge_C,E_CHARGE)

    def test_arrival_tail_is_not_renormalized(self):
        dt=.01
        times=np.arange(-100,101)*dt
        counts=arrival_counts(times,dt,260,'gaussian',sigma_s=1.)
        expected=quad(lambda t: np.exp(-t*t/2)/np.sqrt(2*np.pi),-1.005,1.005)[0]
        self.assertAlmostEqual(counts.sum()/260,expected,places=12)
        self.assertLess(counts.sum()/260,.7)

    def test_packet_charge_conservation(self):
        for kind in ('single','uniform','gaussian','double_gaussian'):
            wave,meta=make_waveform(self.d,kind,ions=260)
            expected=meta['represented_arrivals']*meta['kernel_charge_C']
            self.assertAlmostEqual(wave.total_charge_C/expected,1.,places=12)

    def test_ideal_integration(self):
        wave,_=make_waveform(self.d,'gaussian')
        volts=csa_response(wave,3e-15)
        self.assertAlmostEqual(-volts[-1]*3e-15/wave.total_charge_C,1.,places=11)

    def test_uniform_input_closed_form_and_decay(self):
        dt=1e-4; duration=.01; current=1e-14; cf=3e-15
        q=np.r_[np.full(100,current*dt),np.zeros(200)]
        wave=CurrentWaveform(np.arange(301)*dt,q)
        for tau in (3e-5,.001,.1,10.):
            volts=csa_response(wave,cf,tau)
            amplitude=-current*tau/cf*(-np.expm1(-duration/tau))
            self.assertAlmostEqual(volts[100]/amplitude,1.,places=11)
            expected=amplitude*np.exp(-.02/tau)
            self.assertAlmostEqual(volts[-1],expected,delta=1e-12)

    def test_independent_continuous_kernel_quadrature_convergence(self):
        cf=3e-15; tau=1e-4
        d=self.d; v=d.mobility_m2_V_s*d.gap_field_V_m; a=d.pad_width_m/2
        # Integrate in dimensionless distance, independently of the time bins.
        exact=-E_CHARGE/cf*quad(lambda u:a*square_weight(a*u,a)*np.exp(-a*u/v/tau),
                               0,d.gap_m/a,epsabs=1e-12)[0]
        errors=[]
        for step in (20e-6,10e-6,5e-6):
            computed=csa_response(single_ion_bins(d,step),cf,tau)[-1]
            errors.append(abs(computed/exact-1))
        self.assertLess(errors[1],errors[0]/3)
        self.assertLess(errors[2],errors[1]/3)
        self.assertLess(errors[2],.001)

    def test_linear_charge_and_inverse_capacitance_scaling(self):
        wave,_=make_waveform(self.d,'uniform',ions=1)
        larger=CurrentWaveform(wave.edges_s,260*wave.charge_per_bin_C)
        np.testing.assert_allclose(csa_response(larger,6e-15,.1),
                                   130*csa_response(wave,3e-15,.1),rtol=1e-12,atol=1e-18)

    def test_sampling_phase_is_explicit(self):
        wave,_=make_waveform(self.d,'single',ions=260)
        voltage=csa_response(wave,3e-15,3e-5)
        aligned=waveform_metrics(wave,voltage,3e-15,sample_phase_s=0.)
        delayed=waveform_metrics(wave,voltage,3e-15,sample_phase_s=.005)
        self.assertGreater(aligned['sampled_peak_magnitude_V'],
                           100*delayed['sampled_peak_magnitude_V'])


if __name__=='__main__':
    unittest.main()
