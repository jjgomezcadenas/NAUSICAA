"""Independent physics, circuit, noise and SPICE-export checks."""
from pathlib import Path
from dataclasses import replace
import unittest
import numpy as np
from scipy.integrate import solve_ivp, simpson
from scipy.linalg import solve_continuous_lyapunov
from ion_signal.physics import Detector, E_CHARGE, square_phi
from ion_signal.readout import CurrentWaveform, csa_response
from ion_signal.realistic_csa import (CSA,Noise,response,ions_waveform,
    recovered_charge,frame_charge,sampled_voltage,noise_covariance,noise_metrics,
    estimator_coefficients,pwl_current)

ROOT=Path(__file__).resolve().parents[1]

class RealisticCSATests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=Detector.from_toml(ROOT/'config/detector.toml')
        cls.m=CSA(3e-15,250e-15,1e14,1000.,1000.)

    def test_dc_gain_and_voltage_noise(self):
        hi,he=self.m.transfer(0.)
        a=self.m.open_loop_gain
        self.assertAlmostEqual(hi.real/(-a*self.m.feedback_resistance_ohm/(1+a)),1.)
        self.assertAlmostEqual(he.real,a/(1+a))
        self.assertTrue(np.all(self.m.poles_per_s().real<0))

    def test_state_and_transfer_independently_agree(self):
        m=self.m
        b=np.array([1/(m.input_capacitance_F+m.feedback_capacitance_F),0.])
        for f in (.01,1.,1e3,1e6):
            hi,_=m.transfer(f)
            state=np.linalg.solve(2j*np.pi*f*np.eye(2)-m.state_matrix(),b)
            np.testing.assert_allclose(state[1],hi,rtol=1e-11)

    def test_underdamped_poles_and_frequency_response(self):
        m=replace(self.m,feedback_resistance_ohm=1e8)
        poles=m.poles_per_s()
        self.assertTrue(np.any(abs(poles.imag)>0))
        self.assertTrue(np.all(poles.real<0))
        # A sinusoidal current should give the transfer-function amplitude/phase.
        f=1e5; current=1e-12; duration=30/f; dt=1/(f*400)
        edges=np.linspace(0,duration,12001)
        charges=current*(np.cos(2*np.pi*f*edges[:-1])-np.cos(2*np.pi*f*edges[1:]))/(2*np.pi*f)
        wave=CurrentWaveform(edges,charges)
        voltage=response(wave,m)[:,1]
        hi,_=m.transfer(f)
        expected=current*(hi.real*np.sin(2*np.pi*f*edges)+hi.imag*np.cos(2*np.pi*f*edges))
        self.assertLess(np.max(abs(voltage[8000:]-expected[8000:]))/(current*abs(hi)),3e-5)

    def test_independent_circuit_ode(self):
        m=self.m; current=1e-14; dt=2e-5
        wave=CurrentWaveform(np.arange(101)*dt,np.full(100,current*dt))
        actual=response(wave,m)
        cf=m.feedback_capacitance_F; ci=m.input_capacitance_F
        w=2*np.pi*m.amplifier_pole_Hz; g=1/m.feedback_resistance_ohm
        def kirchhoff(t,state):
            vi,vo=state
            dvo=-w*(vo+m.open_loop_gain*vi)
            return [(current+cf*dvo-g*(vi-vo))/(ci+cf),dvo]
        expected=solve_ivp(kirchhoff,[0,wave.edges_s[-1]],[0,0],method='Radau',
                           t_eval=wave.edges_s,rtol=1e-10,atol=1e-14).y.T
        np.testing.assert_allclose(actual,expected,rtol=2e-7,atol=1e-12)

    def test_high_gain_fast_pole_ideal_limit(self):
        m=replace(self.m,open_loop_gain=1e6,amplifier_pole_Hz=1e4)
        wave=ions_waveform(self.d,[0.],-.002,.005,2e-6)
        voltage=response(wave,m)[:,1]
        ideal=csa_response(wave,m.feedback_capacitance_F,m.nominal_discharge_time_s)
        np.testing.assert_allclose(voltage,ideal,rtol=2e-4,atol=1e-9)

    def test_finite_gain_charge_sharing_limit(self):
        m=replace(self.m,feedback_resistance_ohm=1e30,amplifier_pole_Hz=1e4)
        wave=ions_waveform(self.d,[0.],-.002,.001,2e-6)
        actual=response(wave,m)[-1,1]
        expected=-wave.total_charge_C/(m.feedback_capacitance_F+
            (m.input_capacitance_F+m.feedback_capacitance_F)/m.open_loop_gain)
        self.assertAlmostEqual(actual/expected,1.,places=9)

    def test_exact_k_charge_and_truncated_path(self):
        arrivals=np.linspace(.001,.011,17)
        wave=ions_waveform(self.d,arrivals,-.002,.02,7e-6)
        qg=E_CHARGE*(1-square_phi(self.d.gap_m,self.d.pad_width_m/2))
        self.assertAlmostEqual(wave.total_charge_C/(17*qg),1.,places=12)
        cut=ions_waveform(self.d,arrivals,.005,.006,7e-6)
        self.assertLess(cut.total_charge_C,wave.total_charge_C)

    def test_superposition_and_retained_state(self):
        wave=ions_waveform(self.d,[0.,.01],-.002,.03,5e-6)
        v=response(wave,self.m)
        separate=sum(response(ions_waveform(self.d,[t],-.002,.03,5e-6),self.m)
                     for t in (0.,.01))
        np.testing.assert_allclose(v,separate,rtol=1e-11,atol=1e-16)
        cut=1000
        first=CurrentWaveform(wave.edges_s[:cut+1],wave.charge_per_bin_C[:cut])
        second=CurrentWaveform(wave.edges_s[cut:],wave.charge_per_bin_C[cut:])
        resumed=response(second,self.m,response(first,self.m)[-1])
        np.testing.assert_allclose(v[cut:],resumed,rtol=1e-11,atol=1e-16)

    def test_physical_current_timestep_convergence(self):
        vals=[]
        for dt in (20e-6,10e-6,5e-6,1e-6):
            wave=ions_waveform(self.d,[0.],-.002,0.,dt)
            vals.append(response(wave,self.m)[-1,1])
        errors=np.abs(np.array(vals[:-1])-vals[-1])
        self.assertLess(errors[1],.4*errors[0])
        self.assertLess(errors[2],.4*errors[1])

    def test_sampler_frame_charge_and_estimator(self):
        dt=1e-5; total=.03; current=1e-14
        wave=CurrentWaveform(np.linspace(0,total,3001),np.full(3000,current*dt))
        v=csa_response(wave,3e-15,self.m.nominal_discharge_time_s)
        samples=np.array([0,.01,.02,.03])
        q=recovered_charge(sampled_voltage(wave.edges_s,v,samples),self.m,.01)
        truth=frame_charge(wave,samples)
        np.testing.assert_allclose(q,truth,rtol=1e-10)
        # Exact aperture average of a linear voltage ramp.
        averaged=sampled_voltage(wave.edges_s,2*wave.edges_s,samples[1:],.001)
        np.testing.assert_allclose(averaged,2*(samples[1:]-.0005),rtol=1e-12)

    def test_noise_units_scaling_zero_and_covariance(self):
        n=Noise(5e-9,2.5e-15,1e-19,.1,1e4)
        scaled=replace(n,voltage_white_V_per_sqrtHz=1e-8,
                       voltage_flicker_V2=1e-14,current_white_A_per_sqrtHz=2e-19)
        a=noise_metrics(self.m,n,.01); b=noise_metrics(self.m,scaled,.01)
        self.assertAlmostEqual(b['recovered_charge_rms_C']/a['recovered_charge_rms_C'],2.,places=10)
        self.assertLessEqual(abs(a['adjacent_sample_correlation']),1.)
        z=noise_metrics(self.m,Noise(0,0,0,.1,1e4),.01)
        self.assertEqual(z['recovered_charge_rms_C'],0.)
        hi,he=self.m.transfer(10.)
        expected=abs(hi)**2*1e-38+abs(he)**2*(25e-18+2.5e-15/10)
        self.assertAlmostEqual(n.output_psd(self.m,10.)/expected,1.)

    def test_noise_direct_filter_integral_and_grid_convergence(self):
        n=Noise(5e-9,2.5e-15,1e-19,.1,1e4)
        delta=.01; aperture=1e-5
        result=noise_metrics(self.m,n,delta,aperture)
        alpha,k=estimator_coefficients(self.m,delta)
        direct=[]
        for count in (100001,200001):
            f=np.geomspace(n.low_frequency_Hz,n.high_frequency_Hz,count)
            power=n.output_psd(self.m,f,aperture)*abs(1-alpha*np.exp(-2j*np.pi*f*delta))**2*k*k
            direct.append(simpson(power,x=f))
        self.assertLess(abs(direct[1]/direct[0]-1),1e-7)
        self.assertLess(abs(result['recovered_charge_rms_C']**2/direct[1]-1),2e-6)
        tighter=noise_metrics(self.m,n,delta,aperture,relative_tolerance=1e-9)
        self.assertLess(abs(tighter['recovered_charge_rms_C']/result['recovered_charge_rms_C']-1),1e-6)

    def test_noise_current_white_lyapunov_limit(self):
        m=self.m; n=Noise(0,0,1e-19,1e-5,1e7)
        b=np.array([1/(m.input_capacitance_F+m.feedback_capacitance_F),0.])
        covariance=solve_continuous_lyapunov(m.state_matrix(),-np.outer(b,b)*1e-38/2)
        integrated=noise_covariance(m,n)
        self.assertLess(abs(integrated/covariance[1,1]-1),1e-4)

    def test_pwl_export_charge(self):
        wave=ions_waveform(self.d,[0,.001,.002],-.002,.02,5e-6)
        t,current=pwl_current(wave)
        self.assertEqual(t[0],0.)
        self.assertTrue(np.all(np.diff(t)>0))
        self.assertAlmostEqual(np.trapezoid(current,t)/wave.total_charge_C,1.,places=12)

    def test_invalid_parameters(self):
        with self.assertRaises(ValueError):
            replace(self.m,feedback_resistance_ohm=0)
        with self.assertRaises(ValueError):
            Noise(1e-9,1e-15,0,0,1e4)
        with self.assertRaises(ValueError):
            ions_waveform(self.d,[],0,.01,1e-6)



class UniformMeanResponseTests(unittest.TestCase):
    """Uniform arrival density is independent of an arbitrary discrete count."""
    def test_uniform_charge_conservation_and_scaling(self):
        from ion_signal.realistic_csa import uniform_ion_waveform
        d=Detector.from_toml(ROOT/'config/detector.toml')
        one=uniform_ion_waveform(d,0.,.01,1.,-.002,.012,7e-6)
        many=uniform_ion_waveform(d,0.,.01,68.,-.002,.012,7e-6)
        expected=E_CHARGE*(1-square_phi(d.gap_m,d.pad_width_m/2))
        self.assertAlmostEqual(one.total_charge_C/expected,1.,places=10)
        np.testing.assert_allclose(many.charge_per_bin_C,68*one.charge_per_bin_C,rtol=1e-12,atol=1e-34)

    def test_uniform_bins_against_independent_single_ion_average(self):
        from scipy.integrate import quad
        from ion_signal.realistic_csa import uniform_ion_waveform
        d=Detector.from_toml(ROOT/'config/detector.toml')
        duration=.01
        w=uniform_ion_waveform(d,0.,duration,1.,-.002,.012,2e-5)
        v=d.mobility_m2_V_s*d.gap_field_V_m
        for t in (-.001,0.,.004,.0099):
            index=np.searchsorted(w.edges_s,t)
            lo,hi=w.edges_s[index:index+2]
            # Integrate each ion's exact bin charge over its collection time.
            def integrand(u):
                return float(square_phi(np.clip(v*(u-hi),0,d.gap_m),d.pad_width_m/2)
                             -square_phi(np.clip(v*(u-lo),0,d.gap_m),d.pad_width_m/2))
            points=sorted(x for x in (lo,hi,lo+d.gap_m/v,hi+d.gap_m/v) if 0<x<duration)
            expected=E_CHARGE/duration*quad(integrand,0.,duration,points=points,epsabs=1e-15,epsrel=1e-10)[0]
            self.assertAlmostEqual(w.charge_per_bin_C[index]/expected,1.,places=8)


if __name__=='__main__':
    unittest.main()
