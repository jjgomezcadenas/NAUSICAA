import unittest
import numpy as np
from ion_signal.shaping import box_average,trapezoid,trapezoid_transfer,filtered_rms,averaged_frame_transfer

class ShapingTests(unittest.TestCase):
    def test_average_exact_linear_ramp(self):
        t=np.array([0.,.1,.4,.6,1.]);v=2+3*t
        end=np.array([.33,.55,.94]);T=.21
        np.testing.assert_allclose(box_average(t,v,end,T),2+3*(end-T/2),rtol=1e-14)

    def test_step_produces_trapezoid(self):
        t=np.linspace(-4,5,90001);v=(t>=0).astype(float);x=np.linspace(-.5,3.5,400)
        expected=np.clip(x,0,1)-np.clip(x-2,0,1)
        np.testing.assert_allclose(trapezoid(t,v,x,1.,1.),expected,atol=1.1e-4)

    def test_reject_constant_baseline(self):
        t=np.linspace(-3,5,8001);v=np.full_like(t,1.17);x=np.linspace(0,4,400)
        np.testing.assert_allclose(trapezoid(t,v,x,1.,1.),0,atol=1e-12)

    def test_frequency_response_matches_time_average(self):
        T=.00018;G=.00018;f=2345.
        t=np.linspace(0,.003,30001);v=np.cos(2*np.pi*f*t)
        x=np.linspace(.001,.002,100)
        expected=np.real(trapezoid_transfer(f,T,G)*np.exp(2j*np.pi*f*x))
        np.testing.assert_allclose(trapezoid(t,v,x,T,G),expected,atol=1e-6)

    def test_noise_matches_analytic_time_domain_covariance(self):
        # PSD below has covariance exp(-abs(lag)/tau). Integrating two
        # rectangular averages in time gives this independent exact variance.
        T=.00018;G=.00018;tau=.00005;L=T+G
        f=np.geomspace(.001,1e6,20000);asd=np.sqrt(4*tau/(1+(2*np.pi*f*tau)**2))
        var=2*tau/T-2*(tau/T)**2*(-np.expm1(-T/tau))
        cov=(tau/T)**2*np.exp(-(L-T)/tau)*(-np.expm1(-T/tau))**2
        exact=np.sqrt(2*var-2*cov)
        measured=filtered_rms(f,asd,lambda f:trapezoid_transfer(f,T,G),2*T+G)
        self.assertLess(abs(measured/exact-1),2e-5)

    def test_averaged_frame_noise_against_time_domain_covariance(self):
        T=.0002;dt=.01;tau=.003;alpha=np.exp(-dt/.3)
        f=np.geomspace(.00001,1e6,20000)
        asd=np.sqrt(4*tau/(1+(2*np.pi*f*tau)**2))
        var=2*tau/T-2*(tau/T)**2*(-np.expm1(-T/tau))
        cov=(tau/T)**2*np.exp(-(dt-T)/tau)*(-np.expm1(-T/tau))**2
        exact=np.sqrt((1+alpha**2)*var-2*alpha*cov)
        got=filtered_rms(f,asd,lambda f:averaged_frame_transfer(f,T,dt,alpha),dt+T)
        self.assertLess(abs(got/exact-1),2e-5)

if __name__=='__main__':unittest.main()
