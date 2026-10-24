"""Independent analytic and covariance checks of direct noise integration."""
from dataclasses import dataclass
import unittest
import numpy as np
from ion_signal.realistic_csa import CSA, Noise, estimator_coefficients, noise_metrics
from ion_signal.noise_study import direct_charge_variance


@dataclass
class FlatOutputSpectrum:
    low_frequency_Hz: float = .1
    high_frequency_Hz: float = 12345.
    density_V2_per_Hz: float = 1e-12

    def output_psd(self, model, frequency_Hz, aperture_s):
        return np.full_like(frequency_Hz, self.density_V2_per_Hz)


class NoiseStudyTests(unittest.TestCase):
    def setUp(self):
        self.model = CSA(3e-15, 250e-15, 1e14, 1000., 1000.)

    def test_flat_output_spectrum_exact_integral(self):
        noise = FlatOutputSpectrum()
        dt = .01
        alpha, k = estimator_coefficients(self.model, dt)
        w = 2*np.pi*dt
        lo, hi = noise.low_frequency_Hz, noise.high_frequency_Hz
        expected = k*k*noise.density_V2_per_Hz*((1+alpha*alpha)*(hi-lo)
                    -2*alpha*(np.sin(w*hi)-np.sin(w*lo))/w)
        actual = direct_charge_variance(self.model, noise, dt)
        self.assertAlmostEqual(actual/expected, 1., places=11)

    def test_actual_noise_matches_covariance_for_both_apertures(self):
        noise = Noise(5e-9, 2.5e-15, 0., .1, 1e6, 293.15, True)
        for aperture in (0., 1e-5):
            with self.subTest(aperture_s=aperture):
                metrics = noise_metrics(self.model, noise, .01, aperture, 1e-9)
                direct = direct_charge_variance(self.model, noise, .01, aperture, order=16)
                self.assertAlmostEqual(direct/metrics['recovered_charge_rms_C']**2, 1., places=7)


if __name__ == '__main__':
    unittest.main()
