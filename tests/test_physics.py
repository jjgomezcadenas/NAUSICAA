"""Independent checks of identities, dimensions, limits and numerical stability."""

from pathlib import Path
import unittest
import numpy as np
from scipy.integrate import quad, dblquad
from ion_signal.physics import (Detector, E_CHARGE, K_B, square_phi, square_weight,
                               square_z_at_phi, disk_phi, disk_weight,
                               single_ion_current, thermal_spread)

ROOT = Path(__file__).resolve().parents[1]


class PhysicsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = Detector.from_toml(ROOT/'config/detector.toml')
        cls.a = cls.d.pad_width_m/2

    def test_geometry_and_voltage_gradients(self):
        d = self.d
        self.assertAlmostEqual(d.hole_field_V_m*d.ifg_thickness_m, 5000)
        self.assertAlmostEqual(d.gap_field_V_m*d.gap_m, 3000)
        self.assertGreater(d.open_fraction*d.hole_field_V_m, d.bulk_field_V_m)

    def test_transport_continuity_of_mobility_assumption(self):
        d = self.d
        self.assertAlmostEqual(d.mobility_m2_V_s*d.bulk_field_V_m, d.bulk_speed_m_s)
        for _, length, field in d.regions():
            t = length/(d.mobility_m2_V_s*field)
            self.assertAlmostEqual(np.sqrt(2*d.thermal_D_m2_s*t),
                                   thermal_spread(length, field, d.temperature_K))

    def test_thermal_first_passage_time_variance(self):
        d = self.d
        L, v = d.bulk_reference_length_m, d.bulk_speed_m_s
        sigma_t = thermal_spread(L, d.bulk_field_V_m, d.temperature_K)/v
        self.assertAlmostEqual(sigma_t**2, 2*d.thermal_D_m2_s*L/v**3)

    def test_boundary_values(self):
        for phi in (square_phi, disk_phi):
            self.assertEqual(float(phi(0, self.a)), 1.)
            self.assertLess(float(phi(1e5*self.a, self.a)), 1e-9)

    def test_square_far_field_solid_angle(self):
        z = 1000*self.a
        self.assertAlmostEqual(square_phi(z, self.a)/(2*self.a**2/(np.pi*z*z)),
                               1., delta=2e-6)

    def test_monotonic_and_positive(self):
        z = np.geomspace(self.a*1e-5, self.a*1e4, 300)
        self.assertTrue(np.all(np.diff(square_phi(z, self.a)) < 0))
        self.assertTrue(np.all(square_weight(z, self.a) > 0))

    def test_field_is_negative_gradient(self):
        z = np.geomspace(self.a*.02, self.a*20, 70)
        h = self.a*1e-5
        numerical = (square_phi(z-h, self.a)-square_phi(z+h, self.a))/(2*h)
        np.testing.assert_allclose(numerical, square_weight(z, self.a), rtol=1e-7)

    def test_square_formula_against_surface_integral(self):
        # Independently integrate dimensionless Poisson kernel over the square.
        for u in (.2, 1., 4.):
            val, _ = dblquad(lambda y,x: u/(2*np.pi*(x*x+y*y+u*u)**1.5),
                             -1, 1, lambda x:-1, lambda x:1, epsabs=1e-10)
            self.assertAlmostEqual(val, square_phi(u*self.a, self.a), places=9)

    def test_complete_charge_integral(self):
        for weight in (square_weight, disk_weight):
            val, _ = quad(lambda u: self.a*weight(u*self.a, self.a), 0, np.inf,
                          epsabs=1e-12, epsrel=1e-12)
            self.assertAlmostEqual(val, 1, places=11)

    def test_finite_path_charge_without_renormalization(self):
        z0 = self.d.gap_m
        val, _ = quad(lambda u:self.a*square_weight(u*self.a, self.a),
                      0, z0/self.a, epsabs=1e-12)
        self.assertAlmostEqual(val, 1-square_phi(z0, self.a), places=11)
        self.assertLess(val, 1.)

    def test_fractional_charge_positions(self):
        p = np.array([.001, .01, .1, .5, .9, 1.])
        np.testing.assert_allclose(square_phi(square_z_at_phi(p, self.a), self.a),
                                   p, atol=1e-14)

    def test_collection_termination_and_velocity_scaling(self):
        v = self.d.mobility_m2_V_s*self.d.gap_field_V_m
        self.assertEqual(single_ion_current(0, v, self.a), 0)
        self.assertEqual(single_ion_current(.001, v, self.a), 0)
        t = np.linspace(-.001, -.000001, 100)
        np.testing.assert_allclose(single_ion_current(t/2, 2*v, self.a),
                                   2*single_ion_current(t, v, self.a), rtol=1e-13)

    def test_area_with_independent_uniform_sampling_converges(self):
        z0 = self.d.gap_m
        exact = 1-square_phi(z0, self.a)
        errors = []
        for n in (100, 200, 400):
            z = np.linspace(0,z0,n+1)
            val = np.trapezoid(square_weight(z,self.a),z)
            errors.append(abs(val-exact))
        self.assertLess(errors[2], errors[1])
        self.assertLess(errors[1], errors[0])
        self.assertLess(errors[2], 1e-7)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError): square_phi(-1, self.a)
        with self.assertRaises(ValueError): square_z_at_phi(0, self.a)
        with self.assertRaises(ValueError): thermal_spread(-1, 20000, 293.15)
        with self.assertRaises(ValueError): single_ion_current(-1, 0, self.a)


if __name__ == '__main__':
    unittest.main()
