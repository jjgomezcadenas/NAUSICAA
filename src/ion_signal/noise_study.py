"""Independent positive spectral integral for the two-reading charge estimator."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from .realistic_csa import estimator_coefficients


def direct_charge_variance(model, noise, interval_s, aperture_s=0., order=8):
    """Return charge variance in C² from K² ∫ S_V(f)|1-alpha exp(-i2πfΔt)|² df.

    One-sided PSD in Hz; stationary, independent sources. This quadrature does
    not obtain variance by subtracting covariance integrals. Logarithmic panels
    resolve the low-frequency circuit response; panels at most half an
    estimator oscillation wide resolve high frequencies. Compare orders 8 and
    16 to check convergence over the explicitly specified noise band.
    """
    alpha, k = estimator_coefficients(model, interval_s)
    if not isinstance(order, int) or order < 2:
        raise ValueError('Quadrature order must be an integer >= 2')
    lo, hi = noise.low_frequency_Hz, noise.high_frequency_Hz
    edges = np.unique(np.r_[np.geomspace(lo, hi, 513),
                           np.arange(lo, hi, .5/interval_s), hi])
    nodes, weights = leggauss(order)
    half = np.diff(edges)/2
    f = (edges[1:]+edges[:-1])[:, None]/2 + half[:, None]*nodes
    # Stable equivalent of 1+alpha²-2alpha*cos(2πfΔt).
    estimator_power = (1-alpha)**2 + 4*alpha*np.sin(np.pi*f*interval_s)**2
    integrand = noise.output_psd(model, f, aperture_s)*estimator_power
    return float(k*k*np.sum(half*(integrand @ weights)))
