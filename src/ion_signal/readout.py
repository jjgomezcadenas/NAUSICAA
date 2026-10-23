"""Charge-conserving pad-current benchmarks and a linear CSA, in SI units.

The pad is a central square in an infinite grounded plane. Each ion travels
only the finite collection gap at constant speed. Arrival distributions are
specified for one pad, not inferred from a full event or focusing simulation.
The CSA assumes an ideal amplifier with constant Cf and linear Rf; it is not
a transistor model. Positive input current produces negative output voltage.
"""

from dataclasses import dataclass
import numpy as np
from scipy.signal import fftconvolve, lfilter
from scipy.special import ndtr
from .physics import E_CHARGE, square_phi, thermal_spread


@dataclass(frozen=True)
class CurrentWaveform:
    edges_s: np.ndarray
    charge_per_bin_C: np.ndarray

    @property
    def current_A(self):
        return self.charge_per_bin_C / np.diff(self.edges_s)

    @property
    def total_charge_C(self):
        return float(np.sum(self.charge_per_bin_C))


def single_ion_bins(detector, max_step_s=5e-6):
    """Exact charge in each bin; no renormalization of the finite-gap signal."""
    if not np.isfinite(max_step_s) or max_step_s <= 0:
        raise ValueError('max_step_s must be finite and positive')
    v = detector.mobility_m2_V_s * detector.gap_field_V_m
    transit = detector.gap_m / v
    edges = np.linspace(-transit, 0., int(np.ceil(transit/max_step_s))+1)
    phi = square_phi(np.maximum(-v*edges, 0.), detector.pad_width_m/2)
    return CurrentWaveform(edges, E_CHARGE*np.diff(phi))


def arrival_counts(times_s, step_s, ions, kind, sigma_s=0., duration_s=0.,
                   separation_s=0.):
    """Expected arrivals per bin, placed at bin centres for convolution.

    Gaussian bin probabilities are integrated analytically. Tails outside the
    supplied grid remain excluded and are never rescaled. A uniform packet
    is a diagnostic, not a model of the event ensemble.
    """
    if not np.isfinite(ions) or ions <= 0 or step_s <= 0:
        raise ValueError('ions and step_s must be positive')
    t = np.asarray(times_s, dtype=float)
    if kind == 'single':
        out = np.zeros_like(t)
        i = np.argmin(abs(t))
        if abs(t[i]) > step_s*1e-6:
            raise ValueError('The arrival grid must contain zero')
        out[i] = ions
        return out
    lo, hi = t-step_s/2, t+step_s/2
    if kind in ('gaussian', 'double_gaussian'):
        if sigma_s <= 0 or separation_s < 0:
            raise ValueError('sigma_s > 0 and separation_s >= 0 required')
        shifts = [0.] if kind == 'gaussian' else [-separation_s/2, separation_s/2]
        return ions/len(shifts)*sum(ndtr((hi-m)/sigma_s)-ndtr((lo-m)/sigma_s)
                                    for m in shifts)
    if kind == 'uniform':
        if duration_s <= 0:
            raise ValueError('duration_s must be positive')
        overlap = np.maximum(0., np.minimum(hi, duration_s/2)
                             - np.maximum(lo, -duration_s/2))
        return ions*overlap/duration_s
    raise ValueError('Unknown arrival distribution')


def make_waveform(detector, kind, ions=260., drift_length_m=1.5,
                  duration_s=.01, separation_s=.04, max_step_s=5e-6,
                  tail_sigmas=8., post_window_s=.05):
    """Convolve the finite-gap kernel with a labelled one-pad arrival model.

    Gaussian temporal spread is thermal bulk sigma_L / bulk velocity, carried
    through identical final trajectories. Nonuniform focusing and additional
    final-region diffusion are not modelled. Charge bins integrate the kernel
    exactly; placing arrival-bin mass at its centre is a convergent approximation.
    """
    if drift_length_m < 0 or tail_sigmas <= 0 or post_window_s < 0:
        raise ValueError('Invalid drift distance or time-window settings')
    kernel = single_ion_bins(detector, max_step_s)
    dt = float(kernel.edges_s[1]-kernel.edges_s[0])
    sigma = float(thermal_spread(drift_length_m, detector.bulk_field_V_m,
                                detector.temperature_K)/detector.bulk_speed_m_s)
    if kind == 'single':
        extent = 0.
    elif kind == 'uniform':
        extent = duration_s/2
    elif kind in ('gaussian', 'double_gaussian'):
        extent = tail_sigmas*sigma + (separation_s/2 if kind == 'double_gaussian' else 0.)
    else:
        raise ValueError('Unknown arrival distribution')
    m = int(np.ceil(extent/dt))
    times = np.arange(-m, m+1)*dt
    counts = arrival_counts(times, dt, ions, kind, sigma, duration_s, separation_s)
    q = fftconvolve(counts, kernel.charge_per_bin_C)
    q = np.pad(q, (0, int(np.ceil(post_window_s/dt))))
    edges = kernel.edges_s[0]-m*dt+np.arange(len(q)+1)*dt
    waveform = CurrentWaveform(edges, q)
    metadata = dict(kind=kind, requested_ions=ions,
                    represented_arrivals=float(counts.sum()),
                    arrival_tail_ions_excluded=float(ions-counts.sum()),
                    drift_length_m=drift_length_m, thermal_sigma_t_s=sigma,
                    uniform_duration_s=duration_s if kind == 'uniform' else None,
                    packet_separation_s=separation_s if kind == 'double_gaussian' else None,
                    kernel_charge_C=kernel.total_charge_C,
                    finite_gap_charge_fraction=kernel.total_charge_C/E_CHARGE,
                    step_s=dt, input_charge_C=waveform.total_charge_C,
                    assumptions='Central square pad in grounded infinite plane; constant gap speed; '
                    'finite collection-gap path; specified one-pad arrivals; '
                    'thermal bulk diffusion for Gaussian cases; no focusing-map or event simulation')
    return waveform, metadata


def csa_response(waveform, feedback_capacitance_F, discharge_time_s=np.inf):
    """Output at bin edges, starting from zero relative to the DC baseline.

    Exact zero-order-hold solution for each bin's average current:
    dV/dt + V/tau = -I/Cf. Infinite tau is an ideal charge integrator.
    There is no clipping, fitted MOS resistance, noise or assumed phase margin.
    """
    cf, tau = feedback_capacitance_F, discharge_time_s
    if not np.isfinite(cf) or cf <= 0 or np.isnan(tau) or tau <= 0:
        raise ValueError('Cf must be finite and positive; tau must be positive')
    widths = np.diff(waveform.edges_s)
    if len(widths) != len(waveform.charge_per_bin_C) or np.any(widths <= 0):
        raise ValueError('Invalid waveform bins')
    dt = widths[0]
    if not np.allclose(widths, dt, rtol=1e-9, atol=0):
        raise ValueError('Uniform time bins required')
    if np.isinf(tau):
        return np.r_[0., -np.cumsum(waveform.charge_per_bin_C)/cf]
    alpha = np.exp(-dt/tau)
    gain = -tau*(-np.expm1(-dt/tau))/(dt*cf)
    return np.r_[0., lfilter([gain], [1., -alpha], waveform.charge_per_bin_C)]


def waveform_metrics(waveform, voltage_V, cf_F, sample_interval_s=.01,
                     sample_phase_s=0.):
    """Sample times are n*interval + phase, relative to the arrival centre."""
    if sample_interval_s <= 0:
        raise ValueError('sample_interval_s must be positive')
    edges = waveform.edges_s
    k0 = int(np.ceil((edges[0]-sample_phase_s)/sample_interval_s))
    k1 = int(np.floor((edges[-1]-sample_phase_s)/sample_interval_s))
    sample_t = np.arange(k0,k1+1)*sample_interval_s+sample_phase_s
    sample_v = np.interp(sample_t, edges, voltage_V)
    peak = float(-np.min(voltage_V))
    ideal = waveform.total_charge_C/cf_F
    return dict(peak_output_magnitude_V=peak,
                peak_time_s=float(edges[np.argmin(voltage_V)]),
                ideal_output_magnitude_V=ideal,
                peak_fraction_of_ideal=peak/ideal,
                sampled_peak_magnitude_V=float(np.max(-sample_v)) if len(sample_v) else None,
                sample_interval_s=sample_interval_s, sample_phase_s=sample_phase_s,
                residual_output_V=float(voltage_V[-1]))
