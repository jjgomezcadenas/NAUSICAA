"""Finite-gain, one-pole CSA benchmark; SI units and one-sided PSDs per Hz.

The input capacitance excludes Cf. All voltages are small-signal departures
from a chosen operating point. This is not a transistor/foundry model.
"""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import expm
from scipy.integrate import quad
from .physics import E_CHARGE, K_B, square_phi
from .readout import CurrentWaveform


@dataclass(frozen=True)
class CSA:
    feedback_capacitance_F: float
    input_capacitance_F: float
    feedback_resistance_ohm: float
    open_loop_gain: float
    amplifier_pole_Hz: float

    def __post_init__(self):
        for name, value in vars(self).items():
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f'{name} must be finite and positive')

    @property
    def nominal_discharge_time_s(self):
        return self.feedback_resistance_ohm*self.feedback_capacitance_F

    def amplifier_gain(self, frequency_Hz):
        """Complex open-loop gain A(i2pi f); frequency in Hz, gain dimensionless."""
        f=np.asarray(frequency_Hz,dtype=float)
        if np.any(~np.isfinite(f)) or np.any(f<0):
            raise ValueError('Frequency must be finite and nonnegative')
        return self.open_loop_gain/(1+1j*f/self.amplifier_pole_Hz)

    def transfer(self, frequency_Hz):
        """(output V/input A, output V/input voltage-noise V), s=i2pi f."""
        f=np.asarray(frequency_Hz,dtype=float)
        if np.any(~np.isfinite(f)) or np.any(f<0):
            raise ValueError('Frequency must be finite and nonnegative')
        s=2j*np.pi*f
        a=self.amplifier_gain(f)
        yf=1/self.feedback_resistance_ohm+s*self.feedback_capacitance_F
        denominator=s*self.input_capacitance_F+(1+a)*yf
        return -a/denominator, a*(s*self.input_capacitance_F+yf)/denominator

    def state_matrix(self):
        """State [input-node V, output V]; dx/dt=M x+[i/(Cin+Cf),0]."""
        cf=self.feedback_capacitance_F; ct=cf+self.input_capacitance_F
        g=1/self.feedback_resistance_ohm
        w=2*np.pi*self.amplifier_pole_Hz; a=self.open_loop_gain
        return np.array([[-(cf*w*a+g)/ct,(g-cf*w)/ct],[-w*a,-w]])

    def poles_per_s(self):
        return np.linalg.eigvals(self.state_matrix())

    def charge_gain_reduction(self):
        """Fast-amplifier, negligible-discharge charge-sharing limit."""
        cf=self.feedback_capacitance_F
        return cf/(cf+(self.input_capacitance_F+cf)/self.open_loop_gain)


def response(waveform, model, initial_state_V=None):
    """Exact ZOH response to bin-average currents; no reset between frames.

    Returns Nx2 input/output voltages at waveform edges. Matrix exponential
    integration is exact for the piecewise-constant current representation.
    Refining the physical-current bins is required near collection.
    """
    edges=np.asarray(waveform.edges_s); charges=np.asarray(waveform.charge_per_bin_C)
    widths=np.diff(edges)
    if edges.ndim!=1 or len(widths)!=len(charges) or np.any(widths<=0) or np.any(~np.isfinite(edges)) or np.any(~np.isfinite(charges)):
        raise ValueError('Invalid waveform')
    state=np.zeros(2) if initial_state_V is None else np.asarray(initial_state_V,dtype=float)
    if state.shape!=(2,) or np.any(~np.isfinite(state)):
        raise ValueError('initial_state_V must contain two finite voltages')
    out=np.empty((len(edges),2)); out[0]=state
    augmented=np.zeros((3,3)); augmented[:2,:2]=model.state_matrix(); augmented[0,2]=1.
    ct=model.input_capacitance_F+model.feedback_capacitance_F
    # Uniform grids are overwhelmingly common; tolerate floating subtraction.
    uniform=np.allclose(widths,widths[0],rtol=1e-9,atol=0)
    phi=expm(augmented*widths[0])
    for j,(dt,q) in enumerate(zip(widths,charges)):
        if not uniform:
            phi=expm(augmented*dt)
        state=phi[:2,:2]@state+phi[:2,2]*(q/dt/ct)
        out[j+1]=state
    return out


def ions_waveform(detector, arrival_times_s, start_s, end_s, max_step_s):
    """Explicit ion arrivals; exact finite-gap charge in every numerical bin.

    A single central grounded-plane pad is assumed for every ion. Arrivals are
    collection times. No tails are rescaled when the supplied window truncates
    a trajectory. Each ion has charge +e, including Xe_2^+.
    """
    arrivals=np.asarray(arrival_times_s,dtype=float)
    if arrivals.ndim!=1 or len(arrivals)==0 or np.any(~np.isfinite(arrivals)):
        raise ValueError('Supply a nonempty finite arrival-time vector')
    if not np.isfinite(start_s+end_s+max_step_s) or end_s<=start_s or max_step_s<=0:
        raise ValueError('Invalid time window or step')
    edges=np.linspace(start_s,end_s,int(np.ceil((end_s-start_s)/max_step_s))+1)
    v=detector.mobility_m2_V_s*detector.gap_field_V_m
    half=detector.pad_width_m/2; gap=detector.gap_m
    q=np.zeros(len(edges)-1)
    # At collection phi=1, then stays 1: no appended collection impulse.
    for arrival in arrivals:
        heights=np.clip(v*(arrival-edges),0,gap)
        q+=E_CHARGE*np.diff(square_phi(heights,half))
    return CurrentWaveform(edges,q)


def frame_charge(waveform, sample_times_s):
    """Charge crossing each sample interval, for the bin-average current."""
    sample=np.asarray(sample_times_s,dtype=float)
    if sample.ndim!=1 or len(sample)<2 or np.any(np.diff(sample)<=0):
        raise ValueError('Need increasing sample times')
    if sample[0]<waveform.edges_s[0] or sample[-1]>waveform.edges_s[-1]:
        raise ValueError('Samples lie outside waveform')
    cumulative=np.r_[0.,np.cumsum(waveform.charge_per_bin_C)]
    return np.diff(np.interp(sample,waveform.edges_s,cumulative))


def sampled_voltage(edges_s, voltage_V, sample_times_s, aperture_s=0.):
    """Trailing rectangular aperture ending at each sample; zero means point.

    Uses the piecewise-linear interpolation of calculated voltages. Its
    aperture integral is evaluated exactly on that interpolant. Decrease the
    current step to check its accuracy against the continuous response.
    """
    edges=np.asarray(edges_s); volts=np.asarray(voltage_V); times=np.asarray(sample_times_s)
    if aperture_s<0 or not np.isfinite(aperture_s):
        raise ValueError('Aperture must be finite and nonnegative')
    if np.any(times-aperture_s<edges[0]) or np.any(times>edges[-1]):
        raise ValueError('Sampling aperture lies outside response window')
    if aperture_s==0:
        return np.interp(times,edges,volts)
    dt=np.diff(edges); slopes=np.diff(volts)/dt
    integ=np.r_[0.,np.cumsum(dt*(volts[:-1]+volts[1:])/2)]
    def integral_at(t):
        ix=np.clip(np.searchsorted(edges,t,side='right')-1,0,len(dt)-1)
        h=t-edges[ix]
        return integ[ix]+volts[ix]*h+.5*slopes[ix]*h*h
    return (integral_at(times)-integral_at(times-aperture_s))/aperture_s


def estimator_coefficients(model, interval_s):
    """Nominal ideal constant-current estimator; not an unbiased inverse here."""
    if not np.isfinite(interval_s) or interval_s<=0:
        raise ValueError('interval_s must be finite and positive')
    x=interval_s/model.nominal_discharge_time_s
    alpha=np.exp(-x)
    k=model.feedback_capacitance_F*x/(-np.expm1(-x))
    return alpha,k


def recovered_charge(sampled_voltage_V, model, interval_s):
    alpha,k=estimator_coefficients(model,interval_s)
    voltage=np.asarray(sampled_voltage_V)
    return -k*(voltage[1:]-alpha*voltage[:-1])


@dataclass(frozen=True)
class Noise:
    voltage_white_V_per_sqrtHz: float
    voltage_flicker_V2: float
    current_white_A_per_sqrtHz: float
    low_frequency_Hz: float
    high_frequency_Hz: float
    temperature_K: float=293.15
    passive_feedback_thermal: bool=False

    def __post_init__(self):
        for name in ('voltage_white_V_per_sqrtHz','voltage_flicker_V2','current_white_A_per_sqrtHz'):
            v=getattr(self,name)
            if not np.isfinite(v) or v<0:
                raise ValueError(f'{name} must be finite and nonnegative')
        if not 0<self.low_frequency_Hz<self.high_frequency_Hz or not np.isfinite(self.high_frequency_Hz) or not np.isfinite(self.temperature_K) or self.temperature_K<=0:
            raise ValueError('Finite positive frequency bounds and temperature required')

    def output_psd(self, model, frequency_Hz, aperture_s=0.):
        """One-sided V²/Hz after rectangular-aperture filtering.

        Passive resistor 4kT/Rf is optional, never a calibrated MOS noise law.
        No noise sources are correlated in this benchmark. The explicit finite
        frequency band is part of the model, not an inferred physical cutoff.
        """
        f=np.asarray(frequency_Hz,dtype=float)
        if np.any(f<=0) or np.any(~np.isfinite(f)) or aperture_s<0 or not np.isfinite(aperture_s):
            raise ValueError('Positive frequencies and nonnegative aperture required')
        hi,he=model.transfer(f)
        se=self.voltage_white_V_per_sqrtHz**2+self.voltage_flicker_V2/f
        si=self.current_white_A_per_sqrtHz**2
        if self.passive_feedback_thermal:
            si+=4*K_B*self.temperature_K/model.feedback_resistance_ohm
        band=(f>=self.low_frequency_Hz)&(f<=self.high_frequency_Hz)
        return (abs(hi)**2*si+abs(he)**2*se)*np.sinc(f*aperture_s)**2*band


def noise_covariance(model, noise, lag_s=0., aperture_s=0., relative_tolerance=1e-7):
    """Covariance via independent per-decade quadrature, one-sided Hz PSD.

    Oscillatory covariance integrals use QUADPACK's cosine-weighted rule.
    Rescaling each decade avoids applying an absolute 1e-8 tolerance to tiny
    SI noise powers. Relative convergence should be checked for each study.
    """
    if not np.isfinite(lag_s):
        raise ValueError('Lag must be finite')
    lo,hi=noise.low_frequency_Hz,noise.high_frequency_Hz
    breaks=np.unique(np.r_[lo,10.**np.arange(np.ceil(np.log10(lo)),np.floor(np.log10(hi))+1),hi])
    total=0.
    for low,high in zip(breaks[:-1],breaks[1:]):
        scale=float(np.max(noise.output_psd(model,np.geomspace(low,high,31),aperture_s)))
        if scale==0:
            continue
        width=high-low
        def scaled(f):
            return float(noise.output_psd(model,f,aperture_s))/scale
        if lag_s==0:
            # Log-frequency integration resolves low-frequency current noise.
            answer=quad(lambda u:scaled(np.exp(u))*np.exp(u),np.log(low),np.log(high),
                        epsabs=relative_tolerance*width,epsrel=relative_tolerance,limit=200)[0]
        else:
            answer=quad(scaled,low,high,weight='cos',wvar=2*np.pi*abs(lag_s),
                        epsabs=relative_tolerance*width,epsrel=relative_tolerance,
                        limit=200,maxp1=100)[0]
        total+=scale*answer
    return float(total)


def noise_metrics(model,noise,interval_s,aperture_s=0.,relative_tolerance=1e-7):
    c0=noise_covariance(model,noise,0.,aperture_s,relative_tolerance)
    c1=noise_covariance(model,noise,interval_s,aperture_s,relative_tolerance)
    alpha,k=estimator_coefficients(model,interval_s)
    variance=k*k*((1+alpha*alpha)*c0-2*alpha*c1)
    if variance<0:
        raise ArithmeticError('Negative variance; investigate quadrature')
    sigma=np.sqrt(variance)
    return dict(sample_voltage_rms_V=float(np.sqrt(c0)),
                adjacent_sample_covariance_V2=c1, sample_variance_V2=c0,
                adjacent_sample_correlation=c1/c0 if c0 else 0.,
                recovered_charge_rms_C=float(sigma),
                induced_charge_ENC_e=float(sigma/E_CHARGE),
                interval_s=interval_s,aperture_s=aperture_s,
                low_frequency_Hz=noise.low_frequency_Hz,
                high_frequency_Hz=noise.high_frequency_Hz)


def pwl_current(waveform):
    """Centre-based PWL current preserves total bin charge on the full grid.

    Physical time is shifted so the first point is zero. The endpoint half
    bins are flat; no amplitude renormalization is applied. Exact for uniform
    bin widths. Discontinuities are linearly smoothed over one numerical bin.
    """
    widths=np.diff(waveform.edges_s)
    if not np.allclose(widths,widths[0],rtol=1e-9,atol=0):
        raise ValueError('PWL charge guarantee requires uniform bins')
    centres=(waveform.edges_s[:-1]+waveform.edges_s[1:])/2
    t=np.r_[waveform.edges_s[0],centres,waveform.edges_s[-1]]-waveform.edges_s[0]
    current=np.r_[waveform.current_A[0],waveform.current_A,waveform.current_A[-1]]
    q=float(np.trapezoid(current,t))
    if not np.isclose(q,waveform.total_charge_C,rtol=1e-10,atol=1e-30):
        raise ArithmeticError('PWL charge check failed')
    return t,current
