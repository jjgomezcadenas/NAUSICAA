"""Continuous-time box averages, trapezoidal charge shaping and PSD integration.

Rectangular averages refer to the analog output waveform. Their continuous-time
limit is evaluated; ADC sampling, quantisation and peak-search bias are excluded.
A trapezoid is the difference of equal moving averages whose end times are
separated by T+G: an ideal voltage step produces rise T, flat top G, fall T.
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid


def primitive(time_s, voltage_V, query_s):
    """Exact integral of a piecewise-linear waveform, evaluated at query times."""
    t=np.asarray(time_s);v=np.asarray(voltage_V);x=np.asarray(query_s)
    if t.ndim!=1 or t.shape!=v.shape or np.any(np.diff(t)<=0):
        raise ValueError('Expected increasing times and matching voltages')
    if np.any(x<t[0]-1e-14) or np.any(x>t[-1]+1e-14):
        raise ValueError('Average extends outside available waveform')
    x=np.clip(x,t[0],t[-1]);j=np.minimum(np.searchsorted(t,x,side='right')-1,len(t)-2)
    c=cumulative_trapezoid(v,t,initial=0);dx=x-t[j]
    return c[j]+v[j]*dx+.5*(v[j+1]-v[j])/(t[j+1]-t[j])*dx*dx


def box_average(time_s, voltage_V, end_s, duration_s):
    if duration_s<=0: raise ValueError('Averaging duration must be positive')
    end=np.asarray(end_s)
    return (primitive(time_s,voltage_V,end)-primitive(time_s,voltage_V,end-duration_s))/duration_s


def trapezoid(time_s, voltage_V, end_s, rise_s, flat_s):
    if flat_s<0: raise ValueError('Flat top must be nonnegative')
    return (box_average(time_s,voltage_V,end_s,rise_s)
            -box_average(time_s,voltage_V,np.asarray(end_s)-rise_s-flat_s,rise_s))


def box_transfer(frequency_Hz, duration_s):
    f=np.asarray(frequency_Hz)
    return np.sinc(f*duration_s)*np.exp(-1j*np.pi*f*duration_s)


def trapezoid_transfer(frequency_Hz, rise_s, flat_s):
    f=np.asarray(frequency_Hz)
    return box_transfer(f,rise_s)*(-np.expm1(-2j*np.pi*f*(rise_s+flat_s)))


def averaged_frame_transfer(frequency_Hz, duration_s, interval_s, alpha):
    # Overall phase from centering the boxes does not affect the noise variance.
    f=np.asarray(frequency_Hz)
    return box_transfer(f,duration_s)*((1-alpha)-alpha*np.expm1(-2j*np.pi*f*interval_s))


def filtered_rms(frequency_Hz, output_asd_V_sqrtHz, transfer, phase_span_s, order=8):
    """One-sided stationary PSD; resolve oscillations separately from log grid."""
    f=np.asarray(frequency_Hz);a=np.asarray(output_asd_V_sqrtHz)
    if f.shape!=a.shape or np.any(f<=0) or np.any(np.diff(f)<=0) or np.any(a<0):
        raise ValueError('Invalid noise spectrum')
    if not np.any(a):return 0.
    edges=np.unique(np.r_[f,np.arange(f[0],f[-1],.5/phase_span_s)])
    half=np.diff(edges)/2;mid=(edges[1:]+edges[:-1])/2
    x,w=np.polynomial.legendre.leggauss(order);nodes=mid[:,None]+half[:,None]*x
    psd=np.exp(np.interp(np.log(nodes),np.log(f),2*np.log(np.maximum(a,1e-150))))
    var=np.sum(half[:,None]*w*psd*np.abs(transfer(nodes))**2)
    return float(np.sqrt(var))
