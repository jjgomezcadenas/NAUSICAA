"""SI-unit analytical models. See the LaTeX note for boundary conditions.

Current convention: positive equivalent source current delivered into the
selected pad's readout input for a positive ion approaching that pad.
"""

from dataclasses import dataclass
from pathlib import Path
import tomllib
import numpy as np

E_CHARGE = 1.602176634e-19  # exact SI, C
K_B = 1.380649e-23         # exact SI, J/K


@dataclass(frozen=True)
class Detector:
    hole_diameter_m: float
    pitch_m: float
    ifg_thickness_m: float
    gap_m: float
    pad_width_m: float
    bulk_reference_length_m: float
    pad_voltage_V: float
    ic_voltage_V: float
    gate_voltage_V: float
    bulk_field_V_m: float
    pressure_bar: float
    temperature_K: float
    bulk_speed_m_s: float
    charge_e: int = 1

    def __post_init__(self):
        for name in ("hole_diameter_m", "pitch_m", "ifg_thickness_m", "gap_m",
                     "pad_width_m", "bulk_reference_length_m", "bulk_field_V_m",
                     "pressure_bar", "temperature_K", "bulk_speed_m_s"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive")
        if not self.gate_voltage_V > self.ic_voltage_V > self.pad_voltage_V:
            raise ValueError("Positive-ion benchmark requires gate > IC > pad")
        if self.hole_diameter_m >= self.pitch_m or self.pad_width_m >= self.pitch_m:
            raise ValueError("Hole and pad dimensions must be smaller than pitch")
        if self.charge_e != 1:
            raise ValueError("This milestone models singly charged molecular Xe2+")

    @classmethod
    def from_toml(cls, path: str | Path):
        data = tomllib.loads(Path(path).read_text())
        if data['geometry'].get('pad_model') != 'square_in_grounded_plane':
            raise ValueError('Only the documented square-in-grounded-plane benchmark is implemented')
        flat = {k: v for section in ("geometry", "bias", "gas")
                for k, v in data[section].items() if k in cls.__dataclass_fields__}
        return cls(**flat)

    @property
    def open_fraction(self):
        return np.pi*(self.hole_diameter_m/2)**2/self.pitch_m**2

    @property
    def hole_field_V_m(self):
        return (self.gate_voltage_V-self.ic_voltage_V)/self.ifg_thickness_m

    @property
    def gap_field_V_m(self):
        return (self.ic_voltage_V-self.pad_voltage_V)/self.gap_m

    @property
    def mobility_m2_V_s(self):
        return self.bulk_speed_m_s/self.bulk_field_V_m

    @property
    def thermal_D_m2_s(self):
        return self.mobility_m2_V_s*K_B*self.temperature_K/E_CHARGE

    def regions(self):
        """Nominal uniform-field regions; not a perforated-grid field solution."""
        return [("Bulk reference", self.bulk_reference_length_m, self.bulk_field_V_m),
                ("IFG hole", self.ifg_thickness_m, self.hole_field_V_m),
                ("Collection gap", self.gap_m, self.gap_field_V_m)]


def thermal_spread(length_m, field_V_m, temperature_K):
    """One-coordinate RMS sqrt(2 kT L / eE), uniform thermal drift only."""
    length = np.asarray(length_m, dtype=float)
    if np.any(length < 0) or field_V_m <= 0 or temperature_K <= 0:
        raise ValueError("Length >= 0 and field/temperature > 0 required")
    return np.sqrt(2*K_B*temperature_K/E_CHARGE*length/field_V_m)


def _z_and_a(z_m, half_width_m):
    z = np.asarray(z_m, dtype=float)
    if half_width_m <= 0 or np.any(z < 0):
        raise ValueError("Require a > 0 and z >= 0")
    return z, half_width_m


def square_phi(z_m, half_width_m):
    """On-axis weighting potential, square 2a by 2a in a grounded plane."""
    z, a = _z_and_a(z_m, half_width_m)
    return 2/np.pi*np.arctan2(a*a, z*np.sqrt(z*z+2*a*a))


def square_weight(z_m, half_width_m):
    """Positive z-component of E_w = -grad(phi), in 1/m."""
    z, a = _z_and_a(z_m, half_width_m)
    return 4*a*a/(np.pi*(z*z+a*a)*np.sqrt(z*z+2*a*a))


def disk_phi(z_m, radius_m):
    """On-axis weighting potential of a disk in the same grounded plane."""
    z, a = _z_and_a(z_m, radius_m)
    # Rationalized form avoids cancellation in 1-z/sqrt(z^2+a^2).
    r = np.sqrt(z*z+a*a)
    return a*a/(r*(r+z))


def disk_weight(z_m, radius_m):
    z, a = _z_and_a(z_m, radius_m)
    return a*a/(z*z+a*a)**1.5


def square_z_at_phi(potential, half_width_m):
    """Inverse on-axis phi for 0 < phi <= 1."""
    p = np.asarray(potential, dtype=float)
    if np.any((p <= 0) | (p > 1)) or half_width_m <= 0:
        raise ValueError("Require 0 < phi <= 1 and a > 0")
    return half_width_m*np.sqrt(1/np.sin(np.pi*p/2)-1)


def single_ion_current(t_s, speed_m_s, half_width_m, charge_C=E_CHARGE):
    """Central constant-speed ion; collection at t=0, current zero for t>=0.

    The infinite-past tail is the analytical half-space benchmark, not the
    actual current through the upstream IFG. Integral over the full model = q.
    """
    if speed_m_s <= 0 or half_width_m <= 0:
        raise ValueError("Speed and half-width must be positive")
    t = np.asarray(t_s, dtype=float)
    z = np.maximum(-speed_m_s*t, 0)
    return np.where(t < 0, charge_C*speed_m_s*square_weight(z, half_width_m), 0.)
