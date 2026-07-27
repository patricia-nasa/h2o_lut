"""Runtime vectorized LUT interpolation for sbg_vlidort.py."""

from __future__ import annotations
import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator

M_H2O = 0.01801528
M_DRY = 0.0289647
EPSILON = M_H2O / M_DRY


def q_to_h2o_mole_fraction(q):
    q = np.asarray(q, dtype=np.float64)
    r = q / (1.0 - q)
    return r / (EPSILON + r)


def log_pressure_midpoint(pe_pa):
    p1 = np.asarray(pe_pa[:-1], dtype=np.float64)
    p2 = np.asarray(pe_pa[1:], dtype=np.float64)
    return (p1 - p2) / np.log(p1 / p2)


class H2OAbsorptionLUT:
    def __init__(self, path):
        self.ds = xr.open_dataset(path).load()
        self.channels = np.asarray(self.ds.channel_nm.values)
        self.axes = (
            np.asarray(self.ds.pressure_hpa.values),
            np.asarray(self.ds.temperature_k.values),
            np.log(np.asarray(self.ds.h2o_mole_fraction.values)),
        )
        data = np.asarray(self.ds.k_total.values)
        self._interpolators = [
            RegularGridInterpolator(
                self.axes,
                data[ich],
                method="linear",
                bounds_error=False,
                fill_value=None,
            )
            for ich in range(data.shape[0])
        ]

    def nearest_channel_index(self, channel_nm, tolerance_nm=1.0e-6):
        ich = int(np.argmin(np.abs(self.channels - channel_nm)))
        err = abs(self.channels[ich] - channel_nm)
        if err > tolerance_nm:
            raise ValueError(
                f"No LUT channel matches {channel_nm} nm; nearest is "
                f"{self.channels[ich]} nm (difference {err} nm)."
            )
        return ich

    def alpha(self, channel_nm, pe_pa, te_k, ze_m, qv):
        """
        Return dimensionless alpha with shape (nlev, 1, nobs).
        Inputs:
          pe_pa, te_k, ze_m: (nlev+1,nobs)
          qv:                (nlev,nobs), kg/kg
        """
        pe_pa = np.asarray(pe_pa, dtype=np.float64)
        te_k = np.asarray(te_k, dtype=np.float64)
        ze_m = np.asarray(ze_m, dtype=np.float64)
        qv = np.asarray(qv, dtype=np.float64)

        p_hpa = log_pressure_midpoint(pe_pa) / 100.0
        t_k = 0.5 * (te_k[:-1] + te_k[1:])
        dz_cm = np.abs(np.diff(ze_m, axis=0)) * 100.0
        x = q_to_h2o_mole_fraction(qv)

        if not (p_hpa.shape == t_k.shape == dz_cm.shape == x.shape):
            raise ValueError(
                f"Incompatible shapes p={p_hpa.shape}, T={t_k.shape}, "
                f"dz={dz_cm.shape}, q={x.shape}"
            )

        ich = self.nearest_channel_index(float(channel_nm))
        points = np.column_stack([
            p_hpa.ravel(),
            t_k.ravel(),
            np.log(np.clip(x.ravel(), np.exp(self.axes[2][0]), np.exp(self.axes[2][-1]))),
        ])
        # Explicit clipping avoids uncontrolled extrapolation.
        for j, axis in enumerate(self.axes):
            points[:, j] = np.clip(points[:, j], axis[0], axis[-1])

        k = self._interpolators[ich](points).reshape(p_hpa.shape)
        alpha = k * dz_cm
        return alpha[:, None, :]

