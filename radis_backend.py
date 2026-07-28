from __future__ import annotations
import numpy as np
from radis import SpectrumFactory


class RadisLineBackend:
    """Reusable HITRAN line-absorption calculator for one channel."""

    def __init__(
        self,
        channel_nm,
        databank_name,
        isotopes,
        half_window_cm1,
        wstep_cm1,
        truncation_cm1,
        broadening_method="voigt",
    ):
        self.channel_nm = float(channel_nm)
        self.center_cm1 = 1.0e7 / self.channel_nm
        self.sf = SpectrumFactory(
            wavenum_min=self.center_cm1 - float(half_window_cm1),
            wavenum_max=self.center_cm1 + float(half_window_cm1),
            wstep=float(wstep_cm1),
            molecule="H2O",
            isotope=str(isotopes),
            medium="air",
            truncation=float(truncation_cm1),
            neighbour_lines=float(truncation_cm1),
            cutoff=0.0,
            broadening_method=broadening_method,
            optimization=None,
            db_use_cached=True,
            verbose=False,
        )
        self.sf.load_databank(databank_name, load_columns="equilibrium",load_energies=False)

    def absorption_coefficient_cm1(self, pressure_hpa, temperature_k, x_h2o):
        # Use a 1 cm slab: absorbance = k [cm-1] * 1 cm.
        spec = self.sf.eq_spectrum(
            Tgas=float(temperature_k),
            pressure=float(pressure_hpa) / 1000.0,
            mole_fraction=float(x_h2o),
            path_length=1.0,
        )
        nu, k = spec.get("abscoeff", wunit="cm-1", Iunit="cm-1")
        return float(np.interp(self.center_cm1, np.asarray(nu), np.asarray(k)))

