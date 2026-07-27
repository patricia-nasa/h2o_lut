"""Official MT_CKD executable adapter.

The standalone release's output variable names can differ between builds.
This adapter identifies the spectral coordinate and absorption variables using
NetCDF attributes/names, and fails loudly instead of silently applying wrong units.
"""

from __future__ import annotations
from pathlib import Path
import shutil
import subprocess
import tempfile
import numpy as np
import xarray as xr


def _spectral_coord(ds):
    for name in ds.coords:
        low = name.lower()
        units = str(ds[name].attrs.get("units", "")).lower()
        if "wave" in low or "cm-1" in units or "cm^-1" in units:
            return name
    for name in ds.variables:
        if ds[name].ndim == 1 and ("wave" in name.lower() or "wn" in name.lower()):
            return name
    raise RuntimeError(f"Could not identify MT_CKD spectral coordinate: {list(ds.variables)}")


def _absorption_vars(ds, coord):
    candidates = []
    for name, da in ds.data_vars.items():
        if coord not in da.dims:
            continue
        text = " ".join([
            name,
            str(da.attrs.get("long_name", "")),
            str(da.attrs.get("standard_name", "")),
        ]).lower()
        units = str(da.attrs.get("units", "")).lower()
        if "abs" in text or "continuum" in text or "cm-1" in units or "cm^-1" in units:
            candidates.append(name)
    if not candidates:
        raise RuntimeError(
            "Could not identify continuum absorption variables. "
            f"Variables/attrs: {[(n, dict(ds[n].attrs)) for n in ds.data_vars]}"
        )
    return candidates


class MTCKDBackend:
    def __init__(self, executable, coefficient_file, maximum_wavenumber_cm1=20000.0):
        self.executable = Path(executable).expanduser().resolve()
        self.coefficient_file = Path(coefficient_file).expanduser().resolve()
        self.maximum_wavenumber_cm1 = float(maximum_wavenumber_cm1)
        if not self.executable.exists():
            raise FileNotFoundError(self.executable)
        if not self.coefficient_file.exists():
            raise FileNotFoundError(self.coefficient_file)

    def absorption_coefficient_cm1(self, wavenumber_cm1, pressure_hpa, temperature_k, x_h2o):
        wn = float(wavenumber_cm1)
        if wn > self.maximum_wavenumber_cm1:
            return 0.0

        # Request a narrow interval around the channel. The continuum varies slowly.
        wmin, wmax, dw = max(0.0, wn - 1.0), wn + 1.0, 0.1
        config = f"""&mt_ckd_input
 p_atm={float(pressure_hpa):.12g}
 t_atm={float(temperature_k):.12g}
 h2o_frac={float(x_h2o):.12g}
 wv1={wmin:.12g}
 wv2={wmax:.12g}
 dwv={dw:.12g}
/
"""
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            shutil.copy2(self.coefficient_file, work / "absco-ref_wv-mt-ckd.nc")
            cfg = work / "mt_ckd.config"
            cfg.write_text(config)

            # Most builds read mt_ckd.config from stdin; keeping the file in cwd
            # also supports builds that open it by name.
            proc = subprocess.run(
                [str(self.executable)],
                cwd=work,
                input=config,
                text=True,
                capture_output=True,
                check=False,
            )
            if proc.returncode != 0:
                raise RuntimeError(
                    f"MT_CKD failed ({proc.returncode})\n"
                    f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
                )
            outputs = [p for p in work.glob("*.nc")
                       if p.name != "absco-ref_wv-mt-ckd.nc"]
            if not outputs:
                raise RuntimeError("MT_CKD created no output NetCDF file.")

            ds = xr.open_dataset(outputs[0]).load()
            coord = _spectral_coord(ds)
            vars_ = _absorption_vars(ds, coord)
            x = np.asarray(ds[coord].values, dtype=np.float64)

            total = np.zeros_like(x)
            used = []
            for name in vars_:
                da = ds[name]
                units = str(da.attrs.get("units", "")).strip().lower()
                y = np.asarray(da.values, dtype=np.float64).squeeze()
                if y.shape != x.shape:
                    continue
                # Only accept already state-scaled absorption coefficients.
                if units in ("cm-1", "cm^-1", "1/cm", "cm**-1"):
                    total += y
                    used.append(name)
            ds.close()

            if not used:
                raise RuntimeError(
                    "MT_CKD output did not expose an absorption coefficient in cm-1. "
                    "Inspect the output NetCDF and adapt mtckd_backend.py to its exact "
                    "variable definitions before generating the LUT."
                )
            return float(np.interp(wn, x, total))

