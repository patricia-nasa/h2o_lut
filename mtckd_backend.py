"""Official MT_CKD executable adapter.

The standalone release's output variable names can differ between builds.
This adapter identifies the spectral coordinate and absorption variables using
NetCDF attributes/names, and fails loudly instead of silently applying wrong units.
"""

from __future__ import annotations
import os, shutil, subprocess, tempfile
from pathlib import Path
import numpy as np
import xarray as xr

class MTCKDBackend:
    def __init__(self, executable, coefficient_file,
                 maximum_wavenumber_cm1=20000.0,
                 debug_directory=None):
        self.executable = Path(executable).expanduser().resolve()
        self.coefficient_file = Path(coefficient_file).expanduser().resolve()
        self.maximum_wavenumber_cm1 = maximum_wavenumber_cm1
        self.debug_directory = Path(debug_directory).expanduser().resolve() if debug_directory else None

    def _make_namelist(self, pressure_hpa, temperature_k, x_h2o, wv1, wv2, dwv):
        return (
            "&mt_ckd_input\n"
            f" p_atm={pressure_hpa:.12g},\n"
            f" t_atm={temperature_k:.12g},\n"
            f" h2o_frac={x_h2o:.12g},\n"
            f" wv1={wv1:.12g},\n"
            f" wv2={wv2:.12g},\n"
            f" dwv={dwv:.12g},\n"
            "/\n"
        )

    def _env(self):
        env=os.environ.copy()
        if "CONDA_PREFIX" in env:
            lib=str(Path(env["CONDA_PREFIX"])/"lib")
            env["LD_LIBRARY_PATH"]=lib+":"+env.get("LD_LIBRARY_PATH","")
        return env

    def absorption_coefficient_cm1(self,wavenumber_cm1,pressure_hpa,temperature_k,x_h2o):
        wn=float(wavenumber_cm1)
        if wn>self.maximum_wavenumber_cm1:
            return 0.0
        work=Path(tempfile.mkdtemp(prefix="mtckd_"))
        print(f"Temp directory: {work}")
        shutil.copy2(self.coefficient_file, work/"absco-ref_wv-mt-ckd.nc")
        cfg=work/"mt_ckd.config"
        cfg.write_text(self._make_namelist(pressure_hpa,temperature_k,x_h2o,max(0,wn-1),wn+1,0.1))
        with cfg.open() as fin:
            proc=subprocess.run([str(self.executable)],cwd=work,stdin=fin,
                                text=True,capture_output=True,env=self._env())
        if proc.returncode!=0:
            raise RuntimeError(f"MT_CKD failed ({proc.returncode})\nSTDERR:\n{proc.stderr}")
        ds=xr.open_dataset(work/"mt_ckd_h2o_output.nc")
        return float(np.interp(
            wn,
            ds["wavenumbers"].values,
            ds["self_absorption"].values+ds["frgn_absorption"].values
        ))
