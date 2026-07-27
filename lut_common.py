from __future__ import annotations
from pathlib import Path
import numpy as np
import xarray as xr
import yaml


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def expand_axis(spec):
    if "values" in spec:
        return np.asarray(spec["values"], dtype=np.float64)
    if "logspace" in spec:
        s = spec["logspace"]
        return np.logspace(
            float(s["start_exponent"]),
            float(s["stop_exponent"]),
            int(s["number"]),
        )
    start, stop, step = map(float, (spec["start"], spec["stop"], spec["step"]))
    n = int(np.floor((stop - start) / step + 0.5)) + 1
    return start + np.arange(n, dtype=np.float64) * step


def read_channels_nm(spec):
    kind = spec["type"].lower()
    units = spec.get("input_units", "nm").lower()

    if kind == "ames_brdf":
        ds = xr.open_mfdataset(spec["glob"], combine="by_coords")
        values = np.asarray(ds[spec.get("variable", "nwav")].values, dtype=np.float64)
        ds.close()
    elif kind == "text":
        values = np.loadtxt(Path(spec["path"]), dtype=np.float64)
    else:
        raise ValueError(f"Unsupported channel source type: {kind}")

    values = np.ravel(values)
    if units in ("micron", "microns", "um"):
        values = values * 1000.0
    elif units != "nm":
        raise ValueError(f"Unsupported channel units: {units}")

    if np.any((values < 300.0) | (values > 2500.0)):
        raise ValueError("Channel centers must lie within 300–2500 nm.")
    return values


def q_to_mole_fraction(q):
    epsilon = 0.01801528 / 0.0289647
    q = np.asarray(q, dtype=np.float64)
    r = q / (1.0 - q)
    return r / (epsilon + r)


def validate_monotonic(name, x):
    if np.any(np.diff(x) <= 0):
        raise ValueError(f"{name} must be strictly increasing.")

