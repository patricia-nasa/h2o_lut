#!/usr/bin/env python3
"""Compare LUT interpolation with direct RADIS at random states."""

from __future__ import annotations
import argparse
import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator
from lut_common import load_config
from radis_backend import RadisLineBackend


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--samples", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    cfg = load_config(args.config)
    ds = xr.open_dataset(cfg["output_file"]).load()
    rng = np.random.default_rng(args.seed)

    paxis = ds.pressure_hpa.values
    taxis = ds.temperature_k.values
    xaxis = ds.h2o_mole_fraction.values
    axes = (paxis, taxis, np.log(xaxis))

    worst = (0.0, None)
    for ich, channel in enumerate(ds.channel_nm.values):
        interp = RegularGridInterpolator(
            axes, ds.k_total.values[ich], bounds_error=True
        )
        backend = RadisLineBackend(
            channel, cfg["databank_name"], cfg["isotopes"],
            cfg["spectral"]["half_window_cm1"],
            cfg["spectral"]["reference_wstep_cm1"],
            cfg["spectral"]["truncation_cm1"],
            cfg["spectral"].get("broadening_method", "voigt"),
        )
        for _ in range(args.samples):
            p = rng.uniform(paxis[0], paxis[-1])
            t = rng.uniform(taxis[0], taxis[-1])
            x = np.exp(rng.uniform(np.log(xaxis[0]), np.log(xaxis[-1])))
            reference = backend.absorption_coefficient_cm1(p, t, x)
            test = float(interp([[p, t, np.log(x)]])[0])
            scale = max(abs(reference), 1.0e-30)
            rel = abs(test - reference) / scale
            if rel > worst[0]:
                worst = (rel, (ich, channel, p, t, x, reference, test))

    print(f"Maximum relative line-only interpolation error: {worst[0]:.6e}")
    print("Worst state:", worst[1])
    if worst[0] > 1.0e-3:
        raise SystemExit(
            "FAILED 0.1% interpolation target. Refine the p/T/H2O grid."
        )


if __name__ == "__main__":
    main()

