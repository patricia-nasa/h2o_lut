#!/usr/bin/env python3
"""Generate a channel/p/T/H2O absorption-coefficient LUT."""

from __future__ import annotations
import argparse
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import xarray as xr
from tqdm import tqdm

from lut_common import load_config, expand_axis, read_channels_nm, validate_monotonic
from radis_backend import RadisLineBackend
from mtckd_backend import MTCKDBackend


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--channel-start", type=int, default=0)
    ap.add_argument("--channel-stop", type=int)
    args = ap.parse_args()

    cfg = load_config(args.config)
    channels = read_channels_nm(cfg["channel_source"])
    pressure = expand_axis(cfg["grid"]["pressure_hpa"])
    temperature = expand_axis(cfg["grid"]["temperature_k"])
    xh2o = expand_axis(cfg["grid"]["h2o_mole_fraction"])

    for n, a in [("channel_nm", channels), ("pressure_hpa", pressure),
                 ("temperature_k", temperature), ("h2o_mole_fraction", xh2o)]:
        validate_monotonic(n, a)

    stop = len(channels) if args.channel_stop is None else min(args.channel_stop, len(channels))
    sel = np.arange(args.channel_start, stop)
    if sel.size == 0:
        raise ValueError("Empty channel selection.")

    out = Path(cfg["output_file"]).expanduser().resolve()
    part_dir = out.parent / (out.stem + "_parts")
    part_dir.mkdir(parents=True, exist_ok=True)

    cont_cfg = cfg.get("continuum", {})
    continuum = None
    if cont_cfg.get("enabled", False):
        continuum = MTCKDBackend(
            cont_cfg["executable"],
            cont_cfg["coefficient_file"],
            cont_cfg.get("maximum_wavenumber_cm1", 20000.0),
        )

    sp = cfg["spectral"]
    for ich in sel:
        part = part_dir / f"channel_{ich:05d}.nc"
        if part.exists():
            print(f"Skipping completed {part.name}")
            continue

        channel_nm = channels[ich]
        wn = 1.0e7 / channel_nm
        backend = RadisLineBackend(
            channel_nm=channel_nm,
            databank_name=cfg["databank_name"],
            isotopes=cfg["isotopes"],
            half_window_cm1=sp["half_window_cm1"],
            wstep_cm1=sp["wstep_cm1"],
            truncation_cm1=sp["truncation_cm1"],
            broadening_method=sp.get("broadening_method", "voigt"),
        )

        shape = (pressure.size, temperature.size, xh2o.size)
        k_line = np.empty(shape, dtype=np.float64)
        k_cont = np.zeros(shape, dtype=np.float64)

        total_states = int(np.prod(shape))
        bar = tqdm(total=total_states, desc=f"{ich}: {channel_nm:.3f} nm")
        for ip, p in enumerate(pressure):
            for it, t in enumerate(temperature):
                for ix, x in enumerate(xh2o):
                    k_line[ip, it, ix] = backend.absorption_coefficient_cm1(p, t, x)
                    if continuum is not None:
                        k_cont[ip, it, ix] = continuum.absorption_coefficient_cm1(
                            wn, p, t, x
                        )
                    bar.update(1)
        bar.close()

        ds = xr.Dataset(
            data_vars={
                "k_line": (
                    ("pressure_hpa", "temperature_k", "h2o_mole_fraction"),
                    k_line,
                    {"units": "cm-1", "long_name": "HITRAN2024 local-line absorption coefficient"},
                ),
                "k_continuum": (
                    ("pressure_hpa", "temperature_k", "h2o_mole_fraction"),
                    k_cont,
                    {"units": "cm-1", "long_name": "MT_CKD H2O continuum absorption coefficient"},
                ),
                "k_total": (
                    ("pressure_hpa", "temperature_k", "h2o_mole_fraction"),
                    k_line + k_cont,
                    {"units": "cm-1", "long_name": "total H2O absorption coefficient"},
                ),
            },
            coords={
                "pressure_hpa": pressure,
                "temperature_k": temperature,
                "h2o_mole_fraction": xh2o,
            },
            attrs={
                "channel_index": int(ich),
                "channel_nm": float(channel_nm),
                "wavenumber_cm-1": float(wn),
                "HITRAN_edition": "2024",
                "isotopes": str(cfg["isotopes"]),
                "created_utc": datetime.now(timezone.utc).isoformat(),
            },
        )
        enc = {v: {"zlib": True, "complevel": int(cfg.get("compression", {}).get("complevel", 4))}
               for v in ds.data_vars}
        ds.to_netcdf(part, encoding=enc)

    print(f"\nChannel parts written to {part_dir}")
    print("Run merge_lut.py to assemble the final LUT.")


if __name__ == "__main__":
    main()

