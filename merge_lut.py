#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import xarray as xr
from lut_common import load_config, read_channels_nm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    args = ap.parse_args()
    cfg = load_config(args.config)
    channels = read_channels_nm(cfg["channel_source"])
    out = Path(cfg["output_file"]).expanduser().resolve()
    part_dir = out.parent / (out.stem + "_parts")

    parts = sorted(part_dir.glob("channel_*.nc"))
    if len(parts) != len(channels):
        raise RuntimeError(f"Found {len(parts)} parts for {len(channels)} channels.")

    datasets = []
    for p in parts:
        ds = xr.open_dataset(p)
        ich = int(ds.attrs["channel_index"])
        datasets.append(ds.expand_dims(channel=[ich]))

    merged = xr.concat(datasets, dim="channel")
    merged = merged.assign_coords(
        channel_nm=("channel", channels),
        wavenumber_cm1=("channel", 1.0e7 / channels),
    )
    merged.attrs.update({
        "description": "H2O absorption coefficient LUT for PyVLIDORT SBG",
        "runtime_formula": "alpha = k_total * abs(diff(ze_m)) * 100",
        "alpha_units": "dimensionless",
    })
    enc = {v: {"zlib": True, "complevel": int(cfg.get("compression", {}).get("complevel", 4))}
           for v in merged.data_vars}
    merged.to_netcdf(out, encoding=enc)
    for ds in datasets:
        ds.close()
    print(out)


if __name__ == "__main__":
    main()

