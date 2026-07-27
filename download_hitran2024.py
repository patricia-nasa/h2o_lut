#!/usr/bin/env python3
"""Download/cache HITRAN2024 H2O and register it with RADIS."""

from __future__ import annotations
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
import radis
from radis.io.hitran import fetch_hitran


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--databank-name", default="HITRAN2024-H2O")
    ap.add_argument("--isotopes", default="1,2,3,4,5,6,7")
    args = ap.parse_args()

    out = Path(args.output_dir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    # 300–2500 nm = 33333.333–4000 cm-1. Add 25 cm-1 margins for line wings.
    wmin = 3975.0
    wmax = 33358.3333333333

    result = fetch_hitran(
        molecule="H2O",
        local_databases=str(out),
        databank_name=args.databank_name,
        isotope=args.isotopes,
        load_wavenum_min=wmin,
        load_wavenum_max=wmax,
        cache=True,
        return_local_path=True,
        parse_quanta=False,
        parallel=True,
        verbose=True,
    )
    _, local_path = result

    meta = {
        "database_edition": "HITRAN2024",
        "molecule": "H2O",
        "isotopes": args.isotopes,
        "wavenumber_range_cm-1": [wmin, wmax],
        "databank_name": args.databank_name,
        "local_path": str(local_path),
        "radis_version": radis.__version__,
        "downloaded_utc": datetime.now(timezone.utc).isoformat(),
    }
    (out / "download_metadata.json").write_text(json.dumps(meta, indent=2))
    print(f'\nRegistered databank: {args.databank_name}')
    print(f'Use: sf.load_databank("{args.databank_name}")')


if __name__ == "__main__":
    main()

