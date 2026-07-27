#!/usr/bin/env bash
set -euo pipefail

DEST="${1:-$PWD/MT_CKD_H2O}"

git clone https://github.com/AER-RC/MT_CKD_H2O.git "$DEST"
cd "$DEST"
git checkout 4.3

# The exact makefile location/name is distributed with the release.
# Inspect README/user guide if your platform uses a different target.
find . -maxdepth 3 -type f \( -name 'make*' -o -name 'Makefile*' \) -print

echo
echo "Build the GNU double-precision executable following the release README."
echo "Then locate:"
find . -type f -name 'absco-ref_wv-mt-ckd.nc' -print

