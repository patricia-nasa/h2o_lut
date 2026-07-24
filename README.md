# HITRAN2024 + MT_CKD H₂O absorption LUT for PyVLIDORT

This package creates a lookup table and uses it to populate:

```python
alpha[nlev, 1, nobs]
```

in `sbg_vlidort.py`.

## Design

Offline:

```text
HITRAN2024 local lines + MT_CKD 4.3 continuum
        -> k_total(channel, pressure, temperature, H2O mole fraction)
        -> NetCDF LUT
```

Runtime:

```text
GEOS pe, te, ze, QV
        -> layer p, T, xH2O, thickness
        -> interpolate k_total
        -> alpha = k_total * thickness_cm
```

The LUT stores absorption coefficient in cm⁻¹ rather than optical depth, so it
works with any vertical grid.

## 1. Install

```bash
conda env create -f environment.yml
conda activate h2o-lut
```

## 2. Download/cache HITRAN2024

```bash
python download_hitran2024.py \
  --output-dir /discover/nobackup/$USER/HITRAN2024/H2O \
  --databank-name HITRAN2024-H2O
```

RADIS writes the databank-name/path mapping into its user configuration. Later:

```python
sf.load_databank("HITRAN2024-H2O")
```

uses that mapping.

## 3. Install MT_CKD 4.3

```bash
bash install_mtckd.sh /discover/nobackup/$USER/MT_CKD_H2O
```

Follow the release build README for the exact compiler target. Set the built
executable and `absco-ref_wv-mt-ckd.nc` paths in `config_example.yaml`.

MT_CKD H₂O supplies self and foreign continuum coefficients from 0 to
20,000 cm⁻¹, corresponding to wavelengths ≥500 nm. Above 20,000 cm⁻¹ the
generator records zero H₂O continuum, while HITRAN line absorption remains.

## 4. Configure

Copy:

```bash
cp config_example.yaml config.yaml
```

Update:

- AMES BRDF channel glob;
- output file;
- MT_CKD executable and coefficient paths;
- number of workers.

The generator reads the same `nwav` channel centers used by `sbg_vlidort.py`.

## 5. Generate channel parts

For a small test:

```bash
python generate_lut.py config.yaml --channel-start 0 --channel-stop 2
```

For an HPC job array, split channel ranges among jobs. Each channel writes an
independent restartable file under `<output_stem>_parts/`.

Then merge:

```bash
python merge_lut.py config.yaml
```

## 6. Validate

```bash
python validate_lut.py config.yaml --samples 100
```

The validator compares interpolation with direct RADIS line calculations at
random atmospheric states. A failure means the p/T/H₂O axes need refinement.

Continuum validation should additionally compare selected states against the
official MT_CKD output and, ideally, an LBLRTM reference.

## 7. Integrate with PyVLIDORT

See `sbg_vlidort_integration.txt`. Runtime use is:

```python
alpha = self.h2o_lut.alpha(
    channel_nm=self.channels[ich],
    pe_pa=pe,
    te_k=te,
    ze_m=ze,
    qv=qv,
)
```

## Spectral grid recommendation

The 300–2500 nm interval corresponds to 33,333.33–4,000 cm⁻¹. Do not build one
uniform high-resolution spectrum over that entire interval for every state.

For monochromatic SBG channel centers:

- initial RADIS spacing: 0.001 cm⁻¹;
- reference convergence spacing: 0.0005 cm⁻¹;
- line truncation: 25 cm⁻¹ when adding MT_CKD.

For a finite spectral response function, the LUT should ultimately store the
SRF-effective optical depth rather than a channel-center value. That requires
the actual SRF and the nonlinear average of transmission. This package currently
implements channel-center absorption because that is what the present driver
exposes.

## Isotopologues

Default: HITRAN local isotopologues 1–7.

This is conservative for terrestrial natural abundance. Isotopologues 1–4 are
the minimum recommended subset, but dropping 5–7 must be demonstrated to remain
below the 0.1% channel-level threshold.

## Meaning of the 0.1% target

The configuration is an informed starting grid, not a guarantee. Validate:

1. 0.001 versus 0.0005 cm⁻¹ line sampling;
2. LUT interpolation versus direct calculations;
3. isotopologues 1–4 versus 1–7;
4. continuum against official MT_CKD output;
5. final channel transmittance after applying the instrument SRF.

Spectroscopic parameter uncertainties may exceed 0.1% in some bands even when
the numerical implementation converges more tightly than 0.1%.

