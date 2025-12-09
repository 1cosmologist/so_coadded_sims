# Simons Observatory Coadded Simulations

This repository contains scripts and resources for producing coadded CMB + noise + foreground maps for the Simons Observatory (SO).

## Overview

The coadded maps are produced by combining:
1. **CMB signal** - Lensed CMB alms from SO v4 simulations
2. **Instrumental noise** - Map-based noise realizations
3. **Foregrounds** - Gaussian foreground simulations

## Pipeline Steps

### Step 1: Read CMB alms
- Read lensed CMB alms (T, E, B) from SO v4 simulations
- Source: `/global/cfs/cdirs/sobs/v4_sims/mbs/cmb/fullskyLensedUnabberatedCMB_alm_set00_XXXXX.fits`
- HDU 1, 2, 3 correspond to T, E, B alms respectively

### Step 2: Beam Smoothing
- Apply Gaussian beam smoothing using `healpy.almxfl()`
- Beam FWHM values are specified per channel in the instrument parameter file
- Temperature uses spin-0 beam, polarization (E, B) uses spin-2 beam

### Step 3: Convert to Maps
- Convert smoothed alms to HEALPix maps at NSIDE=512 using `healpy.alm2map()`
- Output: T, Q, U Stokes maps

### Step 4: Add Noise Maps
- Read pre-computed noise realizations per channel
- Source: `/pscratch/sd/s/shamikg/so_mapbased_noise/output/{config}/sobs_noise_{channel}_mcXXX_nside0512.fits`

### Step 5: Add Foreground Maps
- Read Gaussian foreground simulations per channel
- Source: `/pscratch/sd/s/shamikg/so_gaussian_fg/output/foreground_sims/gaussian_fg/sobs_gaussfg_{channel}_mcXXX_nside0512.fits`

### Step 6: Apply Mask
- Multiply coadded map by binary mask
- Mask: `/pscratch/sd/s/shamikg/so_mapbased_noise/resources/so_sat_full-binary_C_nside512.fits`

### Step 7: Save Output
- Output path: `/pscratch/sd/s/shamikg/so_coadded_sims/output/AL{Alens}_{config}/`
- Filename format: `sobs_coadd_{channel}_AL{Alens}_mcXXX_nside0512.fits`

## Configuration Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `Alens` | 1.0 | Lensing amplitude |
| `NSIMS` | 100 | Number of simulations (0-99) |
| `NSIDE` | 512 | HEALPix resolution parameter |

## Instrument Parameters (Baseline Pessimistic)

| Channel | Frequency (GHz) | Beam FWHM (arcmin) | Noise (μK-arcmin) | ℓ_knee | α_knee |
|---------|-----------------|--------------------|--------------------|--------|--------|
| LF027 | 27 | 91.0 | 16.0 | 30 | -2.4 |
| LF039 | 39 | 63.0 | 10.0 | 30 | -2.4 |
| MF093 | 93 | 30.0 | 1.7 | 50 | -2.5 |
| MF145 | 145 | 17.0 | 2.1 | 50 | -3.0 |
| HF225 | 225 | 11.0 | 5.9 | 70 | -3.0 |
| HF280 | 280 | 9.0 | 15.0 | 100 | -3.0 |

## Directory Structure

```
so_coadded_sims/
├── README.md
├── code/
│   ├── coadd_cmb_noise_fg.py    # Main coadding script
│   ├── read_cmb_spectra.py       # CMB spectra utilities
│   └── validate_cmb.ipynb        # Validation notebook
├── output/
│   └── AL1.00_baseline_pessimistic/
│       └── sobs_coadd_{channel}_AL1.00_mcXXX_nside0512.fits
└── resources/
    └── instr_params_baseline_pessimistic.yaml
```

## Usage

```bash
cd code
python coadd_cmb_noise_fg.py
```

## Output FITS Header

Each output FITS file contains the following header keywords:
- `SIMIDX`: Simulation index (0-99)
- `CHANNEL`: Frequency channel name
- `FWHM_ARCMIN`: Beam FWHM in arcminutes
- `ALENS`: Lensing amplitude

## Dependencies

- `numpy`
- `healpy`
- `yaml`
- `tqdm`

## References

- Simons Observatory v4 simulations
- SO SAT baseline pessimistic noise model
