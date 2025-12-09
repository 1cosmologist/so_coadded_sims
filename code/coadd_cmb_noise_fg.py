"""
Script to coadd CMB, noise, and foreground maps for Simons Observatory simulations.

This script:
1. Reads CMB alms from SO simulations
2. Smooths CMB alms with Gaussian beam using instrument parameters
3. Converts to maps at NSIDE=512
4. Adds noise maps and Gaussian foreground maps
5. Applies the full binary mask
"""

import numpy as np
import healpy as hp
import yaml
import os
from pathlib import Path
from tqdm import tqdm

# =============================================================================
# ROOT PATH CONFIGURATION
# =============================================================================
Alens = 0.3
NSIMS = 100  # Number of simulations
fgcomplexity = "gaussian"  # Foreground complexity
# Root path for CMB alms
CMB_ALM_ROOT = "/global/cfs/cdirs/sobs/v4_sims/mbs/cmb"

# Root path for noise maps
NOISE_MAP_ROOT = "/pscratch/sd/s/shamikg/so_mapbased_noise/output"

# Root path for foreground maps
FOREGROUND_MAP_ROOT = f"/pscratch/sd/s/shamikg/so_gaussian_fg/output/foreground_sims/{fgcomplexity}_fg"

# Root path for mask files
MASK_ROOT = "/pscratch/sd/s/shamikg/so_mapbased_noise/resources"

# Output directory for coadded maps
OUTPUT_DIR = "/pscratch/sd/s/shamikg/so_coadded_sims/output"

# =============================================================================


def load_instrument_params(yaml_path):
    """Load instrument parameters from YAML file."""
    with open(yaml_path, 'r') as f:
        params = yaml.safe_load(f)
    return params


def get_gaussian_beam_fl(fwhm_arcmin, lmax):
    """
    Get Gaussian beam window function.
    
    Parameters
    ----------
    fwhm_arcmin : float
        Full-width at half-maximum of the beam in arcminutes
    lmax : int
        Maximum multipole
        
    Returns
    -------
    bl : array
        Beam window function
    """
    fwhm_rad = np.radians(fwhm_arcmin / 60.0)
    bl = hp.gauss_beam(fwhm_rad, lmax=lmax, pol=True)
    return bl


def read_cmb_alms(sim_idx):
    """
    Read CMB alms for a given simulation index.
    
    Parameters
    ----------
    sim_idx : int
        Simulation index (0-99)
        
    Returns
    -------
    alms : list of arrays
        T, E, B alms
    """
    alm_path = os.path.join(CMB_ALM_ROOT, f"fullskyLensedUnabberatedCMB_alm_set00_{sim_idx:05d}.fits")
    
    alms = []
    for hdu in range(1, 4):  # HDU 1, 2, 3 for T, E, B
        alm = hp.read_alm(alm_path, hdu=hdu)
        alms.append(alm)
    
    return alms


def smooth_alms_with_beam(alms, beam_fl):
    """
    Smooth alms with beam window function.
    
    Parameters
    ----------
    alms : list of arrays
        T, E, B alms
    beam_fl : array
        Beam window function (shape: lmax+1, 3 for T, E, B)
        
    Returns
    -------
    smoothed_alms : list of arrays
        Smoothed T, E, B alms
    """
    smoothed_alms = []
    for i, alm in enumerate(alms):
        # beam_fl has columns for T, E, B (indices 0, 1, 2)
        # For polarization: E and B both use the spin-2 beam (column 1)
        if i == 0:  # Temperature
            bl = beam_fl[:, 0]
        else:  # E or B polarization
            bl = beam_fl[:, 1]
        
        smoothed_alm = hp.almxfl(alm, bl)
        smoothed_alms.append(smoothed_alm)
    
    return smoothed_alms


def alms_to_map(alms, nside):
    """
    Convert alms to map.
    
    Parameters
    ----------
    alms : list of arrays
        T, E, B alms
    nside : int
        HEALPix NSIDE parameter
        
    Returns
    -------
    maps : array
        T, Q, U maps (shape: 3, npix)
    """
    alms_array = np.array(alms)
    alms_array[2] *= np.sqrt(Alens)  # Scale B-mode alms by sqrt(Alens)
    maps = hp.alm2map(alms_array, nside, pol=True)
    return maps


def read_noise_map(channel, sim_idx, output_subdir):
    """
    Read noise map for a given channel and simulation index.
    
    Parameters
    ----------
    channel : str
        Channel name (e.g., 'LF027', 'MF093')
    sim_idx : int
        Simulation index (0-99)
    output_subdir : str
        Output subdirectory name (e.g., 'baseline_pessimistic')
        
    Returns
    -------
    noise_map : array
        Noise map (T, Q, U)
    """
    noise_path = os.path.join(NOISE_MAP_ROOT, output_subdir, f"sobs_noise_{channel}_mc{sim_idx:03d}_nside0512.fits")
    noise_map = hp.read_map(noise_path, field=[0, 1, 2])
    return noise_map


def read_foreground_map(channel, sim_idx):
    """
    Read Gaussian foreground map for a given channel and simulation index.
    
    Parameters
    ----------
    channel : str
        Channel name (e.g., 'LF027', 'MF093')
    sim_idx : int
        Simulation index (0-99)
        
    Returns
    -------
    fg_map : array
        Foreground map (T, Q, U)
    """
    fg_path = os.path.join(FOREGROUND_MAP_ROOT, f"sobs_gaussfg_{channel}_mc{sim_idx:03d}_nside0512.fits")
    fg_map = hp.read_map(fg_path, field=[0, 1, 2])
    return fg_map


def read_mask(mask_path):
    """
    Read binary mask.
    
    Parameters
    ----------
    mask_path : str
        Path to mask file
        
    Returns
    -------
    mask : array
        Binary mask
    """
    mask = hp.read_map(mask_path)
    return mask


def coadd_maps(cmb_map, noise_map, fg_map, mask):
    """
    Coadd CMB, noise, and foreground maps and apply mask.
    
    Parameters
    ----------
    cmb_map : array
        CMB map (T, Q, U)
    noise_map : array
        Noise map (T, Q, U)
    fg_map : array
        Foreground map (T, Q, U)
    mask : array
        Binary mask
        
    Returns
    -------
    coadded_map : array
        Coadded and masked map (T, Q, U)
    """
    # coadded_map = cmb_map + noise_map + fg_map
    # Apply mask to each Stokes parameter
    # for i in range(3):
    # coadded_map *= mask
    return ( cmb_map + noise_map + fg_map ) * mask


def process_simulation(sim_idx, instr_params, output_subdir, mask, output_dir):
    """
    Process a single simulation.
    
    Parameters
    ----------
    sim_idx : int
        Simulation index (0-99)
    instr_params : dict
        Instrument parameters
    output_subdir : str
        Output subdirectory name for noise maps
    mask : array
        Binary mask
    output_dir : str
        Output directory for coadded maps
    """
    print(f"Processing simulation {sim_idx:03d}...")
    
    # Read CMB alms
    cmb_alms = read_cmb_alms(sim_idx)
    
    # Get lmax from alms
    lmax = hp.Alm.getlmax(len(cmb_alms[0]))
    
    # Process each channel
    for channel, params in instr_params.items():
        print(f"  Processing channel {channel}...")
        
        # Get beam FWHM and compute beam window function
        fwhm_arcmin = params['beam_fwhm_arcmin']
        nside = params['nside']
        freq = params['central_freq_GHz']
        
        beam_fl = get_gaussian_beam_fl(fwhm_arcmin, lmax)
        
        # Smooth CMB alms with beam
        smoothed_alms = smooth_alms_with_beam(cmb_alms, beam_fl)
        
        # Convert to map
        cmb_map = alms_to_map(smoothed_alms, nside)
        
        # Read noise map
        noise_map = read_noise_map(channel, sim_idx, output_subdir)
        
        # Read foreground map
        fg_map = read_foreground_map(channel, sim_idx)
        
        # Coadd maps and apply mask
        coadded_map = coadd_maps(cmb_map, noise_map, fg_map, mask)
        
        # Save coadded map
        output_path = os.path.join(
            output_dir, 
            f"sobs_coadd_{channel}_AL{Alens:.2f}_mc{sim_idx:03d}_nside{nside:04d}.fits"
        )
        
        header = [
                ('UNITS', 'uK_CMB', 'Map units'),
                ('CHANNEL', channel, 'Channel name'),
                ('FREQ', freq, 'Frequency in GHz'),
                ('BEAM', fwhm_arcmin, 'Beam FWHM in arcmin'),
                ('SIMIDX', sim_idx, 'Simulation index'),
                ('ALENS', Alens, 'Lensing amplitude scaling factor')
            ]
        
        hp.write_map(output_path, coadded_map, overwrite=True, dtype=np.float32, extra_header=header)
        print(f"    Saved: {output_path}")


def main():
    """Main function to run the coadding pipeline."""
    
    # Paths
    yaml_path = "/pscratch/sd/s/shamikg/so_coadded_sims/resources/instr_params_baseline_pessimistic.yaml"
    mask_path = os.path.join(MASK_ROOT, "so_sat_full-binary_C_nside512.fits")
    output_dir = OUTPUT_DIR
    
    # Extract output subdirectory from YAML filename
    # e.g., "instr_params_baseline_pessimistic.yaml" -> "baseline_pessimistic"
    yaml_basename = os.path.basename(yaml_path)
    yaml_name_no_ext = os.path.splitext(yaml_basename)[0]
    # Remove "instr_params_" prefix
    output_subdir = yaml_name_no_ext.replace("instr_params_", "")
    
    # Create output subdirectory with Alens and instrument parameter suffix
    # e.g., "AL1.00_baseline_pessimistic"
    output_subdir_name = f"{fgcomplexity}_AL{Alens:.2f}_{output_subdir}"
    output_dir = os.path.join(OUTPUT_DIR, output_subdir_name)
    
    print(f"Using instrument parameters from: {yaml_path}")
    print(f"Noise map subdirectory: {output_subdir}")
    print(f"Output directory: {output_dir}")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Load instrument parameters
    instr_params = load_instrument_params(yaml_path)
    print(f"Loaded parameters for channels: {list(instr_params.keys())}")
    
    # Read mask
    print(f"Loading mask from: {mask_path}")
    mask = read_mask(mask_path)
    
    # Process all simulations (0-99)
    for sim_idx in tqdm(range(NSIMS), desc="Processing simulations"):
        process_simulation(sim_idx, instr_params, output_subdir, mask, output_dir)
    
    print("Done!")


if __name__ == "__main__":
    main()
