#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Fri Jul 25 14:14 2025
@author: Hao Zhang
[MODIFIED FOR PYTHON 3.13 and Photutils 2.2.0]

The method for background subtraction is a more localized backgroud subtraction technique. For each star, it defines an annular region (a ring) around the star's aperture. The median pixel value within this annulus is then calculated using 
astropy.stats.sigma_clipped_stats. This median value, considered the local sky background for that specific star, is then subtracted. This method is particularly robust when the background varies significantly across the field, 
as it determines the background in close proximity to each object of interest.
run in the ./alipy_out directory.
"""

import os
import time
import glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from astropy.time import Time
from astropy.stats import sigma_clipped_stats
from photutils.aperture import CircularAnnulus
from photutils.psf import fit_fwhm
from astropy.io.fits.verify import VerifyWarning
import warnings

warnings.simplefilter('ignore', category=VerifyWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)

# --- Matplotlib settings remain the same ---
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rc('text', usetex=False) # Changed to False to avoid LaTeX dependency issues, set to True if you have LaTeX installed
plt.rc('font', family='serif')
mpl.rcParams['font.sans-serif'] = ['serif']
mpl.rcParams['font.size'] = 16
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.left'] = True
mpl.rcParams['ytick.right'] = True

def rdr_aper():
    '''
    read measure apertures, r, r_in, r_out
    aper file name: 'aper.txt'
    format: r, r_in, r_out
    '''
    aperlst = open('../aper.txt').readline().split()
    aperlst = [float(i) for i in aperlst]
    return tuple(aperlst)


def get_jd(fn):
    '''
    [COMPLETELY REWRITTEN]
    Gets the Julian Date from a FITS file header.
    It first tries to read the 'JD' keyword.
    If that fails, it calculates the JD from 'DATE-OBS' and 'EXPTIME'
    using astropy.time, which replaces the deprecated pyraf dependency.
    '''
    try:
        with pyfits.open(fn) as hdu:
            # First, try to get the JD directly
            val = hdu[0].header['JD']
            return float(val)
    except (KeyError, TypeError):
        # If 'JD' keyword doesn't exist or is invalid, calculate it
        try:
            with pyfits.open(fn) as hdu:
                header = hdu[0].header
                # Use astropy.time.Time to parse the observation date
                # It can handle most standard FITS date formats
                obs_time = Time(header['DATE-OBS'], format='isot', scale='utc')
                
                # Add half of the exposure time to get the mid-point of the observation
                if 'EXPTIME' in header:
                    exposure_time_days = float(header['EXPTIME']) / (24.0 * 3600.0)
                    mid_obs_time = obs_time.jd + exposure_time_days / 2.0
                    return mid_obs_time
                else:
                    # If no exposure time, return JD at the start of observation
                    return obs_time.jd
        except Exception as e:
            print(f"Error calculating JD for {fn}: {e}")
            print("Could not find 'JD' or 'DATE-OBS' in the header.")
            return None
    
def measure_fwhm(images, positions, fit_shape, r_in, r_out):
    """
    Background subtraction now uses a local annulus for each star,
    consistent with the method in `ca_photutils_noFWHM.py`.
    """
    jdlst, fwhm_tg, fwhm_cp = [], [], []
    
    for img_path in images:
        print(f"Processing image: {img_path}")
        
        # get Julian Date
        jd = get_jd(img_path)
        if jd is None:
            print(f"Skipping {img_path} due to JD calculation failure")
            continue
        jdlst.append(jd)
        
        try:
            # read data
            with pyfits.open(img_path) as hdu:
                data = hdu[0].data.astype(float)
            
            # =============================================
            # Local background subtraction for each star
            # =============================================
            
            pos_obj = positions[-2]
            pos_cmp = positions[-1]
            pos_tobe_fitted = [pos_obj, pos_cmp]
            
            # Create a data copy to subtract background from
            data_bkg_sub = np.copy(data)

            # Define the annulus for background estimation
            annulus_apertures = CircularAnnulus(pos_tobe_fitted, r_in=r_in, r_out=r_out)
            
            # Create masks from the annuli
            annulus_masks = annulus_apertures.to_mask(method='center')
            
            bkg_medians = []
            for mask in annulus_masks:
                annulus_data = mask.multiply(data)
                annulus_data_1d = annulus_data[mask.data > 0]
                _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
                bkg_medians.append(median_sigclip)

            # Subtract the local background from the entire image
            # Note: A more refined approach would be to only subtract the background
            # in the vicinity of each star, but for FWHM fitting, subtracting a constant
            # local background value from the entire array is sufficient.
            # We will perform the background subtraction on the data fed to fit_fwhm.
            
            fwhms = []
            
            # Target star FWHM measurement
            data_sub_tg = data - bkg_medians[0]
            fwhm_tg_val = fit_fwhm(data=data_sub_tg, xypos=[pos_obj], fit_shape=fit_shape)
            fwhms.append(fwhm_tg_val[0])
            
            # Comparison star FWHM measurement
            data_sub_cp = data - bkg_medians[1]
            fwhm_cp_val = fit_fwhm(data=data_sub_cp, xypos=[pos_cmp], fit_shape=fit_shape)
            fwhms.append(fwhm_cp_val[0])

            fwhm_tg.append(fwhms[0])
            fwhm_cp.append(fwhms[1])
            print(f"Target FWHM: {fwhms[0]:.2f}, Comparison FWHM: {fwhms[1]:.2f}")
            
        except Exception as e:
            print(f"An error occurred while processing {img_path}: {str(e)}")
            fwhm_tg.append(np.nan)
            fwhm_cp.append(np.nan)
    
    return jdlst, fwhm_tg, fwhm_cp

               
def write_fwhm_log(objnm, jdlst, fwhm_tg, fwhm_cp):
    '''
    Write the FWHM results to a log file.
    '''
    os.makedirs('../lightcurve', exist_ok=True)
    filenm = '../lightcurve/{}_fwhm_photutils.dat'.format(objnm)
    print(f"Writing FWHM results to {filenm}")
    
    with open(filenm, 'w') as f:
        header = '# JD          obj_fwhm    cmp_fwhm\n'
        f.write(header)
        for jd, tg_fwhm, cp_fwhm in zip(jdlst, fwhm_tg, fwhm_cp):
            line = '%.8f  %.8f      %.8f\n' % (jd, tg_fwhm, cp_fwhm)
            f.write(line)
            
            
def plot_fwhm(objnm, jdlst, fwhm_tg, fwhm_cp):
    '''
    Plot the FWHM results.
    '''
    
    PIXEL_SCALE = 0.53 # arcsec/pixel, CAFOS 2.2m telescope
    
    def align_secondary_ticks(ax, scale_factor):
        """ 强制同步双轴刻度位置 """
        # 获取主轴的刻度位置
        yticks = ax.get_yticks()
        # 设置右侧轴的刻度位置完全对应左侧轴的位置
        sec_ax = ax.secondary_yaxis(
            'right', 
            functions=(lambda x: x * scale_factor, lambda x: x / scale_factor)
        )
        sec_ax.set_yticks(yticks * scale_factor)  # 强制刻度位置映射
        return sec_ax
    
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(16, 10))
    fig.suptitle(f'{objnm} FWHM Variation')
    fig.subplots_adjust(hspace=0.05)
    
    axes[0].scatter(jdlst, fwhm_tg, marker='o', label=f'{objnm}')
    mean_fwhm_tg = np.nanmean(fwhm_tg)
    std_fwhm_tg = np.nanstd(fwhm_tg)
    axes[0].axhline(mean_fwhm_tg, color='red', linestyle='--')
    axes[0].axhline(mean_fwhm_tg + std_fwhm_tg, color='red', linestyle=':')
    axes[0].axhline(mean_fwhm_tg - std_fwhm_tg, color='red', linestyle=':')
    axes[0].text(0.05, 0.95, 
                 rf'{mean_fwhm_tg:.2f}$\pm${std_fwhm_tg:.2f} pix,{mean_fwhm_tg*PIXEL_SCALE:.2f}$\pm${std_fwhm_tg*PIXEL_SCALE:.2f} arcsec',
                 transform=axes[0].transAxes,
                 verticalalignment='top', fontsize=16)
    axes[0].text(0.95, 0.95, f'Object: {objnm}', transform=axes[0].transAxes,
                 horizontalalignment='right', verticalalignment='top', fontsize=20, color='black', fontweight='bold')
    axes[0].set_ylabel('FWHM (pixels)')
    # sec_ax0 = axes[0].secondary_yaxis('right', functions=(lambda x: x * PIXEL_SCALE, lambda x: x / PIXEL_SCALE))
    # sec_ax0.set_ylabel('Seeing (arcsec)')
    sec_ax0 = align_secondary_ticks(axes[0], PIXEL_SCALE)
    sec_ax0.set_ylabel('Seeing (arcsec)')
    sec_ax0.tick_params(axis='y')
    
    axes[1].scatter(jdlst, fwhm_cp, marker='o', label='Comparison Star')
    mean_fwhm_cp = np.nanmean(fwhm_cp)
    std_fwhm_cp = np.nanstd(fwhm_cp)
    axes[1].axhline(mean_fwhm_cp, color='red', linestyle='--')
    axes[1].axhline(mean_fwhm_cp + std_fwhm_cp, color='red', linestyle=':')
    axes[1].axhline(mean_fwhm_cp - std_fwhm_cp, color='red', linestyle=':')
    axes[1].text(0.05, 0.95, 
                 rf'{mean_fwhm_cp:.2f}$\pm${std_fwhm_cp:.2f} pix,{mean_fwhm_cp*PIXEL_SCALE:.2f}$\pm${std_fwhm_cp*PIXEL_SCALE:.2f} arcsec',
                 transform=axes[1].transAxes,
                 verticalalignment='top', fontsize=16)
    axes[1].text(0.95, 0.95, 'Spec.Comp', transform=axes[1].transAxes,
                 horizontalalignment='right', verticalalignment='top', fontsize=20, color='black', fontweight='bold')
    
    axes[1].set_xlabel('Julian Date (JD)')
    axes[1].set_ylabel('FWHM (pixels)')
    # sec_ax1 = axes[1].secondary_yaxis('right', functions=(lambda x: x * PIXEL_SCALE, lambda x: x / PIXEL_SCALE))
    # sec_ax1.set_ylabel('Seeing (arcsec)')
    sec_ax1 = align_secondary_ticks(axes[1], PIXEL_SCALE)
    sec_ax1.set_ylabel('Seeing (arcsec)')
    sec_ax1.tick_params(axis='y')
    
    os.makedirs('../lightcurve', exist_ok=True)
    fwhm_picnm = '../lightcurve/{}_fwhm_photutils.png'.format(objnm)
    fwhm_pdfnm = '../lightcurve/{}_fwhm_photutils.pdf'.format(objnm)
    fig.savefig(fwhm_picnm)
    fig.savefig(fwhm_pdfnm)
    print('save FWHM fig --> {} & {}'.format(fwhm_picnm, fwhm_pdfnm))
    
def main():
    start_time = time.time()

    fname = np.loadtxt('objlst', dtype=str)
    fname.sort()
    images = [i.strip() for i in fname]
    
    posname = glob.glob('../*.pos')
    posname = str(posname[0])
    positions = np.loadtxt(posname)
    
    absdir = os.getcwd() # last directory,not the parent directory
    objpath, objnm = os.path.split(absdir)
    objnm = os.path.basename(objpath)
    objnm = objnm.replace('obj_', '')
    
    r, r_in, r_out = rdr_aper()
    fit_shape = int(2 * r_in + 1)
    
    # jdlst, fwhm_tg, fwhm_cp = measure_fwhm(images, positions, fit_shape)
    # write_fwhm_log(objnm, jdlst, fwhm_tg, fwhm_cp)
    # plot_fwhm(objnm, jdlst, fwhm_tg, fwhm_cp)

    jdlst, fwhm_tg, fwhm_cp = measure_fwhm(images, positions, fit_shape, r_in, r_out)
    write_fwhm_log(objnm, jdlst, fwhm_tg, fwhm_cp)
    plot_fwhm(objnm, jdlst, fwhm_tg, fwhm_cp)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Total processing time: {elapsed_time:.2f} seconds")
    
if __name__ == '__main__':
    main()