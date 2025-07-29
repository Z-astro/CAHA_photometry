#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Fri Jul 25 14:14 2025
@author: Hao Zhang
[MODIFIED FOR PYTHON 3.13 and Photutils 2.2.0]

The method for background subtraction is to employ a 2D background map for the entire image using photutils.background.Backgroud2D.
This approach models the background across the whole frame and subtracts this model from the data before performing photometry. This is effective for images with large-scale background variations.

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
from astropy.stats import SigmaClip
from photutils.psf import fit_fwhm
from photutils.background import Background2D, MedianBackground
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
        
def measure_fwhm(images, positions, fit_shape):
    """ 增加背景扣除的核心修改 """
    jdlst, fwhm_tg, fwhm_cp = [], [], []
    
    # Background modeling parameters
    BKG_BOX_SIZE = (30, 30)  # 根据图像中星体密度调整
    BKG_FILTER_SIZE = (5, 5)  # 背景平滑滤波器
    BKG_SIGMA_CLIP = 3.0       # Sigma裁剪阈值
    
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
            # Background modeling and subtraction
            # =============================================
            # Step 1: estimate 2-D background
            sigma_clip = SigmaClip(sigma=BKG_SIGMA_CLIP)
            bkg_estimator = MedianBackground()
            
            bkg = Background2D(data,
                             box_size=BKG_BOX_SIZE,
                             filter_size=BKG_FILTER_SIZE,
                             sigma_clip=sigma_clip,
                             bkg_estimator=bkg_estimator)
            
            # Step 2: subtract background
            data_bkg_sub = data - bkg.background
            
            # Step 3: check for strong negative residuals
            if np.any(data_bkg_sub < -3 * bkg.background_rms):
                warnings.warn(f"Large negative residuals detected in {img_path}. Over-subtraction may have occurred.")
            
            # =============================================
            # measure FWHM using subtracted data
            # =============================================
            pos_obj = positions[-2]
            pos_cmp = positions[-1]
            pos_tobe_fitted = [pos_obj, pos_cmp]
            
            # fit_fwhm expects background-subtracted data
            fwhms = fit_fwhm(data=data_bkg_sub, 
                            xypos=pos_tobe_fitted, 
                            fit_shape=fit_shape)
            
            fwhm_tg.append(fwhms[0])
            fwhm_cp.append(fwhms[1])
            print(f"Target FWHM: {fwhms[0]:.2f}, Comparison FWHM: {fwhms[1]:.2f}")
            
        except Exception as e:
            print(f"处理 {img_path} 时发生错误: {str(e)}")
            fwhm_tg.append(np.nan)
            fwhm_cp.append(np.nan)
    
    return jdlst, fwhm_tg, fwhm_cp
               
def write_fwhm_log(objnm, jdlst, fwhm_tg, fwhm_cp):
    '''
    Write the FWHM results to a log file.
    '''
    os.makedirs('../lightcurve', exist_ok=True)
    filenm = '../lightcurve/{}_fwhm_photutils_bkg2D.dat'.format(objnm)
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
    fwhm_picnm = '../lightcurve/{}_fwhm_photutils_bkg2D.png'.format(objnm)
    fwhm_pdfnm = '../lightcurve/{}_fwhm_photutils_bkg2D.pdf'.format(objnm)
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
    
    jdlst, fwhm_tg, fwhm_cp = measure_fwhm(images, positions, fit_shape)
    write_fwhm_log(objnm, jdlst, fwhm_tg, fwhm_cp)
    plot_fwhm(objnm, jdlst, fwhm_tg, fwhm_cp)

    # jdlst, fwhm_tg, fwhm_cp = measure_fwhm(images, positions, fit_shape, r_in, r_out)
    # write_fwhm_log(objnm, jdlst, fwhm_tg, fwhm_cp)
    # plot_fwhm(objnm, jdlst, fwhm_tg, fwhm_cp)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Total processing time: {elapsed_time:.2f} seconds")
    
if __name__ == '__main__':
    main()