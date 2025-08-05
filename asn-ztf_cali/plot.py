#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from astropy.stats import sigma_clipped_stats
from astropy.time import Time


plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
mpl.rcParams['font.sans-serif'] = ['serif']
mpl.rcParams['font.size'] = 20
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.left'] = True
mpl.rcParams['ytick.right'] = True


def calmeanstd(mag, err):
    meanmag = sum(mag) / float(len(mag))
    sigma = 0.0
    if len(mag) > 1:
        num = float(len(mag))
        for i in range(0, len(mag)):
            sigma += 1.0 / (num - 1.0) * (mag[i] - meanmag)**2
        sigma = sigma**0.5
    else:
        sigma = err[0]
    err_mean = 0.0
    for i in range(0, len(err)):
        err_mean += err[i]**2
    err_mean = err_mean**0.5 / float(len(err))
    return meanmag, sigma, err_mean


def combine_days(jd, diffmag, differr, nbin):
    data_dict = {'jd': jd, 'mag': diffmag, 'err': differr}
    combine_df = pd.DataFrame(data_dict)
    combine_df['cut_group'] = pd.cut(combine_df.jd, bins=np.arange(np.min(jd), np.max(jd), nbin), right=False)
    grouped = combine_df.groupby(by='cut_group')
    new_jd = []
    new_jderr = []
    new_diffmag = []
    new_differr = []
    for tempname, tempgroup in grouped: 
        if len(tempgroup) > 0:
            new_jd.append(np.sum(tempgroup["jd"].values) / float(len(tempgroup)))
            new_jderr.append((np.max(tempgroup["jd"].values) - np.min(tempgroup["jd"].values)) / 2.0)
            tmpdiffmag, tmpsigma, tmperrmean = calmeanstd(tempgroup["mag"].values, tempgroup["err"].values)
            new_diffmag.append(tmpdiffmag)
            new_differr.append(tmpsigma / len(tempgroup)**0.5)
    return np.array(new_jd), np.array(new_jderr), np.array(new_diffmag), np.array(new_differr)


def plot(objnm):
    # caliDir = '/home/yao/astroWORK/CAHA_LC/cali/'
    # caliDir = f'../../asn-ztf/obj_{objnm}/'
    # figDir = '/home/yao/astroWORK/CAHA_LC/figure/'
    # txtnm = os.sep.join((caliDir, objnm + '.txt_cali'))
    txtnm = f"../../asn-ztf/obj_{objnm}/{objnm}.txt_cali"
    # figpngnm = os.sep.join((figDir, objnm + '.png'))
    figpngnm = f"../../asn-ztf/obj_{objnm}/{objnm}.png"
    # figpdfnm = os.sep.join((figDir, objnm + '.pdf'))
    # figpdfnm = f"../../asn-ztf/obj_{objnm}/{objnm}.pdf"

    alltypes = ['ASAS-SN-g', 'ASAS-SN-V', 'ZTF-r', 'ZTF-g', 'ZTF-i']
    colordic = {'ASAS-SN-g': 'k', 'ASAS-SN-V': 'blue', 'ZTF-r': 'red', 'ZTF-g': 'green', 'ZTF-i': 'orange'}
    daybindic = {'ASAS-SN-g': 10, 'ASAS-SN-V': 10, 'ZTF-r': 3, 'ZTF-g': 3, 'ZTF-i': 3}
    nstddic = {'ASAS-SN-g': 2, 'ASAS-SN-V': 2, 'ZTF-r': 100, 'ZTF-g': 100, 'ZTF-i': 100}

    dat = np.loadtxt(txtnm, dtype=str)
    # jd, flux,
    jdo = np.array(dat[:, 0], dtype=np.float64)
    fluxo = np.array(dat[:, 1], dtype=np.float64)
    fluxerro = np.array(dat[:, 2], dtype=np.float64)
    lctype = dat[:, 3]
    lctypeuniqs = np.unique(lctype)
    lctypeuniqs = sorted(lctypeuniqs)

    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(16, 8))
    fig.subplots_adjust(hspace=0)
    ax[0].invert_yaxis()
    ax[1].invert_yaxis()

    datflag = np.array([[], [], [], []])
    for i in range(len(lctypeuniqs)):
        print(lctypeuniqs[i])
        arg = np.where(lctype == lctypeuniqs[i])
        jd, flux, fluxerr = jdo[arg], fluxo[arg], fluxerro[arg]
        mean_flux, median_flux, std_flux = sigma_clipped_stats(flux)
        arg1 = np.where(np.abs(flux - median_flux) < 4 * std_flux)
        jd, flux, fluxerr = jd[arg1], flux[arg1], fluxerr[arg1]
        mag = -2.5 * np.log10(flux * 10**(-25)) - 48.6
        magerr = fluxerr / flux / np.log(10) * 2.5
        serr = 0.1       

        jd, jderr, mag, magerr = combine_days(jd, mag, magerr, daybindic[lctypeuniqs[i]])
        arg2 = np.where(np.array(magerr) < serr)
        jd, jderr, mag, magerr = jd[arg2], jderr[arg2], mag[arg2], magerr[arg2]
        mean_mag, median_mag, std_mag = sigma_clipped_stats(mag)
        arg3 = np.where(np.abs(mag - median_mag) < nstddic[lctypeuniqs[i]] * std_mag)
        jd, jderr, mag, magerr = jd[arg3], jderr[arg3], mag[arg3], magerr[arg3]
        ax[0].errorbar(jd, mag, xerr=jderr, yerr=magerr, fmt='.', label=lctypeuniqs[i], elinewidth=0.2, ms=6, color=colordic[lctypeuniqs[i]])
        ax[0].legend(loc=3, ncol=4, framealpha=0.3, frameon=False)
        dat = np.array([jd, jderr, mag, magerr])
        # print(dat.shape)
        datflag = np.concatenate((datflag, dat), axis=1)
        # print(datflag.shape)
    [jd, jderr, mag, magerr] = datflag
    print(jdo.shape[0])
    print(jd.shape[0])

    ax[0].legend(loc=2, ncol=4, frameon=False, fontsize=14)
    ax[0].set_title(objnm, pad=10)
    ax[0].tick_params(labelbottom=False)                         # 设置不显示底部x轴刻度值
    ax[0].set_ylabel(r'$\rm{Mag}$')
    ax[0].minorticks_on()

    # x轴上端显示年份
    axt = ax[0].twiny()
    UT = Time(jd + 50000, format='mjd').jyear
    axt.scatter(UT, mag, marker='o', s=20, c='goldenrod', alpha=0.001)
    axt.minorticks_on()

    # 数据点再做一次bin画在图2，方便看光变
    jd, jderr, mag, magerr = combine_days(jd, mag, magerr, 10)
    print(jd.shape[0])
    ax[1].errorbar(jd, mag, xerr=jderr, yerr=magerr, fmt='.', elinewidth=0.2, ms=6, color='k')
    ax[1].set_xlabel('MJD - 50000 (days)')
    ax[1].set_ylabel(r'$\rm{Mag}$')
    ax[1].minorticks_on()
    fig.text(0.79, 0.45, '10 days binned')

    # bintxtnm = '%s_bin10days.lc' % objnm
    bintxtnm = f"../../asn-ztf/obj_{objnm}/{objnm}_bin10days.lc"
    bintxt = open(bintxtnm, 'wt')
    for i in range(len(jd)):
        rowfmt = '%.6f  %.6e  %.6e\n' % (jd[i] + 2450000.5, mag[i], magerr[i])
        bintxt.write(rowfmt)
    bintxt.close()
    
    fig.savefig(figpngnm, format='png', dpi=400, bbox_inches='tight', pad_inches=0.05)
    # plt.show()
    plt.close()


def main():
    objnm = sys.argv[1]
    # objnm = 'NGC7603'
    plot(objnm)


if __name__ == "__main__":
    main()
