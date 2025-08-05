#! /usr/bin/env python
# -*- coding: utf-8 -*-


import sys
import numpy as np


def rd_caliResult(objnm):
    # txtdir = '/home/yao/astroWORK/CAHA_LC/cali/' + objnm + '.txt_cali'
    txtdir = f"../../asn-ztf/obj_{objnm}/{objnm}.txt_cali"
    recali = np.loadtxt(txtdir, dtype=str)
    jd = np.array(recali[:, 0], dtype=np.float64)
    flux = np.array(recali[:, 1], dtype=np.float64)
    fluxerr = np.array(recali[:, 2], dtype=np.float64)
    uniques = recali[:, 3]
    mag = -2.5 * np.log10(flux * 10**(-25.0)) - 48.6
    magerr = fluxerr / flux / np.log(10) * 2.5
    return jd, mag, magerr, uniques


def wt_caliResult(objnm):
    jd, mag, magerr, uniques = rd_caliResult(objnm)
    # txtdir = '/home/yao/astroWORK/CAHA_LC/LC/' + objnm + '.lc'
    txtdir = f"../../asn-ztf/obj_{objnm}/{objnm}.lc"
    print('Write to %s' % txtdir)
    txt = open(txtdir, 'wt')
    for i in range(jd.shape[0]):
        rowformat = '%.6f  %.6e  %.6e  %s\n' % (jd[i] + 2450000.5, mag[i], magerr[i], uniques[i])
        txt.write(rowformat)
    txt.close()


def main():
    objnm = sys.argv[1]
    # objnm = 'IZw1'
    wt_caliResult(objnm)


if __name__ == "__main__":
    main()
