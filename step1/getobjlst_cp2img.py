#!/usr/bin/env python
# coding=utf-8


import os
import sys
import glob
import shutil
from astropy.io import fits as pyfits
from astropy.io.fits.verify import VerifyWarning
import warnings
warnings.simplefilter('ignore', category=VerifyWarning)


def main():
    """
    由于检查fits时已经根据'objradec.lst'匹配过坐标，这一步默认坐标正确
    这一步只需按源名称将文件分开
    """
    object_dir = sys.argv[1]

    cw_dir = os.path.basename(os.getcwd())  # 获取当前目录，判断源名后加什么后缀（V/B/R）
    if cw_dir[-1] == 'B':
        filt_suffix = 'B'
    elif cw_dir[-1] == 'R':
        filt_suffix = 'R'
    else:
        filt_suffix = ''

    fitsnmlst = glob.glob('ftb*.fits')
    for fitsnumb, fitsnm in enumerate(fitsnmlst):
        objnm = pyfits.getval(fitsnm, 'OBJECT').split()[0] + filt_suffix
        print(objnm)
        objsrc_dir = os.sep.join((object_dir, 'obj_' + objnm))
        if not os.path.exists(objsrc_dir):
            os.makedirs(objsrc_dir)
            print('==> mkdir {}'.format(objsrc_dir))
        else:
            pass
        objsrcimg_dir = os.sep.join((objsrc_dir, 'images'))
        if not os.path.exists(objsrcimg_dir):
            os.makedirs(objsrcimg_dir)
            print('==> mkdir {}'.format(objsrcimg_dir))
        else:
            pass
        # print(objsrcimg_dir)
        shutil.copy2(fitsnm, objsrcimg_dir)
        print('> copy {} --> {}'.format(fitsnm, objsrcimg_dir))


if __name__ == '__main__':
    main()