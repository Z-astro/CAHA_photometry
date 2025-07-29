#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
mpl.rcParams['font.sans-serif'] = ['serif']
mpl.rcParams['font.size'] = 22
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.left'] = True
mpl.rcParams['ytick.right'] = True

# 修改fits header的关键字：OBJECT和IMAGETYP
# 判断log中的OBJECT和对应fits文件的是否相同，如果相同直接跳过，如果不同则进行修改
# 修改时判断要把fits中的OBJECT修改成log的OBJECT；还是修改成自己指定的
# 判断log中的IMAGETYP和对应fits文件的是否相同，如果相同直接跳过，如果不同则进行修改
# IMAGETYP只有fits头中有，判断要修改成log读取的还是自己指定修改


def renmfits(fitsnm):
    """
    replace the ":" in a file name to "_"
    """
    newfitsnm = fitsnm.replace(':', '_')
    # print('rename: {} --> {}'.format(fitsnm, newfitsnm))
    os.rename(fitsnm, newfitsnm)


def is_equal(a, b):
    if a == b:
        return True
    else:
        return False


def readimagetype(objectnm):
    '''
    根据OBJECT判断IMAGETYPE应该是什么
    '''

    if '[bias]' in objectnm:
        imagetype = 'bias'
    elif '[flat]' in objectnm:
        imagetype = 'flat'
    elif '[arc]' in objectnm:
        imagetype = 'arc'
    else:
        imagetype = 'science'
    return imagetype


def readlog(dir):
    '''
    log文件内容：
    # FILENAME 0 || DATE 1 || OBJECT 2 || RA 3 || DEC 4 ||
    # EQUINOX 5 || AIRMASS 6 || FILTER 7 || GRISM 8 || 
    # SLIT 9 || P.A. 10 || EXPTIME 11 || COMMENTS 12 ||
    '''
    os.chdir(dir)

    logDic = {'FILENAME': 0, 'DATE': 1, 'OBJECT': 2, 'RA': 3, 'DEC': 4, 
              'EQUINOX': 5, 'AIRMASS': 6, 'FILTER': 7, 'GRISM': 8, 
              'SLIT': 9, 'P.A.': 10, 'EXPTIME': 11, 'COMMENTS': 12}  # logDic写明了一行中每个元素的名称
    lognm = glob.glob('*.log')[0]  # 读取log文件名
    log = open(lognm, 'r').readlines()[1:]  # 读取log文件并将每行写入list
    # print(log[0])
    # print(log[0].split('||'))  # 每行按||分割
    # print([l.strip() for l in log[0].split('||')])  # 每行按||分割并去除空格

    # 首先确认目录下fits文件数目与log文件记录的是否相符
    fitsnumb_log = len(log)  # log文件中记录的fits文件数目
    fitslst = glob.glob('*.fits')
    fitsnumb_dir = len(fitslst)
    for i in range(fitsnumb_dir):
        renmfits(fitslst[i])
    print(fitsnumb_dir, fitsnumb_log)

    if is_equal(fitsnumb_dir, fitsnumb_log):
        # for i in range(1):  # for test
        for i in range(fitsnumb_log):
            infos = [line.strip() for line in log[i].split('||')]
            fitsnm = infos[logDic['FILENAME']]  # 读取fits名称
            # print(fitsnm)
            fitsnm = fitsnm.replace(':', '_')
            hdu = fits.open(fitsnm, ignore_missing_end=True)
            # print(hdu[0].data)
            header = hdu[0].header  # 相比fits.getheader，必须用这一行的才能修改header
            # header = fits.getheader(fitsnm, ignore_missing_end=True)

            # 注意：以下读取的4项都可能错误，但fits头中的错误可能会多一些
            # 读取log中的OBJECT
            object_log = infos[logDic['OBJECT']]
            # 读取fits头中的OBJECT
            object_fits = header['OBJECT']
            # 读取log中的IMAGETYPE
            imagetyp_log = readimagetype(object_log)
            # 读取fits头中的IMAGETYPE
            imagetyp_fits = header['IMAGETYP']

            inputnumb = 0
            if is_equal(imagetyp_log, imagetyp_fits) is False:
                print('{} || "{}" - "{}" || "{}" - "{}"'.format(fitsnm, object_log, object_fits, imagetyp_log, imagetyp_fits))
                imgerr_txt = '!!IMAGETYP error!!'
                print('{:>{}} || "{}" - "{}"'.format(imgerr_txt, len(fitsnm), imagetyp_log, imagetyp_fits))
                inputnumb = int(input(' -- IMAGETYP chioces: [1: imagetyp_log] or [2: imagetyp_fits] or [0: both not] or [3: OBJECT error]: '))
                if inputnumb == 0:
                    imagetyp_input = input(' -- Please input the right IMAGETYP: ')
                    imagetyp_new = imagetyp_input
                elif inputnumb == 1:
                    imagetyp_new = imagetyp_log
                else:
                    imagetyp_new = imagetyp_fits

                header['IMAGETYP'] = imagetyp_new
                hdu.writeto(fitsnm, overwrite=True)
                newheader = hdu[0].header
                print(' -- IMAGETYP: "{}" --> "{}"'.format(imagetyp_fits, newheader['IMAGETYP']))
                print('--' * 50)

            if is_equal(object_log, object_fits) is False or inputnumb == 3:
                print('{} || "{}" - "{}" || "{}" - "{}"'.format(fitsnm, object_log, object_fits, imagetyp_log, imagetyp_fits))
                objerr_txt = '!!OBJECT error!!'
                print('{:>{}} || "{}" - "{}"'.format(objerr_txt, len(fitsnm), object_log, object_fits))
                inputnumb = int(input(' -- OBJECT chioces: [1: object_log] or [2: object_fits] or [0: both not]: '))
                if inputnumb == 0:
                    object_input = input(' -- Please input the right OBJECT: ')
                    object_new = object_input
                elif inputnumb == 1:
                    object_new = object_log
                else:
                    object_new = object_fits

                header['OBJECT'] = object_new
                hdu.writeto(fitsnm, overwrite=True)
                newheader = hdu[0].header
                print(' -- OBJECT: "{}" --> "{}"'.format(object_fits, newheader['OBJECT']))
                print('--' * 50)

            hdu.close()

    else:
        print('Please check the number of fits files!!')

    os.chdir('../')  # 最后要回到上一级目录


def main2():
    dirnm = '250626'
    # print(dirnm[:4])
    # sys.stdout = open('changefits_{}.log'.format(dirnm), 'w')  # 将输出信息写入到log文件中

    dirlst = glob.glob('{}_CAFOS/'.format(dirnm))
    dirlst.sort()
    # for i in range(1):  # for test
    for i in range(len(dirlst)):
        # dataDir = '/media/yao/cahayao/data/' + dirlst[i]
        print(dirlst[i])
        readlog(dirlst[i])


def main():
    dirnm = sys.argv[1]
    readlog(dirnm)


if __name__ == "__main__":
    main()
    # main2()
