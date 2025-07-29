#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import pdb
import glob
import shutil
import subprocess
import numpy as np
import pandas as pd
from scipy import stats
from PyAstronomy import pyasl
from astropy.io import fits as pyfits

# 1. 进入到每天的目录（如`230626_CAFOS/`）下，读取当天所有fits，
#    然后提炼出测光观测了哪些源（这里根据上面所说，按fits头中的坐标信息来识别源），
#    更新出一个observed list（这一步之前没有；这里需要额外注意的是如果源有B/R band的观测，就在源名后加上B/R，如Mrk509B，NGC7603R）
# 0 - all.lst -- [bias], [flat] Dome/sky V, [object]
# 1 - zero.lst -- [bias]
# 3 - domeflat.lst -- [flat] Dome V
# 2 - flat.lst -- [flat] sky V 
# 4 - cor_zero.lst -- [object], [flat] Dome/sky V
# 5 - cor_flat.lst -- [object]
# - obj*.lst -- [object]


def check_objnm(objnmf, raf, decf):
    scriptpath = os.path.dirname(os.path.abspath(__file__))
    objradeclst = [i.split() for i in open(scriptpath + os.sep + 'objradec.lst')]
    # print(objradeclst)
    for i, objradec in enumerate(objradeclst):
        # print(i, objradec)
        angdis = pyasl.getAngDist(float(objradec[1]), float(objradec[2]), raf, decf)
        objnm = ''
        if angdis < 0.6:
            objnm = objradec[0]
            break
        # else:
        #     objnm = ''
        #     print('There is no matched target!!')

    if objnm == '':
        print('{:<14}{:<10.6f}{:<+10.6f} not found in objradec.lst!!!'.format(objnmf, raf, decf))
        print('Please check the fits file, edit and save objradec.lst!!!')
        subprocess.run(['xdg-open', scriptpath + os.sep + 'objradec.lst'])
        input('press any key to reboot')
    return objnm


def savelst(listnm, lst):
    '''
    保存list文件
    '''
    f = open(listnm, 'w')
    lst.sort()
    for content in lst:
        f.write(content + '\n')
    f.close()
    print('generate {}'.format(listnm))


def main():
    # 获取当天进行的测光观测list
    # 将fits所对应的fitslst序号也对应标出来
    # 按序号把fits复制到指定目录下
    widthlst, hightlst = [], []
    allphoto_obslst, bias_obslst = [], []
    allflat_obslst, allskyflat_obslst, alldomeflat_obslst = [], [], []
    insfilt_sets = ['JohnsonV', 'JohnsonB', 'BessellR']
    filt_suffix_sets = ['', 'B', 'R']
    filt_suffix_dicts = dict(zip(insfilt_sets, filt_suffix_sets))
    # print(filt_suffix_dicts)
    fitsnmlst = glob.glob('caf*.fits')
    fitsnmlst.sort()
    for fitsnumb, fitsnm in enumerate(fitsnmlst):  # 一次for循环区分出各种不同类型文件
        imgtypf = pyfits.getval(fitsnm, 'IMAGETYP')  # 'IMAGETYP'有4种：bias，flat，science，arc
        # print(fitsnm, imgtypf, width, hight)
        if imgtypf == 'science':
            insfilt = pyfits.getval(fitsnm, 'INSFLNAM')
            if insfilt in insfilt_sets:
                width = pyfits.getval(fitsnm, 'NAXIS1')
                hight = pyfits.getval(fitsnm, 'NAXIS2')
                widthlst.append(width)
                hightlst.append(hight)
                objnmf = pyfits.getval(fitsnm, 'OBJECT')
                raf = pyfits.getval(fitsnm, 'RA')
                decf = pyfits.getval(fitsnm, 'DEC')
                objnm = ''
                while objnm == '':
                    objnm = check_objnm(objnmf, raf, decf)
                # print(objnm)
                objnm = objnm + filt_suffix_dicts[insfilt]
                allphoto_obslst.append([fitsnumb, objnm, insfilt])
            else:
                pass
        elif imgtypf == 'bias':
            bias_obslst.append(fitsnumb)
        elif imgtypf == 'flat':
            insfilt = pyfits.getval(fitsnm, 'INSFLNAM')
            objnmf = pyfits.getval(fitsnm, 'OBJECT')
            if insfilt in insfilt_sets:
                allflat_obslst.append([fitsnumb, objnmf, insfilt])
                if 'sky' in objnmf.lower():
                    allskyflat_obslst.append([fitsnumb, objnmf, insfilt]) 
                elif 'dome' in objnmf.lower():
                    alldomeflat_obslst.append([fitsnumb, objnmf, insfilt])
                else:
                    print('!! Wrong fits header of {}'.format(fitsnm))
                    print('!! Please check fits again')
        elif imgtypf == 'arc':
            pass
        else:
            print('{}: Can not identify "IMAGETYP"'.format(fitsnm))

    # 统计所有science的fits文件，以NAXIS众数为判断依据，检查fits文件是否有误
    # print(widthlst)
    # print(hightlst)
    width = stats.mode(widthlst)
    hight = stats.mode(hightlst)
    # width = np.min(widthlst)
    # hight = np.min(hightlst)
    print(width, hight)
    if width[0].shape[0] * hight[0].shape[0] == 1:
        # print(width[0][0], hight[0][0])
        trimlst = list([str(width[0][0]), str(hight[0][0])])
        # print(trimlst)
    else:
        print('!! Wrong size of fits')
        print('!! Please check fits again')

    # pdb.set_trace()
    allphoto_df = pd.DataFrame(allphoto_obslst, columns=['fitsnumb', 'objnm', 'filt'])
    # allflat_df = pd.DataFrame(allflat_obslst, columns=['fitsnumb', 'objnm', 'filt'])
    allskyflat_df = pd.DataFrame(allskyflat_obslst, columns=['fitsnumb', 'objnm', 'filt'])
    alldomeflat_df = pd.DataFrame(alldomeflat_obslst, columns=['fitsnumb', 'objnm', 'filt'])

    date_dirnm = os.path.basename(os.getcwd())  # 获取当前目录名称
    date = date_dirnm[:6]

    # 判断文件存在与否并复制
    if len(allphoto_df) == 0:
        print('No photometry observations in {}'.format(date))
    else:
        # 区分不同band
        # print(date, allphoto_obslst)
        for i in range(len(insfilt_sets)):
            # print(insfilt_sets[i][-1])
            photo_dir = '{}_photo{}/'.format(date, filt_suffix_dicts[insfilt_sets[i]])
            # print(photo_dir)
            photo_df = allphoto_df[allphoto_df['filt'] == insfilt_sets[i]]
            photo_df = photo_df.reset_index(drop=True)  # 重置序列号
            skyflat_df = allskyflat_df[allskyflat_df['filt'] == insfilt_sets[i]]
            skyflat_df = skyflat_df.reset_index(drop=True)
            domeflat_df = alldomeflat_df[alldomeflat_df['filt'] == insfilt_sets[i]]
            domeflat_df = domeflat_df.reset_index(drop=True)
            # print(photo_df)
            # print(skyflat_df)
            # print(domeflat_df)

            if len(photo_df) == 0:
                print('!! No {} band photometry observations in {}'.format(insfilt_sets[i][-1], date))
            else:
                lstsets = [[] for i in range(6)]
                lstnms = ['all', 'zero', 'flat', 'domeflat', 'cor_zero', 'cor_flat']
                lstdicts = dict(zip(lstnms, lstsets))

                os.makedirs(photo_dir, exist_ok=True)
                print('==> mkdir {}'.format(photo_dir))                   
                # 复制trimlst
                savelst('{}trimlst'.format(photo_dir), trimlst)

                # 复制object
                for j in range(len(photo_df)):
                    old_photo_fitsnm = fitsnmlst[photo_df.loc[j, 'fitsnumb']]
                    # print(old_photo_fitsnm)
                    shutil.copy2(old_photo_fitsnm, photo_dir)
                    print('> copy [phot]: {} --> {}'.format(old_photo_fitsnm, photo_dir))
                    lstsets[0].append(old_photo_fitsnm)
                    lstsets[4].append(old_photo_fitsnm)
                    lstsets[5].append(old_photo_fitsnm)

                # 复制bias
                if len(bias_obslst) == 0:
                    print('!! No [bias] in {}'.format(date))
                    print('!! Please copy [bias] from other day')
                else:                
                    for j in range(len(bias_obslst)):
                        old_bias_fitsnm = fitsnmlst[bias_obslst[j]]
                        shutil.copy2(old_bias_fitsnm, photo_dir)
                        print('> copy [bias]: {} --> {}'.format(old_bias_fitsnm, photo_dir))
                        lstsets[0].append(old_bias_fitsnm)
                        lstsets[1].append(old_bias_fitsnm)

                # 复制flat（如果没有skyflat就复制domeflat，如果domeflat也没有就去复制其他天的）
                # 默认采用sky flat作为平场，如果没有sky flat则采用dome flat
                if len(skyflat_df) == 0:
                    print('!! No [skyflat] in {}'.format(date))
                    print('!! Please copy [domeflat] first')
                    if len(domeflat_df) == 0:
                        print('!! No [flat] in {}'.format(date))
                        print('!! Please copy [flat] from other day')
                    else:                       
                        for j in range(len(domeflat_df)):  
                            old_flat_fitsnm = fitsnmlst[domeflat_df.loc[j, 'fitsnumb']]
                            shutil.copy2(old_flat_fitsnm, photo_dir)
                            print('> copy [flat] dome: {} --> {}'.format(old_flat_fitsnm, photo_dir))
                            lstsets[0].append(old_flat_fitsnm)
                            lstsets[2].append(old_flat_fitsnm)
                            lstsets[3].append(old_flat_fitsnm)
                            lstsets[4].append(old_flat_fitsnm)
                else:
                    for j in range(len(skyflat_df)):  
                        old_flat_fitsnm = fitsnmlst[skyflat_df.loc[j, 'fitsnumb']]
                        shutil.copy2(old_flat_fitsnm, photo_dir)
                        print('> copy [flat] sky: {} --> {}'.format(old_flat_fitsnm, photo_dir)) 
                        lstsets[0].append(old_flat_fitsnm)
                        lstsets[2].append(old_flat_fitsnm)
                        lstsets[4].append(old_flat_fitsnm)
                    if len(domeflat_df) != 0:
                        for j in range(len(domeflat_df)):
                            old_flat_fitsnm = fitsnmlst[domeflat_df.loc[j, 'fitsnumb']]
                            shutil.copy2(old_flat_fitsnm, photo_dir)
                            print('> copy [flat] dome: {} --> {}'.format(old_flat_fitsnm, photo_dir))
                            lstsets[3].append(old_flat_fitsnm)
                    else:
                        print('!! No [domeflat] in {}'.format(date))

                for lstnm in lstnms:
                    new_lstnm = '{}{}.lst'.format(photo_dir, lstnm)
                    savelst(new_lstnm, lstdicts[lstnm])


if __name__ == "__main__":
    main()
    # main2()
