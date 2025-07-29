#!/usr/bin/env python
# coding=utf-8
# running in python==3.7.16, astroconda.

import os
import sys
import glob
from pyraf import iraf
import numpy as np

        
def delete_file(namelst):
    for name in namelst:
        if os.path.isfile(name):
            os.remove(name)
    

def del_wrong_keyword(fn):
    '''
    caha数据需要删除一些关键词才可以进行操作
    delete the keyword CCDSEC, DATASEC, BIASSEC in fits fn, the value of the
    keyword is wrong in caha's data, iraf command ccdproc can not process
    the data.
    '''
    keywordlst = ['CCDSEC', 'DATASEC', 'BIASSEC']
    for keyword in keywordlst:
        iraf.hedit(images=fn, fields=keyword, value='', add='No', addonly='No',
                   delete='Yes', verify='No', show='Yes', update='Yes')
    # iraf.flpr()


def combine_bias(blst):
    # 合并bias生成Zero
    iraf.zerocombine(input='@' + blst, output='Zero', combine='average',
                     reject='minmax', ccdtype='', process=False, delete=False,
                     clobber=False, scale='none', statsec='', nlow=0, nhigh=1,
                     nkeep=1, mclip=True, lsigma=3.0, hsigma=3.0,
                     rdnoise='CCDRON', gain='CCDGAIN', snoise=0.0, pclip=-0.5,
                     blank=0.0)
    iraf.flpr()
    print('<<<<<combine bias successfully>>>>>')


def cor_zero(corblst):
    # 减bias(Zero.fits)
    iraf.ccdproc(images='@' + corblst, output='b@' + corblst,
                 ccdtype='', max_cache=0, noproc=False, fixpix=False,
                 overscan=False, trim=False, zerocor=True, darkcor=False,
                 flatcor=False, illumcor=False, fringecor=False, readcor=False,
                 scancor=False, readaxis='line', fixfile='', zero='Zero',
                 dark='', flat='', illum='', fringe='', minreplace=1.0,
                 scantype='shortscan', nscan=1, interactive=False,
                 function='chebyshev', order=1, sample='*', naverage=1,
                 niterate=1, low_reject=3.0, high_reject=3.0, grow=1.0)
    iraf.flpr()
    print('<<<<<correct bias successfully>>>>>')


def get_trimsec(width, hight):
    # 得出裁减尺寸
    if width == 1100 and hight == 1600:
        return '[1:1100, 430:1600]'
    elif width == 1100 and hight == 1200:
        return '[1:1100, 1:1200]'
    elif width == 1100 and hight == 1100:
        return '[1:1100, 1:1100]'
    elif width == 1024 and hight == 1024:
        # print '[1:1024,1:1024]'
        return '[1:1024, 1:1024]'
    elif width == 1650 and hight == 1650:
        return '[1:1650, 1:1650]'
    elif width == 1400 and hight == 2048:
        return '[150:1250, 440:1610]'
    else:
        print('can not get the trim sec, width = %d, hight = %d,' % (width, hight))
        sys.exit(1)



def trim(tblst, trimsec):
    # 根据尺寸裁减出所需要的文件, 显式指定需要关闭的改正步骤。这在python3里是非常必要的。
    # darkcor='No' illumcor='No' fringecor='No', readcor='No', scancor='No'
    iraf.ccdproc(images='b//@' + tblst, output='tb//@' + tblst, ccdtype='',
                 noproc=False, fixpix=False, overscan=False, trim=True, 
                 zerocor=False, darkcor='No', flatcor=False, illumcor='No', 
                 fringecor='No', readcor='No', scancor='No',
                 readaxis='line', biassec='', trimsec=trimsec,
                 interactive=False, function='legendre', order=1)
    iraf.flpr()
    print('<<<<<trim section successfully>>>>>')


def combine_flat(ftblst):
    # 合并平场
    iraf.flatcombine(input='tb//@' + ftblst, output='flat', combine='average',
                     reject='avsigclip', ccdtype='', process=False, subsets=False,
                     delete=False, clobber=False, scale='mode', statsec='', nlow=1, 
                     nhigh=1, nkeep=1, mclip=True, lsigma=3.0, hsigma=3.0, rdnoise='CCDRON',
                     gain='CCDGAIN', snoise=0.0, pclip=-0.5, blank=1.0)
    iraf.flpr()
    print('<<<<<combine flat successfully>>>>>')


def cor_flat(corftblst):
    # 除平场
    iraf.ccdproc(images='tb//@' + corftblst, output='ftb//@' + corftblst,
                 ccdtype='', max_cache=0, noproc=False, fixpix=False,
                 overscan=False, trim=False, zerocor=False, darkcor=False,
                 flatcor=True, illumcor=False, fringecor=False, readcor=False,
                 scancor=False, readaxis='line', fixfile='', biassec='',
                 trimsec='', zero='', dark='', flat='flat', illum='',
                 fringe='', minreplace=1.0, scantype='shortscan', nscan=1,
                 interactive=False, function='chebyshev', order=1, sample='*',
                 naverage=1, niterate=1, low_reject=3.0, high_reject=3.0,
                 grow=1.0)
    iraf.flpr()
    print('<<<<<correct flat successfully>>>>>')


def main():
    """
    correct zero and flat for photometric fits
    and also trim the fits.
    finaly get a reduced fits ftb*.fits.(f:flat, t:trim, b:bias)
    """
    print('=' * 50)
    print('<<<run reduce_img.py  # trim, bias, flat>>>>>')

    # 删除之前运行的结果
    delete_file(['Zero.fits'])
    pre_lst = glob.glob('b*.fits')   
    delete_file(pre_lst)
    pre_lst = glob.glob('tb*.fits')   
    delete_file(pre_lst)
    delete_file(['flat.fits'])
    pre_lst = glob.glob('ftb*.fits')   
    delete_file(pre_lst)
    print('<<<<<clear all the previous file>>>>>')

    iraf.imred()
    iraf.ccdred(instrument='ccddb$kpno/specphot.dat')

    alllst = [i.strip() for i in open('all.lst')]   
    for name in alllst:
        del_wrong_keyword(name.strip())
    print('<<<<<correct the wrong keyword>>>>>')    
        
    print('<<<<< combine bias >>>>>')
    combine_bias('zero.lst') 
   
    print('<<<<< correct Zero >>>>>')
    cor_zero('cor_zero.lst')    

    print('<<<<< trim_section >>>>>')
    trimlst = np.loadtxt('trimlst')
    width = int(trimlst[0])
    hight = int(trimlst[1])
    trim_section = get_trimsec(width, hight)
    trim('cor_zero.lst', trim_section)
    
    print('<<<<< combine flat >>>>>')
    combine_flat('flat.lst')

    print('<<<<< correct flat >>>>>')
    cor_flat('cor_flat.lst')    

    print('<<<<< run reduce_img.py successfully >>>>>  ')
    print('=' * 50)


if __name__ == '__main__':
    main()