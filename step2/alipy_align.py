#! /usr/bin/env python
# -*- coding: utf-8 -*-

# @Author: zhixiang zhang <zzx>
# @Date:   24-Aug-2017
# @Email:  zhangzx@ihep.ac.cn
# @Filename: demo1.py
# @Last modified by:   zzx
# @Last modified time: 24-Aug-2017
# @Last modified time: 29-Jul-2025, by Hao Zhang

import os
import sys
import glob
import alipy
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


def save_lst(list_name, contents):
    f = open(list_name, 'w')
    contents.sort()
    for content in contents:
        f.write(content + '\n')
    f.close()


namelst = glob.glob("images/*.fits")

if len(namelst) == 0:
    print('!! No image will be aligned')
    sys.exit(1)

if os.path.isfile('updatelst'):
    updatefile = open('updatelst', 'r')
    updatelst = updatefile.readlines()
    updatefile.close()
    save_lst('updatelst', namelst)
    for i in updatelst:
        namelst.remove(i.split()[0])
    images_to_align = sorted(namelst)    
else:
    print('make a updatelst')
    save_lst('updatelst', namelst)
    images_to_align = sorted(namelst)

ref_image = 'ref.fits'

identifications = alipy.ident.run(ref_image, images_to_align, hdu=0, visu=False)
print(len(identifications))
# That's it !
# Put visu=True to get visualizations in form of png files (nice but much slower)
# On multi-extension data, you will want to specify the hdu (see API doc).
# The output is a list of Identification objects, which contain the transforms :
for id in identifications:  # list of the same length as images_to_align.
    if id.ok is True:  # i.e., if it worked
        
        print("%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio))
        # id.trans is a alipy.star.SimpleTransform object. Instead of printing it out as a string,
        # you can directly access its parameters :
        # print id.trans.v # the raw data, [r*cos(theta)  r*sin(theta)  r*shift_x  r*shift_y]
        # print id.trans.matrixform()
        # print id.trans.inverse() # this returns a new SimpleTransform object
    else:
        print("%20s : no transformation found !" % (id.ukn.name))


# Minimal example of how to align images :
outputshape = alipy.align.shape(ref_image)

# Hao Zhang modifies.
# def convert_fits_to_png(fits_files, png_dir=None, stretch='log', vmin=None, vmax=None, cmap='gray'):
#     """
#     将一批 fits 文件转换为 png, 原则上模拟 f2n 原始功能
#     fits_files: fits 文件路径列表
#     png_dir: 输出 png 文件夹, 默认与源文件同文件夹
#     stretch: 拉伸方式, 目前为'linear'或'log'，后续可扩展
#     vmin, vmax: 拉伸范围，若None自动确定(排除极端像元)
#     cmap: matplotlib的色表
#     """
#     def scale(data, stretch, vmin, vmax):
#         if vmin is None: vmin = np.percentile(data, 0.5)
#         if vmax is None: vmax = np.percentile(data, 99.5)
#         scaled = np.clip(data, vmin, vmax)
#         scaled = (scaled - vmin) / (vmax - vmin)
#         if stretch == 'log':
#             scaled = np.log10(1 + 9 * scaled)  # 0~1 -> log
#         elif stretch == 'sqrt':
#             scaled = np.sqrt(scaled)
#         elif stretch == 'linear':
#             pass
#         else:
#             raise ValueError("Unknown stretch: {}".format(stretch))
#         scaled = np.clip(scaled, 0, 1)
#         return scaled

#     for fits_file in fits_files:
#         with fits.open(fits_file) as hdul:
#             data = hdul[0].data.astype(float)
#         if data.ndim > 2:  # 多通道或有额外维
#             data = data[0]   # 只取第一层
#         scaled_data = scale(data, stretch, vmin, vmax)
#         # 输出目录
#         if png_dir is None:
#             out_dir = os.path.dirname(fits_file)
#         else:
#             out_dir = png_dir
#             os.makedirs(out_dir, exist_ok=True)
#         # 输出文件名
#         base = os.path.splitext(os.path.basename(fits_file))[0]
#         out_png = os.path.join(out_dir, base + '.png')
#         plt.imsave(out_png, scaled_data, cmap=cmap)


# This is simply a tuple (width, height)... you could specify any other shape.
for id in identifications:
    if id.ok is True:
        # Variant 1, using only scipy and the simple affine transorm :
        # alipy.align.affineremap(id.ukn.filepath, id.trans, shape=outputshape, makepng=True)

        # Variant 2, using geomap/gregister, correcting also for distortions :
        alipy.align.irafalign(id.ukn.filepath, id.uknmatchstars,
                              id.refmatchstars, shape=outputshape, makepng=True)
        # id.uknmatchstars and id.refmatchstars are simply lists of corresponding Star objects.

        # By default, the aligned images are written into a directory "alipy_out".

# To be continued ...
# fits_list = glob.glob('alipy_out/ftbcaf-20250609-23_34_19-sci-wanj_gregister.fits')
# convert_fits_to_png(fits_list)
sys.exit(0)
