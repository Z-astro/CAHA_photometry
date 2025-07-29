# -*- coding: utf-8 -*-
import os
from pyraf import iraf

# 设置环境变量以确保非交互模式
os.environ['PYRAF_BETA'] = '1'
os.environ['IRAF_IMCUR'] = ''
os.environ['DISPLAY'] = ':0.0'  # 设置显示环境变量
os.environ['TERM'] = 'xterm'    # 设置终端类型

# 图像文件夹和坐标文件
image_dir = 'alipy_out'  # 图像文件夹路径
# image_dir = 'test'  # 图像文件夹路径
pos_file = 'PHL1092.pos'  # 你的星点坐标文件
log_dir = 'psf_logs'  # 日志输出文件夹

# 创建日志文件夹
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

def run_psfmeasure(image_path, pos_file, log_path):
    iraf.noao()
    iraf.obsutil()
    
    # 设置为非交互模式，使用虚拟显示设备
    iraf.set(stdimage='imt1024')  # 使用虚拟显示设备
    iraf.set(imcur='')  # 清空图像光标
    iraf.set(imtype='fits')  # 设置图像类型
    
    # 运行psfmeasure，只使用支持的参数
    try:
        # 首先尝试最简单的调用
        iraf.psfmeasure(image_path, 
                       imagecur=pos_file, 
                       coords='markall',
                       size='FWHM',
                       logfile=log_path)
    except Exception as e:
        print(f"Error processing {image_path}: {e}")
        # 如果还是失败，尝试更基本的调用
        try:
            iraf.psfmeasure(image_path, 
                           imagecur=pos_file,
                           logfile=log_path)
        except Exception as e2:
            print(f"Second attempt failed for {image_path}: {e2}")

def main():
    for fname in os.listdir(image_dir):
        if fname.endswith('.fits'):
            image_path = os.path.join(image_dir, fname)
            log_path = os.path.join(log_dir, fname.replace('.fits', '_psf.txt'))
            print(f"Processing {fname}...")
            run_psfmeasure(image_path, pos_file, log_path)
    print("PSF measurement completed for all images.")

if __name__ == '__main__':
    main()
