# -*- coding: utf-8 -*-
import os
from pyraf import iraf

image_dir = 'alipy_out'
# image_dir = 'test'
pos_file = 'PHL1092.pos'
log_dir = 'psf_logs'

if not os.path.exists(log_dir):
    os.makedirs(log_dir)

def run_psfmeasure(image_path, pos_file, log_path):
    iraf.noao()
    iraf.obsutil()
    # iraf.set(stdgraph='null')
    iraf.psfmeasure(image_path, coords='markall', imagecur=pos_file, 
                    size='FWHM', display='no', logfile=log_path)

def main():
    for fname in os.listdir(image_dir):
        if fname.endswith('.fits'):
            image_path = os.path.join(image_dir, fname)
            log_path = os.path.join(log_dir, fname.replace('.fits', '_psf.txt'))
            run_psfmeasure(image_path, pos_file, log_path)

if __name__ == '__main__':
    main()
