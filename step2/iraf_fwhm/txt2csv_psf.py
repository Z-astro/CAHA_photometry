# -*- coding: utf-8 -*-
import csv
import re
import os
from astropy.io import fits

header = ['id', 'Column', 'Line', 'Mag', 'Radius', 'Ellip', 'PA', 'JD']
log_dir = 'psf_logs'

def txt2csv(input_txt, output_csv):
    rows = []
    with open(input_txt, 'r') as fin:
        for line in fin:
            line = line.rstrip()
            # 跳过非数据行
            if not line or line.startswith('NOAO/IRAF') or 'Column' in line or 'Average' in line:
                continue
            parts = re.split(r'\s+', line)
            # 只保留数据列，不添加image名
            if len(parts) == 7:
                rows.append(parts[1:])  # 去掉Image名
            elif len(parts) == 6:
                rows.append(parts)

    # 添加JD列(读取csv文件名相对应的fits文件，fits文件在alipy_out文件夹中，csv文件psf_logs中，再获取header中的JD值，给每一行都添加一样的JD值)
    # 推断fits文件名（假设txt文件名格式为xxx_psf.txt，对应fits为alipy_out/xxx.fits）
    fits_name = os.path.basename(input_txt).replace('_psf.txt', '.fits')
    fits_path = os.path.join('alipy_out', fits_name)
    jd = ''
    try:
        with fits.open(fits_path) as hdul:
            hdr = hdul[0].header
            jd = float(hdr.get('JD', 0.0))  # 获取JD值，默认为0.0
    except Exception as e:
        jd = 0.0
    # 为每行添加JD，并加id
    rows = [[i+1] + row + [jd] for i, row in enumerate(rows)]

    with open(output_csv, 'w') as fout:
        writer = csv.writer(fout)
        writer.writerow(header)
        writer.writerows(rows)

def main():
    for fname in os.listdir(log_dir):
        if fname.endswith('.txt'):
            input_txt = os.path.join(log_dir, fname)
            output_csv = os.path.join(log_dir, fname.replace('.txt', '.csv'))
            txt2csv(input_txt, output_csv)
            print('转换完成:', output_csv)

if __name__ == '__main__':
    main()
