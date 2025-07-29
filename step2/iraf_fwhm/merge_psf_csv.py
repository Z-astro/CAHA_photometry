# -*- coding: utf-8 -*-
import os
import csv

csv_dir = 'psf_logs'
output_csv = os.path.join(csv_dir, 'PHL1092_psfall.csv')

header = None
all_rows = []

for fname in os.listdir(csv_dir):
    if fname.endswith('.csv'):
        fpath = os.path.join(csv_dir, fname)
        with open(fpath, 'r') as fin:
            reader = csv.reader(fin)
            file_header = next(reader)
            if header is None:
                header = file_header
            for row in reader:
                all_rows.append(row)

with open(output_csv, 'w') as fout:
    writer = csv.writer(fout)
    writer.writerow(header)
    writer.writerows(all_rows)

print('合并完成，输出文件:', output_csv)
