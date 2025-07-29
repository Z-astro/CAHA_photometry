# -*- coding: utf-8 -*-
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

sys.path.append(os.path.expanduser('~/astroWORK/codesyao/plot/'))
# from pyplotsettings import set_mpl_style

# set_mpl_style()
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix' # stix,cm
plt.rcParams['font.size'] = 12
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rc('text', usetex=True)

markersize, elw = 10, 0.75  # 误差棒的线宽

csv_file = os.path.join('psf_logs', 'PHL1092_psfall.csv')
data = pd.read_csv(csv_file)

jd0 = int(data['JD'].min() - 100)  # 最小JD值，用于归一化

ids = sorted(data['id'].unique())
n = len(ids)
fig, axes = plt.subplots(n, 1, figsize=(16, 3*n), sharex=True)
fig.subplots_adjust(wspace=0, hspace=0.06)
fig.align_ylabels(axes)  # 对齐y轴标签
if n == 1:
    axes = [axes]

for i in range(n):
    ax = axes[i]
    group = data[data['id'] == ids[i]]
    ax.errorbar(group['JD'] - jd0, group['Radius'], fmt='.', 
                markersize=markersize, elinewidth=elw, capsize=3, 
                label='seeing (FWHM)', color='C0')
    # label = 'id={}'.format(i + 1)
    if i  == max(ids):
        label = 'Sepc.comp'
    # elif i == max(ids):
    #     label = 'Spec.comp'
    else:
        label = 'id={}'.format(i + 1)
    ax.text(0.98, 0.95, label, transform=ax.transAxes, ha='right', va='top')
    ax.set_ylabel('FWHM (pixels)', fontsize=20)
    # 画出平均值和1sigma虚线
    mean_radius = group['Radius'].mean()
    std_radius = group['Radius'].std()
    ax.axhline(mean_radius, color='C1', linestyle='--', label='Mean Radius')
    ax.axhline(mean_radius + std_radius, color='C1', linestyle='--', linewidth=elw)
    ax.axhline(mean_radius - std_radius, color='C1', linestyle='--', linewidth=elw)
    ax.text(0.15, 0.95, f'{mean_radius:.3f}({std_radius:.3f})', transform=ax.transAxes, 
            color='C1', ha='right', va='top',fontweight='bold')

axes[-1].set_xlabel('JD - {:.0f}'.format(jd0))
# plt.tight_layout()
# pdfnm = csv_file.replace('.csv', '.pdf')
pdfm = pdfnm = csv_file.replace('_psfall.csv', '_allcomp.pdf')
plt.savefig(pdfnm, format='pdf', dpi=1200, 
            orientation='landscape', bbox_inches='tight')
# plt.show()
plt.close()
