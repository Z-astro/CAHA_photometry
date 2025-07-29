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
# plt.rc('text', usetex=True)

markersize, elw = 10, 0.75  # 误差棒的线宽

csv_file = os.path.join('psf_logs', 'PHL1092_psfall.csv')
data = pd.read_csv(csv_file)

jd0 = int(data['JD'].min() - 100)  # 最小JD值，用于归一化

# 获取所有唯一ID并排序
ids = sorted(data['id'].unique())
# 只保留最后两个ID（倒数第二和倒数第一）
selected_ids = ids[-2:]

# 创建2个子图（上下排列）
fig, axes = plt.subplots(2, 1, figsize=(16, 6), sharex=True)
fig.subplots_adjust(wspace=0, hspace=0.06)
fig.align_ylabels(axes)  # 对齐y轴标签

# 定义目标源和比较星的标签
labels = {
    selected_ids[0]: 'PHL1092',  # 倒数第二是目标源
    selected_ids[1]: 'Spec.comp' # 倒数第一是光谱比较星
}

for i, id_val in enumerate(selected_ids):
    ax = axes[i]
    group = data[data['id'] == id_val]
    
    # 绘制数据点
    ax.errorbar(group['JD'] - jd0, group['Radius'], fmt='.', 
                markersize=markersize, elinewidth=elw, capsize=3, 
                color='C0')
    
    # 添加粗体标签文本
    ax.text(0.98, 0.95, labels[id_val], transform=ax.transAxes, 
            ha='right', va='top', fontweight='bold')
    
    ax.set_ylabel('FWHM (pixels)', fontsize=20)
    
    # 计算并绘制平均值和标准差
    mean_radius = group['Radius'].mean()
    std_radius = group['Radius'].std()
    ax.axhline(mean_radius, color='C1', linestyle='--')
    ax.axhline(mean_radius + std_radius, color='C1', linestyle='--', linewidth=elw)
    ax.axhline(mean_radius - std_radius, color='C1', linestyle='--', linewidth=elw)
    
    # 添加粗体统计文本
    ax.text(0.15, 0.95, f'{mean_radius:.3f}({std_radius:.3f})', 
            transform=ax.transAxes, color='C1', ha='right', va='top',
            fontweight='bold')

# 设置公共X轴标签
axes[-1].set_xlabel('JD - {:.0f}'.format(jd0))

# 保存PDF文件
pdfnm = csv_file.replace('.csv', '.pdf')
plt.savefig(pdfnm, format='pdf', dpi=1200, 
            orientation='landscape', bbox_inches='tight')
plt.close()