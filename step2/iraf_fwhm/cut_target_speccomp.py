# -*- coding: utf-8 -*-
import pandas as pd
import os
header = ['JD', 'FWHM_target', 'FWHM_speccomp']

csv_file = os.path.join('psf_logs', 'PHL1092_psfall.csv')
df = pd.read_csv(csv_file)
print(df.columns)

ids = sorted(df['id'].unique())


df_target = df[df['id'] == ids[-2]]
df_comp = df[df['id'] == ids[-1]]
print(len(df_target))

df_target = df_target[['Radius', 'JD']]
df_comp = df_comp[['Radius', 'JD']]

df_target = df_target[['JD', 'Radius']].rename(columns={'Radius': 'FWHM_target'})
df_comp = df_comp[['JD', 'Radius']].rename(columns={'Radius': 'FWHM_speccomp'})

df_combine = pd.merge(df_target, df_comp, on='JD', how='left')
df_combine = df_combine[header]

df_combine = df_combine.sort_values('JD', ascending=True)

df_combine.to_csv(csv_file.replace('_psfall.csv','_comp_FWHM_measure.csv'), index=False)