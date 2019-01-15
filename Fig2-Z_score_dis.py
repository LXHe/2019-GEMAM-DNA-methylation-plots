# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 10:03:02 2018

@author: u0105352
"""
#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

path1 = r'C:\Users\u0105352\Desktop\Lingxiao\PhD_file\2015-GEMAM-MoveAGE\Raw data\Elderly\DNA_methylation\Epigenetic_sample\methyspl_raw_phenotype.xlsx'
path2 = r'C:\Users\u0105352\Desktop\Lingxiao\PhD_file\2015-GEMAM-MoveAGE\Raw data\Elderly\DNA_methylation\Epigenetic_sample\65-80_z_score.xlsx'

sample = pd.read_excel(path1,sheet_name='methyspl', usecols='A,B')
total = pd.read_excel(path2)
#%%
methyl_sar_id = sample['ID'].loc[sample['STAT']=='Sar'].tolist()
methyl_nonsar_id = sample['ID'].loc[sample['STAT']=='N_Sar'].tolist()
tlt_sar_id = total.loc[total['STAT']=='Sar']
tlt_nonsar_id = total.loc[total['STAT']=='N_Sar']
#%%
tlt_sar_id['selection'] = np.where(tlt_sar_id['ID'].isin(methyl_sar_id), 'Selected', 'Not-selected')
tlt_nonsar_id['selection'] = np.where(tlt_nonsar_id['ID'].isin(methyl_nonsar_id), 'Selected', 'Not-selected')
#%%
# =============================================================================
# This can't show proper legend
# =============================================================================
select_list = ['Selected', 'Not-selected']
color_list = ['red', 'grey']
color_map = dict(zip(select_list, color_list))
#%%
tlt_sar_id.sort_values(by=['z_sum'], inplace=True)
tlt_sar_id.plot(y='z_sum', kind='bar', color=tlt_sar_id['selection'].map(color_map))
plt.xticks([])
plt.legend(loc='best')
plt.show()
#%%
tlt_nonsar_id.sort_values(by=['z_sum'], inplace=True)
tlt_nonsar_id.plot(y='z_sum', kind='bar', color=tlt_nonsar_id['selection'].map(color_map))
plt.xticks([])
plt.legend(loc='best')
plt.show()

#%%
# =============================================================================
# Use this method
# =============================================================================
sar_rank = np.linspace(1, 25, num=25, endpoint=True)
nonsar_rank = np.linspace(1, 143, num=143, endpoint=True)

tlt_sar_sum = tlt_sar_id.sort_values(by=['z_sum'])
tlt_sar_grip = tlt_sar_id.sort_values(by=['z_grip'])
tlt_sar_smi = tlt_sar_id.sort_values(by=['z_SMI'])

tlt_sar_sum['rank'] = sar_rank
tlt_sar_grip['rank'] = sar_rank
tlt_sar_smi['rank'] = sar_rank

sar_sum = tlt_sar_sum.loc[tlt_sar_sum['selection']=='Selected'][['z_sum', 'rank']]
sar_grip = tlt_sar_grip.loc[tlt_sar_grip['selection']=='Selected'][['z_grip', 'rank']]
sar_smi = tlt_sar_smi.loc[tlt_sar_smi['selection']=='Selected'][['z_SMI', 'rank']]

tlt_nonsar_sum = tlt_nonsar_id.sort_values(by=['z_sum'])
tlt_nonsar_grip = tlt_nonsar_id.sort_values(by=['z_grip'])
tlt_nonsar_smi = tlt_nonsar_id.sort_values(by=['z_SMI'])

tlt_nonsar_sum['rank'] = nonsar_rank
tlt_nonsar_grip['rank'] = nonsar_rank
tlt_nonsar_smi['rank'] = nonsar_rank

nonsar_sum = tlt_nonsar_sum.loc[tlt_nonsar_sum['selection']=='Selected'][['z_sum', 'rank']]
nonsar_grip = tlt_nonsar_grip.loc[tlt_nonsar_grip['selection']=='Selected'][['z_grip', 'rank']]
nonsar_smi = tlt_nonsar_smi.loc[tlt_nonsar_smi['selection']=='Selected'][['z_SMI', 'rank']]
#%%
fig = plt.figure()

ax1 = fig.add_subplot(321)
ax1.bar(tlt_sar_sum['rank'], tlt_sar_sum['z_sum'], label='Not-selected')
ax1.bar(sar_sum['rank'], sar_sum['z_sum'], label='Selected')
ax1.set_xticks([])
ax1.legend(loc='best')
ax1.set_title('Sarcopenic group')
ax1.set_ylabel('Summed z score')

ax3 = fig.add_subplot(323)
ax3.bar(tlt_sar_grip['rank'], tlt_sar_grip['z_grip'], label='Not-selected')
ax3.bar(sar_grip['rank'], sar_grip['z_grip'], label='Selected')
ax3.set_xticks([])
ax3.legend(loc='best')
ax3.set_ylabel('z score (hand grip)')

ax5 = fig.add_subplot(325)
ax5.bar(tlt_sar_smi['rank'], tlt_sar_smi['z_SMI'], label='Not-selected')
ax5.bar(sar_smi['rank'], sar_smi['z_SMI'], label='Selected')
ax5.set_xticks([])
ax5.legend(loc='best')
ax5.set_ylabel('z score (SMI)')

ax2 = fig.add_subplot(322)
ax2.bar(tlt_nonsar_sum['rank'], tlt_nonsar_sum['z_sum'], label='Not-selected')
ax2.bar(nonsar_sum['rank'], nonsar_sum['z_sum'], label='Selected')
ax2.set_xticks([])
ax2.legend(loc='best')
ax2.set_title('Non-sarcopenic group')
ax2.set_ylabel('Summed z score')

ax4 = fig.add_subplot(324)
ax4.bar(tlt_nonsar_grip['rank'], tlt_nonsar_grip['z_grip'], label='Not-selected')
ax4.bar(nonsar_grip['rank'], nonsar_grip['z_grip'], label='Selected')
ax4.set_xticks([])
ax4.legend(loc='best')
ax4.set_ylabel('z score (hand grip)')

ax6 = fig.add_subplot(326)
ax6.bar(tlt_nonsar_smi['rank'], tlt_nonsar_smi['z_SMI'], label='Not-selected')
ax6.bar(nonsar_smi['rank'], nonsar_smi['z_SMI'], label='Selected')
ax6.set_xticks([])
ax6.legend(loc='best')
ax6.set_ylabel('z score (SMI)')

plt.suptitle("Distribution of z scores for participant screening for DNA methylation measurement")
plt.show()