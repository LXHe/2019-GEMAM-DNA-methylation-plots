# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 11:11:53 2018

@author: u0105352
"""

#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as sp

path1 = r'C:\Users\u0105352\Desktop\Lingxiao\MOVEAGE2015\Initial_RawData\BC_corrected_methylation_beta\beta_noob_1_drpT_trp.csv'
path2 = r'C:\Users\u0105352\Desktop\Lingxiao\MOVEAGE2015\Initial_RawData\BC_corrected_methylation_beta\beta_noob_1_drpT_gene_annotation.csv'
path3 = r'C:\Users\u0105352\Desktop\Lingxiao\MOVEAGE2015\Initial_RawData\Noob_1_drpT_beta_all_mrg_p01.xlsx'
path4 = r'C:\Users\u0105352\Desktop\Lingxiao\MOVEAGE2015\Partek\Noob_1_drpT\4-p01_S_vs_N\sumdata\p01_S_vs_N_S.csv'

raw_data = pd.read_csv(path1)

all_gr = pd.read_csv(path2, usecols = ['ID', 'UCSC_RefGene_Group'])

status =pd.read_excel(
        path3, sheet_name = 'p01_S_vs_N_S_merged', usecols = 'A,C'
        )

sig_gr = pd.read_csv(path4, usecols = ['Probeset ID', 'UCSC_RefGene_Group'])
#%%
all_gr.rename(
    columns = {'ID': 'Probe', 'UCSC_RefGene_Group': 'Gene_region'},
    inplace = True
    )

all_gr.dropna(axis = 0, subset = ['Gene_region'], inplace = True)

sig_gr.rename(
    columns = {'Probeset ID': 'Probe', 'UCSC_RefGene_Group': 'Gene_region'},
    inplace = True
    )

sig_gr.dropna(axis = 0, subset = ['Gene_region'], inplace = True)

raw_data.set_index('ID', inplace = True)
status.set_index('Column ID', inplace = True)
#%%
gr_list = ["TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR"]
gr_all_dict = {}
gr_sig_dict = {}
beta_avg = pd.DataFrame()
all_avg = pd.DataFrame()
sig_avg = pd.DataFrame()
sum_avg = pd.DataFrame()

for region in gr_list:
    gr_all_dict[region] = all_gr['Probe'].loc[
        all_gr['Gene_region'].str.match(region)
        ].values.tolist()
    
    gr_sig_dict[region] = sig_gr['Probe'].loc[
        sig_gr['Gene_region'].str.match(region)    
        ].values.tolist()
    
    beta_avg['all_{}'.format(region)] = raw_data[
        gr_all_dict[region]
        ].sum(axis = 1) / len(gr_all_dict[region])

    beta_avg['sig_{}'.format(region)] = raw_data[
        gr_sig_dict[region]
        ].sum(axis = 1) / len(gr_sig_dict[region]) 
    
    all_avg['average'] = raw_data[
        gr_all_dict[region]
        ].sum(axis = 1) / len(gr_all_dict[region])
    
    all_avg['Gene Region'] = region 
    
    all_avg['Count'] = len(gr_all_dict[region])
    
    all_avg['Type'] = 'Detected CpGs'
    
    sig_avg['average'] = raw_data[
        gr_sig_dict[region]
        ].sum(axis = 1) / len(gr_sig_dict[region])
    
    sig_avg['Gene Region'] = region
    
    sig_avg['Count'] = len(gr_sig_dict[region])
    
    sig_avg['Type'] = 'dmCpGs'
    
    sum_avg = sum_avg.append(all_avg)
    sum_avg = sum_avg.append(sig_avg)
#%% 
status_sum = pd.merge(beta_avg, status, left_index = True, right_index = True)
status_sum_avg = pd.merge(sum_avg, status, left_index = True, right_index = True)
#%%
writer = pd.ExcelWriter(
    r'C:\Users\u0105352\Desktop\mean_beta_sum.xlsx',
    engine = 'xlsxwriter'
    )
status_sum.to_excel(writer, sheet_name = 'average_beta')
status_sum_avg.to_excel(writer, sheet_name = 'mrg_average_beta')
writer.save()
#%%
# ttest function
def ttest(ds1, ds2, p = 0.05):
    rlt_var, p_var = sp.levene(ds1, ds2)
    eq = p_var > p
    
    # If equal_variance is False, then Welch's ttest is performed
    rlt_tt, p_tt = sp.ttest_ind(ds1, ds2, equal_var = eq)
    return p_tt

p_sum = {}
mean_sum = {}
std_sum = {}

status_list = status_sum.columns.values.tolist()
status_list.remove('STAT')

for region in status_list:
    sar_value = status_sum[region].loc[
        status_sum['STAT'] == 'Sar'
        ]

    nsar_value = status_sum[region].loc[
        status_sum['STAT'] == 'N_Sar'
        ]

    p_sum[region] = ttest(sar_value, nsar_value)
    
    mean_sum['sar_{}'.format(region)] = sar_value.mean()
    mean_sum['nsar_{}'.format(region)] = nsar_value.mean()
    
    std_sum['sar_{}'.format(region)] = sar_value.std()
    std_sum['nsar_{}'.format(region)] = nsar_value.std()
#%%
# Detected CpG global beta value
glb = pd.DataFrame()
glb['value'] = raw_data.sum(axis=1) / len(raw_data.columns)
glb_mrg = pd.merge(glb, status, left_index = True, right_index = True)

sar_value = glb_mrg['value'].loc[glb_mrg['STAT'] == 'Sar']

nsar_value =  glb_mrg['value'].loc[glb_mrg['STAT'] == 'N_Sar']
#%%
# dmCpG global beta value
sig_cpg = sig_gr['Probe'].values.tolist()

sig_glb = pd.DataFrame()
sig_glb['value'] = raw_data[sig_cpg].sum(axis=1) / len(sig_cpg)
sig_glb_mrg = pd.merge(sig_glb, status, left_index = True, right_index = True)

sar_value = sig_glb_mrg['value'].loc[sig_glb_mrg['STAT'] == 'Sar']

nsar_value =  sig_glb_mrg['value'].loc[sig_glb_mrg['STAT'] == 'N_Sar']
#%%
status_sum_avg['Group'] = np.where(status_sum_avg['STAT']=='Sar', 'Sarcopenic', 'Non-sarcopenic')
#%%
# plot of 'All detected CpGs' 
fig, ax = plt.subplots()

sns.barplot(
    x = 'Gene Region', y = 'average', hue = 'Group', data = status_sum_avg.loc[
        status_sum_avg['Type'] == 'All detected CpGs'    
        ], 
    ci = 'sd', ax = ax
    )

ax.set_ylabel('Average β value')
plt.title('Average β value of all detected CpGs in different gene regions')
plt.show()
#%%
# plot of 'dmCpGs'
fig, ax = plt.subplots()

sns.barplot(
    x = 'Gene Region', y = 'average', hue = 'Group', data = status_sum_avg.loc[
        status_sum_avg['Type'] == 'dmCpGs'    
        ],
    ci = 'sd', ax = ax    
    )

# Error mark for TSS200
plt.plot(
    [0.23, 0.30], [0.23, 0.23],
    color = 'k', linewidth = 2,
    transform = ax.transAxes
    )

plt.text(0.265, 0.24,
    '*', ha = 'center', va = 'center', transform = ax.transAxes, fontsize = 17
    )

## Error mark for 1stExon
#plt.plot(
#    [0.54, 0.61], [0.25, 0.25],
#    color = 'k', linewidth = 2,
#    transform = ax.transAxes
#    )
#
#plt.text(0.575, 0.26,
#    '*', ha = 'center', va = 'center', transform = ax.transAxes, fontsize = 17
#    )

# Error mark for Body
plt.plot(
    [0.70, 0.77], [0.86, 0.86],
    color = 'k', linewidth = 2,
    transform = ax.transAxes
    )

plt.text(0.735, 0.87,
    '*', ha = 'center', va = 'center', transform = ax.transAxes, fontsize = 17
    )

# Error mark for 3'UTR
plt.plot(
    [0.855, 0.925], [0.97, 0.97],
    color = 'k', linewidth = 2,
    transform = ax.transAxes
    )

plt.text(0.887, 0.98,
    '*', ha = 'center', va = 'center', transform = ax.transAxes, fontsize = 17
    )

ax.set_ylabel('Average β value')
plt.title('Average β value of dmCpGs in different gene regions')
plt.show()  
