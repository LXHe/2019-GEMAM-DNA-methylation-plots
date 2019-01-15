#%%
import pandas as pd
import scipy.stats as sp
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

all = pd.read_csv('C://Users/u0105352/Desktop/Lingxiao/MOVEAGE2015/Partek/Noob_1_drpT/\
1-STAT_NOCOV/1-Descriptive/T_statistic_all_data.csv')

path = r'C:\Users\u0105352\Desktop\Lingxiao\MOVEAGE2015\Initial_RawData\Noob_1_drpT_beta_all_mrg_p01.xlsx'

glb = pd.read_excel(
    path,
    sheet_name='p01_S_vs_N_S_merged',
    usecols='C, CB, CE'
    )
#%%
# Test sample equal variances
glb['a_avg'] = glb['ABS']/788074
glb['s_avg'] = glb['SBS']/6258
glb['Group'] = np.where(glb['STAT'].str.match('Sar'), 'Sarcopenic', 'Non-sarcopenic')
#%%
glb_a_sar = glb['a_avg'].loc[
    glb['STAT'] == 'Sar'
]
glb_a_nsar = glb['a_avg'].loc[
    glb['STAT'] == 'N_Sar'    
]
rlt_a_ev, p_a_ev = sp.levene(glb_a_sar, glb_a_nsar)
# p value shows that the samples have equal variance
# Welch's t-test
rlt_a_tt, p_a_tt = sp.ttest_ind(glb_a_sar, glb_a_nsar, equal_var = True)
#%%
glb_s_sar = glb['s_avg'].loc[
    glb['STAT'] == 'Sar'
]
glb_s_nsar = glb['s_avg'].loc[
    glb['STAT'] == 'N_Sar'    
]
rlt_s_ev, p_s_ev = sp.levene(glb_s_sar, glb_s_nsar)
# p value shows that the samples have equal variance
# Welch's t-test
rlt_s_tt, p_s_tt = sp.ttest_ind(glb_s_sar, glb_s_nsar, equal_var = True)
#%%
glb_1 = glb[['Group', 'a_avg']]
glb_1['Type'] = 'Tested CpGs'
glb_1.rename(columns={'a_avg': 'average'}, inplace = True)
glb_2 = glb[['Group', 's_avg']]
glb_2['Type'] = 'dmCpGs'
glb_2.rename(columns={'s_avg': 'average'}, inplace = True)
glb_mrg = glb_1.append(glb_2)
#%%
data = all['T(N_Sar vs. Sar)'].values
##print ('Max value is: {}'.format(data.max()))
##print ('Min value is: {}'.format(data.min()))

fig, ax = plt.subplots(nrows=1, ncols=3)

N, bins, patches = ax[0].hist(
    data, bins=16, range=(-6,6), edgecolor='black', 
    linewidth=1, fill=None, log=True
)
positive_range = bins[len(bins)//2:]
negative_range = bins[:len(bins)//2]
ax[0].hist(
    data, bins=positive_range, range=(0,6), edgecolor='black',
    linewidth=1, fill=None, log=True, hatch='/'
)
ax[0].hist(
    data, bins=4, range=(3,6), edgecolor='black',
    linewidth=1, log=True, color='#ff7f0e', alpha=0.5, label='dmCpG sites'
)
ax[0].hist(
    data, bins=4, range=(-6, -3), edgecolor='black',
    linewidth=1, log=True, color='#ff7f0e', alpha=0.5
)
ax[0].text(
    -0.05, 1.03, "A.", ha = 'center', va = 'center', transform = ax[0].transAxes
)
ax[0].legend(loc='upper right')
ax[0].set_title(
    "Distribution of T values in total CpG sites\n(Sarcopenic vs Non-sarcopenic)",
    fontsize=12
)
ax[0].set_xlabel("T value", fontsize=12)
ax[0].set_ylabel("No. of CpG sites", fontsize=12)

#glb.boxplot(
#    column = ['avg'], by = ['STAT'], ax = ax[1], grid = False
#    )
#ax[1].set_xticklabels(['Non-sarcopenic', 'Sarcopenic'])

sns.boxplot(
    x = 'Type', y = 'average', hue = 'Group', data = glb_mrg.loc[
        glb_mrg['Type'] == 'Tested CpGs'
    ],
    ax = ax[1], width = 0.5
)
ax[1].text(
    -0.07, 1.03, 'B.', ha = 'center', va = 'center', transform = ax[1].transAxes    
)
    
ax[1].set_title("Boxplot of average β value of \nall analysed CpGs", fontsize=12)
ax[1].set_xlabel("")
ax[1].set_xticklabels("")
ax[1].set_ylabel("Average β value", fontsize=12)

sns.boxplot(
    x = 'Type', y = 'average', hue = 'Group', data = glb_mrg.loc[
        glb_mrg['Type'] == 'dmCpGs'
    ],
    ax = ax[2], width = 0.5
)
ax[2].text(
    -0.07, 1.03, 'C.', ha = 'center', va = 'center', transform = ax[2].transAxes    
)
    
ax[2].set_title("Boxplot of average β value of \ndmCpGs", fontsize=12)
ax[2].set_xlabel("")
ax[2].set_xticklabels("")
ax[2].set_ylabel("Average β value", fontsize=12)

#plt.tight_layout()
plt.subplots_adjust(wspace = 0.4)
plt.suptitle("")
plt.show()
