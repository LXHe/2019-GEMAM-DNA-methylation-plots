import pandas as pd
import scipy.stats as sp
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

path_glb = r'\Initial_RawData\Noob_1_drpT_beta_all_mrg_p01.xlsx'

all = pd.read_csv(r'\1-Descriptive\T_statistic_all_data.csv')
glb = pd.read_excel(
    path_glb,
    sheet_name='p01_S_vs_N_S_merged',
    usecols='C, CB, CE'
)

# Test sample equal variances
# 788074 is the total No. of analysed CpG sites
# 6258 is the No. of dmCpG sites
glb['a_avg'] = glb['ABS']/788074
glb['s_avg'] = glb['SBS']/6258
glb['Group'] = np.where(glb['STAT'].str.match('Sar'), 'Sarcopenic', 'Non-sarcopenic')

glb_1 = glb[['Group', 'a_avg']]
glb_1['Type'] = 'Tested CpGs'
glb_1.rename(columns={'a_avg': 'average'}, inplace = True)
glb_2 = glb[['Group', 's_avg']]
glb_2['Type'] = 'dmCpGs'
glb_2.rename(columns={'s_avg': 'average'}, inplace = True)
glb_mrg = glb_1.append(glb_2)

# Plot colored histogram
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

# Boxplot
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
