#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

##ir = pd.read_csv('C:/Users/u0105352/Desktop/Lingxiao/MOVEAGE2015/Initial_RawData/\
##BC_corrected_methylation_beta/beta_noob_1_drpT_island_dis.csv', index_col='Island_region',header=0)
sig_ir = pd.read_csv('C:/Users/u0105352/Desktop/Lingxiao/MOVEAGE2015/Partek/Noob_1_drpT/\
4-p01_S_vs_N/sumdata/p01_S_vs_N_S_isl_dis.csv', index_col='Island_region', header=0)
methyl = pd.read_csv('C:/Users/u0105352/Desktop/Lingxiao/MOVEAGE2015/Graphs/unadj_sig_1_all_methl.csv', index_col='Island_region', header=0)

cpg_dis = pd.merge(methyl, sig_ir, left_index=True, right_index=True, how='outer')
cpg_dis['Pct_hypo'] = cpg_dis['Hypo']/cpg_dis['Frequency']*100
cpg_dis['Pct_hyper'] = cpg_dis['Hyper']/cpg_dis['Frequency']*100
sig_sum = cpg_dis['Frequency'].sum()
cpg_dis['Pct_sig'] = cpg_dis['Frequency']/sig_sum*100
#%%
# Plot pie chart
labels = cpg_dis.index
sizes_sig = cpg_dis['Pct_sig'].values
colors_sig = ['orange','deepskyblue','tomato','gold','lightgreen','orchid']
# To explode or not
#explode = (0,0,0.04,0,0,0) 

fig, ax1 = plt.subplots()
ax1.pie(sizes_sig, labels=labels, colors=colors_sig, 
        startangle=65, frame=False, autopct='%1.1f%%', 
        pctdistance=0.75, labeldistance=1.05)

# Create a ring shape
# Draw a circle with white facecolor and null linewidth
centre_circle = plt.Circle((0,0),0.5, fc='white',linewidth=0)
### Get a reference to the current figure
##fig = plt.gcf()
# Get the current Axes instance on the current figure
# Add the pre-drawn circle
ax1.add_artist(centre_circle)
ax1.set_title('Distribution of dmCpG sites in CpG regions\n', fontsize='large')
ax1.axis('equal')

plt.show()
#%%
fig, ax2 = plt.subplots()
### Plot dmCpG distribution
##idx_name = cpg_dis.index
##N = len(idx_name)
##ind = np.arange(N)
##ax2.bar(ind, cpg_dis['Pct_hyper'], label='Hyper', width=0.3)
### Create a stacked barChart by specifying the bottom
##ax2.bar(ind, cpg_dis['Pct_hypo'], bottom=cpg_dis['Pct_hyper'], label='Hypo', width=0.3)
##ax2.set_xticks(ind)
##ax2.set_xticklabels(idx_name)
##ax2.legend(loc='best')

# Plot dmCpG distribution with 0
idx_name = cpg_dis.index
N = len(idx_name)/2
ind = np.arange(N, step=0.5)
hyper_bar = ax2.bar(ind, cpg_dis['Pct_hyper'], label='Hypermethylated', width=0.2)
# Create a stacked barChart by specifying the bottom
hypo_bar = ax2.bar(ind, cpg_dis['Pct_hypo']*(-1), label='Hypomethylated', width=0.2)
ax2.set_xticks(ind)
ax2.set_xticklabels(cpg_dis.index)
ax2.set_yticklabels([0,75,50,25,0,25,50,75])
ax2.legend(loc='upper right')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.set_ylabel('Methylation percentage(%)', fontsize='large')
ax2.set_title('Methylation status in CpG regions\nSarcopenic vs Non-sarcopenic', fontsize='large')

# Display values in bar chart
def autolabel(bar):
    """
    Attach a text label above each bar displaying its height
    """
    for a in bar:
        height = a.get_height()
        if height>0:
            ax2.text(a.get_x() + a.get_width()/2., 
                     height+0.5, '%1.1f%%' % float(height), 
                     ha='center', va='bottom')
        elif height<0:
            ax2.text(a.get_x() + a.get_width()/2., 
                     height-4.5, '%1.1f%%' % abs(float(height)), 
                     ha='center', va='bottom')
        
autolabel(hyper_bar)
autolabel(hypo_bar)

##fig.tight_layout()
plt.show()




### Explode section with max value
### Target array ARRAY
### Get the index of the max value from ARRAY (1D)
##idx = np.argmax(ARRAY)
### If ARRAY is 2D, then specify axis will return an array index
### Then idx = np.argmax(ARRAY, axis=0 as row or 1 as column)
### Create explode zero array
##explode = np.zeros(len(ARRAY))
### Specify the index of explode consistent with that in ARRAY
### Set the value to 0.1
##explode[idx] = 0.1
### An explode code is created


