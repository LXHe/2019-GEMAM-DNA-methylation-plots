import pandas as pd
import matplotlib.pyplot as plt

ir = pd.read_csv(r'\BC_corrected_methylation_beta\beta_noob_1_drpT_island_dis.csv')
sig_ir = pd.read_csv(r'\4-p01_S_vs_N\sumdata/p01_S_vs_N_S_isl_dis.csv')

df_ir = ir.merge(sig_ir, on='Island_region', how='outer')

head = ("Region", "Total CpG sites", "dmCpG sites")
df_ir.columns = head

#====================================================================================
# Draw a broken axis plot
# The idea is to set two subplots and hide the spines between them
# Set two subplots
#====================================================================================
##fig, ax = plt.subplots(2, 1, sharex=True)
##df_ir.plot(kind='bar', ax=ax[0])
##df_ir.plot(kind='bar', ax=ax[1])
##
### Set the y axis range
##ax[0].set_ylim(20000, 480000)
##ax[1].set_ylim(0, 18000)
##
### Set maxium number of ticks that will be displayed in y axis
##ax[0].yaxis.set_major_locator(plt.MaxNLocator(5))
##ax[1].yaxis.set_major_locator(plt.MaxNLocator(5))
##
### Hide the legend of the second subplot
##ax[1].legend().set_visible(False)
##
### Hide the spins between the plots
##ax[0].spines['bottom'].set_visible(False)
##ax[1].spines['top'].set_visible(False)
##
### Hide x ticks of the 1st plot
##ax[0].xaxis.set_visible(False)
### If want ticks to be presented on the top:
###ax[0].xaxis.tick_top()
##
### In the 2nd plot, set x tick at the bottom of the bottom spine
##ax[1].xaxis.tick_bottom()
### Set the x tick label of the 2nd plot and rotate it
### If want a rotation, set 'rotation=20'
##ax[1].set_xticklabels(df_ir['Island_region'].values, rotation='horizontal')
### Set x label
##ax[1].set_xlabel('Island Region', fontsize='large')
##
### Set common y label
### Create a text, set vertical and horizontal alignment, rotate it
##fig.text(0.06, 0.5, 'No. of CpG sites', ha='center', va='center', rotation='vertical', fontsize='large')
##
##### Or add a big axes, hide frame
####ax = fig.add_subplot(111, frameon=False)
##### Hide tick and tick label of the big axes
####plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
####plt.grid(False)
####plt.ylabel('No. of CpG sites')
##
### Create broken axises
##d = 0.01
### Set keyword argument which returns as a dictionary
### transform enables to transform the coordinate system of the Axes;
### (0,0) is bottom left of the axes, and (1,1) is top right of the axes.
### A 'False' clip_on enables the drawings (artists) visible outside the axes
##kwargs = dict(transform=ax[0].transAxes, color='k', clip_on=False)
### Draw lines representing broken axises
##ax[0].plot((-d,+d),(-d,+d), **kwargs)
##ax[0].plot((1-d,1+d),(-d,+d), **kwargs)
##kwargs.update(transform=ax[1].transAxes)
##ax[1].plot((-d,+d),(1-d,1+d), **kwargs)
##ax[1].plot((1-d,1+d),(1-d,1+d), **kwargs)
##
### Set plot title
##ax[0].set_title("Distribution of CpG sites in island region", fontsize='large')
##
### Adjust the distance between the two plots
##fig.subplots_adjust(hspace=0.05)
##plt.show()

#====================================================================================
# Draw a logarithmic scaled plot
#====================================================================================
fig, ax = plt.subplots()
df_ir.plot(kind='bar', ax=ax, log=True)
ax.set_xticklabels(df_ir['Region'].values, rotation='horizontal')
ax.set_xlabel('CpG Regions', fontsize='large')
ax.set_ylabel('No. of CpG sites', fontsize='large')
ax.set_title("Proportion of dmCpG sites in CpG regions", fontsize='large')
plt.show()