import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sig_ir = pd.read_csv(r'\sumdata\p01_S_vs_N_S.csv')

ir = pd.read_csv(r'\BC_corrected_methylation_beta\Sar_vs_N_Sar_ttest.csv')

# Methylation proportion by chromosomes and CpG regions
sig_ir.rename(columns={'Relation_to_UCSC_CpG_Island':'Region',
                       'Difference(Sar vs. N_Sar) (Description)':'description'},
              inplace=True
)
sig_ir['Region'].fillna('OpenSea', inplace=True)
sig_ir['Status'] = np.where(
        sig_ir['description'] == 'Sar down vs N_Sar',
        'Hypomethylated',
        'Hypermethylated'
)
count_sig = sig_ir[
        ['Status','CHR','Region']
        ].groupby(
            ['CHR','Region','Status']
        ).Status.count()

# Count by chromosome
count_sig_chr = sig_ir.groupby(['CHR','Status']).CHR.count()
# Select by the second level of index, specify the index(axis=0) location in loc function
hyper = count_sig_chr.loc(axis=0)[:,'Hypermethylated']
hypo = count_sig_chr.loc(axis=0)[:,'Hypomethylated']

# =============================================================================
# Prepare DataSet for absolute count of dmCpGs in chromosomes
# The following block will sort the CHR.
# First need to convert CHR(x chomosome dropped) into numeric
# Then add x chomosome info back
# =============================================================================
hyper_df = pd.DataFrame(hyper)
hyper_df.rename(columns={"CHR":"No_chr"}, inplace=True)
# If only want to reset specific index, call:
# hyper_df.reset_index('CHR', inplace=True)
hyper_df.reset_index(inplace=True)
hyper_df_dpx = hyper_df.loc[hyper_df["CHR"]!="X"]
hyper_df_dpx["CHR"] = hyper_df_dpx["CHR"].astype('int32')
hyper_df_dpx.sort_values(by="CHR", inplace=True)
hyper_df_cvt = hyper_df_dpx.append(hyper_df.loc[hyper_df["CHR"]=="X"])

hypo_df = pd.DataFrame(hypo)
hypo_df.rename(columns={"CHR":"No_chr"}, inplace=True)
hypo_df.reset_index(inplace=True)
hypo_df_dpx = hypo_df.loc[hypo_df["CHR"]!="X"]
hypo_df_dpx["CHR"] = hypo_df_dpx["CHR"].astype('int32')
hypo_df_dpx.sort_values(by="CHR", inplace=True)
hypo_df_cvt = hypo_df_dpx.append(hypo_df.loc[hypo_df["CHR"]=="X"])
# =============================================================================
# Prepare DataSet for proportion of dmCpGs by chromosomes
# =============================================================================
ir.rename(columns={'Island_region':'Region'}, inplace=True)
ir['Status'] = np.where(ir['Describe'] == 'Sar down vs N_Sar', 'Hypomethylated', 'Hypermethylated')
count_all = ir[['Status','CHR','Region']].groupby(
                   ['CHR','Region','Status']
                   ).Status.count()

concate = pd.concat([count_all, count_sig], axis=1)
concate.fillna(0, inplace=True)

# Methylation proportion by chromosome
count_sig_chr = sig_ir.groupby('CHR').CHR.count()
count_all_chr = ir.groupby('CHR').CHR.count()
# rename SERIES type, only need to specify the name. No need to specify columns as that in DataFrame. 
count_sig_chr.rename('CHR_sig', inplace=True)
count_all_chr.rename('CHR_all', inplace=True)

def cvt_sort(DataSet_all, DataSet_sig):
    count_test = pd.DataFrame(DataSet_all)
    count_test.reset_index(inplace=True)
    count_test_sig = pd.DataFrame(DataSet_sig)
    count_test_sig.reset_index(inplace=True)
    
    count_test["CHR_all"].loc[count_test["CHR"]=="X"] = count_test.loc[count_test["CHR"]=="X"]["CHR_all"].values \
    + count_test.loc[count_test["CHR"]=="Y"]["CHR_all"].values
    # Drop Y CHR
    count_all_dpy = count_test.loc[count_test["CHR"]!="Y"]

    concate_test = count_all_dpy.merge(count_test_sig, on="CHR", how="outer")
    concate_test.fillna(0, inplace=True)
    concate_test['Proportion'] = concate_test['CHR_sig']/concate_test['CHR_all']*100

    concate_test_dpx = concate_test.loc[concate_test["CHR"]!="X"]
    concate_test_dpx["CHR"] = concate_test_dpx["CHR"].astype("int32")
    # To perform a descending sort:
    # concate_test.sort_values(by=['Proportion'], ascending=False, inplace=True)
    concate_test_dpx.sort_values(by="CHR", inplace=True)
    concate_test_cvt = concate_test_dpx.append(concate_test.loc[concate_test["CHR"]=="X"])
    
    return concate_test_cvt

concate_chr_cvt = cvt_sort(count_all_chr, count_sig_chr)
# =============================================================================
# Prepare DataSet for propoertion of hyper/hypo methylated dmCpGs by chromosomes
# =============================================================================
count_sig_region_hyper = sig_ir.loc[sig_ir['Status']=='Hypermethylated'].groupby('CHR').Status.count()
count_all_region_hyper = ir.loc[ir['Status']=='Hypermethylated'].groupby('CHR').Status.count()
count_sig_region_hyper.rename('CHR_sig', inplace=True)
count_all_region_hyper.rename('CHR_all', inplace=True)

concate_hyper_cvt = cvt_sort(count_all_region_hyper, count_sig_region_hyper)

count_sig_region_hypo = sig_ir.loc[sig_ir['Status']=='Hypomethylated'].groupby('CHR').Status.count()
count_all_region_hypo = ir.loc[ir['Status']=='Hypomethylated'].groupby('CHR').Status.count()
count_sig_region_hypo.rename('CHR_sig', inplace=True)
count_all_region_hypo.rename('CHR_all', inplace=True)

concate_hypo_cvt = cvt_sort(count_all_region_hypo, count_sig_region_hypo)
# =============================================================================
# Plot distribution of dmCpGs by chromosomes
# =============================================================================
mean_hyper = np.mean(concate_hyper_cvt["Proportion"])
std_hyper = np.std(concate_hyper_cvt["Proportion"], ddof=1)

mean_hypo = np.mean(concate_hypo_cvt["Proportion"])
std_hypo = np.std(concate_hypo_cvt["Proportion"], ddof=1)

mean_chr = np.mean(concate_chr_cvt["Proportion"])
std_chr = np.std(concate_chr_cvt["Proportion"], ddof=1)

x_tk = np.arange(23)
# Get the values of index by sepcific level
x_tk_label = concate_chr_cvt["CHR"].values
fig = plt.figure()
ax = fig.add_subplot(211)
ax.bar(x_tk-0.2, concate_hyper_cvt["CHR_sig"], width=0.2, label="Hypermethylated CpGs")
ax.bar(x_tk, concate_hypo_cvt["CHR_sig"], width=0.2, label="Hypomethylated CpGs")
ax.bar(x_tk+0.2, concate_chr_cvt["CHR_sig"], width=0.2, label="All detected CpGs")
ax.legend(loc='best')
# Need to specify xticks, otherwise the program will not know where to put xticklabels
ax.set_xticks(x_tk)
ax.set_xticklabels(x_tk_label)
ax.set_xlabel('Chromosome')
ax.set_ylabel('No. of dmCpGs')
ax.set_yscale('log', basey=2)
ax.text(-0.05, 1.05, "A.", ha='center', va='center', transform=ax.transAxes)
ax.set_title("Distribution of dmCpGs in chromosomes")

ax2 = fig.add_subplot(212)
ax2.bar(x_tk-0.2, concate_hyper_cvt["Proportion"], width=0.2, label="Hypermethylated CpGs")
ax2.bar(x_tk, concate_hypo_cvt["Proportion"], width=0.2, label="Hypomethylated CpGs")
ax2.bar(x_tk+0.2, concate_chr_cvt["Proportion"], width=0.2, label="All detected CpGs")
ax2.legend(loc='best')
ax2.set_xticks(x_tk)
ax2.set_xticklabels(x_tk_label)
ax2.set_xlabel('Chromosome')
ax2.set_ylabel('Propertion of dmCpGs (%)')
ax2.text(-0.05, 1.05, "B.", ha='center', va='center', transform=ax2.transAxes)
ax2.text(
    0.16, 0.96, 
    "Mean ± SD:",
    ha='left', va='top', transform=ax2.transAxes
)
ax2.text(
    0.22, 0.96, 
    "Hypermethylated: {:.2f} ± {:.2f}\nHypomethylated: {:.2f} ± {:.2f}\nAll detected: {:.2f} ± {:.2f}".format(
        mean_hyper, std_hyper, mean_hypo, std_hypo, mean_chr, std_chr
    ),
    ha='left', va='top', transform=ax2.transAxes
)
ax2.set_title("Distribution of proportion of dmCpGs in chromosomes")
plt.subplots_adjust(hspace=0.3)
plt.show()