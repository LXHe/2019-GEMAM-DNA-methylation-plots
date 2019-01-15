import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer
import matplotlib.pyplot as plt 
from matplotlib_venn import venn2 

sig_hyper = pd.read_excel(r'\sumdata\cg_S_up.xlsx',
                         usecols="B:F", names=["CHR", "CpG Region", "Gene Region", "Gene"],
                         index_col="Probeset ID"
                         )

sig_hypo = pd.read_excel(r'\sumdata\cg_S_down.xlsx',
                         usecols="B:F", names=["CHR", "CpG Region", "Gene Region", "Gene"],
                         index_col="Probeset ID"
                         )

all_path = r'\BC_corrected_methylation_beta\Sar_vs_N_Sar_ttest.csv'
all = pd.read_csv(
    all_path, 
    usecols=['Probe', 'Gene_name', 'Gene_region', 'T_value'],
    index_col='Probe'
    )
all.rename(columns={'Gene_region':'Gene Region', 'Gene_name':'Gene'}, inplace=True)

all = all.dropna(axis=0, subset=['Gene Region'])
all_hyper = all[
    ['Gene', 'Gene Region']
    ].loc[all['T_value'] > 0]
all_hypo = all[
    ['Gene', 'Gene Region']
    ].loc[all['T_value'] < 0]
#%%
sig_hyper_dp = sig_hyper.dropna(axis=0, subset=['Gene Region'])
sig_hyper_dp['CpG Region'].fillna("OpenSea", inplace=True)
sig_hypo_dp = sig_hypo.dropna(axis=0, subset=['Gene Region'])
sig_hypo_dp['CpG Region'].fillna("OpenSea", inplace=True)
#%%
# =============================================================================
# Test section by "Body", don't use
# =============================================================================

sig_hypo_body = sig_hypo_dp[sig_hypo_dp['Gene Region'].str.match("Body", na=False)]
sig_hypo_body_gr = sig_hypo_body['Gene Region'].str.split(";", expand=True)
sig_hypo_body_gn = sig_hypo_body['Gene'].str.split(";", expand=True)
# Check the contents in each column for incomplete name
col_num = sig_hypo_body_gr.columns.values
unq_gr = []
for i in col_num:
    unq_gr.append(sig_hypo_body_gr[i].unique())
a = unq_gr
#%%
# (3,17) has incomplete value of "Bo", 2 observations have incomplte value of "Bod"
# Replace by "Body"
sig_hypo_body_gr[sig_hypo_body_gr == 'Bo'] = "Body"
sig_hypo_body_gr[sig_hypo_body_gr == 'Bod'] = "Body"
#%%
mask = sig_hypo_body_gr == "Body"

sig_hypo_body_gn_msk = sig_hypo_body_gn[mask]
#%%
# Remove duplicate column values
# It works well with DataFrames with only unique values in each columns
## test_list = np.random.randint(0,20,(4,3))
## df = pd.DataFrame(np.hstack([test_list,test_list]), columns=["N0", "N1", "N2", "N3", "N4", "N5"])
## df_dd = df.T.drop_duplicates().T

# If duplicate data exist in multiple columns, some of the columns will still have duplicate values
# Therefore a column-wise drop function is listed below  

trans = sig_hypo_body_gn_msk.T
df_add = pd.DataFrame()
for probe in sig_hypo_body_gn_msk.index.values:
    a = trans[[probe]].drop_duplicates(subset=[probe])
    df_add = df_add.append(a.T)

#%%
df_add['Gn_body'] = df_add.apply(
    lambda x: ",".join(x.dropna()), axis=1
    )
#%%
def my_tok(txt):
    return [x for x in txt.split(",")]

vectorizer = CountVectorizer(tokenizer=my_tok, lowercase=False)
X = vectorizer.fit_transform(df_add['Gn_body'])
gn_body = vectorizer.get_feature_names()

sig_hypo_body_final = pd.merge(sig_hypo_dp[["CHR","CpG Region"]], df_add[['Gn_body']],
                                left_index=True, right_index=True
                                )
#%%
# =============================================================================
# Formal section
# =============================================================================
def my_tok(txt):
    return [x for x in txt.split(",")]

def gn_extract(DataSet, GeneRegion):
    gr_select = DataSet[DataSet['Gene Region'].str.match(GeneRegion, na=False)]
    gr_select_spt = gr_select['Gene Region'].str.split(";", expand=True)
    gn_select_spt = gr_select['Gene'].str.split(";", expand=True)

    mask = gr_select_spt == GeneRegion
    gn_select_msk = gn_select_spt[mask]
    trans = gn_select_msk.T
    gn_dup = pd.DataFrame()
    for probe in gr_select.index.values:
        dup = trans[[probe]].drop_duplicates(subset=[probe])
        gn_dup = gn_dup.append(dup.T)
        
    gn_dup['Gene Name'] = gn_dup.apply(
        lambda x: ",".join(x.dropna()), axis=1
        )
    
    vectorizer = CountVectorizer(tokenizer=my_tok, lowercase=False)
    vectorizer.fit_transform(gn_dup['Gene Name'])
    gn = vectorizer.get_feature_names()

    df_final = pd.merge(
        DataSet[["CHR","CpG Region"]], 
        gn_dup[['Gene Name']],
        left_index=True, right_index=True
        )
    return gn, df_final
#%%
gr_list = ["TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR"]
gn_dic = {}
df_dic = {}
for i in gr_list:
    gn_dic["hyper_{}".format(i)], df_dic["hyper_{}".format(i)] = gn_extract(DataSet=sig_hyper_dp, GeneRegion=i)
    gn_dic["hypo_{}".format(i)], df_dic["hypo_{}".format(i)] = gn_extract(DataSet=sig_hypo_dp, GeneRegion=i)
#%%
# =============================================================================
# Special process. Because the lengths of gene and gene region don't match, need to manually adjust.
# hyper_body, cg15949044, PCDH
# hypo_body, cg10116456, DIS
# =============================================================================
gn_dic['hyper_Body'].remove('PCDH')
gn_dic['hypo_Body'].remove('DIS')
gn_dic['hypo_Body'].remove('')
#%%
writer = pd.ExcelWriter(r'C:\Users\u0105352\Desktop\Gene_name.xlsx', engine='xlsxwriter')

for key in df_dic:
    df_dic[key].to_excel(writer,sheet_name='{}'.format(key))
writer.save()
#%%
# Compare shared gene by gene region
keys = ["hyper_1stExon", "hyper_3'UTR", "hyper_5'UTR", "hyper_Body", "hyper_TSS1500", "hyper_TSS200", 
        "hypo_1stExon", "hypo_3'UTR", "hypo_5'UTR", "hypo_Body", "hypo_TSS1500", "hypo_TSS200"]
gn_shr_dic = {}
for i in range(12):
    for m in range(i+1, 12):
        gn_shr_dic["{} vs {}".format(keys[i], keys[m])] = list(
                set(gn_dic[keys[i]]) & set(gn_dic[keys[m]])
                )
# Compare shared gene by methylation
gn_hyper = []
gn_hypo = []
for i in range(6):
    gn_hyper = list(set(gn_hyper) | set(gn_dic[keys[i]]))
    gn_hypo = list(set(gn_hypo) | set(gn_dic[keys[i+6]]))
gn_all = list(set(gn_hyper) | set(gn_hypo))
gn_shr = list(set(gn_hyper) & set(gn_hypo))
#%%
# Venn plot for gene dis
venn2(subsets = (2143, 1634, 279), set_labels = ('Gene with hypermethylation', 'Gene with hypomethylation'))
plt.title("Distribution of annotated genes", fontsize=20)
plt.show() 
#%%
# This function has a problem with size, need to check. The function above gives perfect solution.
gn_shr_dic = {}
for key in gn_dic:
    for key2 in gn_dic:
        if key != key2:
            gn_shr_dic["{} vs {}".format(key,key2)] = list(
                set(gn_dic[key]) & set(gn_dic[key2])
                )
#%%
gr_list = ["TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR"]
all_hypercount_dic = {}
all_hypocount_dic = {}
sig_hypercount_dic = {}
sig_hypocount_dic = {}
for i in gr_list:
    all_hypercount_dic[i] = all_hyper['Gene Region'].str.contains(i).sum()
    all_hypocount_dic[i] = all_hypo['Gene Region'].str.contains(i).sum()
    sig_hypercount_dic[i] = sig_hyper_dp['Gene Region'].str.contains(i).sum()
    sig_hypocount_dic[i] = sig_hypo_dp['Gene Region'].str.contains(i).sum()
