import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer
import matplotlib.pyplot as plt 
from matplotlib_venn import venn2 

sig_hyper = pd.read_excel(
    r'\sumdata\cg_S_up.xlsx',
    usecols="B:F", names=["CHR", "CpG Region", "Gene Region", "Gene"],
    index_col="Probeset ID"
)

sig_hypo = pd.read_excel(
    r'\sumdata\cg_S_down.xlsx',
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
all_hyper = all[['Gene', 'Gene Region']].loc[all['T_value'] > 0]
all_hypo = all[['Gene', 'Gene Region']].loc[all['T_value'] < 0]

sig_hyper_dp = sig_hyper.dropna(axis=0, subset=['Gene Region'])
sig_hyper_dp['CpG Region'].fillna("OpenSea", inplace=True)
sig_hypo_dp = sig_hypo.dropna(axis=0, subset=['Gene Region'])
sig_hypo_dp['CpG Region'].fillna("OpenSea", inplace=True)

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

gr_list = ["TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR"]
gn_dic = {}
df_dic = {}
for i in gr_list:
    gn_dic["hyper_{}".format(i)], df_dic["hyper_{}".format(i)] = gn_extract(DataSet=sig_hyper_dp, GeneRegion=i)
    gn_dic["hypo_{}".format(i)], df_dic["hypo_{}".format(i)] = gn_extract(DataSet=sig_hypo_dp, GeneRegion=i)

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

# Venn plot for gene dis
venn2(subsets = (2143, 1634, 279), set_labels = ('Gene with hypermethylation', 'Gene with hypomethylation'))
plt.title("Distribution of annotated genes", fontsize=20)
plt.show() 
