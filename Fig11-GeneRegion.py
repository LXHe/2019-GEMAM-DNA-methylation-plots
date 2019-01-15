import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

path = r'\2018-GEMAM-1st writing\methyl_plt.xlsx'
gr = pd.read_excel(path,index_col=[0,1])

gr_usk = gr.unstack(level='Methylation')
gr_usk_dp = gr_usk.columns.droplevel(0)
gr_usk.columns = gr_usk_dp

# Display values in bar chart
def autolabel(bar):
    for a in bar:
        height = a.get_height()
        ax.text(
                a.get_x()+a.get_width()/2, height+0.5,
                '{}'.format(height), ha='center', va='bottom'
                )

gr_hyper = gr.loc['Hypermethylated']['Count']
gr_hypo = gr.loc['Hypomethylated']['Count']

tick = np.arange(6)
fig = plt.figure()
ax = fig.add_subplot(111)
hyper_bar = ax.bar(tick-0.2, gr_hyper, width=0.4, label="Hypermethylated")
hypo_bar = ax.bar(tick+0.2, gr_hypo, width=0.4, label="Hypomethylated")
ax.set_yscale('log', basey=10)
ax.set_xticks(tick)
ax.set_xticklabels(gr_hyper.index.values)
ax.legend(loc='best')
ax.set_xlabel("Gene region", fontsize=11)
ax.set_ylabel("Gene number", fontsize=11)
plt.title("Distribution of methylated genes by gene regions", fontsize=16)
autolabel(hyper_bar)
autolabel(hypo_bar)
plt.show()
