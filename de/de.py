from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pandas as pd
from sanbomics.tools import id_map
from sanbomics.plots import volcano
import gseapy as gp
from gseapy.plot import gseaplot
import numpy as np
import seaborn as sns
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import textwrap

counts = pd.read_csv("test.csv")
counts = counts.set_index('Geneid')
counts = counts[counts.sum(axis=1)>0]
counts = counts.T

metadata = pd.DataFrame(zip(counts.index, ['C', 'C', 'C', 'C', 'RS', 'RS', 'RS', 'RS']), columns= ['Sample', 'Condition'])
metadata = metadata.set_index('Sample')

dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors='Condition')
dds.deseq2()
stat_res = DeseqStats(dds, n_cpus=8, contrast=('Condition', 'RS', 'C'))
stat_res.summary()
res = stat_res.results_df
mapper = id_map(species='human')
res['Symbol'] = res.index.map(mapper.mapper)
res = res[res.baseMean >= 10]
cols_res = res.columns.tolist()
cols_res = cols_res[-1:] + cols_res[:-1]
ordered_res = res[cols_res]
ordered_res.to_csv('result_table.tsv', sep="\t", index=False)
sigs = res[(res.padj < 0.05) & (abs(res.log2FoldChange) > 0.5)]
print(sigs)
ordered_sigs = sigs[cols_res]
ordered_sigs.to_csv('significant_de.tsv', sep='\t', index=False)
sc.tl.pca(dds)
sc.pl.pca(dds, color = 'Condition', size = 250)

ranking = res[['Symbol', 'stat']].dropna().sort_values('stat', ascending= False)
ranking = ranking.drop_duplicates('Symbol')

pre_res = gp.prerank(rnk= ranking, gene_sets= ['GO_Biological_Process_2021'],
                     seed=6, permutation_num= 200)
out = []

for term in list(pre_res.results):
    out.append([term,
               pre_res.results[term]['fdr'],
               pre_res.results[term]['es'],
               pre_res.results[term]['nes']])

out_df = pd.DataFrame(out, columns = ['Term','fdr', 'es', 'nes']).sort_values('fdr').reset_index(drop = True)
out_df['Term'] = out_df['Term'].str.replace("GO_Biological_Process_2021__", "")

out_df.to_csv("go.tsv", sep='\t', index=False)
#name = out_df.sort_values('nes').iloc[0].Term
#gseaplot(pre_res.ranking, term = 'GO_Biological_Process_2021__mitotic spindle organization (GO:0007052)', **pre_res.results['GO_Biological_Process_2021__mitotic spindle organization (GO:0007052)'])

print(out_df)

dds.layers['log1p'] = np.log1p(dds.layers['normed_counts'])
dds_sigs = dds[:, sigs.index]
grapher = pd.DataFrame(dds_sigs.layers['log1p'].T,
                       index=dds_sigs.var_names, columns=dds_sigs.obs_names)
plot = sns.clustermap(grapher, z_score=0, cmap = 'Blues')
plt.savefig("plot.svg", format='svg')

volcano(res, symbol='Symbol')
plt.savefig("volcano.svg", format='svg')

out_df = out_df.sort_values('nes')
selected_go = out_df[0:10]
sns.set(style="whitegrid", color_codes=True)
plt.figure(figsize = (16,14))
pal = sns.color_palette("rocket", len(selected_go))

ax = sns.barplot(data = selected_go, x = 'nes', y = 'Term', hue= 'Term', legend=False, palette=pal)
ax.set_ylabel("Term", fontsize = 10)
ax.set_xlabel("Normalized Enrichment", fontsize = 10)
ax.set_yticklabels([textwrap.fill(e, 20) for e in selected_go['Term']])
plt.savefig("go_terms.svg", format='svg')
plt.show()

selected_go = out_df.iloc[-10:]
sns.set(style="whitegrid", color_codes=True)
plt.figure(figsize = (16,14))
pal = sns.color_palette("rocket", len(selected_go))

ax = sns.barplot(data = selected_go, x = 'nes', y = 'Term', hue= 'Term', legend=False, palette=pal)
ax.set_ylabel("Term", fontsize = 10)
ax.set_xlabel("Normalized Enrichment", fontsize = 10)
ax.set_yticklabels([textwrap.fill(e, 20) for e in selected_go['Term']])
plt.savefig("go_terms_2.svg", format='svg')
plt.show()