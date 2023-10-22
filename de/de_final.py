from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pandas as pd
from sanbomics.tools import id_map
from sanbomics.plots import volcano
import numpy as np
import seaborn as sns
import scanpy as sc
import matplotlib.pyplot as plt

def readCounts(path):
    if ".csv" in path:
        counts = pd.read_csv(path)
        counts = counts.set_index('Geneid')
        counts = counts[counts.sum(axis=1) > 0]
        counts = counts.T
        return counts
    else:
        counts = pd.read_table(path)
        counts = counts.set_index('Geneid')
        counts = counts[counts.sum(axis=1) > 0]
        counts = counts.T
        return counts

def createMetadata(counts, list):
    metadata = pd.DataFrame(zip(counts.index, list),
                            columns=['Sample', 'Condition'])
    metadata = metadata.set_index('Sample')
    return metadata

def performDe(counts, metadata):
    dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors='Condition')
    dds.deseq2()
    return dds

def createResultTable(dds, contrast1, contrast2):
    stat_res = DeseqStats(dds, n_cpus=8, contrast=('Condition', contrast1, contrast2))
    stat_res.summary()
    res = stat_res.results_df
    mapper = id_map(species='human')
    res['Gene'] = res.index.map(mapper.mapper)
    res = res[res.baseMean >= 10]
    cols_res = res.columns.tolist()
    cols_res = cols_res[-1:] + cols_res[:-1]
    ordered_res = res[cols_res]
    ordered_res.to_csv('tables/result_table.tsv', sep="\t", index=False)
    return res

def identifySignificant(res):
    cols_res = res.columns.tolist()
    cols_res = cols_res[-1:] + cols_res[:-1]
    sigs = res[(res.padj < 0.05) & (abs(res.log2FoldChange) > 0.5)]
    ordered_sigs = sigs[cols_res]
    ordered_sigs.to_csv('tables/significant_de.tsv', sep='\t', index=False)
    return sigs
def performPCA(dds):
    sc.pp.normalize_total(dds)
    sc.tl.pca(dds)
    sc.pl.pca(dds, color='Condition', size=250, color_map="magma", save="_normalized.svg")

def heatmapDE(dds, sigs):
    dds.layers['log1p'] = np.log1p(dds.layers['normed_counts'])
    dds_sigs = dds[:, sigs.index]
    print(dds_sigs.obs_names)
    grapher = pd.DataFrame(dds_sigs.layers['log1p'].T,
                           index=dds_sigs.var_names, columns=dds_sigs.obs_names)
    sns.clustermap(grapher, z_score=0, cmap='Blues')
    plt.savefig("figures/heatmap_de.svg", format='svg')

def volcanoPlot(res):
    volcano(res, symbol='Gene', colors=['blue', 'lightgrey', 'purple'], save="figures/volcano", to_label=5, top_right_frame=True)


counts = readCounts("file.tsv")
metadata = createMetadata(counts, ['V', 'V', '0', 'V', '0', '0'])
dds = performDe(counts, metadata)
res = createResultTable(dds, '0', 'V')
sigs = identifySignificant(res)
performPCA(dds)
heatmapDE(dds, sigs)
volcanoPlot(res)