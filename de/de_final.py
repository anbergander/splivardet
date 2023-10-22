import os
import shutil

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

def createResultTable(dds, contrast1, contrast2, path):
    stat_res = DeseqStats(dds, n_cpus=8, contrast=('Condition', contrast1, contrast2))
    stat_res.summary()
    res = stat_res.results_df
    mapper = id_map(species='human')
    res['Gene'] = res.index.map(mapper.mapper)
    res = res[res.baseMean >= 10]
    cols_res = res.columns.tolist()
    cols_res = cols_res[-1:] + cols_res[:-1]
    ordered_res = res[cols_res]
    ordered_res.to_csv(path + '/tables/result_table.tsv', sep="\t", index=False)
    return res

def identifySignificant(res, path):
    cols_res = res.columns.tolist()
    cols_res = cols_res[-1:] + cols_res[:-1]
    sigs = res[(res.padj < 0.05) & (abs(res.log2FoldChange) > 0.5)]
    ordered_sigs = sigs[cols_res]
    ordered_sigs.to_csv(path + '/tables/significant_de.tsv', sep='\t', index=False)
    return sigs
def performPCA(dds, path):
    sc.pp.normalize_total(dds)
    sc.tl.pca(dds)
    sc.pl.pca(dds, color='Condition', size=250, color_map="magma", save="_normalized.svg")
    shutil.move("figures/pca_normalized.svg", path + "/figures/pca_normalized.svg")
    os.rmdir("/home/ubuntu/splivardet/gui/figures")

def heatmapDE(dds, sigs, path):
    dds.layers['log1p'] = np.log1p(dds.layers['normed_counts'])
    dds_sigs = dds[:, sigs.index]
    print(dds_sigs.obs_names)
    grapher = pd.DataFrame(dds_sigs.layers['log1p'].T,
                           index=dds_sigs.var_names, columns=dds_sigs.obs_names)
    sns.clustermap(grapher, z_score=0, cmap='Blues')
    plt.savefig(path + "/figures/heatmap_de.svg", format='svg')

def volcanoPlot(res, path):
    volcano(res, symbol='Gene', colors=['blue', 'lightgrey', 'purple'], save=path + "/figures/volcano", to_label=5, top_right_frame=True)

def createDir(path):
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)
        print("The new directory is created!")
def main(data, metadata, contrast1, contrast2, path):
    createDir(path + "/figures")
    createDir(path + "/tables")
    print("Start DE")
    counts = readCounts(data)
    metadata = createMetadata(counts, metadata)
    dds = performDe(counts, metadata)
    res = createResultTable(dds, contrast1,contrast2, path)
    sigs = identifySignificant(res, path)
    performPCA(dds, path)
    heatmapDE(dds, sigs, path)
    volcanoPlot(res, path)
    print("Finished DE")