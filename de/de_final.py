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

def createMetadata(counts, list, path):
    metadata = pd.DataFrame(zip(counts.index, list),
                            columns=['Sample', 'Condition'])
    metadata = metadata.set_index('Sample')
    metadata.to_csv(path + '/tables/metadata.tsv', sep="\t", index=False)
    return metadata

def readMetadataTable(path):
    metadata = pd.read_table(path)
    metadata = metadata.set_index('Sample')
    print(metadata)
    return metadata

def performDe(counts, metadata, factors):
    dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors=factors)
    dds.deseq2()
    return dds

def createResultTable(dds, factor, contrast1, contrast2, path):
    stat_res = DeseqStats(dds, n_cpus=8, contrast=(factor, contrast1, contrast2))
    stat_res.summary()
    stat_res.plot_MA(save_path=path + "/figures/ma_plot_" + factor + "_" + contrast1 + "_" + contrast2 + ".svg")
    res = stat_res.results_df
    mapper = id_map(species='human')
    res.to_csv(path + '/tables/result_table_ensembleIds' + "_" + factor + "_" + contrast1 + "_" + contrast2 + ".tsv", sep="\t", index=True)
    sigs = res[(res.padj < 0.05) & (abs(res.log2FoldChange) > 0.5)]
    sigs.to_csv(path + '/tables/significant_de_ensembleIds' + "_" + factor + "_" + contrast1 + "_" + contrast2 + ".tsv", sep='\t', index=True)
    res['Gene'] = res.index.map(mapper.mapper)
    res = res[res.baseMean >= 10]
    cols_res = res.columns.tolist()
    cols_res = cols_res[-1:] + cols_res[:-1]
    ordered_res = res[cols_res]
    ordered_res.to_csv(path + '/tables/result_table' + "_" + factor + "_" + contrast1 + "_" + contrast2 + ".tsv", sep="\t", index=False)
    return res

def identifySignificant(res, path, factor, contrast1, contrast2):
    cols_res = res.columns.tolist()
    cols_res = cols_res[-1:] + cols_res[:-1]
    significant = res[(res.padj < 0.05) & (abs(res.log2FoldChange) > 0.5)]
    ordered_significant = significant[cols_res]
    ordered_significant.to_csv(path + '/tables/significant_de' + "_" + factor + "_" + contrast1 + "_" + contrast2 + ".tsv", sep='\t', index=False)
    return significant
def performPCA(dds, path, factor):
    sc.pp.normalize_total(dds)
    sc.tl.pca(dds)
    sc.pl.pca(dds, color=factor, size=250, color_map="magma", save="_" +factor+ "_normalized.svg")
    shutil.move("figures/pca_" +factor+"_normalized.svg", path + "/figures/pca_" +factor+"_normalized.svg")
    #os.rmdir("/home/ubuntu/splivardet/gui/figures")

def heatmapDE(dds, sigs, path, factor, contrast1, contrast2):
    dds.layers['log1p'] = np.log1p(dds.layers['normed_counts'])
    dds_sigs = dds[:, sigs.index]
    grapher = pd.DataFrame(dds_sigs.layers['log1p'].T,
                           index=dds_sigs.var_names, columns=dds_sigs.obs_names)
    sns.clustermap(grapher, z_score=0, cmap='magma_r', cbar_pos=(0.02, 0.83, 0.03, 0.15))
    plt.savefig(path + "/figures/heatmap_de_"+ factor + "_"+contrast1 + "_"+ contrast2+".svg", format='svg')

def volcanoPlot(res, path, genes, factor, contrast1, contrast2):
    volcano(res, symbol='Gene', colors=['blue', 'lightgrey', 'purple'], save=path + "/figures/volcano_" + factor
                                                                             + "_" + contrast1 + "_" + contrast2
            , to_label=genes, top_right_frame=True)

def createDir(path):
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)
        print("The new directory is created!")
def main(data, metadata, path, contrast1, contrast2, factors, genes):
    createDir(path + "/figures")
    createDir(path + "/tables")
    print("Start DE")
    counts = readCounts(data)
    print(counts)
    metadata = createMetadata(counts, metadata, path)
    dds = performDe(counts, metadata, factors)
    print(dds)
    res = createResultTable(dds, contrast1,contrast2, factors,path)
    significant = identifySignificant(res, path, factors, contrast1, contrast2)
    performPCA(dds, path, factors)
    heatmapDE(dds, significant, path, factors, contrast1, contrast2)
    volcanoPlot(res, path, genes, factors, contrast1, contrast2)
    print("Finished DE")

def main2(data, metadatapath, path, genes):
    createDir(path + "/figures")
    createDir(path + "/tables")
    counts = readCounts(data)
    metadata = readMetadataTable(metadatapath)
    factors = []
    for col in metadata:
        factors.append(col)
    dds = performDe(counts,metadata,factors)
    collection = []
    f = open(path + "/summary.txt", "w")
    f.writelines("Summary for DEG")
    f.write("\n")
    for col in metadata:
        performPCA(dds, path, col)
        for elem in metadata[col]:
            collection.append(elem)
        setList = set(collection)
        conditions = list(setList)
        for item in conditions:
            for items in conditions:
                if not items == item:
                    print(col)
                    fac1 = item
                    print(fac1)
                    fac2 = items
                    print(fac2)
                    res = createResultTable(dds, col, fac1, fac2, path)
                    significant = identifySignificant(res, path, col, fac1, fac2)
                    if significant.empty:
                        f.writelines(
                            "No differentially expressed genes have been identified for " + col + " with conditions: " + fac1 + " and " + fac2)
                        f.write("\n")
                    else:
                        heatmapDE(dds, significant, path, col, fac1, fac2)
                        volcanoPlot(res, path, genes, col, fac1, fac2)
                        f.writelines(
                            "All files generated for " + col + " with conditions: " + fac1 + " and " + fac2)
                        f.write("\n")
        collection = []
    f.close()