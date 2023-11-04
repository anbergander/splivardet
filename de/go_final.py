import os

import matplotlib.pyplot as plt
import pandas as pd
import gseapy as gp
import seaborn as sns
import textwrap


def rankingGO(file, path):
    if ".csv" in file:
        res = pd.read_csv(file)
        ranking = res[['Gene', 'stat']].dropna().sort_values('stat', ascending=False)
        ranking = ranking.drop_duplicates('Gene')

        pre_res = gp.prerank(rnk=ranking, gene_sets=['GO_Biological_Process_2021'],
                             seed=6, permutation_num=200)
        out = []
        for term in list(pre_res.results):
            out.append([term,
                        pre_res.results[term]['fdr'],
                        pre_res.results[term]['es'],
                        pre_res.results[term]['nes']])

        out_df = pd.DataFrame(out, columns=['Term', 'fdr', 'es', 'nes']).sort_values('fdr').reset_index(drop=True)
        pre_res.plot(terms=pre_res.res2d.Term[:5],
                     show_ranking=True,
                     figsize=(3, 4), ofname=path+"/figures/gsea_enrichment_plot.svg")
        out_df['Term'] = out_df['Term'].str.replace("GO_Biological_Process_2021__", "")
        out_df.to_csv(path + "/tables/gsea.tsv", sep='\t', index=False)
        graphGsea(out_df, path)
    else:
        res = pd.read_table(file)
        ranking = res[['Gene', 'stat']].dropna().sort_values('stat', ascending=False)
        ranking = ranking.drop_duplicates('Gene')

        pre_res = gp.prerank(rnk=ranking, gene_sets=['GO_Biological_Process_2021'],
                             seed=6, permutation_num=200)
        out = []

        for term in list(pre_res.results):
            out.append([term,
                        pre_res.results[term]['fdr'],
                        pre_res.results[term]['es'],
                        pre_res.results[term]['nes']])

        out_df = pd.DataFrame(out, columns=['Term', 'fdr', 'es', 'nes']).sort_values('fdr').reset_index(drop=True)
        pre_res.plot(terms=pre_res.res2d.Term[:5],
                     show_ranking=True,
                     figsize=(3, 4), ofname=path + "/figures/gsea_enrichment_plot.svg")
        out_df['Term'] = out_df['Term'].str.replace("GO_Biological_Process_2021__", "")
        out_df.to_csv(path + "/tables/gsea.tsv", sep='\t', index=False)
        graphGsea(out_df, path)

def graphGsea(out_df, path):
    out_df = out_df.sort_values('nes')
    selected_go = out_df[0:10]
    sns.set(style="whitegrid", color_codes=True)
    plt.figure(figsize=(18, 14))
    pal = sns.color_palette("magma", len(selected_go))
    ax = sns.barplot(data=selected_go, x='nes', y='Term', hue='Term', legend=False, palette=pal)
    ax.set_ylabel("Term", fontsize=10)
    ax.set_xlabel("Normalized Enrichment", fontsize=10)
    ax.set_title("Most downregulated terms")
    ax.set_yticklabels([textwrap.fill(e, 20) for e in selected_go['Term']])
    plt.savefig(path + "/figures/gsea_terms_decreased.svg", format='svg')
    plt.show()

    selected_go = out_df.iloc[-10:]
    sns.set(style="whitegrid", color_codes=True)
    plt.figure(figsize=(18, 14))
    pal = sns.color_palette("magma", len(selected_go))
    ax = sns.barplot(data=selected_go, x='nes', y='Term', hue='Term', legend=False, palette=pal)
    ax.set_ylabel("Term", fontsize=10)
    ax.set_xlabel("Normalized Enrichment", fontsize=10)
    ax.set_title("Most upregulated terms")
    ax.set_yticklabels([textwrap.fill(e, 20) for e in selected_go['Term']])
    plt.savefig(path + "/figures/gsea_terms_enriched.svg", format='svg')
    plt.show()

def createDir(path):
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)
        print("The new directory is created!")
def main(path, data):
    createDir(path + "/figures")
    createDir(path + "/tables")
    print("Start GO")
    rankingGO(data, path)
    print("Finished GO")

