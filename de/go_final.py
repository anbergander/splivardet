import matplotlib.pyplot as plt
import pandas as pd
import gseapy as gp
from gseapy.plot import gseaplot
import seaborn as sns
import textwrap

from human_genes import GENEID2NT as GeneID2nt_human
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS


def rankingGO(file):
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
                     figsize=(3, 4), ofname="figures/go_gsea.svg")
        out_df['Term'] = out_df['Term'].str.replace("GO_Biological_Process_2021__", "")
        out_df.to_csv("tables/go_gsea.tsv", sep='\t', index=False)
        graphGsea(out_df)
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
                     figsize=(3, 4), ofname="figures/go_gsea.svg")
        out_df['Term'] = out_df['Term'].str.replace("GO_Biological_Process_2021__", "")
        out_df.to_csv("tables/go_gsea.tsv", sep='\t', index=False)
        graphGsea(out_df)

def graphGsea(out_df):
    out_df = out_df.sort_values('nes')
    selected_go = out_df[0:10]
    sns.set(style="whitegrid", color_codes=True)
    plt.figure(figsize=(18, 14))
    pal = sns.color_palette("rocket", len(selected_go))
    ax = sns.barplot(data=selected_go, x='nes', y='Term', hue='Term', legend=False, palette=pal)
    ax.set_ylabel("Term", fontsize=10)
    ax.set_xlabel("Normalized Enrichment", fontsize=10)
    ax.set_title("Most downregulated terms")
    ax.set_yticklabels([textwrap.fill(e, 20) for e in selected_go['Term']])
    plt.savefig("figures/go_terms_decreased.svg", format='svg')
    plt.show()

    selected_go = out_df.iloc[-10:]
    sns.set(style="whitegrid", color_codes=True)
    plt.figure(figsize=(18, 14))
    pal = sns.color_palette("rocket", len(selected_go))
    ax = sns.barplot(data=selected_go, x='nes', y='Term', hue='Term', legend=False, palette=pal)
    ax.set_ylabel("Term", fontsize=10)
    ax.set_xlabel("Normalized Enrichment", fontsize=10)
    ax.set_title("Most upregulated terms")
    ax.set_yticklabels([textwrap.fill(e, 20) for e in selected_go['Term']])
    plt.savefig("figures/go_terms_enriched.svg", format='svg')
    plt.show()

def go_it(test_genes, mapper, goeaobj, inv_map, GO_items):
    mapped_genes = []
    for gene in test_genes:
        try:
                mapped_genes.append(mapper[gene])
        except:
             pass
    goea_results_all = goeaobj.run_study(mapped_genes)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
    GO = pd.DataFrame(list(map(lambda x: [x.GO, x.goterm.name, x.goterm.namespace, x.p_uncorrected, x.p_fdr_bh,
                                           x.ratio_in_study[0], x.ratio_in_study[1], GO_items.count(x.GO),
                                          list(map(lambda y: inv_map[y], x.study_items)),
                                            ], goea_results_sig)),
                       columns=['GO', 'term', 'class', 'p', 'p_corr', 'n_genes',
                                'n_study', 'n_go', 'study_genes'])
    GO = GO[GO.n_genes > 1]
    return GO

def goaGO(file):
    if ".csv" in file:
        df = pd.read_csv(file)
        df = df[(df.padj < 0.05) & (df.log2FoldChange > 0.5)]
        obo_fname = download_go_basic_obo()
        fin_gene2go = download_ncbi_associations()
        obodag = GODag("go-basic.obo")
        mapper = {}
        for key in GeneID2nt_human:
            mapper[GeneID2nt_human[key].Symbol] = GeneID2nt_human[key].GeneID
        inv_map = {v: k for k, v in mapper.items()}
        objanno = Gene2GoReader(fin_gene2go, taxids=[9606])
        ns2assoc = objanno.get_ns2assc()
        goeaobj = GOEnrichmentStudyNS(
            GeneID2nt_human.keys(),
            ns2assoc,
            obodag,
            propagate_counts=False,
            alpha=0.05,
            methods=['fdr_bh'])
        GO_items = []
        temp = goeaobj.ns2objgoea['BP'].assoc
        for item in temp:
            GO_items += temp[item]
        temp = goeaobj.ns2objgoea['CC'].assoc
        for item in temp:
            GO_items += temp[item]
        temp = goeaobj.ns2objgoea['MF'].assoc
        for item in temp:
            GO_items += temp[item]
        return go_it(df.Gene.values, mapper, goeaobj, inv_map, GO_items)
    else:
        df = pd.read_table(file)
        df = df[(df.padj < 0.05) & (df.log2FoldChange > 0.5)]
        obo_fname = download_go_basic_obo()
        fin_gene2go = download_ncbi_associations()
        obodag = GODag("go-basic.obo")
        mapper ={}
        for key in GeneID2nt_human:
            mapper[GeneID2nt_human[key].Symbol] = GeneID2nt_human[key].GeneID
        inv_map = {v: k for k, v in mapper.items()}
        objanno = Gene2GoReader(fin_gene2go, taxids=[9606])
        ns2assoc = objanno.get_ns2assc()
        goeaobj = GOEnrichmentStudyNS(
            GeneID2nt_human.keys(),
            ns2assoc,
            obodag,
            propagate_counts=False,
            alpha=0.05,
            methods=['fdr_bh'])
        GO_items = []
        temp = goeaobj.ns2objgoea['BP'].assoc
        for item in temp:
            GO_items += temp[item]
        temp = goeaobj.ns2objgoea['CC'].assoc
        for item in temp:
            GO_items += temp[item]
        temp = goeaobj.ns2objgoea['MF'].assoc
        for item in temp:
            GO_items += temp[item]
        return go_it(df.Gene.values, mapper, goeaobj, inv_map, GO_items)

def goatable(res):
    res['per'] = res.n_genes / res.n_go
    res.to_csv("tables/go_goa.tsv", sep='\t', index=False)

def goagraph(res):
    res = res.sort_values('p_corr')
    out = res[0:10]
    sns.set(style="whitegrid", color_codes=True)
    plt.figure(figsize=(16, 14))
    pal = sns.color_palette("rocket", len(out))
    ax = sns.barplot(data=out, x='per', y='term', hue='term', legend=False, palette=pal)
    ax.set_title("10 Most upregulated GO terms", fontsize=14)
    ax.set_ylabel("Term", fontsize=10)
    ax.set_xlabel("Genes/GO", fontsize=10)
    ax.set_yticklabels([textwrap.fill(e, 15) for e in out['term']])
    plt.savefig("figures/go_terms_goa.svg", format='svg')
    plt.show()


rankingGO("tables/result_table.tsv")
res = goaGO("tables/result_table.tsv")
goatable(res)
goagraph(res)

