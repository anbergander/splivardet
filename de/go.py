import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import seaborn as sns
import textwrap

from human_genes import GENEID2NT as GeneID2nt_human
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

df = pd.read_table("result_table.tsv")
df = df[(df.padj < 0.05) & (df.log2FoldChange > 0.5)]
print(df)
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
        GeneID2nt_human.keys(), # List of mouse protein-coding genes
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh'])

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


def go_it(test_genes):
        print(f'input genes: {len(test_genes)}')

        mapped_genes = []
        for gene in test_genes:
                try:
                        mapped_genes.append(mapper[gene])
                except:
                        pass
        print(f'mapped genes: {len(mapped_genes)}')

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

res = go_it(df.Symbol.values)
print(res)
res['per'] = res.n_genes/res.n_go
res.to_csv("go_goa.tsv", sep = '\t', index = False)

out = res[0:10]

print(res)

sns.set(style="whitegrid", color_codes=True)
plt.figure(figsize = (16,14))
pal = sns.color_palette("rocket", len(out))

ax = sns.barplot(data = out, x = 'per', y = 'term', hue= 'term', legend=False, palette=pal)
ax.set_ylabel("Term", fontsize = 10)
ax.set_xlabel("Genes/GO", fontsize = 10)
ax.set_yticklabels([textwrap.fill(e, 15) for e in out['term']])
plt.savefig("go_terms.svg", format='svg')
plt.show()