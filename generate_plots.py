#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 27 11:51:04 2023

@author: denbi
"""
import pickle
import seaborn as sns
from statannotations.Annotator import Annotator
import matplotlib.pyplot as plt
import warnings
import sys
import pandas as pd

warnings.filterwarnings("ignore")

wo = sys.argv[1]+'/'
#wo = '/vol/fastq/splivardet/'


uniprot = pd.read_csv(wo+"resources/uniprot.tsv", sep= "\t")

up = dict(zip(uniprot['Entry'],uniprot['Gene Names']))

file = open(wo+"resources/refseq_ensembl.txt")
annot = file.read()
file.close()

file = open(wo+"resources/transcript_to_gene.txt")
trans_to_gene = file.read()
file.close()

file = open(wo+"resources/ens_gensym.txt")
ens_gensym = file.read()
file.close()




a_d = {}
for n in annot.split("\n"):
    s = n.split("\t")
    if len(s)==2:
        a_d[s[1]] = s[0]
    
t_g = {}
for n in trans_to_gene.split("\n"):
    s = n.split("\t")
    if len(s)==2:
        t_g[s[0]] = s[1]

e_s = {}
s_e = {}
for n in ens_gensym.split("\n"):
    s = n.split("\t")
    if len(s)==2:
        e_s[s[1]] = s[0]
        s_e[s[0]] = s[1]


def gene_name(gene):

    try:
        gn = s_e[gene]
    except: 
        gn = gene

    return gn

def flair(t):
    if t.startswith('ENS'):
        t = t.split('.')[0]
    else:
        t = '0x'+t.replace('-','').upper()

        t = 'FLAIR'+str(int(t, 16))
    return(t)
    
def iso_name(ik):
    i_n = ik.split('_')
    i = i_n[0]
    n = i_n[-1]
    return flair(i)+'\n'+gene_name(n)

file = open('posthoc.pickle','rb')
post_hoc = pickle.load(file)
file.close()
file = open('dataframes_melt.pickle','rb')
data = pickle.load(file)
file.close()

def generate_violin_plot(post_hoc, data, key):

    molten_df = post_hoc.melt(ignore_index=False)
    
    ax = sns.violinplot(data=data, x=data["cohort"], y=data["portion of transcripts"], showfliers=True, palette = 'pastel')
    sns.stripplot(data=data, x=data["cohort"], y=data["portion of transcripts"], color='black', alpha=0.5)
    molten_df.index.name = 'index'
    molten_df = molten_df.reset_index()
    
    pairs = [(i[1]['index'], i[1]["variable"]) for i in molten_df.iterrows()]
    
    p_values = [i[1]["value"] for i in molten_df.iterrows()]
    
    pairs_pvalues = []
    p_gr = []
    p_v =[]
    
    for p in zip(pairs, p_values):
        gr, pv = p
        gr1, gr2 = gr
        if gr1 < gr2:
            gr = (gr2, gr1)
        if gr1 != gr2:
            p_n = (gr,pv)
            pairs_pvalues.append(p_n)
                
    pairs_pvalues = list(set(pairs_pvalues))
                
    for p in pairs_pvalues:
        gr, pv = p
        p_gr.append(gr)
        p_v.append(pv)
                    
    pairs = p_gr
    p_values = p_v
                    

    annotator = Annotator.get_empty_annotator()
    annotator = Annotator(
            ax, pairs, data=data, x="cohort", y="portion of transcripts")
    annotator.configure(text_format="star", loc="inside")
    annotator.set_pvalues_and_annotate(p_values)
    #plt.gca().set_ylim(bottom=0)
    plt.title(iso_name(key), fontsize=9)
    plt.ylabel('transformed portion of respective transcripts')
    #plt.tight_layout()
    fn = iso_name(key).replace('\n', '_')
    print(fn)
    plt.savefig(fn+".violin.pdf", format="pdf")
    plt.clf()

with open('candidates.txt', 'a+') as f:
    for key in post_hoc.keys():
        data_df = data[key]
        post_hoc_df = post_hoc[key]
        data_df.rename(columns={'group':'cohort', 'value':'portion of transcripts'}, inplace=True)
        generate_violin_plot(post_hoc_df, data_df, key)
        f.write(key+'\n')
