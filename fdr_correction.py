#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri Aug 25 17:22:47 2023

@author: denbi
"""
import pickle
import pandas as pd
from statsmodels.stats.multitest import multipletests
import numpy as np
import sys



wo = sys.argv[1]+'/'

#wo = "/vol/fastq/splivardet/"

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
        
        if gene.startswith('ENST'):
            gn = t_g[gene.split('.')[0]]
        elif ':' in gene:
            gn = gene
        else:
            gn = e_s[str(up[gene]).split(' ')[0]]
    except: 
        gn = 'nan'

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
    return flair(i)+' '+gene_name(n)


#file = open('dataframes_melt_raw.pickle','rb')
#dataframes_melt_raw = pickle.load(file)
file.close()
file = open('dataframes_melt.pickle','rb')
dataframes_melt = pickle.load(file)
file.close()
posthoc_df = pd.read_csv('posthoc.csv')

posthoc_df[['p.value']] = posthoc_df[['p.value']].apply(pd.to_numeric)

#print(len(posthoc_df.loc[posthoc_df['p.value'] < 0.05]))
p_value = posthoc_df['p.value']
isoform_keys = np.array(posthoc_df['isoform_key'])

key_list = list(dataframes_melt.keys())

#print(key)
#print(dataframes_melt[key])
#print(dataframes_melt_raw[key])

a = multipletests(p_value, method='fdr_by', alpha=(0.05))
candidate_list = isoform_keys[a[0]]

if len(candidate_list)==0:
    min_p = min(posthoc_df['p.value'])
    candidate_list.append(posthoc_df[posthoc_df['p.value'] == min_p].isoform_key)



ph_df = pd.read_csv("posthoc.csv")



#create unique list of names
isoform_keys = ph_df.isoform_key.unique()

#create a data frame dictionary to store your data frames
DataFrameDict = {elem : pd.DataFrame() for elem in isoform_keys}
grouplist = []
gl_does_not_exist = True
for key in candidate_list:
   
    DataFrameDict[key] = ph_df[:][ph_df.isoform_key == key].drop("Unnamed: 0",axis=1)
    DataFrameDict[key] = DataFrameDict[key].reset_index().drop(columns=['estimate', 'SE', 'df', 'z.ratio', 'index'])
    DataFrameDict[key] = DataFrameDict[key].drop_duplicates()
    if gl_does_not_exist:
        for i in DataFrameDict[key]['contrast'].tolist():
            grouplist.extend(i.split(' - '))
            gl_does_not_exist = False

        
grouplist = set(grouplist)
candidate_list =list(set(candidate_list))
for key in candidate_list:
    dummy_values = [1]*len(grouplist)
    temp_df = pd.DataFrame(data={i:dummy_values for i in grouplist})
    temp_df['gl'] = list(grouplist)
    temp_df.set_index('gl', inplace = True)
    temp_df.index.name = None

    for i in list(DataFrameDict[key][['contrast', 'p.value']].itertuples(index=False, name=None)):
        gg, pv = i
        groups = gg.split(" - ")
        temp_df.at[groups[0], groups[1]] = pv
        temp_df.at[groups[1], groups[0]] = pv
    DataFrameDict[key] = temp_df
keys = list(DataFrameDict.keys())
for key in keys:
    if key not in candidate_list:
        DataFrameDict.pop(key)
with open('posthoc.pickle', 'wb') as f:
    pickle.dump(DataFrameDict, f)

#print(DataFrameDict[isoform_keys[1]])
posthoc_result = pd.DataFrame()
posthoc_result['contrast'] = posthoc_df['contrast']
posthoc_result['pvalue'] = posthoc_df['p.value']
a = list(map(iso_name, list(posthoc_df['isoform_key'])))
posthoc_result['transcript'] = a
posthoc_result.set_index('transcript', inplace=True)
posthoc_result.to_csv('posthoc_result.tsv',sep = '\t')
