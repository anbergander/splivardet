#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 10:07:09 2023

@author: denbi
"""
import numpy as np
import csv
import pandas as pd
from gtfparse import read_gtf
import sys

import warnings
warnings.filterwarnings("ignore")


isofile = sys.argv[1]
#isofile = '/vol/fastq/bc2/barcode02.isoforms.gtf'

pjd = sys.argv[2]+'/'
indp = sys.argv[3]
cand_f = sys.argv[4]

 
#pjd='/vol/fastq/splivardet/' 
#indp='/vol/fastq/index/'
#cand_f='/vol/fastq/bc2/barcode02.candidates.txt'

gtf = pd.read_csv(isofile, sep= "\t")
gtfp = read_gtf(isofile)
ens_annot = pd.read_csv(indp+"hg38.knownGene.gtf", sep= "\t", header=None)

gl = list(gtfp['gene_id'])

uniprot = [ p for p in gl if not(p.startswith("chr") or p.startswith("ENST"))]
with open(pjd+'resources/all_genes.txt', 'w+') as fp:
    fp.write('\n'.join(list(set(uniprot))))


uniprot = pd.read_csv(pjd+"resources/uniprot.tsv", sep= "\t")

up = dict(zip(uniprot['Entry'],uniprot['Gene Names']))

file = open(pjd+"resources/refseq_ensembl.txt")
annot = file.read()
file.close()

file = open(pjd+"resources/transcript_to_gene.txt")
trans_to_gene = file.read()
file.close()

file = open(pjd+"resources/ens_gensym.txt")
ens_gensym = file.read()
file.close()


counts = pd.read_csv("norm_counts.tsv", sep = "\t")

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


unknown = {}
formatted = []

for c in list(set(gtfp['gene_id'])):
    c = c.split("_")

    if c[0].startswith("ENS"):
        n = c[0].split('.')[0]
    elif c[0].startswith("chr") or '-' in c[0]:
        pass
    else:
        n = c[0]
    if n != '':
        formatted.append(n)
    

f="|".join(list(set(formatted)))



gtf.columns = ['0','1','2','3','4','5','6','7','8']
sel_gtf = gtf[gtf['8'].str.contains (f)]

a = sel_gtf['8'].str.replace('gene_id', 'gene_name')

attr_list = []
ts = []


flair = {}

for i in a:
    
    gte = str(i).split(';')
    gte[0] = gte[0].split('.')[0]

    if gte[1].startswith(' transcript_id "ENST'):
        
        t = gte[1].split('"')[1].split('.')[0]
        
        try:
            gene_id = 'gene_id "'+t_g[t]+'"'
        except:
            attr_list.append(np.nan)
            continue
            
        t =' transcript_id "'+ t +'"'
        
    else:
        t = gte[1].split('"')[1].split('.')[0]
        unk_t = t
        
        t = '0x'+t.replace('-','').upper()

        s = 'FLAIR'+str(int(t, 16))
        t = ' transcript_id "FLAIR'+str(int(t, 16))+'"'


        flair[unk_t] = s


        
        gensym = gte[0].split('gene_name "')[1][:-1]
        if gensym.startswith('chr'):
            ensg = gensym.replace(':', '.')
        else:
            transc_gene = gte[0].split('.')[0].split('"')[1]
            if 'ENSG' in transc_gene:
                ensg = transc_gene.split(' ')[1]
            if 'ENST' in transc_gene:
                try:
                    ensg = t_g[transc_gene]
                except:
                    attr_list.append(np.nan)
                    continue
            else: 
                try:
                    ensg = e_s[up[transc_gene].split(' ')[0]]
                except:
                    ensg = []
                if ensg == []:
                    attr_list.append(np.nan)
                    continue

        unk_t = ""
            
        gensym = ""
        gene_id = 'gene_id "'+ensg+'"'


    up_id = gte[0].split('gene_name "')[1][:-1]
    if up_id.startswith('ENST'):
        try:
            gene = s_e[t_g[up_id]]
        except:
            attr_list.append(np.nan)
            continue
    elif not up_id.startswith('chr'):
        try:
            gene = up[gte[0].split('gene_name "')[1][:-1]]
        except:
            attr_list.append(np.nan)
            continue
        gene = str(gene).split(' ')[0]
    else:
        gene = up[up_id]
    gene_name = ' gene_name "'+ gene +'"'
    gene_name = gene_name.replace(';', '')
    attribute = ';'.join([gene_id, t, gene_name, gte[2]]) +';'
    attribute = attribute.replace(';;', ';')
    attr_list.append(attribute)
    ts.append(t.split('"')[1])
sel_gtf['8'] = attr_list
sel_gtf = sel_gtf.dropna()
sel_gtf.to_csv("gtf.gtf", sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)

gene_ID = []
transcript_ID = []
annot_gene_id = []
annot_transcript_id = []
annot_gene_name = []
annot_transcript_name = []
n_exons = []
length = []
gene_novel = []
transcript_novel = []
ISM_subtype = []



counts = counts[~counts['transcript'].str.startswith('unassigned')]
transcripts = counts['transcript']
ens_trans = []
for t in transcripts:
    n = t.split('_')
    gid = t.split('_')[-1]
    tr = n[0]
    if tr.startswith('ENST'):
        s = tr
            
    else:

        tr = t.split('_')[0]

        tr = '0x'+tr.replace('-','').upper()
        sq = 'FLAIR'+str(int(tr, 16))
        s = sq
        

    if tr.startswith('ENS'):
        transcript_novelty = 'Known'
    else:
        transcript_novelty = 'None'
    if gid.startswith('ENSG'):
        gene_novelty = 'Known'#
        gname = s_e[gid]
    else:
        gene_novelty = 'None'
    
    
    gene_ID.append(10)
    transcript_ID.append(10)
    annot_gene_id.append(gid)
    annot_transcript_id.append(s)
    annot_gene_name.append(gname)
    annot_transcript_name.append(s)
    n_exons.append(10)
    length.append(10)
    gene_novel.append(gene_novelty)
    transcript_novel.append(transcript_novelty)
    ISM_subtype.append('None')



talon_header = pd.DataFrame()
talon_header['gene_ID'] = gene_ID
talon_header['transcript_ID'] = transcript_ID
talon_header['annot_gene_id'] = annot_gene_id
talon_header['annot_transcript_id'] = annot_transcript_id
talon_header['annot_gene_name'] = annot_gene_name
talon_header['annot_transcript_name'] = annot_transcript_name
talon_header['n_exons'] = n_exons
talon_header['length'] = length
talon_header['gene_novelty'] = gene_novel
talon_header['transcript_novelty'] = transcript_novel
talon_header['ISM_subtype'] = ISM_subtype
counts = pd.concat([talon_header, counts], axis=1)


counts = counts.drop('transcript', axis = 1)
counts = counts.dropna()
counts = counts[counts['annot_gene_id'].str.startswith('ENSG')]


counts.to_csv("abundance.tsv", sep='\t', header=True, index=False, quoting=csv.QUOTE_NONE)


ens_attr = ens_annot[8]
attr_list = []

for i in ens_attr:

    
    gte = str(i).split(';')
    if gte[1].startswith(' transcript_id "ENST'):
        t = gte[1].split('"')[1].split('.')[0]
        try:
            gene = t_g[t]
            gene_id = 'gene_id "'+ gene +'"'

        except:
            gene_id = 'Not in database'

        t =' transcript_id "'+ t +'"'
        
        if gene_id != 'Not in database':
            try:
                gensym = s_e[gene]
                if gensym == []:
                    gensym = 'Not in database'
            except:
                gensym = 'Not in database'
        else:
            gensym = 'Not in database'
        if gensym == '':
            gensym = gene
        gene_name = ' gene_name "' + gensym + '"'
        attribute = ';'.join([gene_id,t, gene_name])  + ';' + ';'.join(gte[3:-2])
        
    else:
        attribute = 'Not in database'
    attribute = attribute.replace(';;', ';')
    attr_list.append(attribute)
ens_annot[8] = attr_list
ens_annot = ens_annot[~ens_annot[8].str.contains('Not in database')]
ens_annot.to_csv("swan_annot.gtf", sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)

    #up_id = gte[0].split('gene_name "')[1][:-1]
    #if not up_id.startswith('chr'):
    #    gene = up[gte[0].split('gene_name "')[1][:-1]]
    #    gene = str(gene).split(' ')[0]
    #else:
    #    gene = up_id
    #gene_name = ' gene_name "'+ gene +'"'
    #attribute = ';'.join([gene_id, t, gene_name, gte[2]])
    #attr_list.append(attribute)
cand_file = open(cand_f)
cand = cand_file.read()
cand_file.close()

cands = cand.split('\n')

t_g['ENST00000664891'] = 'ENSG00000160808'
hr_cand = []
for c in cands:
    c = c.split('_')[-1]
    if c.startswith('ENS'):
        c = c.split('.')[0]
        hr_cand.append(s_e[c])

            
            
with open('candidates.txt', 'w') as f:
    for line in list(set(hr_cand)):
        f.write(f"{line}\n")