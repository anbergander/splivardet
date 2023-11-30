#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 10:12:37 2023

@author: denbi
"""

from gtfparse import read_gtf
import sys


  

with open(sys.argv[2]) as f:
    ex = f.readlines()

# returns GTF with essential columns such as "feature", "seqname", "start", "end"
# alongside the names of any optional keys which appeared in the attribute column
gtf = read_gtf(sys.argv[1])

# filter DataFrame to gene entries on chrY
#df_genes = df[df["feature"] == "gene"]
#df_genes_chrY = df_genes[df_genes["seqname"] == "Y"]
gene_list = []
transcript_list = []

for l in ex:

    c = l.split("_")
    ident = "_".join(c[0:-1])
    gene_bc = c[-1]
    gene, bc = gene_bc.split(" ")[0],gene_bc.split(" ")[1]
    gene_list.append(gene)
    transcript_list.append(ident)
    
    bc = bc.strip()

    #print("id:",ident,"gene:", gene_bc)
gtf_select = gtf[gtf["transcript_id"].isin(transcript_list)]  
gtf_select['sample'] = bc
#gtf_select['attributes'] = 
gtf_select.to_csv(bc+".sign.gtf",sep="\t", index=False)

